# $Id$
#
# BioPerl module for Bio::Tools::Run::RemoteBlast
#
# Cared for by Jason Stajich, Mat Wiepert
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::RemoteBlast - Object for remote execution of the NCBI Blast
via HTTP

=head1 SYNOPSIS

  #Remote-blast "factory object" creation and blast-parameter initialization

  use Bio::Tools::Run::RemoteBlast;
  use strict;
  my $prog = 'blastp';
  my $db   = 'swissprot';
  my $e_val= '1e-10';

  my @params = ( '-prog' => $prog,
         '-data' => $db,
         '-expect' => $e_val,
         '-readmethod' => 'SearchIO' );

  my $factory = Bio::Tools::Run::RemoteBlast->new(@params);

  #change a paramter
  $Bio::Tools::Run::RemoteBlast::HEADER{'ENTREZ_QUERY'} = 'Homo sapiens [ORGN]';

  #remove a parameter
  delete $Bio::Tools::Run::RemoteBlast::HEADER{'FILTER'};

  my $v = 1;
  #$v is just to turn on and off the messages
  
  my $str = Bio::SeqIO->new(-file=>'amino.fa' , '-format' => 'fasta' );

  while (my $input = $str->next_seq()){
    #Blast a sequence against a database:

    #Alternatively, you could  pass in a file with many
    #sequences rather than loop through sequence one at a time
    #Remove the loop starting 'while (my $input = $str->next_seq())'
    #and swap the two lines below for an example of that.
    my $r = $factory->submit_blast($input);
    #my $r = $factory->submit_blast('amino.fa');

    print STDERR "waiting..." if( $v > 0 );
    while ( my @rids = $factory->each_rid ) {
      foreach my $rid ( @rids ) {
        my $rc = $factory->retrieve_blast($rid);
        if( !ref($rc) ) {
          if( $rc < 0 ) {
            $factory->remove_rid($rid);
          }
          print STDERR "." if ( $v > 0 );
          sleep 5;
        } else {
          my $result = $rc->next_result();
          #save the output
          my $filename = $result->query_name()."\.out";
          $factory->save_output($filename);
          $factory->remove_rid($rid);
          print "\nQuery Name: ", $result->query_name(), "\n";
          while ( my $hit = $result->next_hit ) {
            next unless ( $v > 0);
            print "\thit name is ", $hit->name, "\n";
            while( my $hsp = $hit->next_hsp ) {
              print "\t\tscore is ", $hsp->score, "\n";
            }
          }
        }
      }
    }
  }

  # This example shows how to change a CGI parameter:
  $Bio::Tools::Run::RemoteBlast::HEADER{'MATRIX_NAME'} = 'BLOSUM25';

  # And this is how to delete a CGI parameter:
  delete $Bio::Tools::Run::RemoteBlast::HEADER{'FILTER'};


=head1 DESCRIPTION

Class for remote execution of the NCBI Blast via HTTP.

For a description of the many CGI parameters see:
http://www.ncbi.nlm.nih.gov/BLAST/Doc/urlapi.html

Various additional options and input formats are available.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org               - General discussion
  http://bio.perl.org/MailList.html   - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR -  Jason Stajich

Email jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Tools::Run::RemoteBlast;

use vars qw($AUTOLOAD @ISA %BLAST_PARAMS $URLBASE %HEADER %RETRIEVALHEADER
	    $RIDLINE);
use strict;

use Bio::Root::Root;
use Bio::Root::IO;
use Bio::SeqIO;
use IO::String;
use Bio::Tools::BPlite;
use Bio::SearchIO;
use LWP;
use HTTP::Request::Common;
BEGIN {
    $URLBASE = 'http://www.ncbi.nlm.nih.gov/blast/Blast.cgi';
    %HEADER = ('CMD'                          => 'Put',
	       'PROGRAM'                      => '',
	       'DATABASE'                     => '',
	       'FILTER'                       => 'L',
	       'EXPECT'                       => '',
	       'QUERY'                        =>  '',
	       'CDD_SEARCH'                   => 'off',
	       'COMPOSITION_BASED_STATISTICS' => 'off',
	       'FORMAT_OBJECT'                => 'Alignment',
	       'SERVICE'                      => 'plain',
	       );

    %RETRIEVALHEADER = ('CMD'            => 'Get',
			'RID'            => '',
			'ALIGNMENT_VIEW' => 'Pairwise',
			'DESCRIPTIONS'   => 100,
			'ALIGNMENTS'     => 50,
			'FORMAT_TYPE'    => 'Text',
			);

    $RIDLINE = 'RID\s+=\s+(\d+-\d+-\d+)';

    %BLAST_PARAMS = ( 'prog' => 'blastp',
		       'data' => 'nr',
		       'expect' => '1e-3',
		       'readmethod' => 'SearchIO'
		       );

}

@ISA = qw(Bio::Root::Root Bio::Root::IO);

sub new {
    my ($caller, @args) = @_;
    # chained new
    my $self = $caller->SUPER::new(@args);
    # so that tempfiles are cleaned up
    $self->_initialize_io();
    my ($prog, $data, $expect,
	$readmethod) = $self->_rearrange([qw(PROG DATA
					     EXPECT
					     READMETHOD)],
					 @args);

    $readmethod = $BLAST_PARAMS{'readmethod'} unless defined $readmethod;
    $prog = $BLAST_PARAMS{'prog'}     unless defined $prog;
    $data = $BLAST_PARAMS{'data'}     unless defined $data;
    $expect = $BLAST_PARAMS{'expect'} unless defined $expect;
    $self->readmethod($readmethod);
    $self->program($prog);
    $self->database($data);
    $self->expect($expect);

    return $self;
}

=head2 header

 Title   : header
 Usage   : my $header = $self->header
 Function: Get/Set HTTP header for blast query
 Returns : string
 Args    : none

=cut

sub header {
    my ($self) = @_;
    my %h = %HEADER;
    $h{'PROGRAM'} = $self->program;
    $h{'DATABASE'} = $self->database;
    $h{'EXPECT'}  = $self->expect;
    return %h;
}

=head2 readmethod

 Title   : readmethod
 Usage   : my $readmethod = $self->readmethod
 Function: Get/Set the method to read the blast report
 Returns : string
 Args    : string [ Blast, BPlite ]

=cut

sub readmethod {
    my ($self, $val) = @_;
    if( defined $val ) {
	$self->{'_readmethod'} = $val;
    }
    return $self->{'_readmethod'};
}


=head2 program

 Title   : program
 Usage   : my $prog = $self->program
 Function: Get/Set the program to run
 Returns : string
 Args    : string [ blastp, blastn, blastx, tblastn, tblastx ]

=cut

sub program {
    my ($self, $val) = @_;
    if( defined $val ) {
	$val = lc $val;
	if( $val !~ /t?blast[pnx]/ ) {
	    $self->warn("trying to set program to an invalid program name ($val) -- defaulting to blastp");
	    $val = 'blastp';
	}
#	$self->{'_program'} = $val;
	$HEADER{'PROGRAM'} = $val;
    }
    return $HEADER{'PROGRAM'};
}


=head2 database

 Title   : database
 Usage   : my $db = $self->database
 Function: Get/Set the database to search
 Returns : string
 Args    : string [ swissprot, nr, nt, etc... ]

=cut

sub database {
    my ($self, $val) = @_;
    if( defined $val ) {
#	$self->{'_database'} = $val;
 	$HEADER{'DATABASE'} = $val;
    }
    return $HEADER{'DATABASE'};
}


=head2 expect

 Title   : expect
 Usage   : my $expect = $self->expect
 Function: Get/Set the E value cutoff
 Returns : string
 Args    : string [ '1e-4' ]

=cut

sub expect {
    my ($self, $val) = @_;
    if( defined $val ) {
#	$self->{'_expect'} = $val;
 	$HEADER{'EXPECT'} = $val;
    }
    return $HEADER{'EXPECT'};
}

=head2 ua

 Title   : ua
 Usage   : my $ua = $self->ua or
           $self->ua($ua)
 Function: Get/Set a LWP::UserAgent for use
 Returns : reference to LWP::UserAgent Object
 Args    : none
 Comments: Will create a UserAgent if none has been requested before.

=cut

sub ua {
    my ($self, $value) = @_;
    if( ! defined $self->{'_ua'} ) {
	$self->{'_ua'} = new LWP::UserAgent;
    }
    return $self->{'_ua'};
}

=head2 proxy

 Title   : proxy
 Usage   : $httpproxy = $db->proxy('http')  or
           $db->proxy(['http','ftp'], 'http://myproxy' )
 Function: Get/Set a proxy for use of proxy
 Returns : a string indicating the proxy
 Args    : $protocol : an array ref of the protocol(s) to set/get
           $proxyurl : url of the proxy to use for the specified protocol

=cut

sub proxy {
    my ($self,$protocol,$proxy) = @_;
    return undef if ( !defined $self->ua || !defined $protocol
		      || !defined $proxy );
    return $self->ua->proxy($protocol,$proxy);
}

sub add_rid {
    my ($self, @vals) = @_;
    foreach ( @vals ) {
	$self->{'_rids'}->{$_} = 1;
    }
    return scalar keys %{$self->{'_rids'}};
}

sub remove_rid {
    my ($self, @vals) = @_;
    foreach ( @vals ) {
	delete $self->{'_rids'}->{$_};
    }
    return scalar keys %{$self->{'_rids'}};
}

sub each_rid {
    my ($self) = @_;
    return keys %{$self->{'_rids'}};
}

=head2 submit_blast

 Title   : submit_blast
 Usage   : $self->submit_blast([$seq1,$seq2]);
 Function: Submit blast jobs to ncbi blast queue on sequence(s)
 Returns : Blast report object as defined by $self->readmethod
 Args    : input can be:
           * sequence object
           * array ref of sequence objects
           * filename of file containing fasta formatted sequences

=cut

sub submit_blast {
    my ($self, $input) = @_;
    my @seqs = $self->_load_input($input);
    return 0 unless ( @seqs );
    my $tcount = 0;
    my %header = $self->header;
    foreach my $seq ( @seqs ) {
	#If query has a fasta header, the output has the query line.
	$header{'QUERY'} = ">".(defined $seq->display_id() ? $seq->display_id() : "").
		" ".(defined $seq->desc() ? $seq->desc() : "")."\n".$seq->seq();
	my $request = POST $URLBASE, [%header];
	$self->warn($request->as_string) if ( $self->verbose > 0);
	my $response = $self->ua->request( $request);

	if( $response->is_success ) {
	    if( $self->verbose > 0 ) {
		my ($tempfh) = $self->tempfile();
		# Hmm, what exactly are we trying to do here?
		print $tempfh $response->content;
		close($tempfh);
		undef $tempfh;
	    }
	    my @subdata = split(/\n/, $response->content );
	    my $count = 0;
	    foreach ( @subdata ) {
		if( /$RIDLINE/ ) {
		    $count++;
		    print STDERR $_ if( $self->verbose > 0);
		    $self->add_rid($1);		
		    last;
		}	
	    }
	    if( $count == 0 ) {
		$self->warn("req was ". $request->as_string() . "\n");
		$self->warn(join('', @subdata));
	    }    	
	    $tcount += $count;
	} else {
	    # should try and be a little more verbose here
	    $self->warn("req was ". $request->as_string() . "\n" .
			$response->error_as_HTML);
	    $tcount = -1;
	}
    }
    return $tcount;
}

=head2 retrieve_blast

 Title   : retrieve_blast
 Usage   : my $blastreport = $blastfactory->retrieve_blast($rid);
 Function: Attempts to retrieve a blast report from remote blast queue
 Returns : -1 on error,
           0 on 'job not finished',
           Bio::Tools::BPlite or Bio::Tools::Blast object
           (depending on how object was initialized) on success
 Args    : Remote Blast ID (RID)

=cut

sub retrieve_blast {
    my($self, $rid) = @_;
    my (undef,$tempfile) = $self->tempfile();
    my %hdr = %RETRIEVALHEADER;
    $hdr{'RID'} = $rid;
    my $req = POST $URLBASE, [%hdr];
    if( $self->verbose > 0 ) {
	$self->warn("retrieve request is " . $req->as_string());
    }
    my $response = $self->ua->request($req, $tempfile);
    if( $self->verbose > 0 ) {
	open(TMP, $tempfile) or $self->throw("cannot open $tempfile");
	while(<TMP>) { print $_; }
	close TMP;
    }
    if( $response->is_success ) {	
	my $size = -s $tempfile;
	if( $size > 1000 ) {
	    my $blastobj;
	    if( $self->readmethod =~ /BPlite/ ) {
		$blastobj = new Bio::Tools::BPlite(-file => $tempfile);
	    } else {
		$blastobj = new Bio::SearchIO(-file => $tempfile,
					      -format => 'blast');
	    }
	    #save tempfile
	    $self->file($tempfile);
	    return $blastobj;
	} elsif( $size < 500 ) { # search had a problem
	    open(ERR, "<$tempfile") or $self->throw("cannot open file $tempfile");
	    $self->warn(join("", <ERR>));
	    close ERR;
	    return -1;
	} else { # still working
	    return 0;
	}
    } else {
	$self->warn($response->error_as_HTML);
	return -1;
    }
}

=head2 save_output

 Title   : saveoutput
 Usage   : my $saveoutput = $self->save_output($filename)
 Function: Method to save the blast report
 Returns : 1 (throws error otherwise)
 Args    : string [rid, filename]

=cut

sub save_output {
    my ($self, $filename) = @_;
    if( ! defined $filename ) {
        	$self->throw("Can't save blast output.  You must specify a filename to save to.");
    }
    #should be set when retrieving blast
   	my $blastfile = $self->file;
   	#open temp file and output file, have to filter out some HTML
	open(TMP, $blastfile) or $self->throw("cannot open $blastfile");
	open(SAVEOUT, ">$filename") or $self->throw("cannot open $filename");
	my $seentop=0;
	while(<TMP>) {
		next if (/<pre>/);	
		if( /^(?:[T]?BLAST[NPX])\s*.+$/i ||
	   		/^RPS-BLAST\s*.+$/i ) {
	   		$seentop=1;
	   	}
	   	next if !$seentop;
	    if( $seentop ) {
			print SAVEOUT;
		}
	}
	close SAVEOUT;
	close TMP;
	return 1;	
}

sub _load_input {
    my ($self, $input) = @_;

    if( ! defined $input ) {
	$self->throw("Calling remote blast with no input");	
    }
    my @seqs;
    if( ! ref $input ) {
	if( -e $input ) {
	    my $seqio = new Bio::SeqIO(-format => 'fasta', -file => $input);
	    while( my $seq = $seqio->next_seq ) {
		push @seqs, $seq;
	    }
	} else {
	    $self->throw("Input $input was not a valid filename");
	}	
    } elsif( ref($input) =~ /ARRAY/i ) {
	foreach ( @$input ) {
	    if( ref($_) && $_->isa('Bio::PrimarySeqI') ) {
		push @seqs, $_;
	    } else {
		$self->warn("Trying to add a " . ref($_) .
			    " but expected a Bio::PrimarySeqI");
	    }
	}
	if( ! @seqs) {
	    $self->throw("Did not pass in valid input -- no sequence objects found");
	}
    } elsif( $input->isa('Bio::PrimarySeqI') ) {
	push @seqs, $input;
    }
    return @seqs;
}
1;
__END__
