# $Id$
#
# BioPerl module for Bio::Tools::Run::RemoteBlast
#
# FORMERLY Cared for by Jason Stajich, Mat Wiepert
#
# Somewhat cared for by Roger Hall, Chris Fields (when they have time)
#
# Copyright Jason Stajich, Bioperl
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

  #change a query paramter
  $Bio::Tools::Run::RemoteBlast::HEADER{'ENTREZ_QUERY'} = 'Homo sapiens [ORGN]';

  #change a retrieval parameter
  $Bio::Tools::Run::RemoteBlast::RETRIEVALHEADER{'DESCRIPTIONS'} = 1000;

  #remove a parameter
  delete $Bio::Tools::Run::RemoteBlast::HEADER{'FILTER'};

  #$v is just to turn on and off the messages
  my $v = 1;

  my $str = Bio::SeqIO->new(-file=>'amino.fa' , -format => 'fasta' );

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
  $Bio::Tools::Run::RemoteBlast::HEADER{'MATRIX_NAME'} = 'BLOSUM45';
  $Bio::Tools::Run::RemoteBlast::HEADER{'GAPCOSTS'} = '15 2';

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

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.bioperl.org

=head1 AUTHOR 

Please do NOT contact Jason directly about this module.  Please post to
the bioperl mailing list (L<FEEDBACK>). If you would like to be the
official maintainer of this module, please volunteer on the list and
we will make it official in this POD.

First written by Jason Stajich, many others have helped keep it running.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Tools::Run::RemoteBlast;

use vars qw($AUTOLOAD $URLBASE %HEADER %RETRIEVALHEADER
	    $RIDLINE $MODVERSION %PUTPARAMS %GETPARAMS);
use strict;

use Bio::SeqIO;
use IO::String;
use Bio::Tools::BPlite;
use Bio::SearchIO;
use LWP;
use HTTP::Request::Common;

use base qw(Bio::Root::Root Bio::Root::IO);

BEGIN {
    $MODVERSION = $Bio::Root::Version::VERSION;
    $URLBASE = 'http://blast.ncbi.nlm.nih.gov/Blast.cgi';

    # In GET/PUTPARAMS the values are regexes which validate the input.
    %PUTPARAMS = (
	'AUTO_FORMAT' 	=> '(Off|(Semi|Full)auto)',	# Off, Semiauto, Fullauto
	'COMPOSITION_BASED_STATISTICS'	=> '(yes|no)',	# yes, no
	'DATABASE' 	=>  '.*',
	'DB_GENETIC_CODE' => '([1-9]|1[1-6]|2(1|2))',   # 1..16,21,22
	'ENDPOINTS'	=> '(yes|no)',			# yes,no
	'ENTREZ_QUERY'	=> '.*',
	'EXPECT'	=> '\d+(\.\d+)?([eE]-\d+)?',	# Positive double
	'FILTER'	=> '[LRm]',			# L or R or m
	'GAPCOSTS'	=> '-?\d+(\.\d+)\s+i-?\d+(\.\d+)',
					# Two space separated float values
	'GENETIC_CODE'	=> '([1-9]|1[1-6]|2(1|2))',	# 1..16,21,22
	'HITLIST_SIZE'	=> '\d+',			# Positive integer
	'I_THRESH'	=> '-?\d+(\.\d+)([eE]-\d+)?',	# float
	'LAYOUT'	=> '(One|Two)Windows?',		# onewindow, twowindows
	'LCASE_MASK'	=> '(yes|no)',			# yes, no
	'MATRIX_NAME'	=> '.*',
	'NUCL_PENALTY'	=> '-\d+',			# Negative integer
	'NUCL_REWARD'	=> '-?\d+',			# Integer
	'OTHER_ADVANCED' => '.*',
	'PERC_IDENT'	=> '\d\d+',			# Integer, 0-99 inclusive
	'PHI_PATTERN'	=> '.*',
	'PROGRAM'	=> 't?blast[pnx]',
					# tblastp, tblastn, tblastx, blastp, blastn, blastx
	'QUERY'		=> '.*',
	'QUERY_FILE'	=> '.*',
	'QUERY_BELIEVE_DEFLINE'	=> '(yes|no)',		# yes, no
	'QUERY_FROM'	=> '\d+',			# Positive integer
	'QUERY_TO'	=> '\d+',			# Positive integer
	'SEARCHSP_EFF'	=> '\d+',			# Positive integer
	'SERVICE'	=> '(plain|p[sh]i|(rps|mega)blast)',
					# plain,psi,phi,rpsblast,megablast
	'THRESHOLD'	=> '-?\d+',			# Integer
	'UNGAPPED_ALIGNMENT' => '(yes|no)',		# yes, no
	'WORD_SIZE'	=> '\d+'			# Positive integer
					  );
    %GETPARAMS = (
   'ALIGNMENTS'	=> '\d+',			# Positive integer
	'ALIGNMENT_VIEW' =>
		  '(Pairwise|(Flat)?QueryAnchored(NoIdentities)?|Tabular)',
	 # Pairwise, QueryAnchored, QueryAnchoredNoIdentities, 
  	 # FlatQueryAnchored, FlatQueryAnchoredNoIdentities, Tabular
	 'DESCRIPTIONS'	=> '\d+',			# Positive integer
	 'ENTREZ_LINKS_NEW_WINDOW' => '(yes|no)',	# yes, no
	 'EXPECT_LOW'	=> '\d+(\.\d+)?([eE]-\d+)?',	# Positive double
	 'EXPECT_HIGH'	=> '\d+(\.\d+)?([eE]-\d+)?',	# Positive double
	 'FORMAT_ENTREZ_QUERY' => '',
	 'FORMAT_OBJECT'	=> 
    '(Alignment|Neighbors|PSSM|SearchInfo|TaxBlast(Parent|MultiFrame)?)',
					# Alignment, Neighbors, PSSM,  SearchInfo 
					# TaxBlast, TaxblastParent, TaxBlastMultiFrame 
	 'FORMAT_TYPE'	=> '((HT|X)ML|ASN\.1|Text)',
					# HTML, Text, ASN.1, XML
	 'NCBI_GI'	=> '(yes|no)',			# yes, no
	 'RID' 		=>  '.*',
	 'RESULTS_FILE' 	=>  '(yes|no)',			# yes, no
	 'SERVICE' 	=>  '(plain|p[sh]i|(rps|mega)blast)',
					# plain,psi,phi,rpsblast,megablast
	 'SHOW_OVERVIEW' =>  '(yes|no)'			# yes, no
					  );

    # Default values go in here for PUT
    %HEADER = (
	       'CMD'                          => 'Put',
	       'FORMAT_OBJECT'                => 'Alignment',
	       'COMPOSITION_BASED_STATISTICS' => 'off', 
	       'DATABASE'	    	      => 'nr',
	       'EXPECT'			      => '1e-3', 
	       'FILTER'			      => 'L', 
	       'PROGRAM'		      => 'blastp', 
	       'SERVICE'		      => 'plain' 
	       );
    
    # Default values go in here for GET
    %RETRIEVALHEADER = (
			'CMD'            => 'Get',
			'ALIGNMENTS'	 => '50',
			'ALIGNMENT_VIEW' => 'Pairwise',
			'DESCRIPTIONS'	 => '100',
			'FORMAT_TYPE'	 => 'Text',
			);
    
    $RIDLINE = 'RID\s+=\s+(\S+)';
}

sub new {
	my ($caller, @args) = @_;
	# chained new
	my $self = $caller->SUPER::new(@args);
	# so that tempfiles are cleaned up
	$self->_initialize_io();
	my ($prog, $data, $readmethod, $url_base) =
        $self->_rearrange([qw(PROG DATA READMETHOD URL_BASE)],
					 @args);
	# Use these two parameters for backward-compatibility. 
	# Overridden by PROGRAM and DATABASE if supplied.
	$self->submit_parameter('PROGRAM',$prog) if $prog;
	$self->submit_parameter('DATABASE',$data) if $data;

	$readmethod = 'SearchIO' unless defined $readmethod;
	$self->readmethod($readmethod);

	# Now read the rest of the parameters and set them all

	# PUT parameters first
	my @putValues = $self->_rearrange([keys %PUTPARAMS],@args);
	my %putNames;
	@putNames{keys %PUTPARAMS} = @putValues;
	foreach my $putName (keys %putNames) {
		$self->submit_parameter($putName,$putNames{$putName});
	}
	# GET parameters second
	my @getValues = $self->_rearrange([keys %GETPARAMS],@args);
	my %getNames;
	@getNames{keys %GETPARAMS} = @getValues;
	foreach my $getName (keys %getNames) {
		$self->retrieve_parameter($getName,$getNames{$getName});
	}
        # private variable to keep track of total rids
    $self->{'_total_rids'} = 0;
    $url_base ||= $URLBASE;  # default to regular NCBI BLAST URL
    $self->set_url_base($url_base);
	return $self;
}

=head2 retrieve_parameter

 Title   : retrieve_parameter
 Usage   : my $db = $self->retrieve_parameter
 Function: Get/Set the named parameter for the retrieve_blast operation.
 Returns : string
 Args    : $name : name of GET parameter
	 $val : optional value to set the parameter to

=cut

sub retrieve_parameter {
	my ($self, $name, $val) = @_;
	$name = uc($name);
	$self->throw($name." is not a valid GET parameter.") unless
	  exists $GETPARAMS{$name};
	if (defined $val) {
    	my $regex = $GETPARAMS{$name};
    	$val =~ m/^$regex$/i or 
		  $self->throw("Value ".$val." for GET parameter ".$name." does not match expression ".$regex.". Rejecting.");
		$RETRIEVALHEADER{$name} = $val;
	}
	return $RETRIEVALHEADER{$name};
}

=head2 submit_parameter

 Title   : submit_parameter
 Usage   : my $db = $self->submit_parameter
 Function: Get/Set the named parameter for the submit_blast operation.
 Returns : string
 Args    : $name : name of PUT parameter
    $val : optional value to set the parameter to

=cut

sub submit_parameter {
    my ($self, $name, $val) = @_;
    $name = uc($name);
    $self->throw($name." is not a valid PUT parameter.") unless
	exists $PUTPARAMS{$name};
    if (defined $val) {
    	my $regex = $PUTPARAMS{$name};
    	$val =~ m/^$regex$/i or 
		$self->throw("Value ".$val." for PUT parameter ".$name." does not match expression ".$regex.". Rejecting.");
	$HEADER{$name} = $val;
    }
    return $HEADER{$name};
}

=head2 header

 Title   : header
 Usage   : my $header = $self->header
 Function: Get HTTP header for blast query
 Returns : string
 Args    : none

=cut

sub header {
    my ($self) = @_;
    return %HEADER;
}

=head2 readmethod

 Title   : readmethod
 Usage   : my $readmethod = $self->readmethod
 Function: Get/Set the method to read the blast report
 Returns : string
 Args    : string [ Blast, BPlite, blasttable, xml ]

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
 Function: Get/Set the program to run. Retained for backwards-compatibility.
 Returns : string
 Args    : string [ blastp, blastn, blastx, tblastn, tblastx ]

=cut

sub program {
    my ($self, $val) = @_;
    return $self->submit_parameter('PROGRAM',$val);
}


=head2 database

 Title   : database
 Usage   : my $db = $self->database
 Function: Get/Set the database to search. Retained for backwards-compatibility.
 Returns : string
 Args    : string [ swissprot, nr, nt, etc... ]

=cut

sub database {
    my ($self, $val) = @_;
    return $self->submit_parameter('DATABASE',$val);
}


=head2 expect

 Title   : expect
 Usage   : my $expect = $self->expect
 Function: Get/Set the E value cutoff. Retained for backwards-compatibility.
 Returns : string
 Args    : string [ '1e-4' ]

=cut

sub expect {
    my ($self, $val) = @_;
    return $self->submit_parameter('EXPECT',$val);
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
	$self->{'_ua'} = LWP::UserAgent->new(env_proxy => 1, parse_head => 0);
	my $nm = ref($self);
	$nm =~ s/::/_/g;
	$self->{'_ua'}->agent("bioperl-$nm/$MODVERSION");
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
    return if ( !defined $self->ua || !defined $protocol
		      || !defined $proxy );
    return $self->ua->proxy($protocol,$proxy);
}

sub add_rid {
    my ($self, @vals) = @_;
    foreach ( @vals ) {
	$self->{'_rids'}->{$_} = $self->{'_total_rids'};
        $self->{'_total_rids'}++; 
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
    # sort on key value, a little tricky...
    my @sort_rids = sort {$self->{'_rids'}->{$a} <=> $self->{'_rids'}->{$b}} keys %{$self->{'_rids'}};
    return @sort_rids;
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
    my $url_base = $self->get_url_base;
    return 0 unless ( @seqs );
    my $tcount = 0;
    my %header = $self->header;
    $header{$_} ||= $RETRIEVALHEADER{$_} foreach (keys %RETRIEVALHEADER);    
    foreach my $seq ( @seqs ) {
	#If query has a fasta header, the output has the query line.
	$header{'QUERY'} = ">".(defined $seq->display_id() ? $seq->display_id() : "").
		" ".(defined $seq->desc() ? $seq->desc() : "")."\n".$seq->seq();
	my $request = POST $url_base, [%header];
	$self->warn($request->as_string) if ( $self->verbose > 0);
	my $response = $self->ua->request( $request);

	if( $response->is_success ) {
	    my @subdata = split(/\n/, $response->content );
	    my $count = 0;
	    foreach ( @subdata ) {
			if( /$RIDLINE/ ) {
		    	$count++;
		    	$self->debug("RID: $1\n");
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
           Bio::SearchIO object
 Args    : Remote Blast ID (RID)

=cut

sub retrieve_blast {
    my($self, $rid) = @_;
    my ($fh,$tempfile) = $self->tempfile();
    close $fh;			#explicit close
    my $url_base = $self->get_url_base;
    my %hdr = %RETRIEVALHEADER;
    $hdr{'RID'} = $rid;
    my $req = POST $url_base, [%hdr];
    $self->debug("retrieve request is " . $req->as_string());
    my $response = $self->ua->request($req, $tempfile);
    if( $response->is_success ) {
    	if( $self->verbose > 0 ) {
            #print content of reply if verbose > 1
            open(my $TMP, $tempfile) or $self->throw("cannot open $tempfile");
            while(<$TMP>) { print $_; }
            
    	}
        ## if proper reply 
        open(my $TMP, $tempfile) || $self->throw("Error opening $tempfile");
        my $waiting = 1;
        my $s = 0;
        my $got_content = 0;
        while(<$TMP>) {
            if (/./) {
                $got_content = 1;
            }
            if( /<\?xml version=/ ) { # xml time
                $waiting = 0;
                last;
            }
            if( /QBlastInfoBegin/i ) {
                $s = 1;
            } elsif( $s ) {
                if( /Status=(WAITING|ERROR|READY)/i ) {
                    if( $1 eq 'WAITING' ) {
                        $waiting = 1;
                    } elsif( $1 eq 'ERROR' ) {
                        close($TMP);
                        open(my $ERR, "<$tempfile") or $self->throw("cannot open file $tempfile");
                        $self->warn(join("", <$ERR>));
                        close $ERR;
                        return -1;
                    } elsif( $1 eq 'READY' ) {
                        $waiting = 0;
                        last;
                    } else {
                        $self->warn("Unknown status $1:\n");
                        last;
                    }
                }
            } elsif( /^(?:#\s)?[\w-]*?BLAST\w+/ ) {
                $waiting = 0;
                last;
            } elsif ( /ERROR/i ) {
                close($TMP);
                open(my $ERR, "<$tempfile") or $self->throw("cannot open file $tempfile");
                $self->warn(join("", <$ERR>));
                close $ERR;
                return -1;
            }
            
        }
        close($TMP);
        if( ! $waiting ) {
            my $blastobj;
            my $mthd = $self->readmethod;
            if( $mthd =~ /blasttable/i ) {
            # pre-process
            my ($fh2,$tempfile2) = $self->tempfile();
            open(my $TMP,$tempfile) || $self->throw($!);
            my $s = 0;
            while(<$TMP>) {
                if(/\<PRE\>/i ) {
                $s = 1;
                } elsif( /\<\/PRE\>/i ) {
                $s = 0;
                last;
                } elsif( $s ) {
                print $fh2 $_;
                }
            } 
            close($fh2);
            $blastobj = Bio::SearchIO->new( -file => $tempfile2,
                               -format => 'blasttable');
            } elsif( $mthd =~ /xml/ ) {
            $blastobj = Bio::SearchIO->new( -file => $tempfile,
                               -format => 'blastxml');
            } else {
            $blastobj = Bio::SearchIO->new( -file => $tempfile,
                               -format => 'blast');
            } 
            
            ## store filename in object ##
            $self->file($tempfile);
            return $blastobj;
        } elsif (!$got_content) {
            # server returned no content, can't be good
            $self->warn("Server failed to return any data");
            return -1
        } else {		# still working
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
	my $blastfile = $self->file;
	#open temp file and output file, have to filter out some HTML
	open(my $TMP, $blastfile) or $self->throw("cannot open $blastfile");

	open(my $SAVEOUT, ">", $filename) or $self->throw("cannot open $filename");
	my $seentop = 0;
	while(<$TMP>) {
		next if (/<pre>/);
		if(/^(?:[T]?BLAST[NPX])\s*.+$/i ||
           /^RPS-BLAST\s*.+$/i ||
           /<\?xml\sversion=/ ||
           /^#\s+(?:[T]?BLAST[NPX])\s*.+$/) {
			$seentop=1;
		} 
        next if !$seentop;
		if( $seentop ) {
			print $SAVEOUT $_;
		}
	}
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
			my $seqio = Bio::SeqIO->new(-format => 'fasta',
												-file => $input);
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

=head2 set_url_base

 Title   : set_url_base
 Usage   : $self->set_url_base($url)
 Function: Method to override the default NCBI BLAST database
 Returns : None
 Args    : string (database url like
 NOTE    : This is highly experimental; we cannot maintain support on
           databases other than the default NCBI database at this time

=cut

sub set_url_base {
    my $self = shift;
    $self->{'_urlbase'} = shift if @_;
}

=head2 get_url_base

 Title   : get_url_base
 Usage   : my $url = $self->set_url_base
 Function: Get the current URL for BLAST database searching
 Returns : string (URL used for remote blast searches)
 Args    : None

=cut

sub get_url_base {
    my $self = shift;
    return $self->{'_urlbase'};
}

__END__
