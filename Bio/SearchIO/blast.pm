# $Id$
#
# BioPerl module for Bio::SearchIO::blast
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::blast - Event generator for event based parsing of blast reports 

=head1 SYNOPSIS

#Do not use this object directly - it is used as part of the Bio::SearchIO system.

    use Bio::SearchIO;
    my $searchio = new Bio::SearchIO(-format => 'blast',
				     -file   => 'file1.bls');
    while( my $report = $searchio->next_report ) {
    
    }

=head1 DESCRIPTION

This object encapsulated the necessary methods for generating events suitable for building Bio::Search objects from a BLAST report file. 

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...



package Bio::SearchIO::blast;
use strict;
use vars qw(@ISA %MAPPING %MODEMAP);
use Bio::SearchIO;

@ISA = qw(Bio::SearchIO );

BEGIN { 
    # mapping of NCBI Blast terms to Bioperl hash keys
    %MODEMAP = ('BlastOutput' => 'report',
		'Hit'         => 'subject',
		'Hsp'         => 'hsp'
		);

    # This should really be done more intelligently, like with
    # XSLT

    %MAPPING = ( 
		 'Hsp_bit-score' => 'bits',
		 'Hsp_score'     => 'score',
		 'Hsp_evalue'    => 'evalue',
		 'Hsp_query-from'=> 'querystart',
		 'Hsp_query-to'  => 'queryend',
		 'Hsp_hit-from'  => 'subjectstart',
		 'Hsp_hit-to'    => 'subjectend',
		 'Hsp_positive'  => 'positive',
		 'Hsp_identity'  => 'match',
		 'Hsp_gaps'      => 'gaps',
		 'Hsp_qseq'      => 'queryseq',
		 'Hsp_hseq'      => 'subjectseq',
		 'Hsp_midline'   => 'homolseq',
		 'Hsp_align-len' => 'hsplen',
		 'Hsp_query-frame'=> 'queryframe',
		 'Hsp_hit-frame'  => 'subjectframe',

		 'Hit_id'        => 'subjectname',
		 'Hit_len'       => 'subjectlen',
		 'Hit_accession' => 'subjectacc',
		 'Hit_def'       => 'subjectdesc',
		 
		 'BlastOutput_program'  => 'programname',
		 'BlastOutput_version'  => 'programver',
		 'BlastOutput_query-def'=> 'queryname',
		 'BlastOutput_query-len'=> 'querylen',
		 'BlastOutput_db'       => 'dbname',
		 'BlastOutput_db-len'   => 'dbsize',
		 'Iteration_iter-num'   => 'iternum',
		 'Parameters_matrix'    => { 'param' => 'matrix'},
		 'Parameters_expect'    => { 'param' => 'expect'},
		 'Parameters_include'   => { 'param' => 'include'},
		 'Parameters_sc-match'  => { 'param' => 'match'},
		 'Parameters_sc-mismatch' => { 'param' => 'mismatch'},
		 'Parameters_gap-open'  => { 'param' => 'gapopen'},
		 'Parameters_gap-extend'=> { 'param' => 'gapext'},
		 'Parameters_filter'    => {'param' => 'filter'},
		 'Statistics_db-num'    => { 'stat' => 'dbnum'},
		 'Statistics_db-len'    => { 'stat' => 'dblength'},
		 'Statistics_hsp-len'   => { 'stat' => 'hsplength'},
		 'Statistics_eff-space' => { 'stat' => 'effectivespace'},
		 'Statistics_kappa'     => { 'stat' => 'kappa' },
		 'Statistics_lambda'    => { 'stat' => 'lambda' },
		 'Statistics_entropy'   => { 'stat' => 'entropy'},
		 );
}


=head2 new

 Title   : new
 Usage   : my $obj = new Bio::SearchIO::blast();
 Function: Builds a new Bio::SearchIO::blast object 
 Returns : Bio::SearchIO::blast
 Args    :

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
}


=head2 next_report

 Title   : next_report
 Usage   : my $subject = $searchio->next_report;
 Function: Returns the next Report from a search
 Returns : Bio::Search::ReportI object
 Args    : none

=cut

sub next_report{
   my ($self) = @_;
   
   my $data = '';
   my $seentop = 0;
   my $reporttype;
   $self->start_document();
   while( defined ($_ = $self->_readline )) {
       next if( /^\s+$/); # skip empty lines
       if( /^([T]?BLAST[NPX])\s*(\S+)/i ) {
	   if( $seentop ) {
	       $self->_pushback($_);
	       $self->end_element({ 'Name' => 'BlastOutput'});
	       return $self->end_document();
	   }
	   $self->start_element({ 'Name' => 'BlastOutput' } );
	   $seentop = 1;
	   $reporttype = $1;
	   $self->element({ 'Name' => 'BlastOutput_program',
			    'Data' => $reporttype});
	   
	   $self->element({ 'Name' => 'BlastOutput_version',
			    'Data' => $2});
       } elsif ( /^Query=\s*(.+)$/ ) {
	   my $q = $1;
	   my $size = 0;      
	   $_ = $self->_readline;
	   while( defined ($_) && $_ !~ /^\s+$/ ) {	       
	       if( /([\d,]+)\s+letters/ ) {		   
		   $size = $1;
		   $size =~ s/,//g;
		   last;
	       } else { 
		   $q .= $_;
	       }
	       $_ = $self->_readline;
	   }
	   chomp($q);
	   $self->element({ 'Name' => 'BlastOutput_query-def',
			    'Data' => $q});
	   $self->element({ 'Name' => 'BlastOutput_query-len', 
			    'Data' => $size});

       } elsif ( /^Database:\s*(.+)$/ ) {
	   $self->element({'Name' => 'BlastOutput_db',
			   'Data' => $1});
	   if( $self->_readline =~ /([\d,]+) total letters/ ) {
	       my $s = $1;
	       $s =~ s/,//g;
	       $self->element({'Name' => 'BlastOutput_db-len',
			       'Data' => $s});	
	   }
       } elsif( /^>(\S+)\s*(.*)?/ ) {
	   if( $self->in_element('hsp') ) {
	       $self->end_element({ 'Name' => 'Hsp'});
	   }
	   if( $self->in_element('subject') ) {
	       $self->end_element({ 'Name' => 'Hit'});
	   }
	   $self->start_element({ 'Name' => 'Hit'});
	   my $id = $1;
	   $self->element({ 'Name' => 'Hit_id',
			    'Data' => $id});
	   my @pieces = split(/\|/,$id);
	   my $acc = pop @pieces;
	   $acc =~ s/\.\d+$//;
	   $self->element({ 'Name' =>  'Hit_accession',
			    'Data'  => $acc});	   
	   $self->element({ 'Name' => 'Hit_def',
			    'Data' => $2});	   
       } elsif( $self->_mode eq 'subject' && /Length\s*=\s*(\d+)/ ) {
	   $self->element({ 'Name' => 'Hit_len',
			    'Data' => $1 });
       } elsif( /Score = (\d+) \((\S+) bits\), Expect = ([^,\s]+),/ ) {
	   if( $self->in_element('hsp') ) {
	       $self->end_element({'Name' => 'Hsp'});
	   }
	   $self->start_element({'Name' => 'Hsp'});
       	   $self->element( { 'Name' => 'Hsp_score',
			     'Data' => $1});
	   $self->element( { 'Name' => 'Hsp_bit-score',
			     'Data' => $2});
	   $self->element( { 'Name' => 'Hsp_evalue',
			     'Data' => $3});
       } elsif( /Score = (\S+) bits \((\d+)\), Expect = (\S+)/	) {
	   if( $self->in_element('hsp') ) {
	       $self->end_element({ 'Name' => 'Hsp'});
	   }
	   $self->start_element({'Name' => 'Hsp'});
	   $self->element( { 'Name' => 'Hsp_score',
			     'Data' => $2});
	   $self->element( { 'Name' => 'Hsp_bit-score',
			     'Data' => $1});
	   $self->element( { 'Name' => 'Hsp_evalue',
			     'Data' => $3});
       } elsif( $self->in_element('hsp') &&
		/Identities = (\d+)\/(\d+).+, Positives = (\d+)\/(\d+).+(\, Gaps = (\d+)\/(\d+))?/ ) {
	   
	   $self->element( { 'Name' => 'Hsp_identity',
			     'Data' => $1});
	   $self->element( {'Name' => 'Hsp_align-len',
			    'Data' => $2});
	   $self->element( { 'Name' => 'Hsp_positive',
			     'Data' => $3});
	   if( defined $6 ) { 
	       $self->element( { 'Name' => 'Hsp_gaps',
				 'Data' => $6});	   
	   }
	   $self->{'_Query'} = {'begin' => 0, 'end' => 0};
	   $self->{'_Sbjct'} = { 'begin' => 0, 'end' => 0};
       } elsif( $self->in_element('hsp') &&
		/Frame = ([\+\-][1-3])\s*(\/\s*([\+\-][1-3]))/ ){
	   
	   my ($queryframe,$hitframe);
	   if( $reporttype eq 'TBLASTX' ) {
	       ($queryframe,$hitframe) = ($1,$2);
	       $hitframe =~ s/\/\s*//g;
	   } elsif( $reporttype eq 'TBLASTN' ) {
	       ($hitframe,$queryframe) = ($1,0);
	   } else { 
	       ($queryframe,$hitframe) = ($1,0);
	   }
	   $self->element({'Name' => 'Hsp_query-frame',
			   'Data' => $queryframe});
	   	   
	   $self->element({'Name' => 'Hsp_hit-frame',
			   'Data' => $hitframe});
       } elsif(  /^Parameters:/ || /\s+Database:/ ) {
	   $self->end_element({'Name' => 'Hsp'});
	   $self->end_element({'Name' => 'Hit'});
	   my $blast = ( /Parameters/ ) ? 'wublast'  : 'ncbi'; 
	   my $last = '';
	   while( defined ($_ = $self->_readline ) ) {
	       if( /^([T]?BLAST[NPX])\s*([\d\.]+)/i ) {
		   $self->_pushback($_);
		   # let's handle this in the loop
		   last;
	       }
	       # here is where difference between wublast and ncbiblast
	       # is better handled by different logic
	       if( /Number of Sequences:\s+(\d+)/ ||
			/of sequences in database:\s+(\d+)/) {
		   $self->element({'Name' => 'Statistics_db-num',
				   'Data' => $1});
	       } elsif ( /letters in database:\s+([\d,]+)/) {	   
		   my $s = $1;
		   $s =~ s/,//g;
		   $self->element({'Name' => 'Statistics_db-len',
				   'Data' => $s});
	       } elsif( $blast eq 'wublast' ) {
		   if( /E=(\S+)/ ) {
		       $self->element({'Name' => 'Parameters_expect',
				       'Data' => $1});
		   } elsif( $last =~ /Frame\s+MatID\s+Matrix name/ ) {
		       s/^\s+//;
		       my (undef,undef,$matrix,$lambda,$kappa,$entropy) = split;
		       $self->element({'Name' => 'Parameters_matrix',
				       'Data' => $matrix});
		       $self->element({'Name' => 'Statistics_lambda',
				       'Data' => $lambda});
		       $self->element({'Name' => 'Statistics_kappa',
				       'Data' => $kappa});
		       $self->element({'Name' => 'Statistics_entropy',
				       'Data' => $entropy});
		   } 
	       } elsif ( $blast eq 'ncbi' ) {
		   if( /^Matrix:\s+(\S+)/ ) {
		       $self->element({'Name' => 'Parameters_matrix',
				       'Data' => $1});		       
		   } elsif( /Lambda/ ) {
		       $_ = $self->_readline;
		       s/^\s+//;
		       my ($lambda, $kappa, $entropy) = split;
		       $self->element({'Name' => 'Statistics_lambda',
				       'Data' => $lambda});
		       $self->element({'Name' => 'Statistics_kappa',
				       'Data' => $kappa});
		       $self->element({'Name' => 'Statistics_entropy',
				       'Data' => $entropy});
		   } elsif( /effective search space: (\d+)/ ) {
		       $self->element({'Name' => 'Statistics_eff-space',
				       'Data' => $1});
		   } elsif( /Gap Penalties:\s+Existence:\s+(\d+),\s+Extension:\s+(\d+)/) {
		       $self->element({'Name' => 'Parameters_gap-open',
				       'Data' => $1});
		       $self->element({'Name' => 'Parameters_gap-extend',
				       'Data' => $2});
		   } elsif( /effective HSP length:\s+(\d+)/ ) {
		        $self->element({'Name' => 'Statistics_hsp-len',
					'Data' => $1});
		   }
	       }
	       $last = $_;
	   }
       } elsif( $self->in_element('hsp') ) {
	   
           # let's read 3 lines at a time;
	   my %data = ( 'Query' => '',
			'Mid' => '',
			'Hit' => '' );
	   my $len;
	   for( my $i = 0; 
		defined($_) && $i < 3; 
		$i++ ){
	       chomp;
	       if( /((Query|Sbjct):\s+(\d+)\s+)(\S+)\s+(\d+)/ ) {
		   $data{$2} = $4;
		   $len = length($1);
		   $self->{"\_$2"}->{'begin'} = $3 unless $self->{"_$2"}->{'begin'};
		   $self->{"\_$2"}->{'end'} = $5;
	       } else { 
		   $data{'Mid'} = substr($_,$len);
	       }
	       $_ = $self->_readline();	       
	   }
	   $self->characters({'Name' => 'Hsp_qseq',
			      'Data' => $data{'Query'} });
	   $self->characters({'Name' => 'Hsp_hseq',
			      'Data' => $data{'Sbjct'}});
	   $self->characters({'Name' => 'Hsp_midline',
			      'Data' => $data{'Mid'} });
       } else { 
	  # print "unrecognized line $_";
       }
   } 
   $self->end_element({'Name' => 'BlastOutput'});
   return $self->end_document();
}

=head2 start_element

 Title   : start_element
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub start_element{
   my ($self,$data) = @_;
    # we currently don't care about attributes
    my $nm = $data->{'Name'};    
    if( my $type = $MODEMAP{$nm} ) {
	$self->_mode($type);
	if( $self->_eventHandler->will_handle($type) ) {
	    my $func = sprintf("start_%s",lc $type);
	    $self->_eventHandler->$func($data->{'Attributes'});
	}						 
	unshift @{$self->{'_elements'}}, $type;
    }
    if($nm eq 'BlastOutput') {
	$self->{'_values'} = {};
	$self->{'_report'}= undef;
	$self->{'_mode'} = '';
    }

}

=head2 end_element

 Title   : end_element
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub end_element {
    my ($self,$data) = @_;
    my $nm = $data->{'Name'};
    my $rc;
    if($nm eq 'BlastOutput_program' &&
       $self->{'_last_data'} =~ /(t?blast[npx])/i ) {
	$self->{'_reporttype'} = uc $1; 	    
    }   
    # Hsp are sort of weird, in that they end when another
    # object begins so have to detect this in end_element for now
    if( $nm eq 'Hsp' ) {
	foreach ( qw(Hsp_qseq Hsp_midline Hsp_hseq) ) {
	    $self->element({'Name' => $_,
			    'Data' => $self->{'_last_hspdata'}->{$_}});
	}
	$self->{'_last_hspdata'} = {};
	$self->element({'Name' => 'Hsp_query-from',
			'Data' => $self->{'_Query'}->{'begin'}});
	$self->element({'Name' => 'Hsp_query-to',
			'Data' => $self->{'_Query'}->{'end'}});
	
	$self->element({'Name' => 'Hsp_hit-from',
			'Data' => $self->{'_Sbjct'}->{'begin'}});
	$self->element({'Name' => 'Hsp_hit-to',
			'Data' => $self->{'_Sbjct'}->{'end'}});
    }
    if( my $type = $MODEMAP{$nm} ) {
	if( $self->_eventHandler->will_handle($type) ) {
	    my $func = sprintf("end_%s",lc $type);
	    $rc = $self->_eventHandler->$func($self->{'_reporttype'},
					      $self->{'_values'});	    
	}
	shift @{$self->{'_elements'}};

    } elsif( $MAPPING{$nm} ) { 	
	
	if ( ref($MAPPING{$nm}) =~ /hash/i ) {
	    my $key = (keys %{$MAPPING{$nm}})[0];	    
	    $self->{'_values'}->{$key}->{$MAPPING{$nm}->{$key}} = $self->{'_last_data'};
	} else {
	    $self->{'_values'}->{$MAPPING{$nm}} = $self->{'_last_data'};
	}
    } else { 
	print "unknown nm $nm, ignoring\n";
    }
    $self->{'_last_data'} = ''; # remove read data if we are at 
				# end of an element
    $self->{'_report'} = $rc if( $nm eq 'BlastOutput' );
    return $rc;

}

=head2 element

 Title   : element
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub element{
   my ($self,$data) = @_;
   $self->start_element($data);
   $self->characters($data);
   $self->end_element($data);
}


=head2 characters

 Title   : characters
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub characters{
   my ($self,$data) = @_;   

   if( $self->in_element('hsp') && 
       $data->{'Name'} =~ /Hsp\_(qseq|hseq|midline)/ ) {
       $self->{'_last_hspdata'}->{$data->{'Name'}} .= $data->{'Data'};
   }  
   return unless ( defined $data->{'Data'} && $data->{'Data'} !~ /^\s+$/ );
   
   $self->{'_last_data'} = $data->{'Data'}; 
}

=head2 _mode

 Title   : _mode
 Usage   : $obj->_mode($newval)
 Function: 
 Example : 
 Returns : value of _mode
 Args    : newvalue (optional)


=cut

sub _mode{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'_mode'} = $value;
    }
    return $self->{'_mode'};
}

=head2 in_element

 Title   : in_element
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub in_element{
   my ($self,$name) = @_;  
   return ( $self->{'_elements'}->[0] eq $name)
}

=head2 start_document

 Title   : start_document
 Usage   : 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub start_document{
    my ($self) = @_;
    $self->{'_lasttype'} = '';
    $self->{'_values'} = {};
    $self->{'_report'}= undef;
    $self->{'_mode'} = '';
    $self->{'_elements'} = [];
}

=head2 end_document

 Title   : end_document
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub end_document{
   my ($self,@args) = @_;
   return $self->{'_report'};
}

1;
