# $Id$
#
# BioPerl module for Bio::SearchIO::fasta
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::fasta - A SearchIO parser for FASTA results

=head1 SYNOPSIS

  # Do not use this object directly, use it through the SearchIO system
   use Bio::SearchIO;
   my $searchio = new Bio::SearchIO(-format => 'fasta',
				    -file   => 'report.FASTA');
   while( my $result = $searchio->next_result ) {
	# ....
   }

=head1 DESCRIPTION

This object contains the event based parsing code for FASTA format reports.

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


package Bio::SearchIO::fasta;
use vars qw(@ISA %MODEMAP %MAPPING $IDLENGTH);
use strict;

# Object preamble - inherits from Bio::Root::RootI

use Bio::SearchIO;
use POSIX;

BEGIN { 
    # Set IDLENGTH to a new value if you have
    # compile FASTA with a different ID length
    # (actually newest FASTA allows the setting of this
    #  via -C parameter, default is 6)
    $IDLENGTH = 6;

    # mapping of NCBI Blast terms to Bioperl hash keys
    %MODEMAP = ('FastaOutput' => 'result',
		'Hit'         => 'hit',
		'Hsp'         => 'hsp'
		);

    # This should really be done more intelligently, like with
    # XSLT

    %MAPPING = 
	( 
	  'Hsp_bit-score' => 'HSP-bits',
	  'Hsp_score'     => 'HSP-score',
	  'Hsp_sw-score'  => 'HSP-swscore',
	  'Hsp_evalue'    => 'HSP-evalue',
	  'Hsp_query-from'=> 'HSP-query_start',
	  'Hsp_query-to'  => 'HSP-query_end',
	  'Hsp_hit-from'  => 'HSP-hit_start',
	  'Hsp_hit-to'    => 'HSP-hit_end',
	  'Hsp_positive'  => 'HSP-conserved',
	  'Hsp_identity'  => 'HSP-identical',
	  'Hsp_gaps'      => 'HSP-hsp_gaps',
	  'Hsp_hitgaps'   => 'HSP-hit_gaps',
	  'Hsp_querygaps' => 'HSP-query_gaps',
	  'Hsp_qseq'      => 'HSP-query_seq',
	  'Hsp_hseq'      =>  'HSP-hit_seq',
	  'Hsp_midline'   =>  'HSP-homology_seq',
	  'Hsp_align-len' =>  'HSP-hsp_length',
	  'Hsp_query-frame'=> 'HSP-query_frame',
	  'Hsp_hit-frame'  => 'HSP-hit_frame',

	  'Hit_id'        => 'HIT-name',
	  'Hit_len'       => 'HIT-length',
	  'Hit_accession' => 'HIT-accession',
	  'Hit_def'       => 'HIT-description',
	  'Hit_signif'    => 'HIT-significance',
	  'Hit_score'     => 'HIT-score',

	  'FastaOutput_program'  => 'RESULT-algorithm_name',
	  'FastaOutput_version'  => 'RESULT-algorithm_version',
	  'FastaOutput_query-def'=> 'RESULT-query_name',
	  'FastaOutput_query-len'=> 'RESULT-query_length',
	  'FastaOutput_db'       => 'RESULT-database_name',
	  'FastaOutput_db-len'   => 'RESULT-database_entries',
	  'FastaOutput_db-let'   => 'RESULT-database_letters',

	  'Parameters_matrix'    => { 'RESULT-parameters' => 'matrix'},
	  'Parameters_expect'    => { 'RESULT-parameters' => 'expect'},
	  'Parameters_include'   => { 'RESULT-parameters' => 'include'},
	  'Parameters_sc-match'  => { 'RESULT-parameters' => 'match'},
	  'Parameters_sc-mismatch' => { 'RESULT-parameters' => 'mismatch'},
	  'Parameters_gap-open'  => { 'RESULT-parameters' => 'gapopen'},
	  'Parameters_gap-ext'   => { 'RESULT-parameters' => 'gapext'},
	  'Parameters_word-size' => { 'RESULT-parameters' => 'wordsize'},
	  'Parameters_ktup'      => { 'RESULT-parameters' => 'ktup'},
	  'Parameters_filter'    => {'RESULT-parameters' => 'filter'},
	  'Statistics_db-num'    => { 'RESULT-statistics' => 'dbentries'},
	  'Statistics_db-len'    => { 'RESULT-statistics' => 'dbletters'},
	  'Statistics_hsp-len'   => { 'RESULT-statistics' => 'hsplength'},
	  'Statistics_eff-space' => { 'RESULT-statistics' => 'effectivespace'},
	  'Statistics_kappa'     => { 'RESULT-statistics' => 'kappa' },
	  'Statistics_lambda'    => { 'RESULT-statistics' => 'lambda' },
	  'Statistics_entropy'   => { 'RESULT-statistics' => 'entropy'},
	  );
}


@ISA = qw(Bio::SearchIO );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::SearchIO::fasta();
 Function: Builds a new Bio::SearchIO::fasta object 
 Returns : Bio::SearchIO::fasta
 Args    : -idlength - set ID length to something other 
                       than the default (7), this is only
                       necessary if you have compiled FASTA
                       with a new default id length to display
                       in the HSP alignment blocks

=cut

sub _initialize {
  my($self,@args) = @_;
  $self->SUPER::_initialize(@args);
  return unless @args;
  my ($idlength) = $self->_rearrange([qw(IDLENGTH)],@args);
  $self->idlength($idlength || $IDLENGTH);
  return 1;
}

=head2 next_result

 Title   : next_result
 Usage   : my $hit = $searchio->next_result;
 Function: Returns the next Result from a search
 Returns : Bio::Search::Result::ResultI object
 Args    : none

=cut

sub next_result{
   my ($self) = @_;
   
   my $data = '';
   my $seentop = 0;
   my $current_hsp;
   $self->start_document();
   my @hit_signifs;
   while( defined ($_ = $self->_readline )) {       
       next if( ! $self->in_element('hsp')  &&
		/^\s+$/); # skip empty lines
       if( /(\S+)\s+searches\s+a\s+((protein\s+or\s+DNA\s+sequence)|(sequence\s+database))/i || /(\S+) compares a/ ) {
	   if( $seentop ) {
	       $self->_pushback($_);
	       $self->end_element({ 'Name' => 'FastaOutput'});
	       return $self->end_document();
	   }
	   $self->{'_reporttype'} = $1;
	   $self->start_element({ 'Name' => 'FastaOutput' } );
	   $seentop = 1;
	   
	   $self->element({ 'Name' => 'FastaOutput_program',
			    'Data' => $self->{'_reporttype'}});
	   $_ = $self->_readline();
	   my ($version) = (/version\s+(\S+)/);
	   $version = '' unless defined $version;
	   $self->element({ 'Name' => 'FastaOutput_version',
			    'Data' => $version});
	   my ($last);
	   while( defined($_ = $self->_readline()) ) {
	       if( /\s+>(.+)/ || /^\s*vs\s+/ ) {
		   my $querydef = $1;
		   
		   if( $last =~ /(\S+)[:,]\s*(\d+)\s+(aa|nt)/ ) {
		       if( $self->{'_reporttype'} &&
			   $self->{'_reporttype'} eq 'FASTA' ) {
			   if( $3 eq 'nt') {
			       $self->{'_reporttype'} = 'FASTN' ;
			   } elsif( $3 eq 'aa' ) {
			       $self->{'_reporttype'} = 'FASTP' ;
			   }
		       }
		       
		       $self->element({'Name' => 'FastaOutput_query-def',
				       'Data' => $querydef || $1});
		       $self->element({'Name' => 'FastaOutput_query-len',
				       'Data' => $2});
		   } else {
		       $self->element({'Name' => 'FastaOutput_query-def',
				       'Data' => $querydef });
		       $self->warn("unable to find and set query length");
		   }
		   last;
	       } 
	       $last = $_;
	   }

	   if( $last =~ /^\s*vs\s+(\S+)/ ||	       	       
	       (defined $_ && /^\s*vs\s+(\S+)/) ||
	       (defined ($_ = $self->_readline()) && /^\s*vs\s+(\S+)/)
	       ) {
	       $self->element({'Name' => 'FastaOutput_db',
			       'Data' => $1});
	   }
       } elsif( /(\d+) residues in\s+(\d+)\s+sequences/ ) {
	   $self->element({'Name' => 'FastaOutput_db-let',
			   'Data' => $1});
	   $self->element({'Name' => 'FastaOutput_db-len',
			   'Data' => $2});
	   $self->element({'Name' => 'Statistics_db-len',
			   'Data' => $1});
	   $self->element({'Name' => 'Statistics_db-num',
			   'Data' => $2});	   
       } elsif( /Lambda=\s+(\S+)/ ) {
	   $self->element({'Name' => 'Statistics_lambda',
			   'Data' => $1});	  
       } elsif( /^\s*(Smith-Waterman).+(\S+)\s*matrix/ ) {	   
	   $self->element({'Name' => 'Parameters_matrix',
			   'Data' => $2});
	   $self->{'_reporttype'} = $1;

	   $self->element({ 'Name' => 'FastaOutput_program',
			    'Data' => $self->{'_reporttype'}});
	   
       } elsif( /The best scores are:/ ) {
	   while( defined ($_ = $self->_readline() ) && 
		  ! /^\s+$/ ) {	       
	       my @line = split;
	       # 0 is signif, 1 is raw score
	       push @hit_signifs, [ pop @line, pop @line];
	   }
	   
       } elsif( /^\s*([T]?FAST[XYAF]).+,\s*(\S+)\s*matrix.+ktup:\s*(\d+)/ ) {
	   $self->element({'Name' => 'Parameters_matrix',
			   'Data' => $2});
	   $self->element({'Name' => 'Parameters_ktup',
			   'Data' => $3});
	   $self->{'_reporttype'} = $1 if( $self->{'_reporttype'} !~ /FAST[PN]/i ) ;

	   $self->element({ 'Name' => 'FastaOutput_program',
			    'Data' => $self->{'_reporttype'}});
	   
       } elsif( /gap\-pen:\s+([\-\+]?\d+)\/\s+([\-\+]?\d+).+width:\s+(\d+)/ ) {
	   $self->element({'Name' => 'Parameters_gap-open',
			   'Data' => $1});
	   $self->element({'Name' => 'Parameters_gap-ext',
			   'Data' => $2});
	   $self->element({'Name' => 'Parameters_word-size',
			   'Data' => $3});
       } elsif( /^>>(.+) \((\d+)\s*(aa|nt)\)$/ ) {
	   if( $self->in_element('hsp') ) {
	       $self->end_element({ 'Name' => 'Hsp'});
	   }
	   if( $self->in_element('hit') ) {
	       $self->end_element({ 'Name' => 'Hit'});
	   }
	   
	   $self->start_element({'Name' => 'Hit'});
	   $self->element({ 'Name' => 'Hit_len',
			    'Data' => $2});  
	   my ($id,$desc) = split(/\s+/,$1,2);
	   $self->element({ 'Name' => 'Hit_id',
			    'Data' => $id}); 	   
	   my $v = shift @hit_signifs;
	   if( defined $v ) {
	       $self->element({'Name' => 'Hit_signif',
			       'Data' => $v->[0]});
	       $self->element({'Name' => 'Hit_score',
			       'Data' => $v->[1]});
	   }
	   my @pieces = split(/\|/,$id);
	   my $acc = pop @pieces;
	   $acc =~ s/\.\d+$//;
	   $self->element({ 'Name' =>  'Hit_accession',
			    'Data'  => $acc});	
	   $self->element({ 'Name' => 'Hit_def',
			    'Data' => $desc});	   
	   $self->start_element({'Name' => 'Hsp'});
	   $_ = $self->_readline();
	   my ($score,$bits,$e) = ( /Z-score:\s*(\S+)\s*bits:\s*(\S+)\s+E\(\):\s*(\S+)/ );
	   $self->element({'Name' => 'Hsp_score',
			   'Data' => $score});
	   $self->element({'Name' => 'Hsp_evalue',
			   'Data' => $e});
	   $self->element({'Name' => 'Hsp_bit-score',
			   'Data' => $bits});
	   $_ = $self->_readline();
	   if( /Smith-Waterman score:\s*(\d+)/ ) {
	       $self->element({'Name' => 'Hsp_sw-score',
			       'Data' => $1});
	   }
	   if( /(\d+(\.\d+)?)\%\s*identity(\s*\((\d+(\.\d+)?)\%\s*ungapped\))?\s*in\s*(\d+)\s+(aa|nt)\s+overlap\s*\((\d+)\-(\d+):(\d+)\-(\d+)\)/ ) {
	       my ($identper,$gapper,$len,$querystart,
		   $queryend,$hitstart,$hitend) = ($1,$4,$6,$8,$9,$10,$11);
	       my $ident = POSIX::ceil(($identper/100) * $len);
	       my $gaps = ( defined $gapper ) ? POSIX::ceil ( ($gapper/100) * $len) : undef;
	       
	       $self->element({'Name' => 'Hsp_gaps',
			       'Data' => $gaps});
	       $self->element({'Name' => 'Hsp_identity',
			       'Data' => $ident});
	       $self->element({'Name' => 'Hsp_positive',
			       'Data' => $ident});
	       $self->element({'Name' => 'Hsp_align-len',
			       'Data' => $len});
	       
	       $self->element({'Name' => 'Hsp_query-from',
			       'Data' => $querystart});
	       $self->element({'Name' => 'Hsp_query-to',
			       'Data' => $queryend});
	       $self->element({'Name' => 'Hsp_hit-from',
			       'Data' => $hitstart});
	       $self->element({'Name' => 'Hsp_hit-to',
			       'Data' => $hitend});
	   } else {
	       $self->warn( "unable to parse FASTA score line: $_");
	   }
       } elsif( /\d+\s*residues\s*in\s*\d+\s*query\s*sequences/ ) {
	   if( $self->in_element('hsp') ) {
	       $self->end_element({'Name' => 'Hsp'});
	   } 
	   if( $self->in_element('hit') ) {
	       $self->end_element({'Name' => 'Hit'});
	   }
	   
#	   $_ = $self->_readline();
#	   my ( $liblen,$libsize) = /(\d+)\s+residues\s*in(\d+)\s*library/;
	   # fast forward to the end of the file as there is 
	   # nothing else left to do with this file and want to be sure and
	   # reset it
	   while(defined($_ = $self->_readline() ) ) { 
	       last if( /^Function used was/);
	       if( /(\S+)\s+searches\s+a\s+(protein\s+or\s+DNA\s+sequence)|(sequence\s+database)/ ) { 
		   $self->_pushback($_);
	       }
	   }
	   $self->end_element({ 'Name' => 'FastaOutput'});
	   return $self->end_document();
       } elsif( $self->in_element('hsp' ) ) {
	   
	   my @data = ( '','','');
	   my $count = 0;
	   my $len = $self->idlength + 1;
	   my ($seq1_id);
	   while( defined($_ ) ) {
	       chomp;
	       $self->debug( "$count $_\n");

	       if( /residues in \d+\s+query\s+sequences/) {
		   $self->_pushback($_);
		   last;
	       } elsif( /^>>/ ) {
		   $self->_pushback($_);
		   last;
	       }
	       if( $count == 0 ) { 
		   unless( /^\s+\d+\s+/ || /^\s+$/) {
		       $self->_pushback($_);
		       $count = 2;
		   }
	       } elsif( $count == 1 || $count == 3 ) {
		   if( /^(\S+)\s+/ ) {
		       $len = CORE::length($1) if $len < CORE::length($1);
		       s/\s+$//; # trim trailing spaces,we don't want them 
		       $data[$count-1] = substr($_,$len);
		       
		   } elsif( /^\s+(\d+)\s+/ ) {
		       $self->warn("Unexpected state ($_)");
		   } elsif( /^\s+$/ || length($_) == 0) {
		       $count = 5;  
		       # going to skip these
		   } else {
		       $self->warn("Unrecognized alignment line ($count) '$_'");
		   }
	       } elsif( $count == 2 ) {
		   if( /^\s+\d+\s+/ ) {
		       $self->warn("$_\n");
		       $count = 4;
		   } else {
		       # toss the first IDLENGTH characters of the line
		       if( length($_) >= $len ) {
			   $data[$count-1] = substr($_,$len);
		       }
		   }
	       } 
	       last if( $count++ >= 5);
	       $_ = $self->_readline();	       
	   }
	   if( length($data[0]) > 0 || length($data[2]) > 0 ) {
	       $self->characters({'Name' => 'Hsp_qseq',
				  'Data' => $data[0] });
	       $self->characters({'Name' => 'Hsp_midline',
				  'Data' => $data[1]});
	       $self->characters({'Name' => 'Hsp_hseq',
				  'Data' => $data[2]});
	   }
       } else {
	   if( ! $seentop ) {
	       $self->debug($_);
	       $self->warn("unrecognized FASTA Family report file!");
	       return undef;
	   }
       }
   }
}


=head2 start_element

 Title   : start_element
 Usage   : $eventgenerator->start_element
 Function: Handles a start element event
 Returns : none
 Args    : hashref with at least 2 keys 'Data' and 'Name'


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
    if($nm eq 'FastaOutput') {
	$self->{'_values'} = {};
	$self->{'_result'}= undef;
	$self->{'_mode'} = '';
    }

}

=head2 end_element

 Title   : start_element
 Usage   : $eventgenerator->end_element
 Function: Handles an end element event
 Returns : none
 Args    : hashref with at least 2 keys 'Data' and 'Name'


=cut

sub end_element {
    my ($self,$data) = @_;
    my $nm = $data->{'Name'};
    my $rc;
    # Hsp are sort of weird, in that they end when another
    # object begins so have to detect this in end_element for now
    if( $nm eq 'Hsp' ) {
	foreach ( qw(Hsp_qseq Hsp_midline Hsp_hseq) ) {
	    $self->element({'Name' => $_,
			    'Data' => $self->{'_last_hspdata'}->{$_}});
	}
	$self->{'_last_hspdata'} = {}
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
	$self->warn( "unknown nm $nm, ignoring\n");
    }
    $self->{'_last_data'} = ''; # remove read data if we are at 
				# end of an element
    $self->{'_result'} = $rc if( $nm eq 'FastaOutput' );
    return $rc;

}

=head2 element

 Title   : element
 Usage   : $eventhandler->element({'Name' => $name, 'Data' => $str});
 Function: Convience method that calls start_element, characters, end_element
 Returns : none
 Args    : Hash ref with the keys 'Name' and 'Data'


=cut

sub element{
   my ($self,$data) = @_;
   $self->start_element($data);
   $self->characters($data);
   $self->end_element($data);
}


=head2 characters

 Title   : characters
 Usage   : $eventgenerator->characters($str)
 Function: Send a character events
 Returns : none
 Args    : string


=cut

sub characters{
   my ($self,$data) = @_;   

   return unless ( defined $data->{'Data'} );
   if( $data->{'Data'} =~ /^\s+$/ ) {
       return unless $data->{'Name'} =~ /Hsp\_(midline|qseq|hseq)/;
   }

   if( $self->in_element('hsp') && 
       $data->{'Name'} =~ /Hsp\_(qseq|hseq|midline)/ ) {
       
       $self->{'_last_hspdata'}->{$data->{'Name'}} .= $data->{'Data'};
   }  
   
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

=head2 within_element

 Title   : within_element
 Usage   : if( $eventgenerator->within_element($element) ) {}
 Function: Test if we are within a particular element
           This is different than 'in' because within can be tested
           for a whole block.
 Returns : boolean
 Args    : string element name 


=cut

sub within_element{
   my ($self,$name) = @_;  
   return 0 if ( ! defined $name &&
		 ! defined  $self->{'_elements'} ||
		 scalar @{$self->{'_elements'}} == 0) ;
   foreach (  @{$self->{'_elements'}} ) {
       if( $_ eq $name  ) {
	   return 1;
       } 
   }
   return 0;
}

=head2 in_element

 Title   : in_element
 Usage   : if( $eventgenerator->in_element($element) ) {}
 Function: Test if we are in a particular element
           This is different than 'in' because within can be tested
           for a whole block.
 Returns : boolean
 Args    : string element name 


=cut

sub in_element{
   my ($self,$name) = @_;  
   return 0 if ! defined $self->{'_elements'}->[0];
   return ( $self->{'_elements'}->[0] eq $name)
}


=head2 start_document

 Title   : start_document
 Usage   : $eventgenerator->start_document
 Function: Handles a start document event
 Returns : none
 Args    : none


=cut

sub start_document{
    my ($self) = @_;
    $self->{'_lasttype'} = '';
    $self->{'_values'} = {};
    $self->{'_result'}= undef;
    $self->{'_mode'} = '';
    $self->{'_elements'} = [];
}


=head2 end_document

 Title   : end_document
 Usage   : $eventgenerator->end_document
 Function: Handles an end document event
 Returns : Bio::Search::Result::ResultI object
 Args    : none


=cut

sub end_document{
   my ($self,@args) = @_;
   return $self->{'_result'};
}

=head2 idlength

 Title   : idlength
 Usage   : $obj->idlength($newval)
 Function: Internal storage of the length of the ID desc
           in the HSP alignment blocks.  Defaults to
           $IDLENGTH class variable value
 Returns : value of idlength
 Args    : newvalue (optional)


=cut

sub idlength{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_idlength'} = $value;
    }
    return $self->{'_idlength'} || $IDLENGTH;
}

1;

