#
# BioPerl module for Bio::SearchIO::axt
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::axt - a parser for axt format reports

=head1 SYNOPSIS

  use Bio::SearchIO;
  my $parser = Bio::SearchIO->new(-format => 'axt',
                                 -file   => 't/data/report.blastz');
  while( my $result = $parser->next_result ) {
    while( my $hit = $result->next_hit ) {
      while( my $hsp = $hit->next_hsp ) {
      }
    }
  }

=head1 DESCRIPTION

This is a parser and event-generator for AXT format reports.  BLASTZ
reports (Schwartz et al,(2003) Genome Research, 13:103-107) are normally
in LAV format but are commonly post-processed to AXT format; many precomputed
BLASTZ reports, such as those found in the UCSC Genome
Browser, are in AXT format.   This parser will also parse any
AXT format produced from any lav report and directly out of BLAT.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SearchIO::axt;
use vars qw(%MODEMAP %MAPPING @STATES $GAPCHAR);
use strict;

use Bio::Search::Result::ResultFactory;
use Bio::Search::HSP::HSPFactory;
use base qw(Bio::SearchIO);

use POSIX;

BEGIN { 
    # mapping of NCBI Blast terms to Bioperl hash keys
    %MODEMAP = ('AXTOutput'   => 'result',
		'Hit'         => 'hit',
		'Hsp'         => 'hsp'
		);
    $GAPCHAR = '-';
    %MAPPING = 
	( 
	  'Hsp_score'      => 'HSP-score',
          'Hsp_query-from' => 'HSP-query_start',
          'Hsp_query-to'   => 'HSP-query_end',
          'Hsp_hit-from'   => 'HSP-hit_start',
          'Hsp_hit-to'     => 'HSP-hit_end',
	  'Hsp_positive'   => 'HSP-conserved',
          'Hsp_identity'   => 'HSP-identical',
          'Hsp_gaps'       => 'HSP-hsp_gaps',
          'Hsp_hitgaps'    => 'HSP-hit_gaps',
          'Hsp_querygaps'  => 'HSP-query_gaps',
          'Hsp_qseq'       => 'HSP-query_seq',
          'Hsp_hseq'       => 'HSP-hit_seq',
          'Hsp_midline'    => 'HSP-homology_seq', # ignoring this for now
          'Hsp_align-len'  => 'HSP-hsp_length',
	  
	  'Hit_id'        => 'HIT-name',	  
	  'AXTOutput_query-def'=> 'RESULT-query_name',	  
	  );
}

=head2 new

 Title   : new
 Usage   : my $obj = Bio::SearchIO::axt->new();
 Function: Builds a new Bio::SearchIO::axt object 
 Returns : an instance of Bio::SearchIO::axt
 Args    :


=cut

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  return $self;
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
    local $/ = "\n";
    local $_;

    my ($curquery,$curhit);
    $self->start_document();
    my @hit_signifs;
    while( defined ($_ = $self->_readline )) { 
	next if (/^\s+$/);
	if( m/^(\d+)\s+      # alignment number - we'll throw this away anyways
	    (\S+)\s+         # Query name      
	    (\d+)\s+(\d+)\s+ # Query start Query end (always + strand, 0 based)
	    (\S+)\s+         # Hit name
	    (\d+)\s+(\d+)\s+ # Hit start Hit end (0 based)
	    ([\-\+])\s+      # Hit strand
	    ([\d\.\-]+)\s+   # Score
	    /ox ) {
	    my ($alnnum, $qname,$qstart,$qend, $hname,
		$hstart,$hend,$hstrand, $score) = ($1,$2,$3,$4,$5,
						   $6,$7,$8,$9);
	    $self->{'_reporttype'} = 'AXT';
	    # Jim's code is 0 based
        # yes, but axt format is one-based, see bug 3145 - cjfields 10-11-10
	    #$qstart++;  $qend++; $hstart++; $hend++;
	    if( defined $curquery &&  $curquery ne $qname ) { 
            $self->end_element({'Name' => 'Hit'});
            $self->_pushback($_);
            $self->end_element({'Name' => 'AXTOutput'});
            return $self->end_document();
	    }
	    
	    if( defined $curhit &&
		$curhit ne $hname) {
		# slight duplication here -- keep these in SYNC
		$self->end_element({'Name' => 'Hit'});
		$self->start_element({'Name' => 'Hit'});
		$self->element({'Name' => 'Hit_id',
				'Data' => $hname});
	    } elsif ( ! defined $curquery ) { 
		$self->start_element({'Name' => 'AXTOutput'});
		$self->{'_result_count'}++;
		$self->element({'Name' => 'AXTOutput_query-def',
				'Data' => $qname });
		
		$self->start_element({'Name' => 'Hit'});
		$self->element({'Name' => 'Hit_id',
				'Data' => $hname});
	    }
	    $self->start_element({'Name' => 'Hsp'});
	    my $queryalign = $self->_readline;
	    my $hitalign   = $self->_readline;
	    chomp($queryalign);
	    chomp($hitalign);
	    my $alnlen = length($queryalign);
	    my $qgapnum = ( $queryalign =~ s/\Q$GAPCHAR/$GAPCHAR/g); 
	    my $hgapnum = ( $hitalign =~ s/\Q$GAPCHAR/$GAPCHAR/g); 
	    my $totalgaps = ($qgapnum + $hgapnum);
	    
	    if( $hstrand eq '-' ) { # strand gets inferred by start/end
		($hstart,$hend) = ($hend,$hstart);
	    }
	    $self->element({'Name' => 'Hsp_score',
			    'Data' => $score});
	    $self->element({'Name' => 'Hsp_query-from',
			    'Data' => $qstart});
	    $self->element({'Name' => 'Hsp_query-to',
			    'Data' => $qend});
	    $self->element({'Name' => 'Hsp_hit-from',
			    'Data' => $hstart});
	    $self->element({'Name' => 'Hsp_hit-to',
			    'Data' => $hend});
	    $self->element({'Name' => 'Hsp_gaps',
			    'Data' => $qgapnum + $hgapnum});
	    $self->element({'Name' => 'Hsp_querygaps',
			    'Data' => $qgapnum});
	    $self->element({'Name' => 'Hsp_hitgaps',
			    'Data' => $hgapnum});
	    
	    $self->element({'Name' => 'Hsp_identity',
			    'Data' => $alnlen - $totalgaps});
	    $self->element({'Name' => 'Hsp_positive',
			    'Data' => $alnlen - $totalgaps});
	    $self->element({'Name' => 'Hsp_qseq',
			    'Data' => $queryalign});
	    $self->element({'Name' => 'Hsp_hseq',
			    'Data' => $hitalign});
	    
	    $self->end_element({'Name' => 'Hsp'});	    
	    $curquery = $qname;
	    $curhit   = $hname;	   
	}
    }
    # fence post
    if( defined $curquery  ) {
	$self->end_element({'Name' => 'Hit'});
	$self->end_element({'Name' => 'AXTOutput'});
	return $self->end_document();
    }
    return;
}

sub _initialize {
    my ($self,@args) = @_;
    $self->SUPER::_initialize(@args);
    $self->_eventHandler->register_factory('result', Bio::Search::Result::ResultFactory->new(-type => 'Bio::Search::Result::GenericResult'));

    $self->_eventHandler->register_factory('hsp', Bio::Search::HSP::HSPFactory->new(-type => 'Bio::Search::HSP::GenericHSP'));
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
    if($nm eq 'AXTOutput') {
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
	$self->warn( "unknown nm $nm ignoring\n");
    }
    $self->{'_last_data'} = ''; # remove read data if we are at 
				# end of an element
    $self->{'_result'} = $rc if( $nm eq 'AXTOutput' );
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

=head2 result_count

 Title   : result_count
 Usage   : my $count = $searchio->result_count
 Function: Returns the number of results we have processed
 Returns : integer
 Args    : none


=cut

sub result_count {
    my $self = shift;
    return $self->{'_result_count'};
}

sub report_count { shift->result_count }



1;

