#
# BioPerl module for Bio::SearchIO::waba
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::waba - SearchIO parser for Jim Kent WABA program
alignment output

=head1 SYNOPSIS

    # do not use this object directly, rather through Bio::SearchIO

    use Bio::SearchIO;
    my $in = Bio::SearchIO->new(-format => 'waba',
			       -file   => 'output.wab');
    while( my $result = $in->next_result ) {
	while( my $hit = $result->next_hit ) {
	    while( my $hsp = $result->next_hsp ) {

	    }
	}
    }

=head1 DESCRIPTION

This parser will process the waba output (NOT the human readable format).

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


package Bio::SearchIO::waba;
use vars qw(%MODEMAP %MAPPING @STATES);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Search::Result::ResultFactory;
use Bio::Search::HSP::HSPFactory;

use POSIX;

BEGIN { 
    # mapping of NCBI Blast terms to Bioperl hash keys
    %MODEMAP = ('WABAOutput' => 'result',
		'Hit'         => 'hit',
		'Hsp'         => 'hsp'
		);
    @STATES = qw(Hsp_qseq Hsp_hseq Hsp_stateseq);
    %MAPPING = 
	( 
	  'Hsp_query-from'=>  'HSP-query_start',
	  'Hsp_query-to'  =>  'HSP-query_end',
	  'Hsp_hit-from'  =>  'HSP-hit_start',
	  'Hsp_hit-to'    =>  'HSP-hit_end',
	  'Hsp_qseq'      =>  'HSP-query_seq',
	  'Hsp_hseq'      =>  'HSP-hit_seq',
	  'Hsp_midline'   =>  'HSP-homology_seq',
	  'Hsp_stateseq'  =>  'HSP-hmmstate_seq',
	  'Hsp_align-len' =>  'HSP-hsp_length',
	  
	  'Hit_id'        => 'HIT-name',
	  'Hit_accession' => 'HIT-accession',

	  'WABAOutput_program'  => 'RESULT-algorithm_name',
	  'WABAOutput_version'  => 'RESULT-algorithm_version',
	  'WABAOutput_query-def'=> 'RESULT-query_name',
	  'WABAOutput_query-db' => 'RESULT-query_database',
 	  'WABAOutput_db'       => 'RESULT-database_name',
	  );
}


use base qw(Bio::SearchIO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::SearchIO::waba->new();
 Function: Builds a new Bio::SearchIO::waba object 
 Returns : Bio::SearchIO::waba
 Args    : see Bio::SearchIO

=cut

sub _initialize {
    my ($self,@args) = @_;
    $self->SUPER::_initialize(@args);
    $self->_eventHandler->register_factory('result', Bio::Search::Result::ResultFactory->new(-type => 'Bio::Search::Result::WABAResult'));

    $self->_eventHandler->register_factory('hsp', Bio::Search::HSP::HSPFactory->new(-type => 'Bio::Search::HSP::WABAHSP'));
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
    my $state = -1;
    $self->start_document();
    my @hit_signifs;
    while( defined ($_ = $self->_readline )) { 
	
	if( $state == -1 ) {
	    my ($qid, $qhspid,$qpercent, $junk,
		$alnlen,$qdb,$qacc,$qstart,$qend,$qstrand,
		$hitdb,$hacc,$hstart,$hend,
		$hstrand) =
		    ( /^(\S+)\.(\S+)\s+align\s+ # get the queryid
		      (\d+(\.\d+)?)\%\s+     # get the percentage
		      of\s+(\d+)\s+  # get the length of the alignment
		      (\S+)\s+           # this is the query database
		      (\S+):(\-?\d+)\-(\-?\d+) # The accession:start-end for query
		      \s+([\-\+])        # query strand
		      \s+(\S+)\.         # hit db
		      (\S+):(\-?\d+)\-(\-?\d+) # The accession:start-end for hit
		      \s+([\-\+])\s*$    # hit strand
		      /ox );
	    
	    # Curses.  Jim's code is 0 based, the following is to readjust
	    if( $hstart < 0 ) { $hstart *= -1}
	    if( $hend   < 0 ) { $hend   *= -1}
	    if( $qstart < 0 ) { $qstart *= -1}
	    if( $qend   < 0 ) { $qend   *= -1}
	    $hstart++; $hend++; $qstart++; $qend++;
	    if( ! defined $alnlen ) {
		$self->warn("Unable to parse the rest of the WABA alignment info for: '$_'");
		last;
	    }
	    $self->{'_reporttype'} = 'WABA'; # hardcoded - only 
	                                     # one type of WABA AFAIK	    
	    if( defined $curquery && 
		$curquery ne $qid ) { 
		$self->end_element({'Name' => 'Hit'});
		$self->_pushback($_);
		$self->end_element({'Name' => 'WABAOutput'});
		return $self->end_document();
	    } 
	    
	    if( defined $curhit &&
		$curhit ne $hacc) {
		# slight duplication here -- keep these in SYNC
		$self->end_element({'Name' => 'Hit'});
		$self->start_element({'Name' => 'Hit'});
		$self->element({'Name' => 'Hit_id',
				'Data' => $hacc});
		$self->element({'Name' => 'Hit_accession',
				'Data' => $hacc});

	    } elsif ( ! defined $curquery ) {
		$self->start_element({'Name' => 'WABAOutput'});
		$self->{'_result_count'}++;
		$self->element({'Name' => 'WABAOutput_query-def',
				'Data' => $qid });
		$self->element({'Name' => 'WABAOutput_program',
				'Data' => 'WABA'});
		$self->element({'Name' => 'WABAOutput_query-db',
				'Data' => $qdb});
		$self->element({'Name' => 'WABAOutput_db',
				'Data' => $hitdb});
		
		# slight duplication here -- keep these N'SYNC ;-)
		$self->start_element({'Name' => 'Hit'});
		$self->element({'Name' => 'Hit_id',
				'Data' => $hacc});
		$self->element({'Name' => 'Hit_accession',
				'Data' => $hacc});
	    }

	    
	    # strand is inferred by start,end values
	    # in the Result Builder
	    if( $qstrand eq '-' ) {
		($qstart,$qend) = ($qend,$qstart);
	    }
	    if( $hstrand eq '-' ) {
		($hstart,$hend) = ($hend,$hstart);
	    }

	    $self->start_element({'Name' => 'Hsp'});
	    $self->element({'Name' => 'Hsp_query-from',
			    'Data' => $qstart});
	    $self->element({'Name' => 'Hsp_query-to',
			    'Data' => $qend});
	    $self->element({'Name' => 'Hsp_hit-from',
			    'Data' => $hstart});
	    $self->element({'Name' => 'Hsp_hit-to',
			    'Data' => $hend});
	    $self->element({'Name' => 'Hsp_align-len',
			    'Data' => $alnlen});
	    
	    $curquery = $qid;
	    $curhit   = $hacc;
	    $state = 0;
	} elsif( ! defined $curquery ) {
	    $self->warn("skipping because no Hit begin line was recognized\n$_") if( $_ !~ /^\s+$/ );
	    next;
	} else { 
	    chomp;
	    $self->element({'Name' => $STATES[$state++],
			    'Data' => $_});
	    if( $state >= scalar @STATES ) {
		$state = -1;
		$self->end_element({'Name' => 'Hsp'});
	    }
	}
    }
    if( defined $curquery  ) {
	$self->end_element({'Name' => 'Hit'});
	$self->end_element({'Name' => 'WABAOutput'});
	return $self->end_document();
    }
    return;
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
    if($nm eq 'WABAOutput') {
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
	$self->warn( "unknown nm $nm ignoring\n");
    }
    $self->{'_last_data'} = ''; # remove read data if we are at 
				# end of an element
    $self->{'_result'} = $rc if( $nm eq 'WABAOutput' );
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
