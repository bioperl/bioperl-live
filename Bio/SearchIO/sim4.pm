# $Id $
#
# BioPerl module for Bio::SearchIO::sim4
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::sim4 - parser for Sim4 alignments

=head1 SYNOPSIS

  # do not use this module directly, it is a driver for SearchIO

  use Bio::SearchIO;
  my $searchio = new Bio::SearchIO(-file => 'results.sim4',
                                   -format => 'sim4');

  while( my $r = $searchio->next_result ) {
    print $r->query_name, "\n";
  }

=head1 DESCRIPTION

This is a driver for the SearchIO system for parsing Sim4.
http://globin.cse.psu.edu/html/docs/sim4.html

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

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SearchIO::sim4;
use strict;
use vars qw(@ISA $DEFAULTFORMAT 
	    @STATES %MAPPING %MODEMAP $DEFAULT_WRITER_CLASS $MIN_INTRON);
use Bio::SearchIO;
use Bio::SearchIO::SearchResultEventBuilder;
$DEFAULTFORMAT = 'SIM4';
@ISA = qw(Bio::SearchIO );

use POSIX;

my %align_types = ( 1 => 'Query', 
		    2 => 'Mid', 
		    3 => 'Sbjct');

%MODEMAP = ('Sim4Output' => 'result',
	    'Hit'             => 'hit',
	    'Hsp'             => 'hsp'
	    );
%MAPPING =
    (
    'Hsp_query-from'=>  'HSP-query_start',
    'Hsp_query-to'  =>  'HSP-query_end',
    'Hsp_hit-from'  =>  'HSP-hit_start',
    'Hsp_hit-to'    =>  'HSP-hit_end',
    'Hsp_qseq'      =>  'HSP-query_seq',
    'Hsp_hseq'      =>  'HSP-hit_seq',
    'Hsp_midline'   =>  'HSP-homology_seq',
    'Hsp_score'     =>  'HSP-score',
    'Hsp_qlength'   =>  'HSP-query_length',
    'Hsp_hlength'   =>  'HSP-hit_length',
    'Hsp_align-len' =>  'HSP-hsp_length',
    'Hsp_identity'  =>  'HSP-identical',
    'Hsp_gaps'       => 'HSP-hsp_gaps',
    'Hsp_hitgaps'    => 'HSP-hit_gaps',
    'Hsp_querygaps'  => 'HSP-query_gaps',

    'Hit_id'        => 'HIT-name',
    'Hit_desc'      => 'HIT-description',
    'Hit_len'       => 'HIT-length',
    'Hit_score'     => 'HIT-score',

    'Sim4Output_program'   => 'RESULT-algorithm_name',
    'Sim4Output_query-def' => 'RESULT-query_name',
    'Sim4Output_query-desc'=> 'RESULT-query_description',
    'Sim4Output_query-len' => 'RESULT-query_length',
    );

$DEFAULT_WRITER_CLASS = 'Bio::Search::Writer::HitTableWriter';

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::SearchIO::sim4();
 Function: Builds a new Bio::SearchIO::sim4 object
 Returns : an instance of Bio::SearchIO::sim4
 Args    :


=cut

sub _initialize {
    my ($self,@args) = @_;
    $self->SUPER::_initialize(@args);
    
    my($rpttype,$fmt ) = $self->_rearrange([qw(REPORT_TYPE
					       SIM4_FORMAT)], @args);
    $fmt ||= 0;
    if( $fmt == 2 ) {
	$self->throw("Cannot parse Sim4 format $fmt, only 0, 1, 3, and 4 allowed");
    }
    defined $fmt && ($self->{'_sim4format'} = $fmt);
    defined $rpttype && ($self->{'_reporttype'} = $rpttype);
}

=head2 next_result

 Title   : next_result
 Usage   : my $hit = $searchio->next_result;
 Function: Returns the next Result from a search
 Returns : Bio::Search::Result::ResultI object
 Args    : none

=cut

sub next_result {
    my ($self) = @_;
    $self->{'_last_data'} = '';
    my ($reporttype,$seenquery,$reportline);
    $self->start_document();
    my ($qfull,$seentop);
    $reporttype = $DEFAULTFORMAT;
    # let's parse very basic Sim4

    my $sim4fmt = $self->{'_sim4format'};
    my $hit_direction = 1;
    while( defined($_ = $self->_readline) ) {       
	next if( /^\s+$/);
	chomp;
	if( /^seq1\s*=\s*(.+),\s+(\d+)\s*bp/ ) {
	    if( $seentop ) {
		$self->_pushback($_);
		$self->end_element({'Name' => 'Hit'});
		$self->end_element({'Name' => 'Sim4Output'});
		return $self->end_document();
	    }
	    $seentop = 1;
	    $self->{'_result_count'}++;
	    $self->start_element({'Name' => 'Sim4Output'});
	    $self->element({'Name' => 'Sim4Output_query-def',
			    'Data' => $1 });
	    $self->element({'Name' => 'Sim4Output_query-len',
			    'Data' => $2 });
	} elsif( /^seq2\s*=\s*(.+),\s+(\d+)\s*bp/ ) {
	    my $len = $2;
	    my ($filename,$hitname) = split(/\s+/,$1);

	    $hitname =~ s/[\(\)\>]//g;
	    $self->start_element({'Name' => 'Hit'});
	    $self->element({'Name' => 'Hit_id',
			    'Data' => $hitname});
	    $self->element({'Name' => 'Hit_desc',
			    'Data' => $filename});
	    $self->element({'Name' => 'Hit_len',
			    'Data' => $len});
	} elsif( /^>(\S+)\s+(.+)?/ ) {
	    if( $qfull )  {
		$self->element({'Name' => 'Hit_desc',
				'Data' => $2});
	    } else { 
		$self->element({'Name'  => 'Sim4Output_query-desc',
				'Data'  => $2});
		$qfull = 1;
	    }
	} elsif( /^\(complement\)/ ) {
	    $hit_direction = -1;
	} elsif( $seentop ) { 
	    if( /^\(?(\d+)\-(\d+)\)?\s+\(?(\d+)\-(\d+)\)?\s+(\S+)/ ) {
		next if( $sim4fmt >= 3 );
		if( $self->in_element('hsp') ) {
		    $self->end_element({'Name' => 'Hsp'});
		}
		my ($qs,$qe,$hs,$he,$pid) = ($1,$2,$3,$4,$5);
		$pid =~ s/\%//;
		$self->start_element({'Name' => 'Hsp'});
		$self->element({'Name' => 'Hsp_query-from',
				'Data' => $qs});
		$self->element({'Name' => 'Hsp_query-to',
				'Data' => $qe});

		if( $hit_direction < 0 ) {
		    $self->element({'Name' => 'Hsp_hit-from',
				    'Data' => $he});
		    $self->element({'Name' => 'Hsp_hit-to',
				    'Data' => $hs});
		} else { 
		    $self->element({'Name' => 'Hsp_hit-from',
				    'Data' => $hs});
		    $self->element({'Name' => 'Hsp_hit-to',
				    'Data' => $he});
		}

		my $numiden = int(($pid / 100) * abs($qe - $qs));
		$self->element({'Name' => 'Hsp_identity',
				'Data' => $numiden});
		$self->element({'Name' => 'Hsp_qlength',
				'Data' => abs($qe - $qs)});
		$self->element({'Name' => 'Hsp_hlength',
				'Data' => abs($he - $hs)});
		$self->element({'Name' => 'Hsp_align-len',
				'Data' => abs($qe - $qs)});
		
		$self->end_element({'Name' => 'Hsp'});
	    } elsif( /^\s+(\d+)\s/ ) {
		unless( $sim4fmt ) { $sim4fmt = $self->{'_sim4format'} = 1 }
		unless( $self->in_element('hsp') ) {
		    $self->start_element({'Name' => 'Hsp'});
		}
		my (%data,$len);
		for( my $i = 0; defined($_) && $i < 4; $i++ ) {
		    chomp;
		    $self->debug("i is $i for $_\n");
		    if( $i == 0 && /^\s+$/ ) {
			$self->_pushback($_) if defined $_;
			$self->end_element({'Name' => 'Hsp'});
			last;
		    } elsif( ($i == 1 || $i == 3) && 
			     s/^(\s+(\d+)\s)// ) {
			my $type = $align_types{$i};       
			my $num = $2;
			$len = CORE::length($1) unless $len;
			$data{$type} = $_;
			my $ungapped = $_;
			$ungapped =~ s/[\s\.]//g;
			if( $num ) {
			    $self->{"\_$type"}->{'begin'} ||= $num;
			    $self->{"\_$type"}->{'end'} = 
				$num + CORE::length($ungapped);
			}
		    } elsif( $i == 2 ) {
			$self->throw("no data for midline $_") 
			    unless (defined $_ && defined $len);
			my $mid = substr($_,$len);
			# $mid =~ s/[\>\<\.]//g;
			$data{'Mid'} = $mid;
		    }
		    $_ = $self->_readline();
		}
		$self->characters({'Name' => 'Hsp_qseq',
				   'Data' => $data{'Query'} });
		$self->characters({'Name' => 'Hsp_hseq',
				   'Data' => $data{'Sbjct'}});
		$self->characters({'Name' => 'Hsp_midline',
				   'Data' => $data{'Mid'} });
	    }
	}
    }

    if( $seentop ) {
	$self->end_element({'Name' => 'Hsp'}) if $self->in_element('hsp');
	$self->end_element({'Name' => 'Hit'});
	$self->end_element({'Name' => 'Sim4Output'});
	return $self->end_document();
    } 
    return undef;
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
   my $type = $MODEMAP{$nm};

   if( $type ) {
       if( $self->_will_handle($type) ) {
	   my $func = sprintf("start_%s",lc $type);
	   $self->_eventHandler->$func($data->{'Attributes'});
       }
       unshift @{$self->{'_elements'}}, $type;

       if($type eq 'result') {
	   $self->{'_values'} = {};
	   $self->{'_result'}= undef;
       }
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
    my $type = $MODEMAP{$nm};
    my $rc;

    if( $nm eq 'Hsp' && $self->{'_sim4format'} > 0 ) {	
	foreach ( qw(Hsp_qseq Hsp_midline Hsp_hseq) ) {
            $self->element({'Name' => $_,
                            'Data' => $self->{'_last_hspdata'}->{$_}});
        }
	my $ident = ( $self->{'_last_hspdata'}->{'Hsp_midline'} =~ s/\|//g);
	$self->element({'Name' => 'Hsp_identity',
                        'Data' => $ident});
        $self->{'_last_hspdata'} = {};
        $self->element({'Name' => 'Hsp_query-from',
                        'Data' => $self->{'_Query'}->{'begin'}});
        $self->element({'Name' => 'Hsp_query-to',
                        'Data' => $self->{'_Query'}->{'end'}});
        
        $self->element({'Name' => 'Hsp_hit-from',
                        'Data' => $self->{'_Sbjct'}->{'begin'}});
        $self->element({'Name' => 'Hsp_hit-to',
                        'Data' => $self->{'_Sbjct'}->{'end'}});

#    } elsif( $nm eq 'Iteration' ) {
# Nothing special needs to be done here.
    }
    if( $type = $MODEMAP{$nm} ) {
	if( $self->_will_handle($type) ) {
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
	$self->debug( "unknown nm $nm, ignoring\n");
    }
    $self->{'_last_data'} = ''; # remove read data if we are at
				# end of an element
    $self->{'_result'} = $rc if( defined $type && $type eq 'result' );
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
   return unless ( defined $data->{'Data'} && $data->{'Data'} !~ /^\s+$/ );
   
   if( $self->in_element('hsp') && 
       $data->{'Name'} =~ /Hsp\_(qseq|hseq|midline)/ ) {
       $self->{'_last_hspdata'}->{$data->{'Name'}} .= $data->{'Data'};
   }  

   $self->{'_last_data'} = $data->{'Data'};
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
 Function: Handle a start document event
 Returns : none
 Args    : none


=cut

sub start_document{
    my ($self) = @_;
    $self->{'_lasttype'} = '';
    $self->{'_values'} = {};
    $self->{'_result'}= undef;
    $self->{'_elements'} = [];
    $self->{'_reporttype'} = $DEFAULTFORMAT;
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


sub write_result {
   my ($self, $blast, @args) = @_;

   if( not defined($self->writer) ) {
       $self->warn("Writer not defined. Using a $DEFAULT_WRITER_CLASS");
       $self->writer( $DEFAULT_WRITER_CLASS->new() );
   }
   $self->SUPER::write_result( $blast, @args );
}

sub result_count {
    my $self = shift;
    return $self->{'_result_count'};
}

sub report_count { shift->result_count }

=head2 _will_handle

 Title   : _will_handle
 Usage   : Private method. For internal use only.
              if( $self->_will_handle($type) ) { ... }
 Function: Provides an optimized way to check whether or not an element of a 
           given type is to be handled.
 Returns : Reference to EventHandler object if the element type is to be handled.
           undef if the element type is not to be handled.
 Args    : string containing type of element.

Optimizations:

  1. Using the cached pointer to the EventHandler to minimize repeated lookups.
  2. Caching the will_handle status for each type that is encountered
     so that it only need be checked by calling handler->will_handle($type) once.

This does not lead to a major savings by itself (only 5-10%).
In combination with other optimizations, or for large parse jobs, the
savings good be significant.

To test against the unoptimized version, remove the parentheses from
around the third term in the ternary " ? : " operator and add two
calls to $self-E<gt>_eventHandler().

=cut

sub _will_handle {
    my ($self,$type) = @_;
    my $handler = $self->{'_handler_cache'} ||= $self->_eventHandler;

    my $will_handle = defined($self->{'_will_handle_cache'}->{$type})
                             ? $self->{'_will_handle_cache'}->{$type}
                             : ($self->{'_will_handle_cache'}->{$type} =
                               $handler->will_handle($type));

    return $will_handle ? $handler : undef;
}

1;

