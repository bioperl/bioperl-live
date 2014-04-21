#
# BioPerl module for Bio::SearchIO::sim4
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
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
  my $searchio = Bio::SearchIO->new(-file => 'results.sim4',
                                   -format => 'sim4');

  while ( my $result = $searchio->next_result ) {
      while ( my $hit = $result->next_hit ) {
	  while ( my $hsp = $hit->next_hsp ) {
              # ...
	  }
      }
  }

=head1 DESCRIPTION

This is a driver for the SearchIO system for parsing Sim4.
http://globin.cse.psu.edu/html/docs/sim4.html

Cannot parse LAV or 'exon file' formats (A=2 or A=5)

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

Email jason-at-bioperl-dot-org

=head1 CONTRIBUTORS

Luc Gauthier (lgauthie@hotmail.com)

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SearchIO::sim4;

use strict;
use vars qw($DEFAULTFORMAT %ALIGN_TYPES
            %MAPPING %MODEMAP $DEFAULT_WRITER_CLASS);

use POSIX;
use Bio::SearchIO::SearchResultEventBuilder;

use base qw(Bio::SearchIO);

$DEFAULTFORMAT = 'SIM4';
$DEFAULT_WRITER_CLASS = 'Bio::SearchIO::Writer::HitTableWriter';

%ALIGN_TYPES = (
    0 => 'Ruler',
    1 => 'Query', 
    2 => 'Mid', 
    3 => 'Sbjct'
);

%MODEMAP = (
    'Sim4Output' => 'result',
    'Hit'        => 'hit',
    'Hsp'        => 'hsp'
);

%MAPPING = (
    'Hsp_query-from'=>  'HSP-query_start',
    'Hsp_query-to'  =>  'HSP-query_end',
    'Hsp_qseq'      =>  'HSP-query_seq',
    'Hsp_qlength'   =>  'HSP-query_length',
    'Hsp_querygaps'  => 'HSP-query_gaps',
    'Hsp_hit-from'  =>  'HSP-hit_start',
    'Hsp_hit-to'    =>  'HSP-hit_end',
    'Hsp_hseq'      =>  'HSP-hit_seq',
    'Hsp_hlength'   =>  'HSP-hit_length',
    'Hsp_hitgaps'    => 'HSP-hit_gaps',
    'Hsp_midline'   =>  'HSP-homology_seq',
    'Hsp_score'     =>  'HSP-score',
    'Hsp_align-len' =>  'HSP-hsp_length',
    'Hsp_identity'  =>  'HSP-identical',

    'Hit_id'        => 'HIT-name',
    'Hit_desc'      => 'HIT-description',
    'Hit_len'       => 'HIT-length',

    'Sim4Output_program'   => 'RESULT-algorithm_name',
    'Sim4Output_query-def' => 'RESULT-query_name',
    'Sim4Output_query-desc'=> 'RESULT-query_description',
    'Sim4Output_query-len' => 'RESULT-query_length',
);



=head2 new

 Title   : new
 Usage   : my $obj = Bio::SearchIO::sim4->new();
 Function: Builds a new Bio::SearchIO::sim4 object
 Returns : an instance of Bio::SearchIO::sim4
 Args    :


=cut


=head2 next_result

 Title   : next_result
 Usage   : my $result = $searchio->next_result;
 Function: Returns the next Result from a search
 Returns : Bio::Search::Result::ResultI object
 Args    : none

=cut

sub next_result {
    my ($self) = @_;
    local $/ = "\n";
    local $_;

    # Declare/adjust needed variables
    $self->{'_last_data'} = '';
    my ($seentop, $qfull, @hsps, %alignment, $format);
    my $hit_direction = 1;

    # Start document and main element
    $self->start_document();
    $self->start_element({'Name' => 'Sim4Output'});
    my $lastquery = '';
    # Read output report until EOF
    while( defined($_ = $self->_readline) ) {       
        # Skip empty lines, chomp filled ones
	next if( /^\s+$/); chomp;

        # Make sure sim4 output format is not 2 or 5
        if (!$seentop) {
	    if ( /^\#:lav/ ) { $format = 2; }
            elsif ( /^<|>/ ) { $format = 5; }
            $self->throw("Bio::SearchIO::sim4 module cannot parse 'type $format' outputs.") if $format;
	}

        # This line indicates the start of a new hit
	if( /^seq1\s*=\s*(\S+),\s+(\d+)/ ) {
	    my ($nm,$desc) = ($1,$2);
            # First hit? Adjust some parameters if so
	    if ( ! $seentop ) {
	        $self->element( {'Name' => 'Sim4Output_query-def', 
				 'Data' => $nm} );
	        $self->element( {'Name' => 'Sim4Output_query-len', 
				 'Data' => $desc} );
                $seentop = 1;
	    } elsif( $nm ne $lastquery ) {
		$self->_pushback($_);
		last;
	    }
	    $lastquery = $nm;
            # A previous HSP may need to be ended
            $self->end_element({'Name' => 'Hsp'}) if ( $self->in_element('hsp') );
            # A previous hit exists? End it and reset needed variables
            if ( $self->in_element('hit') ) {
	        foreach (@hsps) {
                    $self->start_element({'Name' => 'Hsp'});
                    while (my ($name, $data) = each %$_) {
                        $self->{'_currentHSP'}{$name} = $data;
		    }
		    $self->end_element({'Name' => 'Hsp'});
                    $self->{'_currentHSP'} = {};
	        }
                $format = 0 if @hsps;
                @hsps = ();
                %alignment = ();
                $qfull = 0;
                $hit_direction = 1;
                $self->end_element({'Name' => 'Hit'});
	    }

        # This line describes the current hit... so let's start it
	} elsif( /^seq2\s*=\s*(\S+)\s+\(>?(\S+)\s*\),\s*(\d+)/ ) {
            $self->start_element({'Name' => 'Hit'});
	    $self->element( {'Name' => 'Hit_id', 'Data' => $2} );
	    $self->element( {'Name' => 'Hit_desc', 'Data' => $1} );
	    $self->element( {'Name' => 'Hit_len', 'Data' => $3} );

        # This line may give additional details about query or subject
	} elsif( /^>(\S+)\s*(.*)?/ ) {
            # Previous line was query details... this time subject details
	    if( $qfull )  {
                $format = 4 if !$format;
		$self->element({'Name' => 'Hit_desc', 'Data' => $2});
            # First line of this type is always query details for a given hit
	    } else { 
		$self->element({'Name' => 'Sim4Output_query-desc', 'Data' => $2});
		$qfull = 1;
	    }

        # This line indicates that subject is on reverse strand
	} elsif( /^\(complement\)/ ) {
	    $hit_direction = -1;

        # This line describes the current HSP... so add it to @hsps array
	} elsif( /^\(?(\d+)\-(\d+)\)?\s+\(?(\d+)\-(\d+)\)?\s+(\d+)/ ) {
		my ($qs,$qe,$hs,$he,$pid) = ($1,$2,$3,$4,$5);
                push @hsps, {
                    'Hsp_query-from' => $qs,
                    'Hsp_query-to' => $qe,
                    'Hsp_hit-from' => $hit_direction >= 0 ? $hs : $he,
                    'Hsp_hit-to' => $hit_direction >= 0 ? $he : $hs,
                    'Hsp_identity' => 0, #can't determine correctly from raw pct
                    'Hsp_qlength' => abs($qe - $qs) + 1,
                    'Hsp_hlength' => abs($he - $hs) + 1,
                    'Hsp_align-len' => abs($qe - $qs) + 1,
	        };

        # This line indicates the start of an alignment block
        } elsif( /^\s+(\d+)\s/ ) {
            # Store the current alignment block in a hash
	    for( my $i = 0; defined($_) && $i < 4; $i++ ) {
                my ($start, $string) = /^\s+(\d*)\s(.*)/;
                $alignment{$ALIGN_TYPES{$i}} = { start => $start, string => $i != 2
                    ? $string
                    : (' ' x (length($alignment{$ALIGN_TYPES{$i-1}}{string}) - length($string))) . $string
                };
                $_ = $self->_readline();
	    }

            # 'Ruler' line indicates the start of a new HSP
            if ($alignment{Ruler}{start} == 0) {
                $format = @hsps ? 3 : 1 if !$format;
                # A previous HSP may need to be ended
                $self->end_element({'Name' => 'Hsp'}) if ( $self->in_element('hsp') );
                # Start the new HSP and fill the '_currentHSP' property with available details
     	        $self->start_element({'Name' => 'Hsp'});
                $self->{'_currentHSP'} = @hsps ? shift @hsps : {
                    'Hsp_query-from' => $alignment{Query}{start},
                    'Hsp_hit-from' => $alignment{Sbjct}{start},
     	        }
	    }

            # Midline indicates a boundary between two HSPs
	    if ( $alignment{Mid}{string} =~ /<|>/g ) {
                my ($hsp_start, $hsp_end);
                # Are we currently in an open HSP?
    	        if ( $self->in_element('hsp') ) {
                    # Find end pos, adjust 'gaps', 'seq' and 'midline' properties... then close HSP
                    $hsp_end = (pos $alignment{Mid}{string}) - 1;
                    $self->{'_currentHSP'}{'Hsp_querygaps'} +=
                        ($self->{'_currentHSP'}{'Hsp_qseq'} .= substr($alignment{Query}{string}, 0, $hsp_end)) =~ s/ /-/g;
                    $self->{'_currentHSP'}{'Hsp_hitgaps'} +=
                        ($self->{'_currentHSP'}{'Hsp_hseq'} .= substr($alignment{Sbjct}{string}, 0, $hsp_end)) =~ s/ /-/g;
                    ($self->{'_currentHSP'}{'Hsp_midline'} .= substr($alignment{Mid}{string}, 0, $hsp_end)) =~ s/-/ /g;
                    $self->end_element({'Name' => 'Hsp'});

                    # Does a new HSP start in the current alignment block?
                    if ( $alignment{Mid}{string} =~ /\|/g ) {
                        # Find start pos, start new HSP and fill it with available details
                        $hsp_start = (pos $alignment{Mid}{string}) - 1;
                        $self->start_element({'Name' => 'Hsp'});
                        $self->{'_currentHSP'} = @hsps ? shift @hsps : {};
                        $self->{'_currentHSP'}{'Hsp_querygaps'} +=
                            ($self->{'_currentHSP'}{'Hsp_qseq'} = substr($alignment{Query}{string}, $hsp_start)) =~ s/ /-/g;
                        $self->{'_currentHSP'}{'Hsp_hitgaps'} +=
                            ($self->{'_currentHSP'}{'Hsp_hseq'} = substr($alignment{Sbjct}{string}, $hsp_start)) =~ s/ /-/g;
                        ($self->{'_currentHSP'}{'Hsp_midline'} = substr($alignment{Mid}{string}, $hsp_start)) =~ s/-/ /g;
		    }
		}
                # No HSP is currently open...
                else {
                    # Find start pos, start new HSP and fill it with available
                    # details then skip to next alignment block
		    $hsp_start = index($alignment{Mid}{string}, '|');
	            $self->start_element({'Name' => 'Hsp'});
                    $self->{'_currentHSP'} = @hsps ? shift @hsps : {
                        'Hsp_query-from' => $alignment{Query}{start},
    	            };
                    $self->{'_currentHSP'}{'Hsp_querygaps'} +=
                        ($self->{'_currentHSP'}{'Hsp_qseq'} = substr($alignment{Query}{string}, $hsp_start)) =~ s/ /-/g;
                    $self->{'_currentHSP'}{'Hsp_hitgaps'} +=
                        ($self->{'_currentHSP'}{'Hsp_hseq'} = substr($alignment{Sbjct}{string}, $hsp_start)) =~ s/ /-/g;
                    ($self->{'_currentHSP'}{'Hsp_midline'} = substr($alignment{Mid}{string}, $hsp_start)) =~ s/-/ /g;
                    next;
		}
	    }
            # Current alignment block does not contain HSPs boundary
            else {
                # Start a new HSP if none is currently open
	        # (Happens if last boundary finished at the very end of previous block)
	        if ( !$self->in_element('hsp') ) {
     	            $self->start_element({'Name' => 'Hsp'});
                    $self->{'_currentHSP'} = @hsps ? shift @hsps : {
                        'Hsp_query-from' => $alignment{Query}{start},
                        'Hsp_hit-from' => $alignment{Sbjct}{start},
     	            }
		}
                # Adjust details of the current HSP
                $self->{'_currentHSP'}{'Hsp_query-from'} ||= 
		    $alignment{Query}{start} - 
		    length($self->{'_currentHSP'}{'Hsp_qseq'} || '');
                $self->{'_currentHSP'}{'Hsp_hit-from'} ||= 
		    $alignment{Sbjct}{start} - 
		    length($self->{'_currentHSP'}{'Hsp_hseq'} || '');
                $self->{'_currentHSP'}{'Hsp_querygaps'} +=
                    ($self->{'_currentHSP'}{'Hsp_qseq'} .= 
		     $alignment{Query}{string}) =~ s/ /-/g;
                $self->{'_currentHSP'}{'Hsp_hitgaps'} +=
                    ($self->{'_currentHSP'}{'Hsp_hseq'} .= 
		     $alignment{Sbjct}{string}) =~ s/ /-/g;
                ($self->{'_currentHSP'}{'Hsp_midline'} .= 
		 $alignment{Mid}{string}) =~ s/-/ /g;
	    }
	}
    }

    # We are done reading the sim4 report, end everything and return
    if( $seentop ) {
        # end HSP if needed
        $self->end_element({'Name' => 'Hsp'}) if ( $self->in_element('hsp') );
        # end Hit if needed
        if ( $self->in_element('hit') ) {
            foreach (@hsps) {
                $self->start_element({'Name' => 'Hsp'});
                while (my ($name, $data) = each %$_) {
                    $self->{'_currentHSP'}{$name} = $data;
    	        }
    	        $self->end_element({'Name' => 'Hsp'});
            }
            $self->end_element({'Name' => 'Hit'});
	}
        # adjust result's algorithm name, end output and return
        $self->element({'Name' => 'Sim4Output_program',
                        'Data' => $DEFAULTFORMAT . ' (A=' . (defined $format ? $format : '?') . ')'});
	$self->end_element({'Name' => 'Sim4Output'});
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
    
    if( $nm eq 'Hsp' ) {
        $self->{'_currentHSP'}{'Hsp_midline'} ||= '';
	$self->{'_currentHSP'}{'Hsp_query-to'} ||=
            $self->{'_currentHSP'}{'Hsp_query-from'} + length($self->{'_currentHSP'}{'Hsp_qseq'}) - 1 - $self->{'_currentHSP'}{'Hsp_querygaps'};
        $self->{'_currentHSP'}{'Hsp_hit-to'} ||=
            $self->{'_currentHSP'}{'Hsp_hit-from'} + length($self->{'_currentHSP'}{'Hsp_hseq'}) - 1 - $self->{'_currentHSP'}{'Hsp_hitgaps'};
        $self->{'_currentHSP'}{'Hsp_identity'} ||= 
	    ($self->{'_currentHSP'}{'Hsp_midline'} =~ tr/\|//);
        $self->{'_currentHSP'}{'Hsp_qlength'} ||= abs($self->{'_currentHSP'}{'Hsp_query-to'} - $self->{'_currentHSP'}{'Hsp_query-from'}) + 1;
        $self->{'_currentHSP'}{'Hsp_hlength'} ||= abs($self->{'_currentHSP'}{'Hsp_hit-to'} - $self->{'_currentHSP'}{'Hsp_hit-from'}) + 1;
        $self->{'_currentHSP'}{'Hsp_align-len'} ||= abs($self->{'_currentHSP'}{'Hsp_query-to'} - $self->{'_currentHSP'}{'Hsp_query-from'}) + 1;
        $self->{'_currentHSP'}{'Hsp_score'} ||= int(100 * ($self->{'_currentHSP'}{'Hsp_identity'} / $self->{'_currentHSP'}{'Hsp_align-len'}));
        foreach (keys %{$self->{'_currentHSP'}}) {
            $self->element({'Name' => $_, 'Data' => delete ${$self->{'_currentHSP'}}{$_}});
	}
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
    return 1; # can a sim4 report contain more than one result?
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
