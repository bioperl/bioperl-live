#
# BioPerl module for Bio::SearchIO::exonerate
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

Bio::SearchIO::exonerate - parser for Exonerate

=head1 SYNOPSIS

  # do not use this module directly, it is a driver for SearchIO

  use Bio::SearchIO;
  my $searchio = Bio::SearchIO->new(-file => 'file.exonerate',
                                   -format => 'exonerate');


  while( my $r = $searchio->next_result ) {
    print $r->query_name, "\n";
  }

=head1 DESCRIPTION

This is a driver for the SearchIO system for parsing Exonerate (Guy
Slater) output.  You can get Exonerate at
http://www.ebi.ac.uk/~guy/exonerate/
[until Guy puts up a Web reference,publication for it.]).

An optional parameter -min_intron is supported by the L<new>
initialization method.  This is if you run Exonerate with a different
minimum intron length (default is 30) the parser will be able to
detect the difference between standard deletions and an intron.  Still
some room to play with there that might cause this to get
misinterpreted that has not been fully tested or explored.

The VULGAR and CIGAR formats should be parsed okay now creating HSPs
where appropriate (so merging match states where appropriate rather
than breaking an HSP at each indel as it may have done in the past).
The GFF that comes from exonerate is still probably a better way to go
if you are doing protein2genome or est2genome mapping.
For example you can see this script:

### TODO: Jason, this link is dead, do we have an updated one?
http://fungal.genome.duke.edu/~jes12/software/scripts/process_exonerate_gff3.perl.txt

If your report contains both CIGAR and VULGAR lines only the first one
will processed for a given Query/Target pair.  If you preferentially
want to use VULGAR or CIGAR add one of these options when initializing
the SearchIO object.

    -cigar  => 1
OR
    -vulgar => 1

Or set them via these methods.

    $parser->cigar(1)
OR
    $parser->vulgar(1)



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


package Bio::SearchIO::exonerate;
use strict;
use vars qw(@STATES %MAPPING %MODEMAP $DEFAULT_WRITER_CLASS $MIN_INTRON);

use base qw(Bio::SearchIO);

%MODEMAP = ( 'ExonerateOutput' => 'result',
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

    'ExonerateOutput_program'   => 'RESULT-algorithm_name',
    'ExonerateOutput_query-def' => 'RESULT-query_name',
    'ExonerateOutput_query-desc'=> 'RESULT-query_description',
    'ExonerateOutput_query-len' => 'RESULT-query_length',
    );

$DEFAULT_WRITER_CLASS = 'Bio::SearchIO::Writer::HitTableWriter';

$MIN_INTRON=30; # This is the minimum intron size

=head2 new

 Title   : new
 Usage   : my $obj = Bio::SearchIO::exonerate->new();
 Function: Builds a new Bio::SearchIO::exonerate object
 Returns : an instance of Bio::SearchIO::exonerate
 Args    : -min_intron => somewhat obselete option, how to determine if a
                          an indel is an intron or a local gap.  Use VULGAR
                          rather than CIGAR to avoid this heuristic,default 30.
           -cigar       => 1   set this to 1 if you want to parse
                               CIGAR exclusively.
           -vulgar      => 1   set this to 1 if you want to parse VULGAR
                               exclusively, setting both to 1 will revert
                               to the default behavior of just parsing the
                               first line that it sees.

=cut

sub new {
    my ($class) = shift;
    my $self = $class->SUPER::new(@_);

    my ($min_intron,$cigar,
	$vulgar) = $self->_rearrange([qw(MIN_INTRON
					 CIGAR
					 VULGAR)], @_);
    if( $min_intron ) {
	$MIN_INTRON = $min_intron;
    }
    if( $cigar && $vulgar ) {
	$self->warn("cannot get HSPs from both CIGAR and VULGAR lines, will just choose whichever comes first (same as if you had chosen neither");
	$cigar = 0; $vulgar=0;
    }
    $self->cigar($cigar);
    $self->vulgar($vulgar);
    $self;
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

   $self->{'_last_data'} = '';
   my ($reporttype,$seenquery,$reportline);
   $self->start_document();
   my @hit_signifs;
   my $seentop;
   my (@q_ex, @m_ex, @h_ex); ## gc addition
   while( defined($_ = $self->_readline) ) {
       # warn( "Reading $_");
       if( /^\s*Query:\s+(\S+)\s*(.+)?/ ) {
	   if( $seentop ) {
	       $self->end_element({'Name' => 'ExonerateOutput'});
	       $self->_pushback($_);
	       return $self->end_document();
	   }
	   $seentop = 1;
	   my ($nm,$desc) = ($1,$2);
	   chomp($desc) if defined $desc;
	   $self->{'_result_count'}++;
	   $self->start_element({'Name' => 'ExonerateOutput'});
	   $self->element({'Name' => 'ExonerateOutput_query-def',
			   'Data' => $nm });
	   $self->element({'Name' => 'ExonerateOutput_query-desc',
			   'Data' => $desc });
	   $self->element({'Name' => 'ExonerateOutput_program',
			   'Data' => 'Exonerate' });
	   $self->{'_seencigar'} = 0;
	   $self->{'_vulgar'}    = 0;

       } elsif ( /^Target:\s+(\S+)\s*(.+)?/ ) {
	   my ($nm,$desc) = ($1,$2);
	   chomp($desc) if defined $desc;
	   $self->start_element({'Name' => 'Hit'});
	   $self->element({'Name' => 'Hit_id',
			   'Data' => $nm});
	   $self->element({'Name' => 'Hit_desc',
			   'Data' => $desc});
	   $self->{'_seencigar'} = 0;
	   $self->{'_vulgar'}    = 0;
       } elsif(  s/^vulgar:\s+(\S+)\s+         # query sequence id
		 (\d+)\s+(\d+)\s+([\-\+\.])\s+ # query start-end-strand
		 (\S+)\s+                      # target sequence id
		 (\d+)\s+(\d+)\s+([\-\+])\s+   # target start-end-strand
		 (-?\d+)\s+                    # score
		 //ox ) {
	   next if( $self->cigar || $self->{'_seencigar'});
	   $self->{'_vulgar'}++;
	   #
	   # Note from Ewan. This is ugly - copy and paste from
	   # cigar line parsing. Should unify somehow...
	   #
	   if( ! $self->within_element('result') ) {
	       $self->start_element({'Name' => 'ExonerateOutput'});
	       $self->element({'Name' => 'ExonerateOutput_query-def',
			       'Data' => $1 });
	   }
	   if( ! $self->within_element('hit') ) {
	       $self->start_element({'Name' => 'Hit'});
	       $self->element({'Name' => 'Hit_id',
			       'Data' => $5});
	   }

	   ## gc note:
	   ## $qe and $he are no longer used for calculating the ends,
	   ## just the $qs and $hs values and the alignment and insert lenghts
	   my ($qs,$qe,$qstrand) = ($2,$3,$4);
	   my ($hs,$he,$hstrand) = ($6,$7,$8);
	   my $score = $9;
#	   $self->element({'Name' => 'ExonerateOutput_query-len',
#			   'Data' => $qe});
#	   $self->element({'Name' => 'Hit_len',
#			   'Data' => $he});

	   ## gc note:
	   ## add one because these values are zero-based
	   ## this calculation was originally done lower in the code,
	   ## but it's clearer to do it just once at the start
	   my @rest = split;
	   my ($qbegin,$qend) = ('query-from', 'query-to');

	   if( $qstrand eq '-' ) {
	       $qstrand = -1; $qe++;
	   } else {
	       $qstrand = 1;
	       $qs++;
	   }
	   my ($hbegin,$hend) = ('hit-from', 'hit-to');

	   if( $hstrand eq '-' ) {
	       $hstrand = -1;
	       $he++;
	   } else {
	       $hstrand = 1;
	       $hs++;
	   }
	   # okay let's do this right and generate a set of HSPs
	   # from the cigar line/home/bio1/jes12/bin/exonerate  --model est2genome --bestn 1 t/data/exonerate_cdna.fa t/data/exonerate_genomic_rev.fa

	   my ($aln_len,$inserts,$deletes) = (0,0,0);
	   my ($laststate,@events,$gaps) =( '' );
	   while( @rest >= 3 ) {
	       my ($state,$len1,$len2) = (shift @rest, shift @rest, shift @rest);
	       #
	       # HSPs are only the Match cases; otherwise we just
	       # move the coordinates on by the correct amount
	       #

	       if( $state eq 'M' ) {
		   if( $laststate eq 'G' ) {
		       # merge gaps across Match states so the HSP
		       # goes across
		       $events[-1]->{$qend} = $qs + $len1*$qstrand - $qstrand;
		       $events[-1]->{$hend}   = $hs + $len2*$hstrand - $hstrand;
		       $events[-1]->{'gaps'} = $gaps;
		   } else {
		       push @events,
		       { 'score'     => $score,
			 'align-len' => $len1,
			 $qbegin => $qs,
			 $qend  => ($qs + $len1*$qstrand - $qstrand),
			 $hbegin => $hs,
			 $hend   => ($hs + $len2*$hstrand - $hstrand),
		     };
		   }
		   $gaps = 0;
	       } else {
		   $gaps = $len1 + $len2 if $state eq 'G';
	       }
	       $qs += $len1*$qstrand;
	       $hs += $len2*$hstrand;
	       $laststate= $state;
	   }
	   for my $event ( @events ) {
	       $self->start_element({'Name' => 'Hsp'});
	       while( my ($key,$val) = each %$event ) {
		   $self->element({'Name' => "Hsp_$key",
				   'Data' => $val});
	       }
	       $self->element({'Name' => 'Hsp_identity',
			       'Data' => 0});
	       $self->end_element({'Name' => 'Hsp'});
	   }

	   # end of hit
	   $self->element({'Name' => 'Hit_score',
			   'Data' => $score});
	   # issued end...
	   $self->end_element({'Name' => 'Hit'});
	   $self->end_element({'Name' => 'ExonerateOutput'});

	   return $self->end_document();

       } elsif(  s/^cigar:\s+(\S+)\s+          # query sequence id
		 (\d+)\s+(\d+)\s+([\-\+])\s+   # query start-end-strand
		 (\S+)\s+                      # target sequence id
		 (\d+)\s+(\d+)\s+([\-\+])\s+   # target start-end-strand
		 (-?\d+)\s+                    # score
		 //ox ) {
	   next if( $self->vulgar || $self->{'_seenvulgar'});
	   $self->{'_cigar'}++;

	   if( ! $self->within_element('result') ) {
	       $self->start_element({'Name' => 'ExonerateOutput'});
	       $self->element({'Name' => 'ExonerateOutput_query-def',
			       'Data' => $1 });
	   }
	   if( ! $self->within_element('hit') ) {
	       $self->start_element({'Name' => 'Hit'});
	       $self->element({'Name' => 'Hit_id',
			       'Data' => $5});
	   }
	   ## gc note:
	   ## $qe and $he are no longer used for calculating the ends,
	   ## just the $qs and $hs values and the alignment and insert lenghts
	   my ($qs,$qe,$qstrand) = ($2,$3,$4);
	   my ($hs,$he,$hstrand) = ($6,$7,$8);
	   my $score = $9;
#	   $self->element({'Name' => 'ExonerateOutput_query-len',
#			   'Data' => $qe});
#	   $self->element({'Name' => 'Hit_len',
#			   'Data' => $he});

	   my @rest = split;
	   if( $qstrand eq '-' ) {
	       $qstrand = -1;
	       ($qs,$qe) = ($qe,$qs); # flip-flop if we're on opp strand
	       $qs--; $qe++;
	   } else { $qstrand = 1; }
	   if( $hstrand eq '-' ) {
	       $hstrand = -1;
	       ($hs,$he) = ($he,$hs); # flip-flop if we're on opp strand
	       $hs--; $he++;
	   } else { $hstrand = 1; }
	   # okay let's do this right and generate a set of HSPs
	   # from the cigar line

	   ## gc note:
	   ## add one because these values are zero-based
	   ## this calculation was originally done lower in the code,
	   ## but it's clearer to do it just once at the start
	   $qs++; $hs++;

	   my ($aln_len,$inserts,$deletes) = (0,0,0);
	   while( @rest >= 2 ) {
	       my ($state,$len) = (shift @rest, shift @rest);
	       if( $state eq 'I' ) {
		   $inserts+=$len;
	       } elsif( $state eq 'D' ) {
		   if( $len >= $MIN_INTRON ) {
		       $self->start_element({'Name' => 'Hsp'});

		       $self->element({'Name' => 'Hsp_score',
				       'Data' => $score});
		       $self->element({'Name' => 'Hsp_align-len',
				       'Data' => $aln_len});
		       $self->element({'Name' => 'Hsp_identity',
				       'Data' => $aln_len -
					   ($inserts + $deletes)});

		       # HSP ends where the other begins
		       $self->element({'Name' => 'Hsp_query-from',
				       'Data' => $qs});
		       ## gc note:
		       ## $qs is now the start of the next hsp
		       ## the end of this hsp is 1 before this position
		       ## (or 1 after in case of reverse strand)
		       $qs += $aln_len*$qstrand;
		       $self->element({'Name' => 'Hsp_query-to',
				       'Data' => $qs - ($qstrand*1)});

		       $hs += $deletes*$hstrand;
		       $self->element({'Name' => 'Hsp_hit-from',
				       'Data' => $hs});
		       $hs += $aln_len*$hstrand;
		       $self->element({'Name' => 'Hsp_hit-to',
				       'Data' => $hs-($hstrand*1)});

		       $self->element({'Name' => 'Hsp_align-len',
				       'Data' => $aln_len + $inserts
					   + $deletes});
		       $self->element({'Name' => 'Hsp_identity',
				       'Data' => $aln_len });

		       $self->element({'Name' => 'Hsp_gaps',
				       'Data' => $inserts + $deletes});
		       $self->element({'Name' => 'Hsp_querygaps',
				       'Data' => $inserts});
		       $self->element({'Name' => 'Hsp_hitgaps',
				       'Data' => $deletes});

## gc addition start

		       $self->element({'Name' => 'Hsp_qseq',
				       'Data' => shift @q_ex,
				   });
		       $self->element({'Name' => 'Hsp_hseq',
				       'Data' => shift @h_ex,
				   });
		       $self->element({'Name' => 'Hsp_midline',
				       'Data' => shift @m_ex,
				   });
## gc addition end
		       $self->end_element({'Name' => 'Hsp'});

		       $aln_len = $inserts = $deletes = 0;
		   }
		   $deletes+=$len;
	       } else {
		   $aln_len += $len;
	       }
	   }
	   $self->start_element({'Name' => 'Hsp'});

## gc addition start

	   $self->element({'Name' => 'Hsp_qseq',
			   'Data' => shift @q_ex,
		       });
	   $self->element({'Name' => 'Hsp_hseq',
			   'Data' => shift @h_ex,
		       });
	   $self->element({'Name' => 'Hsp_midline',
			   'Data' => shift @m_ex,
		       });
## gc addition end

	   $self->element({'Name' => 'Hsp_score',
			   'Data' => $score});

	   $self->element({'Name' => 'Hsp_query-from',
			   'Data' => $qs});

	   $qs += $aln_len*$qstrand;
	   $self->element({'Name' => 'Hsp_query-to',
			   'Data' => $qs - ($qstrand*1)});

	   $hs += $deletes*$hstrand;
	   $self->element({'Name' => 'Hsp_hit-from',
			   'Data' => $hs});
	   $hs += $aln_len*$hstrand;
	   $self->element({'Name' => 'Hsp_hit-to',
			   'Data' => $hs -($hstrand*1)});

	   $self->element({'Name' => 'Hsp_align-len',
			   'Data' => $aln_len});

	   $self->element({'Name' => 'Hsp_identity',
			   'Data' => $aln_len - ($inserts + $deletes)});

	   $self->element({'Name' => 'Hsp_gaps',
			   'Data' => $inserts + $deletes});

	   $self->element({'Name' => 'Hsp_querygaps',
			   'Data' => $inserts});
	   $self->element({'Name' => 'Hsp_hitgaps',
			   'Data' => $deletes});
	   $self->end_element({'Name' => 'Hsp'});

	   $self->element({'Name' => 'Hit_score',
			   'Data' => $score});

	   $self->end_element({'Name' => 'Hit'});
	   $self->end_element({'Name' => 'ExonerateOutput'});

	   return $self->end_document();
       } else {
	   # skipping this line
       }
   }
   return $self->end_document() if( $seentop );
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
       if( $self->_eventHandler->will_handle($type) ) {
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

    if( $type = $MODEMAP{$nm} ) {
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
    $self->{'_reporttype'} = 'exonerate';
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

=head2 vulgar

 Title   : vulgar
 Usage   : $obj->vulgar($newval)
 Function: Get/Set flag, do you want to build HSPs from VULGAR string?
 Returns : value of vulgar (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub vulgar{
    my $self = shift;
    my $x = shift if @_;
    if( @_ ) {
	if( $_[0] && $self->{'_cigar'} ) {
	    $self->warn("Trying to set vulgar and cigar both to 1, must be either or");
	    $self->{'_cigar'}  = 0;
	    return $self->{'_vulgar'} = 0;
	}
    }
    return $self->{'_vulgar'};
}

=head2 cigar

 Title   : cigar
 Usage   : $obj->cigar($newval)
 Function: Get/Set boolean flag do you want to build HSPs from CIGAR strings?
 Returns : value of cigar (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub cigar{
    my $self = shift;
    my $x = shift if @_;
    if( @_ ) {
	if( $_[0] && $self->{'_vulgar'} ) {
	    $self->warn("Trying to set vulgar and cigar both to 1, must be either or");
	    $self->{'_vulgar'}  = 0;
	    return $self->{'_cigar'} = 0;
	}
    }
    return $self->{'_cigar'};
}

1;

