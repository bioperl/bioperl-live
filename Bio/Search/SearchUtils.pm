=head1 NAME

Bio::Search::SearchUtils - Utility functions for Bio::Search:: objects

=head1 SYNOPSIS

  # This module is just a collection of subroutines, not an object.

=head1 DESCRIPTION

The SearchUtils.pm module is a collection of subroutines used
primarily by Bio::Search::Hit::HitI objects for some of the additional
functionality, such as HSP tiling. Right now, the SearchUtils is just
a collection of methods, not an object.

=head1 AUTHOR 

Steve Chervitz E<lt>sac@bioperl.orgE<gt>

=head1 CONTRIBUTORS

Sendu Bala, bix@sendu.me.uk

=cut

package Bio::Search::SearchUtils;
use Bio::Root::Version;

use strict;

=head2 tile_hsps

 Usage     : tile_hsps( $sbjct );
           : This is called automatically by methods in Bio::Search::Hit::GenericHit 
           : that rely on having tiled data.
           :
           : If you are interested in getting data about the constructed HSP contigs:
           : my ($qcontigs, $scontigs) = Bio::Search::SearchUtils::tile_hsps($hit);
           : if (ref $qcontigs) {
           :    print STDERR "Query contigs:\n";
           :    foreach (@{$qcontigs}) {
           :         print "contig start is $_->{'start'}\n";
           :         print "contig stop is $_->{'stop'}\n";
           :    }
           : }
           : See below for more information about the contig data structure.
           :
 Purpose   : Collect statistics about the aligned sequences in a set of HSPs.
           : Calculates the following data across all HSPs: 
           :    -- total alignment length 
           :    -- total identical residues 
           :    -- total conserved residues
 Returns   : If there was only a single HSP (so no tiling was necessary)
               tile_hsps() returns a list of two non-zero integers.
             If there were multiple HSP, 
               tile_hsps() returns a list of two array references containin HSP contig data.
             The first array ref contains a list of HSP contigs on the query sequence.
             The second array ref contains a list of HSP contigs on the subject sequence.
             Each contig is a hash reference with the following data fields:
               'start' => start coordinate of the contig
               'stop'  => start coordinate of the contig
               'iden'  => number of identical residues in the contig
               'cons'  => number of conserved residues in the contig
               'strand'=> strand of the contig
               'frame' => frame of the contig
 Argument  : A Bio::Search::Hit::HitI object 
 Throws    : n/a
 Comments  :
           : This method performs more careful summing of data across
           : all HSPs in the Sbjct object. Only HSPs that are in the same strand 
           : and frame are tiled. Simply summing the data from all HSPs
           : in the same strand and frame will overestimate the actual 
           : length of the alignment if there is overlap between different HSPs 
           : (often the case).
           :
           : The strategy is to tile the HSPs and sum over the
           : contigs, collecting data separately from overlapping and
           : non-overlapping regions of each HSP. To facilitate this, the
           : HSP.pm object now permits extraction of data from sub-sections
           : of an HSP.
           : 
           : Additional useful information is collected from the results
           : of the tiling. It is possible that sub-sequences in
           : different HSPs will overlap significantly. In this case, it
           : is impossible to create a single unambiguous alignment by
           : concatenating the HSPs. The ambiguity may indicate the
           : presence of multiple, similar domains in one or both of the
           : aligned sequences. This ambiguity is recorded using the
           : ambiguous_aln() method.
           : 
           : This method does not attempt to discern biologically
           : significant vs. insignificant overlaps. The allowable amount of 
           : overlap can be set with the overlap() method or with the -OVERLAP
           : parameter used when constructing the Hit object.
           : 
           : For a given hit, both the query and the sbjct sequences are
           : tiled independently.
           : 
           :    -- If only query sequence HSPs overlap, 
           :          this may suggest multiple domains in the sbjct.
           :    -- If only sbjct sequence HSPs overlap, 
           :          this may suggest multiple domains in the query.
           :    -- If both query & sbjct sequence HSPs overlap, 
           :          this suggests multiple domains in both.
           :    -- If neither query & sbjct sequence HSPs overlap, 
           :          this suggests either no multiple domains in either
           :          sequence OR that both sequences have the same
           :          distribution of multiple similar domains.
           : 
           : This method can deal with the special case of when multiple
           : HSPs exactly overlap.
           : 
           : Efficiency concerns:
           :  Speed will be an issue for sequences with numerous HSPs.
           : 
 Bugs      : Currently, tile_hsps() does not properly account for
           : the number of non-tiled but overlapping HSPs, which becomes a problem
           : as overlap() grows. Large values overlap() may thus lead to 
           : incorrect statistics for some hits. For best results, keep overlap()
           : below 5 (DEFAULT IS 2). For more about this, see the "HSP Tiling and
           : Ambiguous Alignments" section in L<Bio::Search::Hit::GenericHit>.

See Also   : L<_adjust_contigs>(), L<Bio::Search::Hit::GenericHit|Bio::Search::Hit::GenericHit>

=cut

#--------------
sub tile_hsps {
#--------------
    my $sbjct = shift;

    #print STDERR "Calling tile_hsps(): $sbjct\n";
    #$sbjct->verbose(1);  # to activate debugging
    $sbjct->tiled_hsps(1);

    # changed to not rely on n() (which is unreliable here) --cjfields 4/6/10
    if( $sbjct->num_hsps == 0) { 		
	#print STDERR "_tile_hsps(): no hsps, nothing to tile! (", $sbjct->num_hsps, ")\n";
        _warn_about_no_hsps($sbjct);
        return (undef, undef);

    } elsif($sbjct->num_hsps == 1) {
        ## Simple summation scheme. Valid if there is only one HSP.
        #print STDERR "_tile_hsps(): single HSP, easy stats.\n";
        my $hsp = $sbjct->hsp;
        $sbjct->length_aln('query', $hsp->length('query'));
        $sbjct->length_aln('hit', $hsp->length('sbjct'));
        $sbjct->length_aln('total', $hsp->length('total'));
        $sbjct->matches( $hsp->matches() );
        $sbjct->gaps('query', $hsp->gaps('query'));
        $sbjct->gaps('sbjct', $hsp->gaps('sbjct'));

        _adjust_length_aln($sbjct);
        return (1, 1);
    } else {
        #print STDERR "Sbjct: _tile_hsps: summing multiple HSPs\n";
        $sbjct->length_aln('query', 0);
        $sbjct->length_aln('sbjct', 0);
        $sbjct->length_aln('total', 0); 
        $sbjct->matches( 0, 0);
        $sbjct->gaps('query', 0);
        $sbjct->gaps('hit', 0);
    }

    ## More than one HSP. Must tile HSPs.
#    print "\nTiling HSPs for $sbjct\n";
    my($hsp, $qstart, $qstop, $sstart, $sstop);
    my($frame, $strand, $qstrand, $sstrand);
    my(@qcontigs, @scontigs);
    my $qoverlap = 0;
    my $soverlap = 0;
    my $max_overlap = $sbjct->overlap;
    my $hit_qgaps = 0;
    my $hit_sgaps = 0;
    my $hit_len_aln = 0;
    my %start_stop;
    my $v = $sbjct->verbose;
    foreach $hsp ( $sbjct->hsps() ) {
	#$sbjct->debug( sprintf("  HSP: %s %d..%d\n",$hsp->query->seq_id, $hsp->query->start, $hsp->hit->end)) if $v > 0; #$hsp->str('query');
#	printf "  Length = %d; Identical = %d; Conserved = %d; Conserved(1-10): %d",$hsp->length, $hsp->length(-TYPE=>'iden'), 
#	$hsp->length(-TYPE=>'cons'),
#	$hsp->length(-TYPE=>'cons',
#		     -START=>0,-STOP=>10); 

	($qstart, $qstop) = $hsp->range('query');
	($sstart, $sstop) = $hsp->range('sbjct');
	$frame = $hsp->frame('hit');
	$frame = -1 unless defined $frame;
	
	($qstrand, $sstrand) = ($hsp->query->strand,
				$hsp->hit->strand);

        # Note: No correction for overlap.
	
	my ($qgaps, $sgaps)  = ($hsp->gaps('query'), $hsp->gaps('hit'));
	$hit_qgaps += $qgaps;
	$hit_sgaps += $sgaps;
	$hit_len_aln += $hsp->length;

	## Collect contigs in the query sequence.
 	$qoverlap += &_adjust_contigs('query', $hsp, $qstart, $qstop, 
				     \@qcontigs, $max_overlap, $frame, 
				     $qstrand);

	## Collect contigs in the sbjct sequence 
	#  (needed for domain data and gapped Blast).
	$soverlap += &_adjust_contigs('sbjct', $hsp, $sstart, $sstop, 
				     \@scontigs, $max_overlap, $frame, 
				     $sstrand);

	## Collect overall start and stop data for query and 
	#  sbjct over all HSPs.
	unless ( defined $start_stop{'qstart'} ) {
            $start_stop{'qstart'}  = $qstart;
            $start_stop{'qstop'}   = $qstop;
            $start_stop{'sstart'}  = $sstart;
            $start_stop{'sstop'}   = $sstop;
	} else {
	    $start_stop{'qstart'} = ($qstart < $start_stop{'qstart'} ? 
				     $qstart : $start_stop{'qstart'} );
	    $start_stop{'qstop'}  = ($qstop  > $start_stop{'qstop'} ? 
				     $qstop  : $start_stop{'qstop'} );
	    $start_stop{'sstart'} = ($sstart < $start_stop{'sstart'} ? 
				     $sstart : $start_stop{'sstart'} );
	    $start_stop{'sstop'}  = ($sstop  > $start_stop{'sstop'} ? 
				     $sstop  : $start_stop{'sstop'} );
	}	    
    }

    # Store the collected data in the Hit object
    $sbjct->gaps('query', $hit_qgaps);
    $sbjct->gaps('hit', $hit_sgaps);
    $sbjct->length_aln('total', $hit_len_aln);
    
    $sbjct->start('query',$start_stop{'qstart'});
    $sbjct->end('query', $start_stop{'qstop'});
    $sbjct->start('hit', $start_stop{'sstart'});
    $sbjct->end('hit', $start_stop{'sstop'});
    ## Collect data across the collected contigs.

    #$sbjct->debug( "\nQUERY CONTIGS:\n"."  gaps = $sbjct->{'_gaps_query'}\n");

    # Account for strand/frame.
    # Strategy: collect data on a per strand+frame basis and 
    #           save the most significant one.
    my (%qctg_dat);
    foreach (@qcontigs) {
        ($frame, $strand) = ($_->{'frame'}, $_->{'strand'});
        
        if( $v > 0 ) {
            #$sbjct->debug(sprintf( "$frame/$strand len is getting %d for %d..%d\n", 
            #			   ($_->{'stop'} - $_->{'start'} + 1), $_->{'start'}, $_->{'stop'}));
        }
        
        $qctg_dat{ "$frame$strand" }->{'length_aln_query'} += $_->{'stop'} - $_->{'start'} + 1;
        $qctg_dat{ "$frame$strand" }->{'totalIdentical'}   += $_->{'iden'};
        $qctg_dat{ "$frame$strand" }->{'totalConserved'}   += $_->{'cons'};
        $qctg_dat{ "$frame$strand" }->{'qstrand'}   = $strand;
    }

    # Find longest contig.
    my @sortedkeys = sort { $qctg_dat{$b}->{'length_aln_query'} 
			    <=> $qctg_dat{$a}->{'length_aln_query'} }
    keys %qctg_dat;

    # Save the largest to the sbjct:
    my $longest = $sortedkeys[0];
    #$sbjct->debug( "longest is ". $qctg_dat{ $longest }->{'length_aln_query'}. "\n");
    $sbjct->length_aln('query', $qctg_dat{ $longest }->{'length_aln_query'});
    $sbjct->matches($qctg_dat{ $longest }->{'totalIdentical'},
                    $qctg_dat{ $longest }->{'totalConserved'});
    $sbjct->strand('query', $qctg_dat{ $longest }->{'qstrand'});

    ## Collect data for sbjct contigs. Important for gapped Blast.
    ## The totalIdentical and totalConserved numbers will be the same
    ## as determined for the query contigs.

    #$sbjct->debug( "\nSBJCT CONTIGS:\n"."  gaps = ". $sbjct->gaps('sbjct'). "\n");
    my (%sctg_dat);
    foreach(@scontigs) {
	#$sbjct->debug("  sbjct contig: $_->{'start'} - $_->{'stop'}\n".
	#	     "         iden = $_->{'iden'}; cons = $_->{'cons'}\n");
	($frame, $strand) = ($_->{'frame'}, $_->{'strand'});
	$sctg_dat{ "$frame$strand" }->{'length_aln_sbjct'}   += $_->{'stop'} - $_->{'start'} + 1;
	$sctg_dat{ "$frame$strand" }->{'frame'}  = $frame;
	$sctg_dat{ "$frame$strand" }->{'sstrand'}  = $strand;
    }
    
    @sortedkeys = sort { $sctg_dat{ $b }->{'length_aln_sbjct'} 
			 <=> $sctg_dat{ $a }->{'length_aln_sbjct'} 
		     } keys %sctg_dat;

    # Save the largest to the sbjct:
    $longest = $sortedkeys[0];

    $sbjct->length_aln('sbjct', $sctg_dat{ $longest }->{'length_aln_sbjct'});
    $sbjct->frame( $sctg_dat{ $longest }->{'frame'} );
    $sbjct->strand('hit', $sctg_dat{ $longest }->{'sstrand'});

    if($qoverlap) {
	if($soverlap) { $sbjct->ambiguous_aln('qs'); 
			#$sbjct->debug("\n*** AMBIGUOUS ALIGNMENT: Query and Sbjct\n\n");
		    }
	else { $sbjct->ambiguous_aln('q');
	       #$sbjct->debug( "\n*** AMBIGUOUS ALIGNMENT: Query\n\n");
	   }
    } elsif($soverlap) { 
	$sbjct->ambiguous_aln('s'); 
	#$sbjct->debug( "\n*** AMBIGUOUS ALIGNMENT: Sbjct\n\n");
    }

    _adjust_length_aln($sbjct);

    return ( [@qcontigs], [@scontigs] );
}



# Title    : _adjust_length_aln  
# Usage    : n/a; internal use only; called by tile_hsps.
# Purpose  : Adjust length of aligment based on BLAST flavor.
# Comments : See comments in logica_length()
sub _adjust_length_aln {
    my $sbjct = shift;
    my $algo = $sbjct->algorithm;
    my $hlen = $sbjct->length_aln('sbjct');
    my $qlen = $sbjct->length_aln('query');

    $sbjct->length_aln('sbjct', logical_length($algo, 'sbjct', $hlen));
    $sbjct->length_aln('query', logical_length($algo, 'query', $qlen));
}

=head2 logical_length

 Usage     : logical_length( $alg_name, $seq_type, $length );
 Purpose   : Determine the logical length of an aligned sequence based on 
           : algorithm name and sequence type.
 Returns   : integer representing the logical aligned length.
 Argument  : $alg_name = name of algorigthm (e.g., blastx, tblastn)
           : $seq_type = type of sequence (e.g., query or hit)
           : $length = physical length of the sequence in the alignment.
 Throws    : n/a
 Comments  : This function is used to account for the fact that number of identities 
             and conserved residues is reported in peptide space while the query 
             length (in the case of BLASTX and TBLASTX) and/or the hit length 
             (in the case of TBLASTN and TBLASTX) are in nucleotide space.
             The adjustment affects the values reported by the various frac_XXX 
             methods in GenericHit and GenericHSP.

=cut

sub logical_length {
    my ($algo, $type, $len) = @_;
    my $logical = $len;
    if($algo =~ /^(?:PSI)?T(?:BLASTN|FAST(?:X|Y|XY))/oi ) {
        $logical = $len/3 if $type =~ /sbjct|hit|tot/i;
    } elsif($algo =~ /^(?:BLASTX|FAST(?:X|Y|XY))/oi ) {
        $logical = $len/3 if $type =~ /query|tot/i;
    } elsif($algo =~ /^TBLASTX/oi ) {
        $logical = $len/3;
    }
    return $logical;
}


#=head2 _adjust_contigs
#
# Usage     : n/a; internal function called by tile_hsps
# Purpose   : Builds HSP contigs for a given BLAST hit.
#           : Utility method called by _tile_hsps()
# Returns   : 
# Argument  : 
# Throws    : Exceptions propagated from Bio::Search::Hit::BlastHSP::matches()
#           : for invalid sub-sequence ranges.
# Status    : Experimental
# Comments  : This method supports gapped alignments through a patch by maj
#           : to B:S:HSP:HSPI::matches().
#           : It does not keep track of the number of HSPs that
#           : overlap within the amount specified by overlap().
#           : This will lead to significant tracking errors for large
#           : overlap values.
#
#See Also   : L<tile_hsps>(), L<Bio::Search::Hit::BlastHSP::matches|Bio::Search::Hit::BlastHSP>
#
#=cut

sub _adjust_contigs {
    my ($seqType, $hsp, $start, $stop, $contigs_ref, 
	$max_overlap, $frame, $strand) = @_;
    my $overlap = 0;
    my ($numID, $numCons);
    
    foreach (@$contigs_ref) {
        # Don't merge things unless they have matching strand/frame.
        next unless ($_->{'frame'} == $frame && $_->{'strand'} == $strand);
        
        # Test special case of a nested HSP. Skip it.
        if ($start >= $_->{'start'} && $stop <= $_->{'stop'}) {
            $overlap = 1; 
            next;
        }
        
        # Test for overlap at beginning of contig, or precedes consecutively
        if ($start < $_->{'start'} && $stop >= ($_->{'start'} + $max_overlap - 1)) {
            eval {
                ($numID, $numCons) = $hsp->matches(-SEQ   =>$seqType, 
                               -START => $start, 
                               -STOP  => $_->{'start'} - 1); 
		if ($numID eq '') {
		    $hsp->warn("\$hsp->matches() returned '' for number identical; setting to 0");
		    $numID = 0;
		}
		if ($numCons eq '') {
		    $hsp->warn("\$hsp->matches() returned '' for number conserved; setting to 0");
		    $numCons = 0;
		}
            };
            if($@) { warn "\a\n$@\n"; }
            else {
                $_->{'start'} = $start;	# Assign a new start coordinate to the contig
                $_->{'iden'} += $numID;	# and add new data to #identical, #conserved.
                $_->{'cons'} += $numCons;
                push(@{$_->{hsps}}, $hsp);
                $overlap     = 1; 
            }
        }
        
        # Test for overlap at end of contig, or follows consecutively
        if ($stop > $_->{'stop'} and $start <= ($_->{'stop'} - $max_overlap + 1)) {
            eval {
                ($numID,$numCons) = $hsp->matches(-SEQ   =>$seqType, 
                              -START => $_->{'stop'} + 1, 
                              -STOP  => $stop); 
		if ($numID eq '') {
		    $hsp->warn("\$hsp->matches() returned '' for number identical; setting to 0");
		    $numID = 0;
		}
		if ($numCons eq '') {
		    $hsp->warn("\$hsp->matches() returned '' for number conserved; setting to 0");
		    $numCons = 0;
		}
	    };
            if($@) { warn "\a\n$@\n"; }
            else {
                $_->{'stop'}  = $stop; # Assign a new stop coordinate to the contig
                $_->{'iden'} += $numID;	# and add new data to #identical, #conserved.
                $_->{'cons'} += $numCons;
                push(@{$_->{hsps}}, $hsp);
                $overlap    = 1; 
            }
        }
        
        last if $overlap;
    }
    
    if ($overlap && @$contigs_ref > 1) {
        ## Merge any contigs that now overlap
        my $max = $#{$contigs_ref};
        for my $i (0..$max) {
            ${$contigs_ref}[$i] || next;
            my ($i_start, $i_stop) = (${$contigs_ref}[$i]->{start}, ${$contigs_ref}[$i]->{stop});
            
            for my $u ($i+1..$max) {
                ${$contigs_ref}[$u] || next;
                my ($u_start, $u_stop) = (${$contigs_ref}[$u]->{start}, ${$contigs_ref}[$u]->{stop});
                
                if ($u_start < $i_start && $u_stop >= ($i_start + $max_overlap - 1)) {
                    # find the hsps within the contig that have sequence
                    # extending before $i_start
                    my ($ids, $cons) = (0, 0);
                    my $use_start = $i_start;
                    foreach my $hsp (sort { $b->end($seqType) <=> $a->end($seqType) } @{${$contigs_ref}[$u]->{hsps}}) {
                        my $hsp_start = $hsp->start($seqType);
                        $hsp_start < $use_start || next;
                        
                        my ($these_ids, $these_cons);
                        eval {
                            ($these_ids, $these_cons) = $hsp->matches(-SEQ => $seqType, -START => $hsp_start, -STOP => $use_start - 1);
			    if ($these_ids eq '') {
				$hsp->warn("\$hsp->matches() returned '' for number identical; setting to 0");
				$these_ids = 0;
			    }
			    if ($these_cons eq '') {
				$hsp->warn("\$hsp->matches() returned '' for number conserved; setting to 0");
				$these_cons = 0;
			    }
			};
                        if($@) { warn "\a\n$@\n"; }
                        else {
                            $ids  += $these_ids;
                            $cons += $these_cons;
                        }
                        
                        last if $hsp_start == $u_start;
                        $use_start = $hsp_start;
                    }
                    ${$contigs_ref}[$i]->{start} = $u_start;
                    ${$contigs_ref}[$i]->{'iden'} += $ids;
                    ${$contigs_ref}[$i]->{'cons'} += $cons;
                    push(@{${$contigs_ref}[$i]->{hsps}}, @{${$contigs_ref}[$u]->{hsps}});
                    
                    ${$contigs_ref}[$u] = undef;
                }
                elsif ($u_stop > $i_stop && $u_start <= ($i_stop - $max_overlap + 1)) {
                    # find the hsps within the contig that have sequence
                    # extending beyond $i_stop
                    my ($ids, $cons) = (0, 0);
                    my $use_stop = $i_stop;
                    foreach my $hsp (sort { $a->start($seqType) <=> $b->start($seqType) } @{${$contigs_ref}[$u]->{hsps}}) {
                        my $hsp_end = $hsp->end($seqType);
                        $hsp_end > $use_stop || next;
                        
                        my ($these_ids, $these_cons);
                        eval {
                            ($these_ids, $these_cons) = $hsp->matches(-SEQ => $seqType, -START => $use_stop + 1, -STOP => $hsp_end);
			    if ($these_ids eq '') {
				$hsp->warn("\$hsp->matches() returned '' for number identical; setting to 0");
				$these_ids = 0;
			    }
			    if ($these_cons eq '') {
				$hsp->warn("\$hsp->matches() returned '' for number conserved; setting to 0");
				$these_cons = 0;
			    }
                        };
                        if($@) { warn "\a\n$@\n"; }
                        else {
                            $ids  += $these_ids;
                            $cons += $these_cons;
                        }
                        
                        last if $hsp_end == $u_stop;
                        $use_stop = $hsp_end;
                    }
                    ${$contigs_ref}[$i]->{'stop'}  = $u_stop;
                    ${$contigs_ref}[$i]->{'iden'} += $ids;
                    ${$contigs_ref}[$i]->{'cons'} += $cons;
                    push(@{${$contigs_ref}[$i]->{hsps}}, @{${$contigs_ref}[$u]->{hsps}});
                    
                    ${$contigs_ref}[$u] = undef;
                }
                elsif ($u_start >= $i_start && $u_stop <= $i_stop) {
                    # nested, drop this contig
                    #*** ideally we might do some magic to keep the stats of the
                    #    better hsp...
                    ${$contigs_ref}[$u] = undef;
                }
            }
        }
        
        my @merged;
        foreach (@$contigs_ref) {
            push(@merged, $_ || next);
        }
        @{$contigs_ref} = @merged;
    }
    elsif (! $overlap) {
        ## If there is no overlap, add the complete HSP data.
        ($numID,$numCons) = $hsp->matches(-SEQ=>$seqType);
	if ($numID eq '') {
	    $hsp->warn("\$hsp->matches() returned '' for number identical; setting to 0");
	    $numID = 0;
	}
	if ($numCons eq '') {
	    $hsp->warn("\$hsp->matches() returned '' for number conserved; setting to 0");
	    $numCons = 0;
	}

        push @$contigs_ref, {'start' =>$start, 'stop' =>$stop,
			     'iden'  =>$numID, 'cons' =>$numCons,
			     'strand'=>$strand,'frame'=>$frame,'hsps'=>[$hsp]};
    }
    
    return $overlap;
}

=head2 get_exponent

 Usage     : &get_exponent( number );
 Purpose   : Determines the power of 10 exponent of an integer, float, 
           : or scientific notation number.
 Example   : &get_exponent("4.0e-206");
           : &get_exponent("0.00032");
           : &get_exponent("10.");
           : &get_exponent("1000.0");
           : &get_exponent("e+83");
 Argument  : Float, Integer, or scientific notation number
 Returns   : Integer representing the exponent part of the number (+ or -).
           : If argument == 0 (zero), return value is "-999".
 Comments  : Exponents are rounded up (less negative) if the mantissa is >= 5.
           : Exponents are rounded down (more negative) if the mantissa is <= -5.

=cut

sub get_exponent {
    my $data = shift;

    my($num, $exp) = split /[eE]/, $data;

    if( defined $exp) { 
	$num = 1 if not $num;
	$num >= 5 and $exp++;
	$num <= -5 and $exp--;
    } elsif( $num == 0) {
	$exp = -999;
    } elsif( not $num =~ /\./) {
	$exp = CORE::length($num) -1;
    } else {
	$exp = 0;
	$num .= '0' if $num =~ /\.$/;
	my ($c);
	my $rev = 0;
	if($num !~ /^0/) {
	    $num = reverse($num);
	    $rev = 1;
	}
	do { $c = chop($num);
	     $c == 0 && $exp++; 
	 } while( $c ne '.');

	$exp = -$exp if $num == 0 and not $rev;
	$exp -= 1 if $rev;
    }
    return $exp;
}

=head2 collapse_nums

 Usage     : @cnums = collapse_nums( @numbers );
 Purpose   : Collapses a list of numbers into a set of ranges of consecutive terms:
           : Useful for condensing long lists of consecutive numbers.
           :  EXPANDED:
           :     1 2 3 4 5 6 10 12 13 14 15 17 18 20 21 22 24 26 30 31 32
           :  COLLAPSED:
           :     1-6 10 12-15 17 18 20-22 24 26 30-32
 Argument  : List of numbers sorted numerically.
 Returns   : List of numbers mixed with ranges of numbers (see above).
 Throws    : n/a

See Also   : L<Bio::Search::Hit::BlastHit::seq_inds()|Bio::Search::Hit::BlastHit>

=cut

sub collapse_nums {
# This is probably not the slickest connectivity algorithm, but will do for now.
    my @a = @_;
    my ($from, $to, $i, @ca, $consec);
    
    $consec = 0;
    for($i=0; $i < @a; $i++) {
	not $from and do{ $from = $a[$i]; next; };
    # pass repeated positions (gap inserts)
    next if $a[$i] == $a[$i-1];
	if($a[$i] == $a[$i-1]+1) {
	    $to = $a[$i];
	    $consec++;
	} else {
	    if($consec == 1) { $from .= ",$to"; }
	    else { $from .= $consec>1 ? "\-$to" : ""; }
	    push @ca, split(',', $from);
	    $from =  $a[$i];
	    $consec = 0;
	    $to = undef;
	}
    }
    if(defined $to) {
	if($consec == 1) { $from .= ",$to"; }
	else { $from .= $consec>1 ? "\-$to" : ""; }
    }
    push @ca, split(',', $from) if $from;

    @ca;
}


=head2 strip_blast_html

 Usage     : $boolean = &strip_blast_html( string_ref );
           : This method is exported.
 Purpose   : Removes HTML formatting from a supplied string.
           : Attempts to restore the Blast report to enable
           : parsing by Bio::SearchIO::blast.pm
 Returns   : Boolean: true if string was stripped, false if not.
 Argument  : string_ref = reference to a string containing the whole Blast
           :              report containing HTML formatting.
 Throws    : Croaks if the argument is not a scalar reference.
 Comments  : Based on code originally written by Alex Dong Li
           : (ali@genet.sickkids.on.ca).
           : This method does some Blast-specific stripping 
           : (adds back a '>' character in front of each HSP 
           : alignment listing).
           :   
           : THIS METHOD IS VERY SENSITIVE TO BLAST FORMATTING CHANGES!
           :
           : Removal of the HTML tags and accurate reconstitution of the
           : non-HTML-formatted report is highly dependent on structure of
           : the HTML-formatted version. For example, it assumes that first 
           : line of each alignment section (HSP listing) starts with a
           : <a name=..> anchor tag. This permits the reconstruction of the 
           : original report in which these lines begin with a ">".
           : This is required for parsing.
           :
           : If the structure of the Blast report itself is not intended to
           : be a standard, the structure of the HTML-formatted version
           : is even less so. Therefore, the use of this method to
           : reconstitute parsable Blast reports from HTML-format versions
           : should be considered a temorary solution.

=cut

sub strip_blast_html {
      # This may not best way to remove html tags. However, it is simple.
      # it won't work under following conditions:
      #    1) if quoted > appears in a tag  (does this ever happen?)
      #    2) if a tag is split over multiple lines and this method is
      #       used to process one line at a time.
      
    my ($string_ref) = shift;

    ref $string_ref eq 'SCALAR' or 
	croak ("Can't strip HTML: ".
	       "Argument is should be a SCALAR reference not a ${\ref $string_ref}\n");

    my $str = $$string_ref;
    my $stripped = 0;

    # Removing "<a name =...>" and adding the '>' character for 
    # HSP alignment listings.
    $str =~ s/(\A|\n)<a name ?=[^>]+> ?/>/sgi and $stripped = 1;

    # Removing all "<>" tags. 
    $str =~ s/<[^>]+>|&nbsp//sgi and $stripped = 1;

    # Re-uniting any lone '>' characters.
    $str =~ s/(\A|\n)>\s+/\n\n>/sgi and $stripped = 1;

    $$string_ref = $str;
    $stripped;
}

=head2 result2hash

 Title    : result2hash
 Usage    : my %data = &Bio::Search::SearchUtils($result)
 Function : converts ResultI data to simple hash
 Returns  : hash
 Args     : ResultI
 Note     : used mainly as a utility for running SearchIO tests

=cut

sub result2hash {
    my ($result) = @_;
    my %hash;
    $hash{'query_name'} = $result->query_name;
    my $hitcount = 1;
    my $hspcount = 1;
    foreach my $hit ( $result->hits ) {
	$hash{"hit$hitcount\_name"}   =  $hit->name;
	# only going to test order of magnitude
	# too hard as these don't always match
#	$hash{"hit$hitcount\_signif"} =  
#	    ( sprintf("%.0e",$hit->significance) =~ /e\-?(\d+)/ );
	$hash{"hit$hitcount\_bits"}   =  sprintf("%d",$hit->bits);

	foreach my $hsp ( $hit->hsps ) {
	    $hash{"hsp$hspcount\_bits"}   = sprintf("%d",$hsp->bits);
	    # only going to test order of magnitude
 	    # too hard as these don't always match
#	    $hash{"hsp$hspcount\_evalue"} =  
#		( sprintf("%.0e",$hsp->evalue) =~ /e\-?(\d+)/ );
	    $hash{"hsp$hspcount\_qs"}     = $hsp->query->start;
	    $hash{"hsp$hspcount\_qe"}     = $hsp->query->end;
	    $hash{"hsp$hspcount\_qstr"}   = $hsp->query->strand;
	    $hash{"hsp$hspcount\_hs"}     = $hsp->hit->start;
	    $hash{"hsp$hspcount\_he"}     = $hsp->hit->end;
	    $hash{"hsp$hspcount\_hstr"}   = $hsp->hit->strand;

	    #$hash{"hsp$hspcount\_pid"}     = sprintf("%d",$hsp->percent_identity);
	    #$hash{"hsp$hspcount\_fid"}     = sprintf("%.2f",$hsp->frac_identical);
	    $hash{"hsp$hspcount\_gaps"}    = $hsp->gaps('total');
	    $hspcount++;
	}
	$hitcount++;
    }
    return %hash;
}

sub _warn_about_no_hsps {
    my $hit = shift;
    my $prev_func=(caller(1))[3];
    $hit->warn("There is no HSP data for hit '".$hit->name."'.\n".
               "You have called a method ($prev_func)\n".
               "that requires HSP data and there was no HSP data for this hit,\n".
               "most likely because it was absent from the BLAST report.\n".
               "Note that by default, BLAST lists alignments for the first 250 hits,\n".
               "but it lists descriptions for 500 hits. If this is the case,\n".
               "and you care about these hits, you should re-run BLAST using the\n".
               "-b option (or equivalent if not using blastall) to increase the number\n".
               "of alignments.\n"
              );
}

1;
