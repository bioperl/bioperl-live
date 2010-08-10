=head1 NAME

Bio::Search::BlastUtils - Utility functions for Bio::Search:: BLAST objects

=head1 SYNOPSIS

 # This module is just a collection of subroutines, not an object.

See L<Bio::Search::Hit::BlastHit>.

=head1 DESCRIPTION

The BlastUtils.pm module is a collection of subroutines used primarily by
Bio::Search::Hit::BlastHit objects for some of the additional
functionality, such as HSP tiling. Right now, the BlastUtils is just a
collection of methods, not an object, and it's tightly coupled to
Bio::Search::Hit::BlastHit. A goal for the future is to generalize it
to work based on the Bio::Search interfaces, then it can work with any
objects that implements them.

=head1 AUTHOR 

Steve Chervitz E<lt>sac@bioperl.orgE<gt>

=cut

#'

package Bio::Search::BlastUtils;
use Bio::Root::Version;


=head2 tile_hsps

 Usage     : tile_hsps( $sbjct );
           : This is called automatically by Bio::Search::Hit::BlastHit 
           : during object construction or
           : as needed by methods that rely on having tiled data.
 Purpose   : Collect statistics about the aligned sequences in a set of HSPs.
           : Calculates the following data across all HSPs: 
           :    -- total alignment length 
           :    -- total identical residues 
           :    -- total conserved residues
 Returns   : n/a
 Argument  : A Bio::Search::Hit::BlastHit object 
 Throws    : n/a
 Comments  :
 	   : This method is *strongly* coupled to Bio::Search::Hit::BlastHit
 	   : (it accesses BlastHit data members directly).
 	   : TODO: Re-write this to the Bio::Search::Hit::HitI interface.
 	   :
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
           : parameter used when constructing the Blast & Sbjct objects. 
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
           : Ambiguous Alignments" section in L<Bio::Search::Hit::BlastHit>.

See Also   : L<_adjust_contigs>(), L<Bio::Search::Hit::BlastHit|Bio::Search::Hit::BlastHit>

=cut

#--------------
sub tile_hsps {
#--------------
    my $sbjct = shift;

    $sbjct->{'_tile_hsps'} = 1;
    $sbjct->{'_gaps_query'} = 0;
    $sbjct->{'_gaps_sbjct'} = 0;

    ## Simple summation scheme. Valid if there is only one HSP.
    if((defined($sbjct->{'_n'}) and $sbjct->{'_n'} == 1) or $sbjct->num_hsps == 1) {
	my $hsp = $sbjct->hsp;
	$sbjct->{'_length_aln_query'} = $hsp->length('query');
	$sbjct->{'_length_aln_sbjct'} = $hsp->length('sbjct');
	$sbjct->{'_length_aln_total'} = $hsp->length('total');
	($sbjct->{'_totalIdentical'},$sbjct->{'_totalConserved'}) = $hsp->matches();
	$sbjct->{'_gaps_query'} = $hsp->gaps('query');
	$sbjct->{'_gaps_sbjct'} = $hsp->gaps('sbjct');

#	print "_tile_hsps(): single HSP, easy stats.\n";
	return;
    } else {
#	print STDERR "Sbjct: _tile_hsps: summing multiple HSPs\n";
	$sbjct->{'_length_aln_query'} = 0;
	$sbjct->{'_length_aln_sbjct'} = 0;
	$sbjct->{'_length_aln_total'} = 0;
	$sbjct->{'_totalIdentical'}   = 0;
	$sbjct->{'_totalConserved'}   = 0;
    }

    ## More than one HSP. Must tile HSPs.
#    print "\nTiling HSPs for $sbjct\n";
    my($hsp, $qstart, $qstop, $sstart, $sstop);
    my($frame, $strand, $qstrand, $sstrand);
    my(@qcontigs, @scontigs);
    my $qoverlap = 0;
    my $soverlap = 0;
    my $max_overlap = $sbjct->{'_overlap'};

    foreach $hsp ($sbjct->hsps()) {
#	printf "  HSP: %s\n%s\n",$hsp->name, $hsp->str('query');
#	printf "  Length = %d; Identical = %d; Conserved = %d; Conserved(1-10): %d",$hsp->length, $hsp->length(-TYPE=>'iden'), $hsp->length(-TYPE=>'cons'), $hsp->length(-TYPE=>'cons',-START=>0,-STOP=>10); 
	($qstart, $qstop) = $hsp->range('query');
	($sstart, $sstop) = $hsp->range('sbjct');
	$frame = $hsp->frame('hit');
	$frame = -1 unless defined $frame;
	($qstrand, $sstrand) = $hsp->strand;

	my ($qgaps, $sgaps)  = $hsp->gaps();
	$sbjct->{'_gaps_query'} += $qgaps;
	$sbjct->{'_gaps_sbjct'} += $sgaps;

	$sbjct->{'_length_aln_total'} += $hsp->length;
	## Collect contigs in the query sequence.
	$qoverlap = &_adjust_contigs('query', $hsp, $qstart, $qstop, \@qcontigs, $max_overlap, $frame, $qstrand);

	## Collect contigs in the sbjct sequence (needed for domain data and gapped Blast).
	$soverlap = &_adjust_contigs('sbjct', $hsp, $sstart, $sstop, \@scontigs, $max_overlap, $frame, $sstrand);

	## Collect overall start and stop data for query and sbjct over all HSPs.
	if(not defined $sbjct->{'_queryStart'}) {
	    $sbjct->{'_queryStart'} = $qstart;
	    $sbjct->{'_queryStop'}  = $qstop;
	    $sbjct->{'_sbjctStart'} = $sstart;
	    $sbjct->{'_sbjctStop'}  = $sstop;
	} else {
	    $sbjct->{'_queryStart'} = ($qstart < $sbjct->{'_queryStart'} ? $qstart : $sbjct->{'_queryStart'});
	    $sbjct->{'_queryStop'}  = ($qstop  > $sbjct->{'_queryStop'}  ? $qstop  : $sbjct->{'_queryStop'});
	    $sbjct->{'_sbjctStart'} = ($sstart < $sbjct->{'_sbjctStart'} ? $sstart : $sbjct->{'_sbjctStart'});
	    $sbjct->{'_sbjctStop'}  = ($sstop  > $sbjct->{'_sbjctStop'}  ? $sstop  : $sbjct->{'_sbjctStop'});
	}	    
    }

    ## Collect data across the collected contigs.

#    print "\nQUERY CONTIGS:\n";
#    print "  gaps = $sbjct->{'_gaps_query'}\n";

    # TODO: Account for strand/frame issue!
    # Strategy: collect data on a per strand+frame basis and save the most significant one.
    my (%qctg_dat);
    foreach(@qcontigs) {
#	print "  query contig: $_->{'start'} - $_->{'stop'}\n";
#	print "         iden = $_->{'iden'}; cons = $_->{'cons'}\n";
	($frame, $strand) = ($_->{'frame'}, $_->{'strand'});
	$qctg_dat{ "$frame$strand" }->{'length_aln_query'} += $_->{'stop'} - $_->{'start'} + 1;
	$qctg_dat{ "$frame$strand" }->{'totalIdentical'}   += $_->{'iden'};
	$qctg_dat{ "$frame$strand" }->{'totalConserved'}   += $_->{'cons'};
	$qctg_dat{ "$frame$strand" }->{'qstrand'}   = $strand;
    }

    # Find longest contig.
    my @sortedkeys = reverse sort { $qctg_dat{ $a }->{'length_aln_query'} <=> $qctg_dat{ $b }->{'length_aln_query'} } keys %qctg_dat;

    # Save the largest to the sbjct:
    my $longest = $sortedkeys[0];
    $sbjct->{'_length_aln_query'} = $qctg_dat{ $longest }->{'length_aln_query'};
    $sbjct->{'_totalIdentical'}   = $qctg_dat{ $longest }->{'totalIdentical'};
    $sbjct->{'_totalConserved'}   = $qctg_dat{ $longest }->{'totalConserved'};
    $sbjct->{'_qstrand'} = $qctg_dat{ $longest }->{'qstrand'};

    ## Collect data for sbjct contigs. Important for gapped Blast.
    ## The totalIdentical and totalConserved numbers will be the same
    ## as determined for the query contigs.

#    print "\nSBJCT CONTIGS:\n";
#    print "  gaps = $sbjct->{'_gaps_sbjct'}\n";

    my (%sctg_dat);
    foreach(@scontigs) {
#	print "  sbjct contig: $_->{'start'} - $_->{'stop'}\n";
#	print "         iden = $_->{'iden'}; cons = $_->{'cons'}\n";
	($frame, $strand) = ($_->{'frame'}, $_->{'strand'});
	$sctg_dat{ "$frame$strand" }->{'length_aln_sbjct'}   += $_->{'stop'} - $_->{'start'} + 1;
	$sctg_dat{ "$frame$strand" }->{'frame'}  = $frame;
	$sctg_dat{ "$frame$strand" }->{'sstrand'}  = $strand;
    }

    @sortedkeys = reverse sort { $sctg_dat{ $a }->{'length_aln_sbjct'} <=> $sctg_dat{ $b }->{'length_aln_sbjct'} } keys %sctg_dat;

    # Save the largest to the sbjct:
    $longest = $sortedkeys[0];

    $sbjct->{'_length_aln_sbjct'} = $sctg_dat{ $longest }->{'length_aln_sbjct'};
    $sbjct->{'_frame'} = $sctg_dat{ $longest }->{'frame'};
    $sbjct->{'_sstrand'} = $sctg_dat{ $longest }->{'sstrand'};

    if($qoverlap) {
	if($soverlap) { $sbjct->ambiguous_aln('qs'); 
#			print "\n*** AMBIGUOUS ALIGNMENT: Query and Sbjct\n\n";
		      }
	else { $sbjct->ambiguous_aln('q');
#	       print "\n*** AMBIGUOUS ALIGNMENT: Query\n\n";
	   }
    } elsif($soverlap) { 
	$sbjct->ambiguous_aln('s'); 
#	print "\n*** AMBIGUOUS ALIGNMENT: Sbjct\n\n";
    }

    # Adjust length based on BLAST flavor.
    my $prog = $sbjct->algorithm;
    if($prog eq 'TBLASTN') {
	$sbjct->{'_length_aln_sbjct'} /= 3;
    } elsif($prog eq 'BLASTX' ) {
	$sbjct->{'_length_aln_query'} /= 3;
    } elsif($prog eq 'TBLASTX') {
	$sbjct->{'_length_aln_query'} /= 3;
	$sbjct->{'_length_aln_sbjct'} /= 3;
    }
}



=head2 _adjust_contigs

 Usage     : n/a; called automatically during object construction.
 Purpose   : Builds HSP contigs for a given BLAST hit.
           : Utility method called by _tile_hsps()
 Returns   : 
 Argument  : 
 Throws    : Exceptions propagated from Bio::Search::Hit::BlastHSP::matches()
           : for invalid sub-sequence ranges.
 Status    : Experimental
 Comments  : This method does not currently support gapped alignments.
           : Also, it does not keep track of the number of HSPs that
           : overlap within the amount specified by overlap().
           : This will lead to significant tracking errors for large
           : overlap values.

See Also   : L<tile_hsps>(), L<Bio::Search::Hit::BlastHSP::matches|Bio::Search::Hit::BlastHSP>

=cut

#-------------------
sub _adjust_contigs {
#-------------------
    my ($seqType, $hsp, $start, $stop, $contigs_ref, $max_overlap, $frame, $strand) = @_;

    my $overlap = 0;
    my ($numID, $numCons);

#    print STDERR "Testing $seqType data: HSP (${\$hsp->name});  $start, $stop, strand=$strand, frame=$frame\n"; 
    foreach(@$contigs_ref) {
#	print STDERR "  Contig: $_->{'start'} - $_->{'stop'}, strand=$_->{'strand'}, frame=$_->{'frame'}, iden= $_->{'iden'}, cons= $_->{'cons'}\n";

	# Don't merge things unless they have matching strand/frame.
	next unless ($_->{'frame'} == $frame and $_->{'strand'} == $strand);

	## Test special case of a nested HSP. Skip it.
	if($start >= $_->{'start'} and $stop <= $_->{'stop'}) { 
#	    print STDERR "----> Nested HSP. Skipping.\n";
	    $overlap = 1; 
	    next;
	}

	## Test for overlap at beginning of contig.
	if($start < $_->{'start'} and $stop > ($_->{'start'} + $max_overlap)) { 
#	    print STDERR "----> Overlaps beg: existing beg,end: $_->{'start'},$_->{'stop'}, new beg,end: $start,$stop\n";
	    # Collect stats over the non-overlapping region.
	    eval {
		($numID, $numCons) = $hsp->matches(-SEQ   =>$seqType, 
						   -START =>$start, 
						   -STOP  =>$_->{'start'}-1); 
	    };
	    if($@) { warn "\a\n$@\n"; }
	    else {
		$_->{'start'} = $start; # Assign a new start coordinate to the contig
		$_->{'iden'} += $numID; # and add new data to #identical, #conserved.
		$_->{'cons'} += $numCons;
		$overlap     = 1; 
	    }
	}

	## Test for overlap at end of contig.
	if($stop > $_->{'stop'} and $start < ($_->{'stop'} - $max_overlap)) { 
#	    print STDERR "----> Overlaps end: existing beg,end: $_->{'start'},$_->{'stop'}, new beg,end: $start,$stop\n";
	    # Collect stats over the non-overlapping region.
	    eval {
		($numID,$numCons) = $hsp->matches(-SEQ   =>$seqType, 
						  -START =>$_->{'stop'}, 
						  -STOP  =>$stop); 
	    };
	    if($@) { warn "\a\n$@\n"; }
	    else {
		$_->{'stop'}  = $stop;  # Assign a new stop coordinate to the contig
		$_->{'iden'} += $numID; # and add new data to #identical, #conserved.
		$_->{'cons'} += $numCons;
		$overlap    = 1; 
	    }
	}
	$overlap && do {
#		print STDERR " New Contig data:\n";
#		print STDERR "  Contig: $_->{'start'} - $_->{'stop'}, iden= $_->{'iden'}, cons= $_->{'cons'}\n";
		last;
	    };
    }
    ## If there is no overlap, add the complete HSP data.
    !$overlap && do {
#	print STDERR "No overlap. Adding new contig.\n";
	($numID,$numCons) = $hsp->matches(-SEQ=>$seqType); 
	push @$contigs_ref, {'start'=>$start, 'stop'=>$stop,
			     'iden'=>$numID,  'cons'=>$numCons,
			     'strand'=>$strand, 'frame'=>$frame};
    };
    $overlap;
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

#------------------
sub get_exponent {
#------------------
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

#------------------
sub collapse_nums {
#------------------
# This is probably not the slickest connectivity algorithm, but will do for now.
    my @a = @_;
    my ($from, $to, $i, @ca, $consec);
    
    $consec = 0;
    for($i=0; $i < @a; $i++) {
	not $from and do{ $from = $a[$i]; next; };
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

#--------------------
sub strip_blast_html {
#--------------------
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


1;


