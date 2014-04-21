#
# BioPerl module for Bio::Assembly::ContigAnalysis
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Robson francisco de Souza <rfsouza@citri.iq.usp.br>
#
# Copyright Robson Francisco de Souza
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Assembly::ContigAnalysis - 
    Perform analysis on sequence assembly contigs.

=head1 SYNOPSIS

    # Module loading
    use Bio::Assembly::ContigAnalysis;

    # Assembly loading methods
    my $ca = Bio::Assembly::ContigAnalysis->new( -contig=>$contigOBJ );

    my @lcq = $ca->low_consensus_quality;
    my @hqd = $ca->high_quality_discrepancies;
    my @ss  = $ca->single_strand_regions;

=head1 DESCRIPTION

A contig is as a set of sequences, locally aligned to each other, when
the sequences in a pair may be aligned. It may also include a
consensus sequence. Bio::Assembly::ContigAnalysis is a module
holding a collection of methods to analyze contig objects. It was
developed around the Bio::Assembly::Contig implementation of contigs and
can not work with another contig interface.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Robson Francisco de Souza

Email: rfsouza@citri.iq.usp.br

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Assembly::ContigAnalysis;

use strict;

use base qw(Bio::Root::Root);

=head1 Object creator

=head2 new

 Title     : new
 Usage     : my $contig = Bio::Assembly::ContigAnalysis->new(-contig=>$contigOBJ);
 Function  : Creates a new contig analysis object
 Returns   : Bio::Assembly::ContigAnalysis
 Args      :
             -contig : a Bio::Assembly::Contig object

=cut

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($contigOBJ) = $self->_rearrange([qw(CONTIG)],@args);
    unless ($contigOBJ->isa("Bio::Assembly::Contig")) {
	$self->throw("ContigAnal works only on Bio::Assembly::Contig objects\n");
    }

    $self->{'_objref'} = $contigOBJ;
    return $self;
}

=head1 Analysis methods

=head2 high_quality_discrepancies

 Title     : high_quality_discrepancies
 Usage     : my $sfc = $ContigAnal->high_quality_discrepancies();
 Function  : 

             Locates all high quality discrepancies among aligned
             sequences and the consensus sequence.

             Note: see Bio::Assembly::Contig POD documentation,
             section "Coordinate System", for a definition of
             available types. Default coordinate system type is
             "gapped consensus", i.e. consensus sequence (with gaps)
             coordinates. If limits are not specified, the entire
             alignment is analyzed.

 Returns   : Bio::SeqFeature::Collection
 Args      : optional arguments are
             -threshold : cutoff value for low quality (minimum high quality)
                          Default: 40
             -ignore    : number of bases that will not be analysed at
                          both ends of contig aligned elements
                          Default: 5
             -start     : start of interval that will be analyzed
             -end       : start of interval that will be analyzed
             -type      : coordinate system type for interval

=cut

sub high_quality_discrepancies {
    my ($self,@args) = shift; # Package reference

    my ($threshold,$ignore,$start,$end,$type) = 
	$self->_rearrange([qw(THRESHOLD IGNORE START END TYPE)],@args);

    # Defining default threhold and HQD_ignore
    $threshold  = 40 unless (defined($threshold));
    $ignore = 5  unless (defined($ignore));
    $type = 'gapped consensus' unless (defined($type));

    # Changing input coordinates system (if needed)
    if (defined $start && $type ne 'gapped consensus') {
	$start = $self->{'_objref'}->change_coord($type,'gapped consensus',$start);
    } elsif (!defined $start) {
	$start = 1;
    }
    if (defined $end && $type ne 'gapped consensus') {
	$end = $self->{'_objref'}->change_coord($type,'gapped consensus',$end);
    } elsif (!defined $end) {
	$end = $self->{'_objref'}->get_consensus_length();
    }

    # Scanning each read sequence and the contig sequence and
    # adding discrepancies to Bio::SeqFeature::Collection
    my @seqIDs = $self->{'_objref'}->get_seq_ids(-start=>$start,
						 -end=>$end,
						 -type=>$type);
    my $consensus = $self->{'_objref'}->get_consensus_sequence()->seq;

    my @HQD = ();
    foreach my $seqID (@seqIDs) {
	# Setting aligned read sub-sequence limits and loading data
	my $seq  = $self->{'_objref'}->get_seq_by_name($seqID);
	my $qual = $self->{'_objref'}->get_qual_by_name($seqID);
	unless (defined $qual) {
	    $self->warn("Can't correctly evaluate HQD without aligned sequence qualities for $seqID");
	    next;
	}
	my $sequence = $seq->seq;
	my @quality  = @{ $qual->qual };

	# Scanning the aligned region of each read
	my $seq_ix = 0;
	my $coord = $self->{'_objref'}->get_seq_feat_by_tag($seq,"_align_clipping:$seqID");
    if (!$coord) {
        $self->warn("Read $seqID has no alignment coordinates; considered low quality.\nSkipping...");
        next;
    }
	my ($astart,$aend) = ($coord->start,$coord->end);
	$astart = $astart + $ignore; # Redefining limits to evaluate HQDs (jump $ignore at start)
	$aend   = $aend   - $ignore; # Redefining limits to evaluate HQDs (stop $ignore before end)

	my ($d_start,$d_end,$i);
	for ($i=$astart-1; $i<=$aend-1; $i++) {
	    # Changing coordinate $i+1 from 'gapped consensus' mode to "aligned $seqID" (coordinate $seq_ix)
	    $seq_ix = $self->{'_objref'}->change_coord('gapped consensus',"aligned $seqID",$i+1);
	    next unless (($i >= $start) && ($i <= $end));

	    my $r_base = uc(substr($sequence,$seq_ix-1,1));
	    my $c_base = uc(substr($consensus,$i,1));

	    # Discrepant region start: store $d_start and $type
	    (!defined($d_start) &&
	     ($r_base ne $c_base) &&
	     ($quality[$seq_ix-1] >= $threshold)) && do {
		 $d_start = $self->{'_objref'}->change_coord('gapped consensus','ungapped consensus',$i+1);
		 #print $seqID," ",$r_base," ",$i+1," ",$c_base," ",$contig_ix-1," ",$quality[$i]," $type\n";
		 next;
	     };

	    # Quality change or end of discrepant region: store limits and undef $d_start
	    if (defined($d_start) &&
		(($quality[$seq_ix-1] < $threshold) ||
		 (uc($r_base) eq uc($c_base)))) {
		$d_end = $self->{'_objref'}->change_coord('gapped consensus','ungapped consensus',$i);
		#print $seqID," ",$r_base," ",$i+1," ",$c_base," ",$contig_ix-1," ",$quality[$i]," $type\n";
		push(@HQD, Bio::SeqFeature::Generic->new(-primary=>"high_quality_discrepancy:$seqID",
							 -start=>$d_start,
							 -end=>$d_end,
							 -strand=>$seq->strand()) );
		$d_start = undef;
		next;
	    }
	} # for ($i=$astart-1; $i<=$aend-1; $i++)

	# Loading discrepancies located at sub-sequence end, if any.
	if (defined($d_start)) {
	    $d_end = $self->{'_objref'}->change_coord('gapped consensus','ungapped consensus',$i);
	    push(@HQD, Bio::SeqFeature::Generic->new(-primary=>"high_quality_discrepancy:$seqID",
						     -start=>$d_start,
						     -end=>$d_end,
						     -strand=>$seq->strand()) );
	}
    } # foreach my $seqID (@seqIDs)

    return @HQD;
}

=head2 low_consensus_quality

 Title     : low_consensus_quality
 Usage     : my $sfc = $ContigAnal->low_consensus_quality();
 Function  : Locates all low quality regions in the consensus
 Returns   : an array of Bio::SeqFeature::Generic objects
 Args      : optional arguments are
             -threshold : cutoff value for low quality (minimum high quality)
                          Default: 25
             -start     : start of interval that will be analyzed
             -end       : start of interval that will be analyzed
             -type      : coordinate system type for interval

=cut

sub low_consensus_quality {
    my ($self,@args) = shift; # Packege reference

    my ($threshold,$start,$end,$type) = 
	$self->_rearrange([qw(THRESHOLD START END TYPE)],@args);

    # Setting default value for threshold
    $threshold = 25 unless (defined($threshold));

    # Loading qualities
    my @quality = @{ $self->{'_objref'}->get_consensus_quality()->qual };

    # Changing coordinates to gap mode noaln (consed: consensus without alignments)
    $start = 1 unless (defined($start));
    if (defined $start && defined $type && ($type ne 'gapped consensus')) {
	$start = $self->{'objref'}->change_coord($type,'gapped consensus',$start);
	$end   = $self->{'objref'}->change_coord($type,'gapped consensus',$end) if (defined($end));
    }
    $end = $self->{'_objref'}->get_consensus_length unless (defined $end);

    # Scanning @quality vector and storing intervals limits with base qualities less then
    # the threshold value
    my ($lcq_start);
    my ($i,@LCQ);
    for ($i=$start-1; $i<=$end-1; $i++) {
#	print $quality[$i],"\t",$i,"\n";
	if (!defined($lcq_start) && (($quality[$i] <= $threshold) || ($quality[$i] == 98))) {
	    $lcq_start = $i+1;
	} elsif (defined($lcq_start) && ($quality[$i] > $threshold)) {
	    $lcq_start  = $self->{'_objref'}->change_coord('gapped consensus','ungapped consensus',$lcq_start);
	    my $lcq_end = $self->{'_objref'}->change_coord('gapped consensus','ungapped consensus',$i);
	    push(@LCQ, Bio::SeqFeature::Generic->new(-start=>$lcq_start,
						     -end=>$lcq_end,
						     -primary=>'low_consensus_quality') );
	    $lcq_start = undef;
	}
    }

    if (defined $lcq_start) {
	$lcq_start  = $self->{'_objref'}->change_coord('gapped consensus','ungapped consensus',$lcq_start);
	my $lcq_end = $self->{'_objref'}->change_coord('gapped consensus','ungapped consensus',$i);
	push(@LCQ, Bio::SeqFeature::Generic->new(-start=>$lcq_start,
						 -end=>$lcq_end,
						 -primary=>'low_consensus_quality') );
    }

    return @LCQ;
}

=head2 not_confirmed_on_both_strands

 Title     : low_quality_consensus
 Usage     : my $sfc = $ContigAnal->low_quality_consensus();
 Function  : 

             Locates all regions whose consensus bases were not
             confirmed by bases from sequences aligned in both
             orientations, i.e., in such regions, no bases in aligned
             sequences of either +1 or -1 strand agree with the
             consensus bases.

 Returns   : an array of Bio::SeqFeature::Generic objects
 Args      : optional arguments are
             -start : start of interval that will be analyzed
             -end   : start of interval that will be analyzed
             -type  : coordinate system type for interval

=cut

sub not_confirmed_on_both_strands {
    my ($self,@args) = shift; # Package reference

    my ($start,$end,$type) = 
	$self->_rearrange([qw(START END TYPE)],@args);

    # Changing coordinates to default system 'align' (contig sequence with alignments)
    $start = 1 unless (defined($start));
    if (defined($type) && ($type ne 'gapped consensus')) {
	$start = $self->{'_objref'}->change_coord($type,'gapped consensus',$start);
	$end   = $self->{'_objref'}->change_coord($type,'gapped consensus',$end) if (defined($end));
    }
    $end = $self->{'_objref'}->get_consensus_length unless (defined($end));

    # Scanning alignment
    my %confirmed = (); # If ($confirmed{$orientation}[$i] > 0) then $i is confirmed in $orientation strand
    my ($i);
    my $consensus = $self->{'_objref'}->get_consensus_sequence()->seq;
    foreach my $seqID ($self->{'_objref'}->get_seq_ids) {
	# Setting aligned read sub-sequence limits and loading data
	my $seq = $self->{'_objref'}->get_seq_by_name($seqID);
	my $sequence = $seq->seq;

	# Scanning the aligned regions of each read and registering confirmed sites
	my $contig_ix = 0;
	my $coord = $self->{'_objref'}->get_seq_feat_by_tag($seq,"_align_clipping:$seqID");
	my ($astart,$aend,$orientation) = ($coord->start,$coord->end,$coord->strand);
	$astart = $self->{'_objref'}->change_coord('gapped consensus',"aligned $seqID",$astart);
	$aend   = $self->{'_objref'}->change_coord('gapped consensus',"aligned $seqID",$aend);
	for ($i=$astart-1; $i<=$aend-1; $i++) {
	    # $i+1 in 'align' mode is $contig_ix
	    $contig_ix = $self->{'_objref'}->change_coord("aligned $seqID",'gapped consensus',$i+1);
	    next unless (($contig_ix >= $start) && ($contig_ix <= $end));
	    my $r_base = uc(substr($sequence,$i,1));
	    my $c_base = uc(substr($consensus,$contig_ix-1,1));
	    if ($c_base eq '-') {
		$confirmed{$orientation}[$contig_ix] = -1;
	    } elsif (uc($r_base) eq uc($c_base)) { # Non discrepant region found
		$confirmed{$orientation}[$contig_ix]++;
	    }
	} # for ($i=$astart-1; $i<=$aend-1; $i++)
    } # foreach $seqID (@reads)

    # Locating non-confirmed aligned regions for each orientation in $confirmed registry
    my ($orientation);
    my @NCBS = ();
    foreach $orientation (keys %confirmed) {
	my ($ncbs_start,$ncbs_end);

	for ($i=$start; $i<=$end; $i++) {
	    if (!defined($ncbs_start) &&
		(!defined($confirmed{$orientation}[$i]) || ($confirmed{$orientation}[$i] == 0))) {
		$ncbs_start = $self->{'_objref'}->change_coord('gapped consensus','ungapped consensus',$i);
	    } elsif (defined($ncbs_start) &&
		     defined($confirmed{$orientation}[$i]) &&
		     ($confirmed{$orientation}[$i] > 0)) {
		$ncbs_end   = $self->{'_objref'}->change_coord('gapped consensus','ungapped consensus',$i-1);
		push(@NCBS, Bio::SeqFeature::Generic->new(-start=>$ncbs_start,
							  -end=>$ncbs_end,
							  -strand=>$orientation,
							  -primary=>"not_confirmed_on_both_strands") );
		$ncbs_start = undef;
	    }
	}

	if (defined($ncbs_start)) { # NCBS at the end of contig
	    $ncbs_end = $self->{'_objref'}->change_coord('gapped consensus','ungapped consensus',$end);
	    push(@NCBS, Bio::SeqFeature::Generic->new(-start=>$ncbs_start,
						      -end=>$ncbs_end,
						      -strand=>$orientation,
						      -primary=>'not_confirmed_on_both_strands') );
	}
    }

    return @NCBS;
}

=head2 single_strand

 Title     : single_strand
 Usage     : my $sfc = $ContigAnal->single_strand();
 Function  : 

             Locates all regions covered by aligned sequences only in
             one of the two strands, i.e., regions for which aligned
             sequence's strand() method returns +1 or -1 for all
             sequences.

 Returns   : an array of Bio::SeqFeature::Generic objects
 Args      : optional arguments are
             -start : start of interval that will be analyzed
             -end   : start of interval that will be analyzed
             -type  : coordinate system type for interval

=cut

#'
sub single_strand {
    my ($self,@args) = shift; # Package reference

    my ($start,$end,$type) = 
	$self->_rearrange([qw(START END TYPE)],@args);

    # Changing coordinates to gap mode align (consed: consensus sequence with alignments)
    $type  = 'gapped consensus' unless(defined($type));
    $start = 1 unless (defined($start));
    if (defined($type) && $type ne 'gapped consensus') {
	$start = $self->{'objref'}->change_coord($type,'gapped consensus',$start);
	$end   = $self->{'objref'}->change_coord($type,'gapped consensus',$end) if (defined($end));
    }
    ($end) = $self->{'_objref'}->get_consensus_length unless (defined($end));

    # Loading complete list of coordinates for aligned sequences
    my $sfc = $self->{'_objref'}->get_features_collection();
    my @forward = grep { $_->primary_tag =~ /^_aligned_coord:/ } 
    $sfc->features_in_range(-start=>$start,
			    -end=>$end,
			    -contain=>0,
			    -strand=>1,
			    -strandmatch=>'strong');
    my @reverse = grep { $_->primary_tag =~ /^_aligned_coord:/ } 
    $sfc->features_in_range(-start=>$start,
			    -end=>$end,
			    -contain=>0,
			    -strand=>-1,
			    -strandmatch=>'strong');
    # Merging overlapping features
    @forward = $self->_merge_overlapping_features(@forward);
    @reverse = $self->_merge_overlapping_features(@reverse);

    # Finding single stranded regions
    my ($length) = $self->{'_objref'}->get_consensus_length;
    $length  = $self->{'_objref'}->change_coord('gapped consensus','ungapped consensus',$length);
    @forward = $self->_complementary_features_list(1,$length,@forward);
    @reverse = $self->_complementary_features_list(1,$length,@reverse);

    my @SS = ();
    foreach my $feat (@forward, @reverse) {
	$feat->primary_tag('single_strand_region');
	push(@SS,$feat);
    }

    return @SS;
}

=head1 Internal Methods

=head2 _merge_overlapping_features

 Title     : _merge_overlapping_features
 Usage     : my @feat = $ContigAnal->_merge_overlapping_features(@features);
 Function  : Merge all overlapping features into features
             that hold original features as sub-features
 Returns   : array of Bio::SeqFeature::Generic objects
 Args      : array of Bio::SeqFeature::Generic objects

=cut

sub _merge_overlapping_features {
    my ($self,@feat) = @_;

    $self->throw_not_implemented();
}

=head2 _complementary_features_list

 Title     : _complementary_features_list
 Usage     : @feat = $ContigAnal->_complementary_features_list($start,$end,@features);
 Function  : Build a list of features for regions
             not covered by features in @features array
 Returns   : array of Bio::SeqFeature::Generic objects
 Args      : 
             $start    : [integer] start of first output feature
             $end      : [integer] end of last output feature
             @features : array of Bio::SeqFeature::Generic objects

=cut

sub _complementary_features_list {
    my ($self,$start,$end,@feat) = @_;

    $self->throw_not_implemented();
}

1;

__END__
