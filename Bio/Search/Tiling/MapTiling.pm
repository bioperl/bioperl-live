#
# BioPerl module for Bio::Search::Tiling::MapTiling
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Mark A. Jensen <maj@fortinbras.us>
#
# Copyright Mark A. Jensen
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Tiling::MapTiling - An implementation of an HSP tiling
algorithm, with methods to obtain frequently-requested statistics

=head1 SYNOPSIS

 # get a BLAST $hit from somewhere, then
 $tiling = Bio::Search::Tiling::MapTiling->new($hit);

 # stats
 $numID = $tiling->identities();
 $numCons = $tiling->conserved();
 $query_length = $tiling->length('query');
 $subject_length = $tiling->length('subject'); # or...
 $subject_length = $tiling->length('hit');

 # get a visual on the coverage map
 print $tiling->coverage_map_as_text('query',$context,'LEGEND');

 # tilings
 $context = $tiling->_context( -type => 'subject', -strand=> 1, -frame=>1);
 @covering_hsps_for_subject = $tiling->next_tiling('subject',$context);
 $context = $tiling->_context( -type => 'query', -strand=> -1, -frame=>0);
 @covering_hsps_for_query   = $tiling->next_tiling('query', $context);

=head1 DESCRIPTION

Frequently, users want to use a set of high-scoring pairs (HSPs)
obtained from a BLAST or other search to assess the overall level of
identity, conservation, or coverage represented by matches between a
subject and a query sequence. Because a set of HSPs frequently
describes multiple overlapping sequence fragments, a simple summation of
statistics over the HSPs will generally overestimate those
statistics. To obtain an accurate estimate of global hit statistics, a
'tiling' of HSPs onto either the subject or the query sequence must be
performed, in order to properly correct for this. 

This module will execute a tiling algorithm on a given hit based on an
interval decomposition I'm calling the "coverage map". Internal object
methods compute the various statistics, which are then stored in
appropriately-named public object attributes. See
L<Bio::Search::Tiling::MapTileUtils> for more info on the algorithm. 

=head2 STRAND/FRAME CONTEXTS

In BLASTX, TBLASTN, and TBLASTX reports, strand and frame information
are reported for the query, subject, or query and subject,
respectively, for each HSP. Tilings for these sequence types are only
meaningful when they include HSPs in the same strand and frame, or 
"context". So, in these situations, the context must be specified
in the method calls or the methods will throw. 

Contexts are specified as strings: C<[ 'all' | [m|p][_|0|1|2] ]>, where
C<all> = all HSPs (will throw if context must be specified), C<m> = minus
strand, C<p> = plus strand, and C<_> = no frame info, C<0,1,2> = respective
(absolute) frame. The L</_make_context_key> method will convert a (strand,
frame) specification to a context string, e.g.:

    $context = $self->_context(-type=>'query', -strand=>-1, -frame=>-2);

returns C<m2>.

The contexts present among the HSPs in a hit are identified and stored
for convenience upon object construction. These are accessed off the
object with the L</contexts> method. If contexts don't apply for the
given report, this returns C<('all')>. 

=head1 TILED ALIGNMENTS

The experimental method L<ALIGNMENTS/get_tiled_alns> will use a tiling
to concatenate tiled hsps into a series of L<Bio::SimpleAlign>
objects:

 @alns = $tiling->get_tiled_alns($type, $context);

Each alignment contains two sequences with ids 'query' and 'subject',
and consists of a concatenation of tiling HSPs which overlap or are
directly adjacent. The alignment are returned in C<$type> sequence
order. When HSPs overlap, the alignment sequence is taken from the HSP
which comes first in the coverage map array.

The sequences in each alignment contain features (even though they are
L<Bio::LocatableSeq> objects) which map the original query/subject
coordinates to the new alignment sequence coordinates. You can
determine the original BLAST fragments this way:

 $aln = ($tiling->get_tiled_alns)[0];
 $qseq = $aln->get_seq_by_id('query');
 $hseq = $aln->get_seq_by_id('subject');
 foreach my $feat ($qseq->get_SeqFeatures) {
    $org_start = ($feat->get_tag_values('query_start'))[0];
    $org_end = ($feat->get_tag_values('query_end'))[0];
    # original fragment as represented in the tiled alignment:
    $org_fragment = $feat->seq;
 }
 foreach my $feat ($hseq->get_SeqFeatures) {
    $org_start = ($feat->get_tag_values('subject_start'))[0];
    $org_end = ($feat->get_tag_values('subject_end'))[0];
    # original fragment as represented in the tiled alignment:
    $org_fragment = $feat->seq;
 }

=head1 DESIGN NOTE

The major calculations are made just-in-time, and then memoized. So,
for example, for a given MapTiling object, a coverage map would
usually be calculated only once (for the query), and at most twice (if
the subject perspective is also desired), and then only when a
statistic is first accessed. Afterward, the map and/or any statistic
is read from storage. So feel free to call the statistic methods
frequently if it suits you.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Mark A. Jensen

Email maj -at- fortinbras -dot- us

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Search::Tiling::MapTiling;
use strict;
use warnings;

# Object preamble - inherits from Bio::Root::Root
#use lib '../../..';

use Bio::Root::Root;
use Bio::Search::Tiling::TilingI;
use Bio::Search::Tiling::MapTileUtils;
use Bio::LocatableSeq;

# use base qw(Bio::Root::Root Bio::Search::Tiling::TilingI);
use base qw(Bio::Root::Root Bio::Search::Tiling::TilingI);

=head1 CONSTRUCTOR

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Search::Tiling::GenericTiling();
 Function: Builds a new Bio::Search::Tiling::GenericTiling object 
 Returns : an instance of Bio::Search::Tiling::GenericTiling
 Args    : -hit    => $a_Bio_Search_Hit_HitI_object
           general filter function:
           -hsp_filter => sub { my $this_hsp = shift; 
                                ...;
                                return 1 if $wanted;
                                return 0; }

=cut

sub new {
    my $class = shift;
    my @args = @_;
    my $self = $class->SUPER::new(@args);
    my($hit, $filter) = $self->_rearrange( [qw( HIT HSP_FILTER)],@args );

    $self->throw("HitI object required") unless $hit;
    $self->throw("Argument must be HitI object") unless ( ref $hit && $hit->isa('Bio::Search::Hit::HitI') );
    $self->{hit} = $hit;
    $self->_set_attributes();
    $self->{"_algorithm"} = $hit->algorithm;

    my @hsps = $hit->hsps;
    # apply filter function if requested
    if ( defined $filter ) {
	if ( ref($filter) eq 'CODE' ) {
	    @hsps = map { $filter->($_) ? $_ : () } @hsps;
	}
	else {
	    $self->warn("-filter is not a coderef; ignoring");
	}
    }
    
    # identify available contexts
    for my $t (qw( query hit )) {
	my %contexts;
	for my $i (0..$#hsps) {
	    my $ctxt = $self->_context(
		-type => $t,
		-strand => $hsps[$i]->strand($t),
		-frame  => $hsps[$i]->frame($t));
	    $contexts{$ctxt} ||= [];
	    push @{$contexts{$ctxt}}, $i;
	}
	$self->{"_contexts_${t}"} = \%contexts;
    }

    $self->warn("No HSPs present in hit after filtering") unless (@hsps);
    $self->hsps(\@hsps);
    return $self;
}

# a tiling is based on the set of hsps contained in a single hit.
# check all the boundaries - zero hsps, one hsp, all disjoint hsps

=head1 TILING ITERATORS

=head2 next_tiling

 Title   : next_tiling
 Usage   : @hsps = $self->next_tiling($type);
 Function: Obtain a tiling: a minimal set of HSPs covering the $type
           ('hit', 'subject', 'query') sequence
 Example :
 Returns : an array of HSPI objects
 Args    : scalar $type: one of 'hit', 'subject', 'query', with
           'subject' an alias for 'hit'

=cut

sub next_tiling{
    my $self = shift;
    my ($type, $context) = @_;
    $self->_check_type_arg(\$type);
    $self->_check_context_arg($type, \$context);
    return $self->_tiling_iterator($type, $context)->();
}

=head2 rewind_tilings

 Title   : rewind_tilings
 Usage   : $self->rewind_tilings($type)
 Function: Reset the next_tilings($type) iterator
 Example :
 Returns : True on success
 Args    : scalar $type: one of 'hit', 'subject', 'query';
           default is 'query'

=cut

sub rewind_tilings{
    my $self = shift;
    my ($type,$context) = @_;
    $self->_check_type_arg(\$type);
    $self->_check_context_arg($type, \$context);
    return $self->_tiling_iterator($type, $context)->('REWIND');
}

=head1 ALIGNMENTS

=head2 get_tiled_alns()

 Title   : get_tiled_alns
 Usage   : @alns = $tiling->get_tiled_alns($type, $context)
 Function: Use a tiling to construct a minimal set of alignment 
           objects covering the region specified by $type/$context
           by splicing adjacent HSP tiles
 Returns : an array of Bio::SimpleAlign objects; see Note below
 Args    : scalar $type: one of 'hit', 'subject', 'query'
           default is 'query'
           scalar $context: strand/frame context string
           Following $type and $context, an array of 
           ordered, tiled HSP objects can be specified; this is 
           the tiling that will directly the alignment construction
           default -- the first tiling provided by a tiling iterator
 Notes   : Each returned alignment is a concatenation of adjacent tiles.
           The set of alignments will cover all regions described by the 
           $type/$context pair in the hit. The pair of sequences in each 
           alignment have ids 'query' and 'subject', and each sequence 
           possesses SeqFeatures that map the original query or subject 
           coordinates to the sequence coordinates in the tiled alignment.
           
=cut

sub get_tiled_alns {
    my $self = shift;
    my ($type, $context) = @_;
    $self->_check_type_arg(\$type);
    $self->_check_context_arg($type, \$context);
    my $t = shift; # first arg after type/context, arrayref to a tiling
    my @tiling;
    if ($t && (ref($t) eq 'ARRAY')) {
	@tiling = @$t;
    }
    else { # otherwise, take the first tiling available

	@tiling = $self->_make_tiling_iterator($type,$context)->(); 
    }
    my @ret;

    my @map = $self->coverage_map($type, $context);
    my @intervals = map {$_->[0]} @map; # disjoint decomp
    # divide into disjoint covering groups
    my @groups = covering_groups(@intervals);

    require Bio::SimpleAlign;
    require Bio::SeqFeature::Generic;
    # cut hsp aligns along the disjoint decomp
    # look for gaps...or add gaps?
    my ($q_start, $h_start, $q_strand, $h_strand);
    # build seqs
    for my $grp (@groups) {
	my $taln = Bio::SimpleAlign->new();
	my (@qfeats, @hfeats);
	my $query_string = '';
	my $hit_string = '';
	my ($qlen,$hlen) = (0,0);
	my ($qinc, $hinc, $qstart, $hstart);
	for my $intvl (@$grp) {
	    # following just chooses the first available hsp containing the
	    # interval -- this is arbitrary, could be smarter.
	    my $aln = ( containing_hsps($intvl, @tiling) )[0]->get_aln;
	    my $qseq = $aln->get_seq_by_pos(1);
	    my $hseq = $aln->get_seq_by_pos(2);
	    $qstart ||= $qseq->start;
	    $hstart ||= $hseq->start;
	    # calculate the slice boundaries
	    my ($beg, $end);
	    for ($type) {
		/query/ && do {
		    $beg = $aln->column_from_residue_number($qseq->id, $intvl->[0]);
		    $end = $aln->column_from_residue_number($qseq->id, $intvl->[1]);
		    last;
		};
		/subject|hit/ && do {
		    $beg = $aln->column_from_residue_number($hseq->id, $intvl->[0]);
		    $end = $aln->column_from_residue_number($hseq->id, $intvl->[1]);
		    last;
		};
	    }
	    $aln = $aln->slice($beg, $end);
	    $qseq = $aln->get_seq_by_pos(1);
	    $hseq = $aln->get_seq_by_pos(2);
	    $qinc = $qseq->length - $qseq->num_gaps($Bio::LocatableSeq::GAP_SYMBOLS);
	    $hinc = $hseq->length - $hseq->num_gaps($Bio::LocatableSeq::GAP_SYMBOLS);

	    push @qfeats, Bio::SeqFeature::Generic->new(
		-start => $qlen+1,
		-end => $qlen+$qinc,
		-strand => $qseq->strand,
		-primary => 'query',
		-source_tag => 'BLAST',
		-display_name => 'query coordinates',
		-tag => { query_id => $qseq->id,
			  query_desc => $qseq->desc,
			  query_start => $qstart + (($qseq->strand && $qseq->strand < 0) ? -1 : 1)*$qlen,
			  query_end => $qstart + (($qseq->strand && $qseq->strand < 0) ? -1 : 1)*($qlen+$qinc-1),
		}
		);
	    push @hfeats, Bio::SeqFeature::Generic->new(
		-start => $hlen+1,
		-end => $hlen+$hinc,
		-strand => $hseq->strand,
		-primary => 'subject/hit',
		-source_tag => 'BLAST',
		-display_name => 'subject/hit coordinates',
		-tag => { subject_id => $hseq->id,
			  subject_desc => $hseq->desc,
			  subject_start => $hstart + (($hseq->strand && $hseq->strand < 0) ? -1 : 1)*$hlen,
			  subject_end => $hstart + (($hseq->strand && $hseq->strand < 0) ? -1 : 1)*($hlen+$hinc-1)
		}
		);
	    $query_string .= $qseq->seq;
	    $hit_string .= $hseq->seq;
	    $qlen += $qinc;
	    $hlen += $hinc;
	}
	# create the LocatableSeqs and add the features to each
	# then add the seqs to the new aln
	# note in MapTileUtils, Bio::FeatureHolderI methods have been
	# mixed into Bio::LocatableSeq
	my $qseq = Bio::LocatableSeq->new( -id => 'query',
					   -seq => $query_string);
	$qseq->add_SeqFeature(@qfeats);
	my $hseq = Bio::LocatableSeq->new( -id => 'subject',
					   -seq => $hit_string );
	$hseq->add_SeqFeature(@hfeats);
	$taln->add_seq($qseq);
	$taln->add_seq($hseq);
	push @ret, $taln;
    }
    
    return @ret;
}

=head1 STATISTICS

=head2 identities

 Title   : identities
 Usage   : $tiling->identities($type, $action, $context)
 Function: Retrieve the calculated number of identities for the invocant
 Example : 
 Returns : value of identities (a scalar)
 Args    : scalar $type: one of 'hit', 'subject', 'query'
           default is 'query'
           option scalar $action: one of 'exact', 'est', 'fast', 'max'
           default is 'exact'
           option scalar $context: strand/frame context string
 Note    : getter only

=cut

sub identities{
    my $self = shift;
    my ($type, $action, $context) = @_;
    $self->_check_type_arg(\$type);
    $self->_check_action_arg(\$action);
    $self->_check_context_arg($type, \$context);
    if (!defined $self->{"identities_${type}_${action}_${context}"}) {
	$self->_calc_stats($type, $action, $context);
    }
    return $self->{"identities_${type}_${action}_${context}"};
}

=head2 conserved

 Title   : conserved
 Usage   : $tiling->conserved($type, $action)
 Function: Retrieve the calculated number of conserved sites for the invocant
 Example : 
 Returns : value of conserved (a scalar)
 Args    : scalar $type: one of 'hit', 'subject', 'query'
           default is 'query'
           option scalar $action: one of 'exact', 'est', 'fast', 'max'
           default is 'exact'
           option scalar $context: strand/frame context string
 Note    : getter only 

=cut

sub conserved{
    my $self = shift;
    my ($type, $action, $context) = @_;
    $self->_check_type_arg(\$type);
    $self->_check_action_arg(\$action);
    $self->_check_context_arg($type, \$context);
    if (!defined $self->{"conserved_${type}_${action}_${context}"}) {
	$self->_calc_stats($type, $action, $context);
    }
    return $self->{"conserved_${type}_${action}_${context}"};
}

=head2 length

 Title   : length
 Usage   : $tiling->length($type, $action, $context)
 Function: Retrieve the total length of aligned residues for 
           the seq $type
 Example : 
 Returns : value of length (a scalar)
 Args    : scalar $type: one of 'hit', 'subject', 'query'
           default is 'query'
           option scalar $action: one of 'exact', 'est', 'fast', 'max'
           default is 'exact'
           option scalar $context: strand/frame context string
 Note    : getter only 

=cut

sub length{
    my $self = shift;
    my ($type,$action,$context) = @_;
    $self->_check_type_arg(\$type);
    $self->_check_action_arg(\$action);
    $self->_check_context_arg($type, \$context);
    if (!defined $self->{"length_${type}_${action}_${context}"}) {
	$self->_calc_stats($type, $action, $context);
    }
    return $self->{"length_${type}_${action}_${context}"};
}

=head2 frac

 Title   : frac
 Usage   : $tiling->frac($type, $denom, $action, $context, $method)
 Function: Return the fraction of sequence length consisting
           of desired kinds of pairs (given by $method), 
           with respect to $denom
 Returns : scalar float
 Args    : -type => one of 'hit', 'subject', 'query'
           -denom => one of 'total', 'aligned'
           -action => one of 'exact', 'est', 'fast', 'max'
           -context => strand/frame context string
           -method => one of 'identical', 'conserved'
 Note    : $denom == 'aligned', return desired_stat/num_aligned
           $denom == 'total', return desired_stat/_reported_length
             (i.e., length of the original input sequences)
 Note    : In keeping with the spirit of Bio::Search::HSP::HSPI, 
           reported lengths of translated dna are reduced by 
           a factor of 3, to provide fractions relative to 
           amino acid coordinates. 

=cut

sub frac {
    my $self = shift;
    my @args = @_;
    my ($type, $denom, $action, $context, $method) = $self->_rearrange([qw(TYPE DENOM ACTION CONTEXT METHOD)],@args);
    $self->_check_type_arg(\$type);
    $self->_check_action_arg(\$action);
    $self->_check_context_arg($type, \$context);
    unless ($method and grep(/^$method$/, qw( identical conserved ))) {
	$self->throw("-method must specified; one of ('identical', 'conserved')");
    }
    $denom ||= 'total';
    unless (grep /^$denom/, qw( total aligned )) {
	$self->throw("Denominator selection must be one of ('total', 'aligned'), not '$denom'");
    }
    my $key = "frac_${method}_${type}_${denom}_${action}_${context}";
    my $stat;
    for ($method) {
	$_ eq 'identical' && do {
	    $stat = $self->identities($type, $action, $context);
	    last;
	};
	$_ eq 'conserved' && do {
	    $stat = $self->conserved($type, $action, $context);
	    last;
	};
	do {
	    $self->throw("What are YOU doing here?");
	};
    }
    if (!defined $self->{$key}) {
	for ($denom) {
	    /total/ && do {
		$self->{$key} =
		    $stat/$self->_reported_length($type); # need fudge fac??
		last;
	    };
	    /aligned/ && do {
		$self->{$key} =
		    $stat/$self->length($type,$action,$context);
		last;
	    };
	    do {
		$self->throw("What are YOU doing here?");
	    };
	}
    }
    return $self->{$key};
}

=head2 frac_identical

 Title   : frac_identical
 Usage   : $tiling->frac_identical($type, $denom, $action, $context)
 Function: Return the fraction of sequence length consisting
           of identical pairs, with respect to $denom
 Returns : scalar float
 Args    : -type => one of 'hit', 'subject', 'query'
           -denom => one of 'total', 'aligned'
           -action => one of 'exact', 'est', 'fast', 'max'
           -context => strand/frame context string
 Note    : $denom == 'aligned', return conserved/num_aligned
           $denom == 'total', return conserved/_reported_length
             (i.e., length of the original input sequences)
 Note    : In keeping with the spirit of Bio::Search::HSP::HSPI, 
           reported lengths of translated dna are reduced by 
           a factor of 3, to provide fractions relative to 
           amino acid coordinates. 
 Note    : This an alias that calls frac()

=cut

sub frac_identical{
    my $self = shift;
    my @args = @_;
    my ($type, $denom, $action,$context) = $self->_rearrange( [qw[ TYPE DENOM ACTION CONTEXT]],@args );
    $self->frac( -type=>$type, -denom=>$denom, -action=>$action, -method=>'identical', -context=>$context);
}

=head2 frac_conserved

 Title   : frac_conserved
 Usage   : $tiling->frac_conserved($type, $denom, $action, $context)
 Function: Return the fraction of sequence length consisting
           of conserved pairs, with respect to $denom
 Returns : scalar float
 Args    : -type => one of 'hit', 'subject', 'query'
           -denom => one of 'total', 'aligned'
           -action => one of 'exact', 'est', 'fast', 'max'
           -context => strand/frame context string
 Note    : $denom == 'aligned', return conserved/num_aligned
           $denom == 'total', return conserved/_reported_length
             (i.e., length of the original input sequences)
 Note    : In keeping with the spirit of Bio::Search::HSP::HSPI, 
           reported lengths of translated dna are reduced by 
           a factor of 3, to provide fractions relative to 
           amino acid coordinates. 
 Note    : This an alias that calls frac()

=cut

sub frac_conserved{
    my $self = shift;
    my @args = @_;
    my ($type, $denom, $action, $context) = $self->_rearrange( [qw[ TYPE DENOM ACTION CONTEXT]],@args );
    $self->frac( -type=>$type, -denom=>$denom, -action=>$action, -context=>$context, -method=>'conserved');
}

=head2 frac_aligned

 Title   : frac_aligned
 Aliases : frac_aligned_query - frac_aligned(-type=>'query',...)
           frac_aligned_hit   - frac_aligned(-type=>'hit',...)
 Usage   : $tiling->frac_aligned(-type=>$type,
                                 -action=>$action,
                                 -context=>$context)
 Function: Return the fraction of input sequence length
           that was aligned by the algorithm
 Returns : scalar float
 Args    : -type => one of 'hit', 'subject', 'query'
           -action => one of 'exact', 'est', 'fast', 'max'
           -context => strand/frame context string

=cut

sub frac_aligned{
    my ($self, @args) = @_;
    my ($type, $action, $context) = $self->_rearrange([qw(TYPE ACTION CONTEXT)],@args);
    $self->_check_type_arg(\$type);
    $self->_check_action_arg(\$action);
    $self->_check_context_arg($type, \$context);
    if (!$self->{"frac_aligned_${type}_${action}_${context}"}) {
	$self->{"frac_aligned_${type}_${action}_${context}"} = $self->num_aligned($type,$action,$context)/$self->_reported_length($type);
    }
    return $self->{"frac_aligned_${type}_${action}_${context}"};
}

sub frac_aligned_query { shift->frac_aligned(-type=>'query', @_) }
sub frac_aligned_hit { shift->frac_aligned(-type=>'hit', @_) }
    

=head2 num_aligned

 Title   : num_aligned
 Usage   : $tiling->num_aligned(-type=>$type)
 Function: Return the number of residues of sequence $type
           that were aligned by the algorithm
 Returns : scalar int
 Args    : -type => one of 'hit', 'subject', 'query'
           -action => one of 'exact', 'est', 'fast', 'max'
           -context => strand/frame context string
 Note    : Since this is calculated from reported coordinates,
           not symbol string counts, it is already in terms of
           "logical length"
 Note    : Aliases length()

=cut

sub num_aligned { shift->length( @_ ) };

=head2 num_unaligned

 Title   : num_unaligned
 Usage   : $tiling->num_unaligned(-type=>$type)
 Function: Return the number of residues of sequence $type
           that were left unaligned by the algorithm
 Returns : scalar int
 Args    : -type => one of 'hit', 'subject', 'query'
           -action => one of 'exact', 'est', 'fast', 'max'
           -context => strand/frame context string
 Note    : Since this is calculated from reported coordinates,
           not symbol string counts, it is already in terms of
           "logical length"

=cut

sub num_unaligned {
    my $self = shift;
    my ($type,$action,$context) = @_;
    my $ret;
    $self->_check_type_arg(\$type);
    $self->_check_action_arg(\$action);
    $self->_check_context_arg($type, \$context);
    if (!defined $self->{"num_unaligned_${type}_${action}_${context}"}) {
	$self->{"num_unaligned_${type}_${action}_${context}"} = $self->_reported_length($type)-$self->num_aligned($type,$action,$context);
    }
    return $self->{"num_unaligned_${type}_${action}_${context}"};
}
	

=head2 range

 Title   : range
 Usage   : $tiling->range(-type=>$type)
 Function: Returns the extent of the longest tiling
           as ($min_coord, $max_coord)
 Returns : array of two scalar integers
 Args    : -type => one of 'hit', 'subject', 'query'
           -context => strand/frame context string

=cut

sub range {
    my ($self, $type, $context) = @_;
    $self->_check_type_arg(\$type);
    $self->_check_context_arg($type, \$context);
    my @a = $self->_contig_intersection($type,$context);
    return ($a[0][0], $a[-1][1]);
}



=head1 ACCESSORS

=head2 coverage_map

 Title   : coverage_map
 Usage   : $map = $tiling->coverage_map($type)
 Function: Property to contain the coverage map calculated
           by _calc_coverage_map() - see that for 
           details
 Example : 
 Returns : value of coverage_map_$type as an array
 Args    : scalar $type: one of 'hit', 'subject', 'query'
           default is 'query'
 Note    : getter 

=cut

sub coverage_map{
    my $self = shift;
    my ($type, $context) = @_;
    $self->_check_type_arg(\$type);
    $self->_check_context_arg($type, \$context);

    if (!defined $self->{"coverage_map_${type}_${context}"}) {
	# following calculates coverage maps in all strands/frames
	# if necessary
	$self->_calc_coverage_map($type, $context);
    }
    # if undef is returned, then there were no hsps for given strand/frame
    if (!defined $self->{"coverage_map_${type}_${context}"}) {
	$self->warn("No HSPS present for type '$type' in context '$context' for this hit");
	return undef;
    }
    return @{$self->{"coverage_map_${type}_${context}"}};
}

=head2 coverage_map_as_text

 Title   : coverage_map_as_text
 Usage   : $tiling->coverage_map_as_text($type, $legend_flag)
 Function: Format a text-graphic representation of the
           coverage map
 Returns : an array of scalar strings, suitable for printing
 Args    : $type: one of 'query', 'hit', 'subject'
           $context: strand/frame context string
           $legend_flag: boolean; add a legend indicating
            the actual interval coordinates for each component
            interval and hsp (in the $type sequence context)
 Example : print $tiling->coverage_map_as_text('query',1);

=cut

sub coverage_map_as_text{
    my $self = shift;
    my ($type, $context, $legend_q) = @_;
    $self->_check_type_arg(\$type);
    $self->_check_context_arg($type, \$context);

    my @map = $self->coverage_map($type, $context);
    my @ret;
    my @hsps = $self->hit->hsps;
    my %hsps_i;
    require Tie::RefHash;
    tie %hsps_i, 'Tie::RefHash';
    @hsps_i{@hsps} = (0..$#hsps);
    my @mx;
    foreach (0..$#map) {
	my @hspx = ('') x @hsps;
	my @these_hsps = @{$map[$_]->[1]};
	@hspx[@hsps_i{@these_hsps}] = ('*') x @these_hsps;
	$mx[$_] = \@hspx;
    }
    untie %hsps_i;

    push @ret, "\tIntvl\n";
    push @ret, "HSPS\t", join ("\t", (0..$#map)), "\n";
    foreach my $h (0..$#hsps) {
	push @ret, join("\t", $h, map { $mx[$_][$h] } (0..$#map)  ),"\n";
    }
    if ($legend_q) {
	push @ret, "Interval legend\n";
	foreach (0..$#map) {
	    push @ret, sprintf("%d\t[%d, %d]\n", $_, @{$map[$_][0]});
	}
	push @ret, "HSP legend\n";
	my @ints = get_intervals_from_hsps($type,@hsps);
	foreach (0..$#hsps) {
	    push @ret, sprintf("%d\t[%d, %d]\n", $_, @{$ints[$_]});
	}
    }
    return @ret;
}

=head2 hit

 Title   : hit
 Usage   : $tiling->hit
 Function: 
 Example : 
 Returns : The HitI object associated with the invocant
 Args    : none
 Note    : getter only 

=cut

sub hit{
    my $self = shift;
    $self->warn("Getter only") if @_;
    return $self->{'hit'};
}

=head2 hsps

 Title   : hsps
 Usage   : $tiling->hsps()
 Function: Container for the HSP objects associated with invocant
 Example : 
 Returns : an array of hsps associated with the hit
 Args    : on set, new value (an arrayref or undef, optional)

=cut

sub hsps{
    my $self = shift;
    return $self->{'hsps'} = shift if @_;
    return @{$self->{'hsps'}};
}

=head2 contexts

 Title   : contexts
 Usage   : @contexts = $tiling->context($type) or
           @indices = $tiling->context($type, $context)
 Function: Retrieve the set of available contexts in the hit,
           or the indices of hsps having the given context
           (integer indices for the array returned by $self->hsps)
 Returns : array of scalar context strings or 
           array of scalar positive integers
           undef if no hsps in given context
 Args    : $type: one of 'query', 'hit', 'subject'
           optional $context: context string

=cut

sub contexts{
    my $self = shift;
    my ($type, $context) = @_;
    $self->_check_type_arg(\$type);
    return keys %{$self->{"_contexts_$type"}} unless defined $context;
    return undef unless $self->{"_contexts_$type"}{$context};
    return @{$self->{"_contexts_$type"}{$context}};
}

=head2 mapping

 Title   : mapping
 Usage   : $tiling->mapping($type)
 Function: Retrieve the mapping coefficient for the sequence type
           based on the underlying algorithm
 Returns : scalar integer (mapping coefficient)
 Args    : $type: one of 'query', 'hit', 'subject'
 Note    : getter only (set in constructor)

=cut

sub mapping{
    my $self = shift;
    my $type = shift;
    $self->_check_type_arg(\$type);
    return $self->{"_mapping_${type}"};
}

=head2 default_context

 Title   : default_context
 Usage   : $tiling->default_context($type)
 Function: Retrieve the default strand/frame context string
           for the sequence type based on the underlying algorithm
 Returns : scalar string (context string)
 Args    : $type: one of 'query', 'hit', 'subject'
 Note    : getter only (set in constructor)

=cut

sub default_context{
    my $self = shift;
    my $type = shift;
    $self->_check_type_arg(\$type);
    return $self->{"_def_context_${type}"};
}

=head2 algorithm

 Title   : algorithm
 Usage   : $tiling->algorithm
 Function: Retrieve the algorithm name associated with the 
           invocant's hit object
 Returns : scalar string 
 Args    : none
 Note    : getter only (set in constructor)

=cut

sub algorithm{
    my $self = shift;
    $self->warn("Getter only") if @_;
    return $self->{"_algorithm"};
}

=head1 "PRIVATE" METHODS

=head2 Calculators

See L<Bio::Search::Tiling::MapTileUtils> for lower level
calculation methods.

=head2 _calc_coverage_map

 Title   : _calc_coverage_map
 Usage   : $tiling->_calc_coverage_map($type)
 Function: Calculates the coverage map for the object's associated
           hit from the perspective of the desired $type (see Args:) 
           and sets the coverage_map() property
 Returns : True on success
 Args    : optional scalar $type: one of 'hit'|'subject'|'query'
           default is 'query'
 Note    : The "coverage map" is an array with the following format:
           ( [ $component_interval => [ @containing_hsps ] ], ... ),
           where $component_interval is a closed interval (see 
           DESCRIPTION) of the form [$a0, $a1] with $a0 <= $a1, and
           @containing_hsps is an array of all HspI objects in the hit 
           which completely contain the $component_interval.
           The set of $component_interval's is a disjoint decomposition
           of the minimum set of minimal intervals that completely
           cover the hit's HSPs (from the perspective of the $type)
 Note    : This calculates the map for all strand/frame contexts available
           in the hit

=cut

sub _calc_coverage_map {
    my $self = shift;
    my ($type) = @_;
    $self->_check_type_arg(\$type);

    # obtain the [start, end] intervals for all hsps in the hit (relative
    # to the type)
    unless ($self->{'hsps'}) {
	$self->warn("No HSPs for this hit");
	return;
    }

    my (@map, @hsps, %filters, @intervals);
    

    # conversion here?
    my $c = $self->mapping($type);
    
    # create the possible maps 
    for my $context ($self->contexts($type)) {
	@map = ();
	@hsps = ($self->hsps)[$self->contexts($type, $context)];
	@intervals = get_intervals_from_hsps( $type, @hsps );
	# the "frame"
	my $f = ($intervals[0]->[0] - 1) % $c;

	# convert interval endpoints...
	for (@intervals) {
	    $$_[0] = ($$_[0] - $f + $c - 1)/$c;
	    $$_[1]  = ($$_[1] - $f)/$c;
	}
	
	# determine the minimal set of disjoint intervals that cover the
	# set of hsp intervals
	my @dj_set = interval_tiling(\@intervals);

	# decompose each disjoint interval into another set of disjoint 
	# intervals, each of which is completely contained within the
	# original hsp intervals with which it overlaps
	my $i=0;
	my @decomp;
	for my $dj_elt (@dj_set) {
	    my ($covering, $indices) = @$dj_elt;
	    my @covering_hsps = @hsps[@$indices];
	    my @coverers = @intervals[@$indices];
	    @decomp = decompose_interval( \@coverers );
	    for (@decomp) {
		my ($component, $container_indices) = @{$_};
		push @map, [ $component, 
			     [@covering_hsps[@$container_indices]] ];
	    }
	    1;
	}
    
	# unconvert the components:
#####
	foreach (@map) {
	    $$_[0][0] = $c*$$_[0][0] - $c + 1 + $f;
	    $$_[0][1] = $c*$$_[0][1] + $f;
	}
	foreach (@dj_set) {
	    $$_[0][0] = $c*$$_[0][0] - $c + 1 + $f;
	    $$_[0][1] = $c*$$_[0][1] + $f;
	}	    

	# sort the map on the interval left-ends
	@map = sort { $a->[0][0]<=>$b->[0][0] } @map;
	$self->{"coverage_map_${type}_${context}"} = [@map];
	# set the _contig_intersection attribute here (side effect)
	$self->{"_contig_intersection_${type}_${context}"} = [map { $$_[0] } @map];
    }

    return 1; # success
}

=head2 _calc_stats

 Title   : _calc_stats
 Usage   : $tiling->_calc_stats($type, $action, $context)
 Function: Calculates [estimated] tiling statistics (identities, conserved sites
           length) and sets the public accessors
 Returns : True on success
 Args    : scalar $type: one of 'hit', 'subject', 'query'
           default is 'query'
           optional scalar $action: requests calculation method
            currently one of 'exact', 'est', 'fast', 'max'
           option scalar $context: strand/frame context string
 Note    : Action: The statistics are calculated by summing quantities
           over the disjoint component intervals, taking into account
           coverage of those intervals by multiple HSPs. The action
           tells the algorithm how to obtain those quantities--
           'exact' will use Bio::Search::HSP::HSPI::matches
            to count the appropriate segment of the homology string;
           'est' will estimate the statistics by multiplying the 
            fraction of the HSP overlapped by the component interval
            (see MapTileUtils) by the BLAST-reported identities/postives
            (this may be convenient for BLAST summary report formats)
           * Both exact and est take the average over the number of HSPs
             that overlap the component interval.
           'max' uses the exact method to calculate the statistics, 
            and returns only the maximum identites/positives over 
            overlapping HSP for the component interval. No averaging
            is involved here.
           'fast' doesn't involve tiling at all (hence the name),
            but it seems like a very good estimate, and uses only
            reported values, and so does not require sequence data. It
            calculates an average of reported identities, conserved
            sites, and lengths, over unmodified hsps in the hit,
            weighted by the length of the hsps.  

=cut

sub _calc_stats {
    my $self = shift;
    my ($type, $action, $context) = @_;
    # need to check args here, in case method is called internally.
    $self->_check_type_arg(\$type);
    $self->_check_action_arg(\$action);
    $self->_check_context_arg($type, \$context);
    my ($ident, $cons, $length) = (0,0,0);

    # fast : avoid coverage map altogether, get a pretty damn
    # good estimate with a weighted average of reported hsp
    # statistics

    ($action eq 'fast') && do {
	my @hsps = $self->hit->hsps;
	@hsps = @hsps[$self->contexts($type, $context)];
	# weights for averages
	my @wt = map {$_->length($type)} @hsps;
	my $sum = eval( join('+',@wt) );
	$_ /= $sum for (@wt);
	for (@hsps) { 
	    my $wt = shift @wt;
	    $ident  += $wt*$_->matches_MT($type,'identities');
	    $cons   += $wt*$_->matches_MT($type,'conserved');
	    $length += $wt*$_->length($type);
	}
    };

    # or, do tiling

    # calculate identities/conserved sites in tiling
    # estimate based on the fraction of the component interval covered
    # and ident/cons reported by the HSPs
    ($action ne 'fast') && do {
	foreach ($self->coverage_map($type, $context)) {
	    my ($intvl, $hsps) = @{$_};
	    my $len = ($$intvl[1]-$$intvl[0]+1);
	    my $ncover = ($action eq 'max') ? 1 : scalar @$hsps;
	    my ($acc_i, $acc_c) = (0,0);
	    foreach my $hsp (@$hsps) {
		for ($action) {
		    ($_ eq 'est') && do {
			my ($inc_i, $inc_c) = $hsp->matches_MT(
			    -type   => $type,
			    -action => 'searchutils',
			    );
			my $frac = $len/$hsp->length($type);
			$acc_i += $inc_i * $frac;
			$acc_c += $inc_c * $frac;
			last;
		    };
		    ($_ eq 'max') && do {
			my ($inc_i, $inc_c) = $hsp->matches_MT(
			    -type   => $type,
			    -action => 'searchutils',
			    -start => $$intvl[0], 
			    -end   => $$intvl[1]
			    );
			$acc_i = ($acc_i > $inc_i) ? $acc_i : $inc_i;
			$acc_c = ($acc_c > $inc_c) ? $acc_c : $inc_c;
			last;
		    };
		    (!$_ || ($_ eq 'exact')) && do {
			my ($inc_i, $inc_c) = $hsp->matches_MT(
			    -type   => $type, 
			    -action => 'searchutils',
			    -start  => $$intvl[0], 
			    -end    => $$intvl[1]
			    );
			$acc_i += $inc_i;
			$acc_c += $inc_c;
			last;
		    };
		}
	    }
	    $ident += ($acc_i/$ncover);
	    $cons  += ($acc_c/$ncover);
	    $length += $len;
	}
    };
    
    $self->{"identities_${type}_${action}_${context}"} = $ident;
    $self->{"conserved_${type}_${action}_${context}"} = $cons;
    $self->{"length_${type}_${action}_${context}"} = $length;
    
    return 1;
}

=head2 Tiling Helper Methods

=cut

# coverage_map is of the form
# ( [ $interval, \@containing_hsps ], ... )

# so, for each interval, pick one of the containing hsps,
# and return the union of all the picks.

# use the combinatorial generating iterator, with 
# the urns containing the @containing_hsps for each
# interval

=head2 _make_tiling_iterator

 Title   : _make_tiling_iterator
 Usage   : $self->_make_tiling_iterator($type)
 Function: Create an iterator code ref that will step through all 
           minimal combinations of HSPs that produce complete coverage
           of the $type ('hit', 'subject', 'query') sequence, 
           and set the correct iterator property of the invocant
 Example :
 Returns : The iterator
 Args    : scalar $type, one of 'hit', 'subject', 'query';
           default is 'query'

=cut

sub _make_tiling_iterator {
    ### create the urns
    my $self = shift;
    my ($type, $context) = @_;
    $self->_check_type_arg(\$type);
    $self->_check_context_arg($type, \$context);

    # initialize the urns
    my @urns = map { [0,  $$_[1]] } $self->coverage_map($type, $context);

    my $FINISHED = 0;
    my $iter = sub {
	# rewind?
	if (my $rewind = shift) {
	    # reinitialize urn indices
	    $$_[0] = 0 for (@urns);
	    $FINISHED = 0;
	    return 1;
	}	    
	# check if done...
        return if $FINISHED;

        my $finished_incrementing = 0;
	# @ret is the collector of urn choices
	my @ret;

	for my $urn (@urns) {
	    my ($n, $hsps) = @$urn;
	    push @ret, $$hsps[$n];
	    unless ($finished_incrementing) {
		if ($n == $#$hsps) { $$urn[0] = 0; }
		else { ($$urn[0])++; $finished_incrementing = 1 }
	    }
	}

	# backstop...
        $FINISHED = 1 unless $finished_incrementing;
	# uniquify @ret
	# $hsp->rank is a unique identifier for an hsp in a hit.
	# preserve order in @ret
	
	my (%order, %uniq);
	@order{(0..$#ret)} = @ret;
	$uniq{$order{$_}->rank} = $_ for (0..$#ret);
	@ret = @order{ sort {$a<=>$b} values %uniq };

        return @ret;
    };
    return $iter;
}

=head2 _tiling_iterator

 Title   : _tiling_iterator
 Usage   : $tiling->_tiling_iterator($type,$context)
 Function: Retrieve the tiling iterator coderef for the requested 
           $type ('hit', 'subject', 'query')
 Example : 
 Returns : coderef to the desired iterator
 Args    : scalar $type, one of 'hit', 'subject', 'query'
           default is 'query'
           option scalar $context: strand/frame context string
 Note    : getter only

=cut

sub _tiling_iterator {
    my $self = shift;
    my ($type, $context) = @_;
    $self->_check_type_arg(\$type);
    $self->_check_context_arg($type, \$context);

    if (!defined $self->{"_tiling_iterator_${type}_${context}"}) {
	$self->{"_tiling_iterator_${type}_${context}"} =
	    $self->_make_tiling_iterator($type,$context);
    }
    return $self->{"_tiling_iterator_${type}_${context}"};
}
=head2 Construction Helper Methods

See also L<Bio::Search::Tiling::MapTileUtils>.

=cut

sub _check_type_arg {
    my $self = shift;
    my $typeref = shift;
    $$typeref ||= 'query';
    $self->throw("Unknown type '$$typeref'") unless grep(/^$$typeref$/, qw( hit query subject ));    
    $$typeref = 'hit' if $$typeref eq 'subject';
    return 1;
}

sub _check_action_arg {
    my $self = shift;
    my $actionref = shift;
    if (!$$actionref) {
	$$actionref = ($self->_has_sequence_data ? 'exact' : 'est');
    }
    else {
	$self->throw("Calc action '$$actionref' unrecognized") unless grep /^$$actionref$/, qw( est exact fast max );
	if ($$actionref ne 'est' and !$self->_has_sequence_data) {
	    $self->warn("Blast file did not possess sequence data; defaulting to 'est' action");
	    $$actionref = 'est';
	}
    }
    return 1;
}

sub _check_context_arg {
    my $self = shift;
    my ($type, $contextref) = @_;
    if (!$$contextref) {
	$self->throw("Type '$type' requires strand/frame context for algorithm ".$self->algorithm) unless ($self->mapping($type) == 1);
	# set default according to default_context attrib
	$$contextref = $self->default_context($type);
    }
    else {
	($$contextref =~ /^[mp]$/) && do { $$contextref .= '_' };
	$self->throw("Context '$$contextref' unrecognized") unless
	    $$contextref =~ /all|[mp][0-2_]/;
    }
	
}

=head2 _make_context_key

 Title   : _make_context_key
 Alias   : _context
 Usage   : $tiling->_make_context_key(-strand => $strand, -frame => $frame)
 Function: create a string indicating strand/frame context; serves as 
           component of memoizing hash keys
 Returns : scalar string
 Args    : -type => one of ('query', 'hit', 'subject')
           -strand => one of (1,0,-1)
           -frame  => one of (-2, 1, 0, 1, -2)
           called w/o args: returns 'all'

=cut

sub _make_context_key {
    my $self = shift;
    my @args = @_;
    my ($type, $strand, $frame) = $self->_rearrange([qw(TYPE STRAND FRAME)], @args);
    _check_type_arg(\$type);
    return 'all' unless (defined $strand or defined $frame);
    if ( defined $strand && $self->_has_strand($type) ) {
	if (defined $frame && $self->_has_frame($type)) {
	    return ($strand >= 0 ? 'p' : 'm').abs($frame);
	}
	else {
	    return ($strand >= 0 ? 'p_' : 'm_');
	}
    }
    else {
	if (defined $frame && $self->_has_frame($type)) {
	    $self->warn("Frame defined without strand; punting with plus strand");
	    return 'p'.abs($frame);
	}
	else {
	    return 'all';
	}
    }
}

=head2 _context

 Title   : _context
 Alias   : _make_context_key
 Usage   : $tiling->_make_context_key(-strand => $strand, -frame => $frame)
 Function: create a string indicating strand/frame context; serves as 
           component of memoizing hash keys
 Returns : scalar string
 Args    : -type => one of ('query', 'hit', 'subject')
           -strand => one of (1,0,-1)
           -frame  => one of (-2, 1, 0, 1, -2)
           called w/o args: returns 'all'

=cut

sub _context { shift->_make_context_key(@_) }

=head2 Predicates

Most based on a reading of the algorithm name with a configuration lookup.

=over

=item _has_sequence_data()

=cut 

sub _has_sequence_data {
    my $self = shift;
    $self->throw("Hit attribute  not yet set") unless defined $self->hit;
    return (($self->hit->hsps)[0]->seq_str('match') ? 1 : 0);
}

=item _has_logical_length()

=cut

sub _has_logical_length {
    my $self = shift;
    my $type = shift;
    $self->_check_type_arg(\$type);
    # based on mapping coeff
    $self->throw("Mapping coefficients not yet set") unless defined $self->mapping($type);
    return ($self->mapping($type) > 1);
}

=item _has_strand()

=cut

sub _has_strand {
    my $self = shift;
    my $type = shift;
    $self->_check_type_arg(\$type);
    return $self->{"_has_strand_${type}"};
}

=item _has_frame()

=cut

sub _has_frame {
    my $self = shift;
    my $type = shift;
    $self->_check_type_arg(\$type);
    return $self->{"_has_frame_${type}"};
}

=back

=head1 Private Accessors

=head2 _contig_intersection

 Title   : _contig_intersection
 Usage   : $tiling->_contig_intersection($type)
 Function: Return the minimal set of $type coordinate intervals
           covered by the invocant's HSPs
 Returns : array of intervals (2-member arrayrefs; see MapTileUtils)
 Args    : scalar $type: one of 'query', 'hit', 'subject'

=cut

sub _contig_intersection {
    my $self = shift;
    my ($type, $context) = @_;
    $self->_check_type_arg(\$type);
    $self->_check_context_arg($type, \$context);
    if (!defined $self->{"_contig_intersection_${type}_${context}"}) {
	$self->_calc_coverage_map($type);
    }
    return @{$self->{"_contig_intersection_${type}_${context}"}};
}

=head2 _reported_length

 Title   : _reported_length
 Usage   : $tiling->_reported_length($type)
 Function: Get the total length of the seq $type
           for the invocant's hit object, as reported
           by (not calculated from) the input data file
 Returns : scalar int
 Args    : scalar $type: one of 'query', 'hit', 'subject'
 Note    : This is kludgy; the hit object does not currently
           maintain accessors for these values, but the 
           hsps possess these attributes. This is a wrapper
           that allows a consistent access method in the 
           MapTiling code.
 Note    : Since this number is based on a reported length,
           it is already a "logical length". 

=cut

sub _reported_length {
    my $self = shift;
    my $type = shift;
    $self->_check_type_arg(\$type);
    my $key = uc( $type."_LENGTH" );
    return ($self->hsps)[0]->{$key};
}

1;
    
