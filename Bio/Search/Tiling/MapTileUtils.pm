#$Id$
package Bio::Search::Tiling::MapTileUtils;
use strict;
use warnings;
use Exporter;

BEGIN {
    our @ISA = qw( Exporter );
    our @EXPORT = qw( get_intervals_from_hsps 
                      interval_tiling 
                      decompose_interval
                      containing_hsps
                      covering_groups
                      _allowable_filters 
                      _set_attributes
                      _mapping_coeff);
}

# tiling trials
# assumed: intervals are [$a0, $a1], with $a0 <= $a1
=head1 NAME

Bio::Search::Tiling::MapTileUtils - utilities for manipulating closed intervals for an HSP tiling algorithm

=head1 SYNOPSIS

Not used directly.

=head1 DESCRIPTION

Not used directly.

=head1 NOTE

An "interval" in this module is defined as an arrayref C<[$a0, $a1]>, where
C<$a0, $a1> are scalar numbers satisfying C<$a0 E<lt>= $a1>.

=head1 AUTHOR

Mark A. Jensen - maj -at- fortinbras -dot- us

=head1 APPENDIX

=head2 interval_tiling    

 Title   : interval_tiling()
 Usage   : @tiling = interval_tiling( \@array_of_intervals )
 Function: Find minimal set of intervals covering the input set
 Returns : array of arrayrefs of the form
  ( [$interval => [ @indices_of_collapsed_input_intervals ]], ...)
 Args    : arrayref of intervals

=cut

sub interval_tiling {
    return unless $_[0]; # no input
    my $n = scalar @{$_[0]};
    my %input;
    @input{(0..$n-1)} = @{$_[0]};
    my @active = (0..$n-1);
    my @hold;
    my @tiled_ints;
    my @ret;
    while (@active) {
	my $tgt = $input{my $tgt_i = shift @active};
	push @tiled_ints, $tgt_i;
	my $tgt_not_disjoint = 1;
	while ($tgt_not_disjoint) {
	    $tgt_not_disjoint = 0;
	    while (my $try_i = shift @active) {
		my $try = $input{$try_i};
		if ( !are_disjoint($tgt, $try) ) {
		    $tgt = min_covering_interval($tgt,$try);
		    push @tiled_ints, $try_i;
		    $tgt_not_disjoint = 1;
		}
		else {
		    push @hold, $try_i;
		}
	    }
	    if (!$tgt_not_disjoint) {
		push @ret, [ $tgt => [@tiled_ints] ];
		@tiled_ints = ();
	    }
	    @active = @hold;
	    @hold = ();
	}
    }
    return @ret;
}

=head2 decompose_interval

 Title   : decompose_interval
 Usage   : @decomposition = decompose_interval( \@overlappers )
 Function: Calculate the disjoint decomposition of a set of
           overlapping intervals, each annotated with a list of
           covering intervals
 Returns : array of arrayrefs of the form
           ( [[@interval] => [@indices_of_coverers]], ... )
 Args    : arrayref of intervals (arrayrefs like [$a0, $a1], with
 Note    : Each returned interval is associated with a list of indices of the
           original intervals that cover that decomposition component
           (scalar size of this list could be called the 'coverage coefficient')
 Note    : Coverage: each component of the decomp is completely contained
           in the input intervals that overlap it, by construction.
 Caveat  : This routine expects the members of @overlappers to overlap,
           but doesn't check this.

=cut

### what if the input intervals don't overlap?? They MUST overlap; that's
### what interval_tiling() is for.

sub decompose_interval {
    return unless $_[0]; # no input
    my @ints = @{$_[0]};
    my (%flat,@flat);
    ### this is ok, but need to handle the case where a lh and rh endpoint
    ### coincide...
    # decomposition --
    # flatten:
    # every lh endpoint generates (lh-1, lh)
    # every rh endpoint generates (rh, rh+)
    foreach (@ints) {
	$flat{$$_[0]-1}++;
	$flat{$$_[0]}++;
	$flat{$$_[1]}++;
	$flat{$$_[1]+1}++;
    }
    # sort, create singletons if nec.
    my @a;
    @a = sort {$a<=>$b} keys %flat;
    # throw out first and last (meeting a boundary condition)
    shift @a; pop @a;
    # look for singletons
    @flat = (shift @a, shift @a);
    if ( $flat[1]-$flat[0] == 1 ) {
	@flat = ($flat[0],$flat[0], $flat[1]);
    }
    while (my $a = shift @a) {
	if ($a-$flat[-2]==2) {
	    push @flat, $flat[-1]; # create singleton interval
	}
	push @flat, $a;
    }
    if ($flat[-1]-$flat[-2]==1 and @flat % 2) {
	push @flat, $flat[-1];
    }
    # component intervals are consecutive pairs
    my @decomp;
    while (my $a = shift @flat) {
	push @decomp, [$a, shift @flat];
    }

    # for each component, return a list of the indices of the input intervals
    # that cover the component.
    my @coverage;
    foreach my $i (0..$#decomp) {
	foreach my $j (0..$#ints) {
	    unless (are_disjoint($decomp[$i], $ints[$j])) {
		if (defined $coverage[$i]) {
		    push @{$coverage[$i]}, $j;
		}
		else {
		    $coverage[$i] = [$j];
		}
	    }
	}
    }
    return map { [$decomp[$_] => $coverage[$_]] } (0..$#decomp);
}    

=head2 are_disjoint

 Title   : are_disjoint
 Usage   : are_disjoint( [$a0, $a1], [$b0, $b1] )
 Function: Determine if two intervals are disjoint
 Returns : True if the intervals are disjoint, false if they overlap
 Args    : array of two intervals

=cut

sub are_disjoint {
    my ($int1, $int2) = @_;
    return 1 if ( $$int1[1] < $$int2[0] ) || ( $$int2[1] < $$int1[0]);
    return 0;
}

=head2 min_covering_interval

 Title   : min_covering_interval 
 Usage   : $interval = min_covering_interval( [$a0,$a1],[$b0,$b1] )
 Function: Determine the minimal covering interval for two intervals
 Returns : an interval
 Args    : two intervals

=cut

sub min_covering_interval {
    my ($int1, $int2) = @_;
    my @a = sort {$a <=> $b} (@$int1, @$int2);
    return [$a[0],$a[-1]];
}

=head2 get_intervals_from_hsps

 Title   : get_intervals_from_hsps
 Usage   : @intervals = get_intervals_from_hsps($type, @hsp_objects)
 Function: Return array of intervals of the form [ $start, $end ],
           from an array of hsp objects
 Returns : an array of intervals
 Args    : scalar $type, array of HSPI objects; where $type is one of 'hit',
           'subject', 'query'

=cut

sub get_intervals_from_hsps {
    my $type = shift;
    my @hsps;
    if (ref($type) =~ /HSP/) {
	push @hsps, $type;
	$type = 'query';
    }
    elsif ( !grep /^$type$/,qw( hit subject query ) ) {
	die "Unknown HSP type '$type'";
    }
    $type = 'hit' if $type eq 'subject';
    push @hsps, @_;
    my @ret;
    foreach (@hsps) {
#	push @ret, [ $_->start($type), $_->end($type)];
	push @ret, [ $_->range($type) ];
    }
    return @ret;
}

# fast, clear, nasty, brutish and short.
# for _allowable_filters(), _set_mapping()
# covers BLAST, FAST families
# FASTA is ambiguous (nt or aa) based on alg name only

my $alg_lookup = {
    'N'  => { 'mapping' => [1,1],
              'def_context' => ['p_','p_'],
	      'has_strand' => [1, 1],
              'has_frame' => [0, 0]},
    'P'  => { 'mapping' => [1,1],
              'def_context' => ['all','all'],
	      'has_strand' => [0, 0],
	      'has_frame' => [0, 0]},
    'X'  => { 'mapping' => [3, 1],
              'def_context' => [undef,'all'],
	      'has_strand' => [1, 0],
	      'has_frame' => [1, 0]}, 
    'Y'  => { 'mapping' => [3, 1],
              'def_context' => [undef,'all'],
	      'has_strand' => [1, 0],
	      'has_frame' => [1, 0]}, 
    'TA' => { 'mapping' => [1, 3],
              'def_context' => ['all',undef],
	      'has_strand' => [0, 1],
	      'has_frame' => [0, 1]}, 
    'TN' => { 'mapping' => [1, 3],
              'def_context' => ['p_',undef],
	      'has_strand' => [1,1],
	      'has_frame' => [0, 1]}, 
    'TX' => { 'mapping' => [3, 3],
              'def_context' => [undef,undef],
	      'has_strand' => [1, 1],
	      'has_frame' => [1, 1]}, 
    'TY' => { 'mapping' => [3, 3],
              'def_context' => [undef,undef],
	      'has_strand' => [1, 1],
	      'has_frame' => [1, 1]}
};
   
=head2 _allowable_filters
    
 Title   : _allowable_filters
 Usage   : _allowable_filters($Bio_Search_Hit_HitI, $type)
 Function: Return the HSP filters (strand, frame) allowed, 
           based on the reported algorithm
 Returns : String encoding allowable filters: 
           s = strand, f = frame
           Empty string if no filters allowed
           undef if algorithm unrecognized
 Args    : A Bio::Search::Hit::HitI object,
           scalar $type, one of 'hit', 'subject', 'query';
           default is 'query'

=cut

sub _allowable_filters {
    my $hit = shift;
    my $type = shift;
    $type ||= 'q';
    unless (grep /^$type$/, qw( h q s ) ) {
	warn("Unknown type '$type'; returning ''");
	return '';
    }
    $type = 'h' if $type eq 's';
    my $alg = $hit->algorithm;

    # pretreat (i.e., kludge it)
    $alg =~ /^RPS/ and ($alg) = ($alg =~ /\(([^)]+)\)/);

    for ($hit->algorithm) {
	/MEGABLAST/i && do {
	    return qr/[s]/;
	};
	/(.?)BLAST(.?)/i && do {
	    return $$alg_lookup{$1.$2}{$type};
	};
	/(.?)FAST(.?)/ && do {
	    return $$alg_lookup{$1.$2}{$type};
	};
	do { # unrecognized
	    last;
	};
    }
    return;
}


=head2 _set_attributes

 Title   : _set_attributes
 Usage   : $tiling->_set_attributes()
 Function: Sets attributes for invocant
           that depend on algorithm name
 Returns : True on success
 Args    : none
 Note    : setting based on the configuration table
           %alg_lookup

=cut

sub _set_attributes {
    my $self = shift;
    my $alg = $self->hit->algorithm;

    # pretreat (i.e., kludge it)
    $alg =~ /^RPS/ and ($alg) = ($alg =~ /\(([^)]+)\)/);
    
    for ($alg) {
	/MEGABLAST/i && do {
	    ($self->{_mapping_query},$self->{_mapping_hit}) = (1,1);
	    ($self->{_def_context_query},$self->{_def_context_hit}) =
		('p_','p_');
	    ($self->{_has_frame_query},$self->{_has_frame_hit}) =
		(0, 0);
	    ($self->{_has_strand_query},$self->{_has_strand_hit}) =
		(1, 1);
	    last;
	};
	/(.?)BLAST(.?)/i && do {
	    ($self->{_mapping_query},$self->{_mapping_hit}) = 
		@{$$alg_lookup{$1.$2}{mapping}};
	    ($self->{_def_context_query},$self->{_def_context_hit}) =
		@{$$alg_lookup{$1.$2}{def_context}};
	    ($self->{_has_frame_query},$self->{_has_frame_hit}) =	    
		@{$$alg_lookup{$1.$2}{has_frame}};
	    ($self->{_has_strand_query},$self->{_has_strand_hit}) =	    
		@{$$alg_lookup{$1.$2}{has_strand}};
	    last;
	};
	/(.?)FAST(.?)/ && do {
	    ($self->{_mapping_query},$self->{_mapping_hit}) = 
		@{$$alg_lookup{$1.$2}{mapping}};
	    ($self->{_def_context_query},$self->{_def_context_hit}) =
		@{$$alg_lookup{$1.$2}{def_context}};
	    ($self->{_has_frame_query},$self->{_has_frame_hit}) =	    
		@{$$alg_lookup{$1.$2}{has_frame}};
	    ($self->{_has_strand_query},$self->{_has_strand_hit}) =	    
		@{$$alg_lookup{$1.$2}{has_strand}};
	    last;
	};
	do { # unrecognized
	    $self->warn("Unrecognized algorithm '$alg'; defaults may not work");
	    ($self->{_mapping_query},$self->{_mapping_hit}) = (1,1);
	    ($self->{_def_context_query},$self->{_def_context_hit}) =
		('all','all');
	    ($self->{_has_frame_query},$self->{_has_frame_hit}) =	    
		(0,0);
	    ($self->{_has_strand_query},$self->{_has_strand_hit}) =	    
		(0,0);
	    return 0;
	    last;
	};
    }
    return 1;
}
           
sub _mapping_coeff {
    my $obj = shift;
    my $type = shift;
    my %type_i = ( 'query' => 0, 'hit' => 1 );
    unless ( ref($obj) && $obj->can('algorithm') ) {
	$obj->warn("Object type unrecognized");
	return undef;
    }
    $type ||= 'query';
    unless ( grep(/^$type$/, qw( query hit subject ) ) ) {
	$obj->warn("Sequence type unrecognized");
	return undef;
    }
    $type = 'hit' if $type eq 'subject';
    my $alg = $obj->algorithm;

    # pretreat (i.e., kludge it)
    $alg =~ /^RPS/ and ($alg) = ($alg =~ /\(([^)]+)\)/);
    
    for ($alg) {
	/MEGABLAST/i && do {
	    return 1;
	};
	/(.?)BLAST(.?)/i && do {
	    return $$alg_lookup{$1.$2}{'mapping'}[$type_i{$type}];
	};
	/(.?)FAST(.?)/ && do {
	    return $$alg_lookup{$1.$2}{'mapping'}[$type_i{$type}];
	};
	do { # unrecognized
	    last;
	};
    }
    return;
}

# a graphical depiction of a set of intervals
sub _ints_as_text {
    my $ints = shift;
    my @ints = @$ints;
    my %pos;
    for (@ints) {
	$pos{$$_[0]}++;
	$pos{$$_[1]}++;
    }
    
    my @pos = sort {$a<=>$b} keys %pos;
    @pos = map {sprintf("%03d",$_)} @pos;
#header
    my $max=0;
    $max = (length > $max) ? length : $max for (@pos);
    for my $j (0..$max-1) {
	my $i = $max-1-$j; 
	my @line = map { substr($_, $j, 1) || '0' } @pos;
	print join('', @line), "\n";
    }
    print '-' x @pos, "\n";
    undef %pos;
    @pos{map {sprintf("%d",$_)} @pos} = (0..@pos);
    foreach (@ints) {
	print ' ' x $pos{$$_[0]}, '[', ' ' x ($pos{$$_[1]}-$pos{$$_[0]}-1), ']', ' ' x (@pos-$pos{$$_[1]}), "\n";
    }
}



=head2 containing_hsps()

 Title   : containing_hsps
 Usage   : @hsps = containing_hsps($interval, @hsps_to_search)
 Function: Return a list of hsps whose coordinates completely contain the
           given $interval
 Returns : Array of HSP objects
 Args    : $interval : [$int1, $int2],
           array of HSP objects

=cut
# could be more efficient if hsps are assumed ordered...
sub containing_hsps {
    my $intvl = shift;
    my @hsps = @_;
    my @ret;
    my ($beg, $end) = @$intvl;
    foreach my $hsp (@hsps) {
	my ($start, $stop) = ($hsp->start, $hsp->end);
	push @ret, $hsp if ( $start <= $beg and $end <= $stop );
    }
    return @ret;
}



=head2 covering_groups()

 Title   : covering_groups
 Usage   : 
 Function: divide a list of **ordered,disjoint** intervals (as from a 
           coverage map) into a set of disjoint covering groups
 Returns : array of arrayrefs, each arrayref a covering set of 
           intervals
 Args    : array of intervals

=cut

sub covering_groups {
    my @intervals = @_;
    return unless @intervals;
    my (@groups, $grp);
    push @{$groups[0]}, shift @intervals;
    $grp = $groups[0];
    for (my $intvl = shift @intervals; @intervals; $intvl = shift @intervals) {
	if ( $intvl->[0] - $grp->[-1][1] == 1 ) { # intervals are direcly adjacent
	    push @$grp, $intvl;
	}
	else {
	    $grp = [$intvl];
	    push @groups, $grp;
	}
    }
    return @groups;
}

1;
# need our own subsequencer for hsps. 

package Bio::Search::HSP::HSPI;

use strict;
use warnings;

=head2 matches_MT

 Title   : matches_MT
 Usage   : $hsp->matches($type, $action, $start, $end)
 Purpose   : Get the total number of identical or conserved matches 
             in the query or sbjct sequence for the given HSP. Optionally can
             report data within a defined interval along the seq.
 Returns   : scalar int 
 Args      : 
 Comments  : Relies on seq_str('match') to get the string of alignment symbols
             between the query and sbjct lines which are used for determining
             the number of identical and conservative matches.
 Note      : Modeled on Bio::Search::HSP::HSPI::matches

=cut

sub matches_MT {
    my( $self, @args ) = @_;
    my($type, $action, $beg, $end) = $self->_rearrange( [qw(TYPE ACTION START END)], @args);
    my @actions = qw( identities conserved searchutils );
    
    # prep $type
    $self->throw("Type not specified") if !defined $type;
    $self->throw("Type '$type' unrecognized") unless grep(/^$type$/,qw(query hit subject));
    $type = 'hit' if $type eq 'subject';

    # prep $action
    $self->throw("Action not specified") if !defined $action;
    $self->throw("Action '$action' unrecognized") unless grep(/^$action$/, @actions);

    my ($len_id, $len_cons);
    my $c = Bio::Search::Tiling::MapTileUtils::_mapping_coeff($self, $type);
    if ((defined $beg && !defined $end) || (!defined $beg && defined $end)) {
	$self->throw("Both start and end are required");
    }
    elsif ( (!defined($beg) && !defined($end)) || !$self->seq_str('match') ) {
        ## Get data for the whole alignment.
	# the reported values x mapping 
	$self->debug("Sequence data not present in report; returning data for entire HSP") unless $self->seq_str('match');
	($len_id, $len_cons) = map { $c*$_ } ($self->num_identical, $self->num_conserved);
	for ($action) {
	    $_ eq 'identities'  && do {
		return $len_id;
	    };
	    $_ eq 'conserved'   && do {
		return $len_cons;
	    };
	    $_ eq 'searchutils' && do {
		return ($len_id, $len_cons);
	    };
	    do {
		$self->throw("What are YOU doing here?");
	    };
	}
    }
    else {
        ## Get the substring representing the desired sub-section of aln.
        my($start,$stop) = $self->range($type);
	if ( $beg < $start or $stop < $end ) {
	    $self->throw("Start/stop out of range [$start, $stop]");
	}

	# handle gaps
	my $match_str = $self->seq_str('match');
	if ($self->gaps) {
	    # strip the homology string of gap positions relative
	    # to the target type
	    $match_str = $self->seq_str('match');
	    my $tgt   = $self->seq_str($type);
	    my $encode = $match_str ^ $tgt;
	    my $zap = '-'^' ';
	    $encode =~ s/$zap//g;
	    $tgt =~ s/-//g;
	    $match_str = $tgt ^ $encode;
	    # match string is now the correct length for substr'ing below,
	    # given that start and end are gapless coordinates in the 
	    # blast report
	}

        my $seq = "";
	$seq = substr( $match_str, 
		       int( ($beg-$start)/Bio::Search::Tiling::MapTileUtils::_mapping_coeff($self, $type) ),
		       int( 1+($end-$beg)/Bio::Search::Tiling::MapTileUtils::_mapping_coeff($self, $type) ) 
	    );

        if(!CORE::length $seq) {
            $self->throw("Undefined sub-sequence ($beg,$end). Valid range = $start - $stop");
        }

        $seq =~ s/ //g;  # remove space (no info).
        $len_cons = (CORE::length $seq)*(Bio::Search::Tiling::MapTileUtils::_mapping_coeff($self,$type));
        $seq =~ s/\+//g;  # remove '+' characters (conservative substitutions)
        $len_id = (CORE::length $seq)*(Bio::Search::Tiling::MapTileUtils::_mapping_coeff($self,$type));
	for ($action) {
	    $_ eq 'identities' && do {
		return $len_id;
	    };
	    $_ eq 'conserved' && do {
		return $len_cons;
	    };
	    $_ eq 'searchutils' && do {
		return ($len_id, $len_cons);
	    };
	    do {
		$self->throw("What are YOU doing here?");
	    };
	}
    }
}

1;

package Bio::LocatableSeq;
use strict;
use warnings;

# mixin the Bio::FeatureHolderI implementation of 
# Bio::Seq -- for get_tiled_aln

=head2 get_SeqFeatures

 Title   : get_SeqFeatures
 Usage   :
 Function: Get the feature objects held by this feature holder.

           Features which are not top-level are subfeatures of one or
           more of the returned feature objects, which means that you
           must traverse the subfeature arrays of each top-level
           feature object in order to traverse all features associated
           with this sequence.

           Top-level features can be obtained by tag, specified in 
           the argument.

           Use get_all_SeqFeatures() if you want the feature tree
           flattened into one single array.

 Example :
 Returns : an array of Bio::SeqFeatureI implementing objects
 Args    : [optional] scalar string (feature tag)


=cut

sub get_SeqFeatures{
   my $self = shift;
   my $tag = shift;

   if( !defined $self->{'_as_feat'} ) {
       $self->{'_as_feat'} = [];
   }
   if ($tag) {
       return map { $_->primary_tag eq $tag ? $_ : () } @{$self->{'_as_feat'}};
   }
   else {
       return @{$self->{'_as_feat'}};
   }
}

=head2 feature_count

 Title   : feature_count
 Usage   : $seq->feature_count()
 Function: Return the number of SeqFeatures attached to a sequence
 Returns : integer representing the number of SeqFeatures
 Args    : None


=cut

sub feature_count {
    my ($self) = @_;
    if (defined($self->{'_as_feat'})) {
	return ($#{$self->{'_as_feat'}} + 1);
    } else {
	return 0;
    }
}

=head2 add_SeqFeature

 Title   : add_SeqFeature
 Usage   : $seq->add_SeqFeature($feat);
           $seq->add_SeqFeature(@feat);
 Function: Adds the given feature object (or each of an array of feature
           objects to the feature array of this
           sequence. The object passed is required to implement the
           Bio::SeqFeatureI interface.
 Returns : 1 on success
 Args    : A Bio::SeqFeatureI implementing object, or an array of such objects.


=cut

sub add_SeqFeature {
   my ($self,@feat) = @_;
   $self->{'_as_feat'} = [] unless $self->{'_as_feat'};
   foreach my $feat ( @feat ) {
       if( !$feat->isa("Bio::SeqFeatureI") ) {
	   $self->throw("$feat is not a SeqFeatureI and that's what we expect...");
       }
       $feat->attach_seq($self);
       push(@{$self->{'_as_feat'}},$feat);
   }
   return 1;
}

=head2 remove_SeqFeatures

 Title   : remove_SeqFeatures
 Usage   : $seq->remove_SeqFeatures();
 Function: Flushes all attached SeqFeatureI objects.

           To remove individual feature objects, delete those from the returned
           array and re-add the rest.
 Example :
 Returns : The array of Bio::SeqFeatureI objects removed from this seq.
 Args    : None


=cut

sub remove_SeqFeatures {
    my $self = shift;

    return () unless $self->{'_as_feat'};
    my @feats = @{$self->{'_as_feat'}};
    $self->{'_as_feat'} = [];
    return @feats;
}

1;
