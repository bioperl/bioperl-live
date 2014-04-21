#
# BioPerl module for Bio::Map::Mappable
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Sendu Bala
# 
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::Mappable - An object representing a generic map element
that can have multiple locations in several maps.

=head1 SYNOPSIS

  # a map element in two different positions on the same map
  $map1 = Bio::Map::SimpleMap->new();
  $position1 = Bio::Map::Position->new(-map => $map1, -value => 100);
  $position2 = Bio::Map::Position->new(-map => $map1, -value => 200);
  $mappable = Bio::Map::Mappable->new(-positions => [$position1, $position2] );

  # add another position on a different map
  $map2 = Bio::Map::SimpleMap->new();
  $position3 = Bio::Map::Position->new(-map => $map2, $value => 50);
  $mappable->add_position($position3);

  # get all the places our map element is found, on a particular map of interest
  foreach $pos ($mappable->get_positions($map1)) {
     print $pos->value, "\n";
  }

=head1 DESCRIPTION

This object handles the notion of a generic map element. Mappables are
entities with one or more positions on one or more maps.

This object is a pure perl implementation of L<Bio::Map::MappableI>. That
interface implements some of its own methods so check the docs there for
those.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing list.  Your participation is much appreciated.

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

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::Mappable;
use strict;
use Bio::Map::Relative;
use Bio::Map::Position;

use base qw(Bio::Root::Root Bio::Map::MappableI);

=head2 new

 Title   : new
 Usage   : my $mappable = Bio::Map::Mappable->new();
 Function: Builds a new Bio::Map::Mappable object
 Returns : Bio::Map::Mappable
 Args    : -name => string : name of the mappable element
           -id   => string : id of the mappable element

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    my ($name, $id) = $self->_rearrange([qw(NAME ID)], @args);
    $self->name($name) if $name;
    $self->id($id) if $id;
    
    return $self;
}

=head2 name

 Title   : name
 Usage   : $mappable->name($new_name);
	       my $name = $mappable->name();
 Function: Get/Set the name for this Mappable
 Returns : A scalar representing the current name of this Mappable
 Args    : none to get
           string to set

=cut

sub name {
    my $self = shift;
    if (@_) { $self->{_name} = shift }
    return $self->{_name} || '';
}

=head2 id

 Title   : id
 Usage   : my $id = $mappable->id();
           $mappable->id($new_id);
 Function: Get/Set the id for this Mappable.
 Returns : A scalar representing the current id of this Mappable
 Args    : none to get
           string to set

=cut

sub id {
    my $self = shift;
    if (@_) { $self->{_id} = shift }
    return $self->{_id} || return;
}

=head2 in_map

 Title   : in_map
 Usage   : if ($mappable->in_map($map)) {...}
 Function: Tests if this mappable is found on a specific map
 Returns : boolean
 Args    : L<Bio::Map::MapI>

=cut

sub in_map {
	my ($self, $query_map) = @_;
	$self->throw("Must supply an argument") unless $query_map;
    $self->throw("This is [$query_map], not an object") unless ref($query_map);
    $self->throw("This is [$query_map], not a Bio::Map::MapI object") unless $query_map->isa('Bio::Map::MapI');
    
    foreach my $map ($self->known_maps) {
        ($map eq $query_map) && return 1;
    }
    
    return 0;
}

=head2 Comparison methods

=cut

=head2 equals

 Title   : equals
 Usage   : if ($mappable->equals($other_mappable)) {...}
           my @equal_positions = $mappable->equals($other_mappable);
 Function: Finds the positions in this mappable that are equal to any
           comparison positions.
 Returns : array of L<Bio::Map::PositionI> objects
 Args    : arg #1 = L<Bio::Map::MappableI> OR L<Bio::Map::PositionI> to compare
                    this one to (mandatory)
           arg #2 = optionally, one or more of the key => value pairs below
		   -map => MapI           : a Bio::Map::MapI to only consider positions
		                            on the given map
		   -relative => RelativeI : a Bio::Map::RelativeI to calculate in terms
                                    of each Position's relative position to the
                                    thing described by that Relative

=cut

sub equals {
    my $self = shift;
    return $self->_compare('equals', @_);
}

=head2 less_than

 Title   : less_than
 Usage   : if ($mappable->less_than($other_mappable)) {...}
           my @lesser_positions = $mappable->less_than($other_mappable);
 Function: Finds the positions in this mappable that are less than all
           comparison positions.
 Returns : array of L<Bio::Map::PositionI> objects
 Args    : arg #1 = L<Bio::Map::MappableI> OR L<Bio::Map::PositionI> to compare
                    this one to (mandatory)
           arg #2 = optionally, one or more of the key => value pairs below
		   -map => MapI           : a Bio::Map::MapI to only consider positions
		                            on the given map
		   -relative => RelativeI : a Bio::Map::RelativeI to calculate in terms
                                    of each Position's relative position to the
                                    thing described by that Relative

=cut

sub less_than {
    my $self = shift;
    return $self->_compare('less_than', @_);
}

=head2 greater_than

 Title   : greater_than
 Usage   : if ($mappable->greater_than($other_mappable)) {...}
           my @greater_positions = $mappable->greater_than($other_mappable);
 Function: Finds the positions in this mappable that are greater than all
           comparison positions.
 Returns : array of L<Bio::Map::PositionI> objects
 Args    : arg #1 = L<Bio::Map::MappableI> OR L<Bio::Map::PositionI> to compare
                    this one to (mandatory)
           arg #2 = optionally, one or more of the key => value pairs below
		   -map => MapI           : a Bio::Map::MapI to only consider positions
		                            on the given map
		   -relative => RelativeI : a Bio::Map::RelativeI to calculate in terms
                                    of each Position's relative position to the
                                    thing described by that Relative

=cut

sub greater_than {
    my $self = shift;
    return $self->_compare('greater_than', @_);
}

=head2 overlaps

 Title   : overlaps
 Usage   : if ($mappable->overlaps($other_mappable)) {...}
           my @overlapping_positions = $mappable->overlaps($other_mappable);
 Function: Finds the positions in this mappable that overlap with any
           comparison positions.
 Returns : array of L<Bio::Map::PositionI> objects
 Args    : arg #1 = L<Bio::Map::MappableI> OR L<Bio::Map::PositionI> to compare
                    this one to (mandatory)
           arg #2 = optionally, one or more of the key => value pairs below
		   -map => MapI           : a Bio::Map::MapI to only consider positions
		                            on the given map
		   -relative => RelativeI : a Bio::Map::RelativeI to calculate in terms
                                    of each Position's relative position to the
                                    thing described by that Relative

=cut

sub overlaps {
    my $self = shift;
    return $self->_compare('overlaps', @_);
}

=head2 contains

 Title   : contains
 Usage   : if ($mappable->contains($other_mappable)) {...}
           my @container_positions = $mappable->contains($other_mappable);
 Function: Finds the positions in this mappable that contain any comparison
           positions.
 Returns : array of L<Bio::Map::PositionI> objects
 Args    : arg #1 = L<Bio::Map::MappableI> OR L<Bio::Map::PositionI> to compare
                    this one to (mandatory)
           arg #2 = optionally, one or more of the key => value pairs below
		   -map => MapI           : a Bio::Map::MapI to only consider positions
		                            on the given map
		   -relative => RelativeI : a Bio::Map::RelativeI to calculate in terms
                                    of each Position's relative position to the
                                    thing described by that Relative

=cut

sub contains {
    my $self = shift;
    return $self->_compare('contains', @_);
}

=head2 overlapping_groups

 Title   : overlapping_groups
 Usage   : my @groups = $mappable->overlapping_groups($other_mappable);
           my @groups = Bio::Map::Mappable->overlapping_groups(\@mappables);
 Function: Look at all the positions of all the supplied mappables and group
           them according to overlap.
 Returns : array of array refs, each ref containing the Bio::Map::PositionI
           objects that overlap with each other
 Args    : arg #1 = L<Bio::Map::MappableI> OR L<Bio::Map::PositionI> to  compare
                    this one to, or an array ref of such objects (mandatory)
           arg #2 = optionally, one or more of the key => value pairs below
		   -map => MapI           : a Bio::Map::MapI to only consider positions
		                            on the given map
		   -relative => RelativeI : a Bio::Map::RelativeI to calculate in terms
                                    of each Position's relative position to the
                                    thing described by that Relative
           -min_pos_num => int    : the minimum number of positions that must
                                    be in a group before it will be returned
                                    [default is 1]
           -min_mappables_num => int        : the minimum number of different
                                              mappables represented by the
                                              positions in a group before it
                                              will be returned [default is 1]
           -min_mappables_percent => number : as above, but the minimum
                                              percentage of input mappables
                                              [default is 0]
           -min_map_num => int              : the minimum number of different
                                              maps represented by the positions
                                              in a group before it will be
                                              returned [default is 1]
           -min_map_percent => number       : as above, but the minimum
                                              percentage of maps known by the
                                              input mappables [default is 0]
           -require_self => 1|0             : require that at least one of the
                                              calling object's positions be in
                                              each group [default is 1, has no
                                              effect when the second usage form
                                              is used]
           -required => \@mappables         : require that at least one position
                                              for each mappable supplied in this
                                              array ref be in each group

=cut

sub overlapping_groups {
    my $self = shift;
    return $self->_compare('overlapping_groups', @_);
}

=head2 disconnected_intersections

 Title   : disconnected_intersections
 Usage   : @positions = $mappable->disconnected_intersections($other_mappable);
           @positions = Bio::Map::Mappable->disconnected_intersections(\@mappables);
 Function: Make the positions that are at the intersection of each group of
           overlapping positions, considering all the positions of the supplied
           mappables.
 Returns : new Bio::Map::Mappable who's positions on maps are the calculated
           disconnected unions
 Args    : arg #1 = L<Bio::Map::MappableI> OR L<Bio::Map::PositionI> to  compare
                    this one to, or an array ref of such objects (mandatory)
           arg #2 = optionally, one or more of the key => value pairs below
		   -map => MapI           : a Bio::Map::MapI to only consider positions
		                            on the given map
		   -relative => RelativeI : a Bio::Map::RelativeI to calculate in terms
                                    of each Position's relative position to the
                                    thing described by that Relative
           -min_pos_num => int    : the minimum number of positions that must
                                    be in a group before the intersection will
                                    be calculated and returned [default is 1]
           -min_mappables_num => int        : the minimum number of different
                                              mappables represented by the
                                              positions in a group before the
                                              intersection will be calculated
                                              and returned [default is 1]
           -min_mappables_percent => number : as above, but the minimum
                                              percentage of input mappables
                                              [default is 0]
           -min_map_num => int              : the minimum number of different
                                              maps represented by the positions
                                              in a group before the intersection
                                              will be calculated and returned
                                              [default is 1]
           -min_map_percent => number       : as above, but the minimum
                                              percentage of maps known by the
                                              input mappables [default is 0]
           -require_self => 1|0             : require that at least one of the
                                              calling object's positions be in
                                              each group [default is 1, has no
                                              effect when the second usage form
                                              is used]
           -required => \@mappables         : require that at least one position
                                              for each mappable supplied in this
                                              array ref be in each group

=cut

sub disconnected_intersections {
    my $self = shift;
    return $self->_compare('intersection', @_);
}

=head2 disconnected_unions

 Title   : disconnected_unions
 Usage   : my @positions = $mappable->disconnected_unions($other_mappable);
           my @positions = Bio::Map::Mappable->disconnected_unions(\@mappables);
 Function: Make the positions that are the union of each group of overlapping
           positions, considering all the positions of the supplied mappables.
 Returns : new Bio::Map::Mappable who's positions on maps are the calculated
           disconnected unions
 Args    : arg #1 = L<Bio::Map::MappableI> OR L<Bio::Map::PositionI> to  compare
                    this one to, or an array ref of such objects (mandatory)
           arg #2 = optionally, one or more of the key => value pairs below
		   -map => MapI           : a Bio::Map::MapI to only consider positions
		                            on the given map
		   -relative => RelativeI : a Bio::Map::RelativeI to calculate in terms
                                    of each Position's relative position to the
                                    thing described by that Relative
           -min_pos_num => int    : the minimum number of positions that must
                                    be in a group before the union will be
                                    calculated and returned [default is 1]
           -min_mappables_num => int        : the minimum number of different
                                              mappables represented by the
                                              positions in a group before the
                                              union will be calculated and
                                              returned [default is 1]
           -min_mappables_percent => number : as above, but the minimum
                                              percentage of input mappables
                                              [default is 0]
           -min_map_num => int              : the minimum number of different
                                              maps represented by the positions
                                              in a group before the union will
                                              be calculated and returned
                                              [default is 1]
           -min_map_percent => number       : as above, but the minimum
                                              percentage of maps known by the
                                              input mappables [default is 0]
           -require_self => 1|0             : require that at least one of the
                                              calling object's positions be in
                                              each group [default is 1, has no
                                              effect when the second usage form
                                              is used]
           -required => \@mappables         : require that at least one position
                                              for each mappable supplied in this
                                              array ref be in each group

=cut

sub disconnected_unions {
    my $self = shift;
    return $self->_compare('union', @_);
}

# do a RangeI-related comparison by calling the corresponding PositionI method
# on all the requested Positions of our Mappables
sub _compare {
    my ($self, $method, $input, @extra_args) = @_;
    $self->throw("Must supply an object or array ref of them") unless ref($input);
    $self->throw("Wrong number of extra args (should be key => value pairs)") unless @extra_args % 2 == 0;
    my @compares = ref($input) eq 'ARRAY' ? @{$input} : ($input);
    
    my %args = (-map => undef, -relative => undef, -min_pos_num => 1,
                -min_mappables_num => 1, -min_mappables_percent => 0,
                -min_map_num => 1, -min_map_percent => 0,
                -require_self => 0, -required => undef, -min_overlap_percent => 0, @extra_args);
    my $map = $args{-map};
    my $rel = $args{-relative};
    my $overlap = $args{-min_overlap_percent};
    my $min_pos_num = $args{-min_pos_num};
    my $min_pables_num = $args{-min_mappables_num};
    if ($args{-min_mappables_percent}) {
        my $mn = (@compares + (ref($self) ? 1 : 0)) / 100 * $args{-min_mappables_percent};
        if ($mn > $min_pables_num) {
            $min_pables_num = $mn;
        }
    }
    my $min_map_num = $args{-min_map_num};
    if ($args{-min_map_percent}) {
        my %known_maps;
        foreach my $pable (@compares, ref($self) ? ($self) : ()) {
            foreach my $known ($pable->known_maps) {
                $known_maps{$known->unique_id} = 1;
            }
        }
        my $mn = scalar(keys %known_maps) / 100 * $args{-min_map_percent};
        if ($mn > $min_map_num) {
            $min_map_num = $mn;
        }
    }
    my %required = map { $_ => 1 } $args{-required} ? @{$args{-required}} : ();
    my (@mine, @yours);
    
    if (ref($self)) {
        @mine = $self->get_positions($map);
        if ($args{-require_self}) {
            @mine > 0 or return;
            $required{$self} = 1;
        }
    }
    my @required = sort keys %required;
    
    foreach my $compare (@compares) {
        if ($compare->isa('Bio::Map::PositionI')) {
            push(@yours, $compare);
        }
        elsif ($compare->isa('Bio::Map::MappableI')) {
            push(@yours, $compare->get_positions($map));
        }
        else {
            $self->throw("This is [$compare], not a Bio::Map::MappableI or Bio::Map::PositionI");
        }
    }
    @yours > 0 or return;
    
    my @ok;
    SWITCH: for ($method) {
        /equals|overlaps|contains/ && do {
            @mine > 0 or return;
            foreach my $my_pos (@mine) {
                foreach my $your_pos (@yours) {
                    if ($my_pos->$method($your_pos, undef, $rel)) {
                        push(@ok, $my_pos);
                        last;
                    }
                }
            }
            last SWITCH;
        };
        /less_than|greater_than/ && do {
            @mine > 0 or return;
            if ($method eq 'greater_than') {
                @mine =  map { $_->[1] }
                         sort { $b->[0] <=> $a->[0] }
                         map { [$_->end($_->absolute_relative), $_] }
                         @mine;
                @yours = map { $_->[1] }
                         sort { $b->[0] <=> $a->[0] }
                         map { [$_->end($_->absolute_relative), $_] }
                         @yours;
            }
            my $test_pos = shift(@yours);
            
            foreach my $my_pos (@mine) {
                if ($my_pos->$method($test_pos, $rel)) {
                    push(@ok, $my_pos);
                }
                else {
                    last;
                }
            }
            
            if ($method eq 'greater_than') {
                @ok = map { $_->[1] }
                      sort { $a->[0] <=> $b->[0] }
                      map { [$_->sortable, $_] }
                      @ok;
            }
            
            last SWITCH;
        };
        /overlapping_groups|intersection|union/ && do {
            my @positions = (@mine, @yours);
            my $start_pos = shift(@positions);
            
            my $dr_able = $start_pos->disconnected_ranges(\@positions, $rel, $overlap) || return;
            my @disconnected_ranges = $dr_able->get_positions;
            
            #print "got ", scalar(@disconnected_ranges), " disconnected_ranges, first has range ", $disconnected_ranges[0]->toString, "\n";
            
            #use Benchmark qw(:all);
            #my $t0 = new Benchmark;
            
            my %all_groups;
            my %done_ranges;
            for my $i (0..$#disconnected_ranges) {
                my $range = $disconnected_ranges[$i];
                my $range_string = $range->toString;
                next if $done_ranges{$range_string};
                $done_ranges{$range_string} = 1;
                
                foreach my $pos ($start_pos, @positions) {
                    if ($pos->overlaps($range, undef, $rel)) {
                        $all_groups{$range_string}->{$pos} = $pos;
                    }
                }
            }
            
            #my $t1 = new Benchmark;
            #my $td = timediff($t1, $t0);
            #print "grouping took: ",timestr($td),"\n";
            
            # purge the temporary working (not $dr_able->purge_positions since
            # that removes the element from each position, but leaves it on
            # the map. *** need complete purge that removes position from
            # memory...
            foreach my $pos (@disconnected_ranges) {
                my $map = $pos->map || next;
                $map->purge_positions($pos);
            }
            
            my @groups;
            GROUPS: foreach my $group_range (sort keys %all_groups) { 
                my $group = $all_groups{$group_range};
                my @group = sort values %{$group};
                #print "* in group $group_range, there are ", scalar(@group), " members\n";
                
                @group >= $min_pos_num or next;
                @group >= $min_pables_num or next; # shortcut before having to work it out properly
                @group >= $min_map_num or next; # shortcut before having to work it out properly
                
                my %mappables;
                foreach my $pos (@group) {
                    my $mappable = $pos->element || next;
                    $mappables{$mappable} = 1;
                }
                keys %mappables >= $min_pables_num || next;
                
                my %maps;
                foreach my $pos (@group) {
                    my $map = $pos->map || next;
                    $maps{$map->unique_id} = 1;
                }
                keys %maps >= $min_map_num || next;
                
                foreach my $required (@required) {
                    exists $mappables{$required} or next GROUPS;
                }
                
                my @sorted = map { $_->[1] }
                             sort { $a->[0] <=> $b->[0] }
                             map { [$_->sortable, $_] }
                             @group;
                push(@groups, \@sorted);
            }
            
            if ($method eq 'overlapping_groups') {
                return @groups;
            }
            else {
                foreach my $group (@groups) {
                    my $start_pos = shift(@{$group});
                    
                    unless (@{$group}) {
                        # we'll consider the 'intersection' or 'union' of just
                        # one position as the position itself
                        push(@ok, Bio::Map::Position->new(-map => $start_pos->map,
                                                          -start => $start_pos->start,
                                                          -end => $start_pos->end));
                    }
                    else {
                        my @rel_arg = $method eq 'intersection' ? (undef, $rel) : ($rel);
                        my $result = $start_pos->$method($group, @rel_arg) || next;
                        push(@ok, $result->get_positions);
                    }
                }
                
                # assign all the positions to a result mappable
                my $result = $self->new();
                $result->add_position(@ok) if @ok; # add_position can actually take a list
                
                return $result;
            }
            
            last SWITCH;
        };
        
        $self->throw("Unknown method '$method'");
    }
    
    return @ok;
}

=head2 tuple

 Title   : tuple
 Usage   : Do Not Use!
 Function: tuple was supposed to be a private method; this method no longer
           does anything
 Returns : warning
 Args    : none
 Status  : deprecated, will be removed in next version

=cut

sub tuple {
    my $self = shift;
    $self->warn("The tuple method was supposed to be a private method, don't call it!");
}

=head2 annotation

 Title   : annotation
 Usage   : $mappable->annotation($an_col);
           my $an_col = $mappable->annotation();
 Function: Get the annotation collection (see Bio::AnnotationCollectionI)
           for this annotatable object.
 Returns : a Bio::AnnotationCollectionI implementing object, or undef
 Args    : none to get, OR
           a Bio::AnnotationCollectionI implementing object to set

=cut

sub annotation {
    my $self = shift;
    if (@_) { $self->{_annotation} = shift }
    return $self->{_annotation} || return;
}

1;
