package Bio::Coordinate::Collection;
use utf8;
use strict;
use warnings;
use Bio::Coordinate::Result;
use Bio::Coordinate::Result::Gap;
use parent qw(Bio::Root::Root Bio::Coordinate::MapperI);

# ABSTRACT: Noncontinuous match between two coordinate sets.
# AUTHOR:   Heikki Lehvaslaiho <heikki@bioperl.org>
# OWNER:    Heikki Lehvaslaiho
# LICENSE:  Perl_5

# CONTRIBUTOR: Ewan Birney <birney@ebi.ac.uk>

=head1 SYNOPSIS

  # create Bio::Coordinate::Pairs or other Bio::Coordinate::MapperIs somehow
  $pair1; $pair2;

  # add them into a Collection
  $collection = Bio::Coordinate::Collection->new;
  $collection->add_mapper($pair1);
  $collection->add_mapper($pair2);

  # create a position and map it
  $pos = Bio::Location::Simple->new (-start => 5, -end => 9 );
  $res = $collection->map($pos);
  $res->match->start == 1;
  $res->match->end == 5;

  # if mapping is many to one (*>1) or many-to-many (*>*)
  # you have to give seq_id not get unrelevant entries
  $pos = Bio::Location::Simple->new
      (-start => 5, -end => 9 -seq_id=>'clone1');

=head1 DESCRIPTION

Generic, context neutral mapper to provide coordinate transforms
between two B<disjoint> coordinate systems. It brings into Bioperl the
functionality from Ewan Birney's Bio::EnsEMBL::Mapper ported into
current bioperl.

This class is aimed for representing mapping between whole chromosomes
and contigs, or between contigs and clones, or between sequencing
reads and assembly. The submaps are automatically sorted, so they can
be added in any order.

To map coordinates to the other direction, you have to swap() the
collection. Keeping track of the direction and ID restrictions
are left to the calling code.

=cut

=head2 new
=cut

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->{'_mappers'} = [];

    my($in, $out, $strict, $mappers, $return_match) =
        $self->_rearrange([qw(IN
                              OUT
                              STRICT
                              MAPPERS
                              RETURN_MATCH
                             )],
                         @args);

    $in  && $self->in($in);
    $out  && $self->out($out);
    $mappers && $self->mappers($mappers);
    $return_match && $self->return_match('return_match');
    return $self; # success - we hope!
}

=head2 add_mapper

 Title   : add_mapper
 Usage   : $obj->add_mapper($mapper)
 Function: Pushes one Bio::Coordinate::MapperI into the list of mappers.
           Sets _is_sorted() to false.
 Example :
 Returns : 1 when succeeds, 0 for failure.
 Args    : mapper object

=cut

sub add_mapper {
  my ($self,$value) = @_;

  $self->throw("Is not a Bio::Coordinate::MapperI but a [$self]")
      unless defined $value && $value->isa('Bio::Coordinate::MapperI');

  # test pair range lengths
  $self->warn("Coordinates in pair [". $value . ":" .
              $value->in->seq_id . "/". $value->out->seq_id .
              "] are not right.")
      unless $value->test;

  $self->_is_sorted(0);
  push(@{$self->{'_mappers'}},$value);
}

=head2 mappers

 Title   : mappers
 Usage   : $obj->mappers();
 Function: Returns or sets a list of mappers.
 Example :
 Returns : array of mappers
 Args    : array of mappers

=cut

sub mappers{
        my ($self,@args) = @_;

        if (@args) {
                if (@args == 1 && ref $args[0] eq 'ARRAY') {
                        @args = @{$args[0]};
                }
                $self->throw("Is not a Bio::Coordinate::MapperI but a [$self]")
                        unless defined $args[0] && $args[0]->isa('Bio::Coordinate::MapperI');
                push(@{$self->{'_mappers'}}, @args);
        }

        return @{$self->{'_mappers'}};
}

=head2 each_mapper

 Title   : each_mapper
 Usage   : $obj->each_mapper();
 Function: Returns a list of mappers.
 Example :
 Returns : list of mappers
 Args    : none

=cut

sub each_mapper{
   my ($self) = @_;
   return @{$self->{'_mappers'}};
}

=head2 mapper_count

 Title   : mapper_count
 Usage   : my $count = $collection->mapper_count;
 Function: Get the count of the number of mappers stored
           in this collection
 Example :
 Returns : integer
 Args    : none

=cut

sub mapper_count{
   my $self = shift;
   return scalar @{$self->{'_mappers'} || []};
}

=head2 swap

 Title   : swap
 Usage   : $obj->swap;
 Function: Swap the direction of mapping;input <-> output
 Example :
 Returns : 1
 Args    :

=cut

sub swap {
   my ($self) = @_;
   use Data::Dumper;

   $self->sort unless $self->_is_sorted;
   map {$_->swap;} @{$self->{'_mappers'}};
   ($self->{'_in_ids'}, $self->{'_out_ids'}) =
       ($self->{'_out_ids'}, $self->{'_in_ids'});
   1;
}

=head2 test

 Title   : test
 Usage   : $obj->test;
 Function: test that both components of all pairs are of the same length.
           Ran automatically.
 Example :
 Returns : boolean
 Args    :

=cut

sub test {
   my ($self) = @_;

   my $res = 1;

   foreach my $mapper ($self->each_mapper) {
       unless( $mapper->test ) {
           $self->warn("Coordinates in pair [". $mapper . ":" .
                       $mapper->in->seq_id . "/". $mapper->out->seq_id .
                       "] are not right.");
           $res = 0;
       }
   }
   $res;
}

=head2 map

 Title   : map
 Usage   : $newpos = $obj->map($pos);
 Function: Map the location from the input coordinate system
           to a new value in the output coordinate system.
 Example :
 Returns : new value in the output coordinate system
 Args    : integer

=cut

sub map {
   my ($self,$value) = @_;

   $self->throw("Need to pass me a value.")
       unless defined $value;
   $self->throw("I need a Bio::Location, not [$value]")
       unless $value->isa('Bio::LocationI');
   $self->throw("No coordinate mappers!")
       unless $self->each_mapper;

   $self->sort unless $self->_is_sorted;

   if ($value->isa("Bio::Location::SplitLocationI")) {

       my $result = Bio::Coordinate::Result->new();
       foreach my $loc ( $value->sub_Location(1) ) {

           my $res = $self->_map($loc);
           map { $result->add_sub_Location($_) } $res->each_Location;

       }
       return $result;

   } else {
       return $self->_map($value);
   }

}

=head2 _map

 Title   : _map
 Usage   : $newpos = $obj->_map($simpleloc);
 Function: Internal method that does the actual mapping. Called multiple times
           by map() if the location  to be mapped is a split location

 Example :
 Returns : new location in the output coordinate system or undef
 Args    : Bio::Location::Simple

=cut

sub _map {
   my ($self,$value) = @_;

   my $result = Bio::Coordinate::Result->new(-is_remote=>1);

IDMATCH: {

       # bail out now we if are forcing the use of an ID
       # and it is not in this collection
       last IDMATCH if defined $value->seq_id &&
           ! $self->{'_in_ids'}->{$value->seq_id};

       foreach my $pair ($self->each_mapper) {

           # if we are limiting input to a certain ID
           next if defined $value->seq_id && $value->seq_id ne $pair->in->seq_id;

           # if we haven't even reached the start, move on
           next if $pair->in->end < $value->start;
           # if we have over run, break
           last if $pair->in->start > $value->end;

           my $subres = $pair->map($value);
           $result->add_result($subres);
       }
   }

   $result->seq_id($result->match->seq_id) if $result->match;
   unless ($result->each_Location) {
       #build one gap;
       my $gap = Bio::Location::Simple->new(-start => $value->start,
                                            -end => $value->end,
                                            -strand => $value->strand,
                                            -location_type => $value->location_type
                                           );
       $gap->seq_id($value->seq_id) if defined $value->seq_id;
       bless $gap, 'Bio::Coordinate::Result::Gap';
       $result->seq_id($value->seq_id) if defined $value->seq_id;
       $result->add_sub_Location($gap);
   }
   return $result;
}

=head2 sort

 Title   : sort
 Usage   : $obj->sort;
 Function: Sort function so that all mappings are sorted by
           input coordinate start
 Example :
 Returns : 1
 Args    :

=cut

sub sort{
   my ($self) = @_;

   @{$self->{'_mappers'}} = map { $_->[0] }
                            sort { $a->[1] <=> $b->[1] }
                            map { [ $_, $_->in->start] }
                            @{$self->{'_mappers'}};

   #create hashes for sequence ids
   $self->{'_in_ids'} = ();
   $self->{'_out_ids'} = ();
   foreach ($self->each_mapper) {
       $self->{'_in_ids'}->{$_->in->seq_id} = 1;
       $self->{'_out_ids'}->{$_->out->seq_id} = 1;
   }

   $self->_is_sorted(1);
}

=head2 _is_sorted

 Title   : _is_sorted
 Usage   : $newpos = $obj->_is_sorted;
 Function: toggle for whether the (internal) coodinate mapper data are sorted
 Example :
 Returns : boolean
 Args    : boolean

=cut

sub _is_sorted{
   my ($self,$value) = @_;

   $self->{'_is_sorted'} = 1 if defined $value && $value;
   return $self->{'_is_sorted'};
}

1;
