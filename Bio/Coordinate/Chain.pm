#
# bioperl module for Bio::Coordinate::Chain
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Coordinate::Chain - Mapping locations through a chain of  coordinate mappers

=head1 SYNOPSIS

  # create Bio::Coordinate::Pairs, or any MapperIs, somehow
  $pair1; $pair2;

  # add them into a Chain
  $collection = Bio::Coordinate::Chain->new;
  $collection->add_mapper($pair1);
  $collection->add_mapper($pair2);

  # create a position and map it
  $pos = Bio::Location::Simple->new (-start => 5, -end => 9 );
  $match = $collection->map($pos);
  if ($match) {
      sprintf "Matches at %d-%d\n", $match->start, $match->end,
  } else {
      print "No match\n";
  }

=head1 DESCRIPTION

This class assumes that you have built several mappers and want to
link them together so that output from the previous mapper is the next
mappers input. This way you can build arbitrarily complex mappers from
simpler components.

Note that Chain does not do any sanity checking on its mappers. You
are solely responsible that input and output coordinate systems,
direction of mapping and parameters internal to mappers make sense
when chained together.

To put it bluntly, the present class is just a glorified foreach loop
over an array of mappers calling the map method.

It would be neat to an internal function that would generate a new
single step mapper from those included in the chain. It should speed
things up considerably. Any volunteers?

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

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

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki-at-bioperl-dot-org

=head1 CONTRIBUTORS

Ewan Birney, birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Coordinate::Chain;
use strict;

# Object preamble - inherits from Bio::Root::Root
use Bio::Root::Root;
use Bio::Coordinate::Result;

use base qw(Bio::Coordinate::Collection Bio::Coordinate::MapperI);


=head2 map

 Title   : map
 Usage   : $newpos = $obj->map($pos);
 Function: Map the location through all the mappers in the chain.
 Example :
 Returns : new Location in the output coordiante system
 Args    : a Bio::Location::Simple object

=cut

sub map {
    my ($self,$value) = @_;

    $self->throw("Need to pass me a value.")
	unless defined $value;
    $self->throw("I need a Bio::Location, not [$value]")
	unless $value->isa('Bio::LocationI');
    $self->throw("No coordinate mappers!")
	unless $self->each_mapper;

    my $res = Bio::Coordinate::Result->new();

    foreach my $mapper ($self->each_mapper) {

	my $res = $mapper->map($value);
	return unless $res->each_match;
	$value = $res->match;
    }

   return $value;
}


=head2 Inherited methods

=cut

=head2 add_mapper

 Title   : add_mapper
 Usage   : $obj->add_mapper($mapper)
 Function: Pushes one Bio::Coodinate::MapperI into the list of mappers.
           Sets _is_sorted() to false.
 Example : 
 Returns : 1 when succeeds, 0 for failure.
 Args    : mapper object

=cut

=head2 mappers

 Title   : mappers
 Usage   : $obj->mappers();
 Function: Returns or sets a list of mappers.
 Example : 
 Returns : array of mappers
 Args    : array of mappers

=cut

=head2 each_mapper

 Title   : each_mapper
 Usage   : $obj->each_mapper();
 Function: Returns a list of mappers.
 Example : 
 Returns : array of mappers
 Args    : none

=cut

=head2 swap

 Title   : swap
 Usage   : $obj->swap;
 Function: Swap the direction of mapping;input <-> output
 Example :
 Returns : 1
 Args    : 

=cut

=head2 test

 Title   : test
 Usage   : $obj->test;
 Function: test that both components of all pairs are of the same length.
           Ran automatically.
 Example :
 Returns : boolean
 Args    :

=cut



sub sort{
   my ($self) = @_;
   $self->warn("You do not really want to sort your chain, do you!\nDoing nothing.");
}

1;

