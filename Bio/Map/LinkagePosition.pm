# BioPerl module for Bio::Map::LinkagePosition
#
# Cared for by Chad Matsalla <bioinformatics1@dieselwurks.com>
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::LinkagePosition - Create a Position for a Marker that will be placed
	on a Bio::Map::LinkageMap

=head1 SYNOPSIS

    use Bio::Map::Position;
    my $position = new Bio::Map::LinkagePosition(-positions => 1,
						 -distance => 22.1 );
    
	    # can get listing of positions
    my @positions = $position->each_position;


=head1 DESCRIPTION

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Chad Matsalla

Email bioinformatics1@dieselwurks.com

=head1 CONTRIBUTORS

Lincoln Stein, lstein@cshl.org
Heikki Lehvaslaiho, heikki@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Map::LinkagePosition;
use vars qw(@ISA);
use strict;
require 'dumpvar.pl';

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::Map::PositionI;

@ISA = qw(Bio::Root::Root Bio::Map::PositionI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Map::LinkagePosition(-positions => $position,
				-distance => $distance );
 Function: Builds a new Bio::Map::LinkagePosition object 
 Returns : Bio::Map::LinkagePosition
 Args    : -positions => the relative order of this marker on a linkage map
	-distance => the centimorgam distance of this marker from the previous
		marker. Can be 0!
 Notes   : In this case, if -position is a list the first element will be used
	as the position for this marker. It is much better to simply provide a
	scalar for this value in the context of a LinkagePosition object.
	( position_s_ is used simply because this object inherits from
	PositionI )
	If -distance = 0 or is omitted it is assumed that the marker here
	is linked to the marker before it.	

=cut

sub new {
  my($class,@args) = @_;
	my $self = $class->SUPER::new(@args);
	my %param = @args;
	my ($positions,$distance) = $self->_rearrange([qw(POSITIONS DISTANCE)], @args);
	if( ref($positions) =~ /array/i ) {
		$self->add_position(pop(@$positions));
	}
	else {
		$self->add_position($positions);
	}
	if ($distance) {
		$self->distance($distance);
	}
	else {
		$self->distance('0.0');
	}
  return bless $self;
}


=head2 Bio::Map::PositionI methods

=head2 each_position

 Title   : each_position
 Usage   : my @positions = $position->each_position();
 Function: Retrieve a list of positions
 Returns : An array of position values containing a single element.
 Args    : 
 Notes   : I would like this to return a scalar but that would
	break the interface.

=cut

sub each_position {
   my ($self) = @_;
   return @{$self->{'_positions'}};
}

=head2 add_position($position)

 Title   : add_position($position)
 Usage   : $position->add_position('100')
 Function: Add a position to the LinkagePosition container
 Returns : Nothing.
 Args    : String or Numeric coding for a position on a map
 Notes   : In the context of a LinkagePosition this is somewhat
	moot but is included for compatibility and fascist interface
	conformation. What should be done if you use this method
	when there already is an element in the _positions array?
	As it is, the element you pass in as an argument will
	_replace_ the current one. _check the source_, Luke.

=cut

sub add_position {
   my ($self,$value) = @_;
   if( ! $value ) { 
       $self->warn("Attempting to add a position with a null value");
       return;
   }
   push @{$self->{'_positions'}}, $value;
   return;
}

=head2 purge()

 Title   : purge()
 Usage   : $position->purge()
 Function: Remove all the position values stored for a position
 Returns : Nothing
 Args    : None

=cut

sub purge_positions {
   my ($self) = @_;
   $self->{'_positions'} = [];
}

=head2 distance($new_distance)

 Title   : distance($new_distance)
 Usage   : $position->distance(new_distance) _or_
	$position->distance()
 Function: get/set the distance of this position from the previous marker
 Returns : A scalar representing the current distance for this position.
 Args    : If $new_distance is provided the distance of this Position will
	be set to $new_distance

=cut

sub distance {
	my ($self,$distance) = @_;
	if ($distance) {
	   $self->{'_distance'} = $distance;
	}
	return $self->{'_distance'};	
}

=head2 position($new_postion)

 Title   : Get/set the position for this LinkagePosition
 Usage   : $o_position->position($new_position) _or_
	$o_position->position()
 Function: get/set the position of this LinkagePosition
 Returns : A scalar representing the current position.
 Args    : If $new_position is provided, the current position of this Position
	will be set to $new_position.
 Notes   : Ha! This is _mine_ so it returns a scalar. Bwahaha.

=cut

sub position {
	my ($self,$position) = @_;
	if ($position) {
		# no point in keeing the old ones
		$self->purge_positions();
		$self->add_position($position);
	}
		# ::dumpValue($self);
	return $self->{'_positions'};	
}



1;
