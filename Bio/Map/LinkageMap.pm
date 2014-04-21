# BioPerl module for Bio::Map::LinkageMap
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::LinkageMap - A representation of a genetic linkage map.

=head1 SYNOPSIS

    use Bio::Map::LinkageMap;
	# create a new map
    my $map = Bio::Map::LinkageMap->new(-name => 'Chads Superterriffic Map',
                                      -type => 'Linkage',
                                      -units=> 'cM');
	# create the location of a marker for that map
    my $position = Bio::Map::LinkagePosition->new( -positions => 1,
                -distance => "22.3");
	# create a marker and place it at that position
    my $marker = Bio::Map::Marker::Microsatellite->new(
			-name => 'SuuuperMarker',
			-position => $position);
	# place that marker on that map
    $map->add_element($marker);

	# done!

=head1 DESCRIPTION

This object describes the basic functionality of a genetic linkage map in
Bioperl. Each 'position' can have one or more markers that map some number of
units from the markers at the previous position.

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

=head1 AUTHOR - Chad Matsalla

Email bioinformatics1@dieselwurks.com

=head1 CONTRIBUTORS

Lincoln Stein       lstein@cshl.org
Heikki Lehvaslaiho  heikki-at-bioperl-dot-org
Jason Stajich       jason@bioperl.org
Sendu Bala          bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::LinkageMap;
use strict;

use base qw(Bio::Map::SimpleMap);

=head2 new

 Title   : new
 Usage   : my $linkage_map = Bio::Map::LinkageMap->new();
 Function: Builds a new Bio::Map::LinkageMap object
 Returns : Bio::Map::LinkageMap
 Args    : -name    => the name of the map (string) [optional]
	       -type    => the type of this map (string, defaults to Linkage) [optional]
           -species => species for this map (Bio::Species) [optional]
           -units   => the map units (string, defaults to cM) [optional]
           -elements=> elements to initialize with
                       (arrayref of Bio::Map::MappableI objects) [optional]
           -uid      => Unique ID of this map

=cut

=head2 length

 Title   : length
 Usage   : my $length = $map->length();
 Function: Retrieves the length of the map. In the case of a LinkageMap, the
	       length is the sum of all marker distances.
 Returns : An integer representing the length of this LinkageMap. Will return
	       0 if length is not calculateable
 Args    : None.


=cut

sub length {
    my ($self) = @_;
    $self->throw("Not yet implemented correctly");
    
    my $total_distance;
    foreach my $element (@{$self->get_elements}) {
        #*** there is no such method ->each_position_value!
        $total_distance += ($element->position->each_position_value($self))[0];
    }
    return $total_distance;
}

=head2 add_element($marker)

 Title   : add_element($marker)
 Usage   : $map->add_element($marker)
 Function: Add a Bio::Map::MappableI object to the Map
 Returns : none
 Args    : Bio::Map::MappableI object
 Notes   : It is strongly recommended that you use a
	   Bio::Map::LinkagePosition as the position in any
	   Bio::Map::Mappable that you create to place on this
	   map. Using some other Bio::Map::Position might work but might
	   be unpredictable.
           N.B. I've added Bio::Map::OrderedPosition which should achieve
                similar things from LinkagePosition and will work for
                RH markers too.
=cut

#*** what is this? what calls it? note that it seems to be private
sub _add_element_will_be_deleted {
    my ($self,$marker) = @_;

    my $o_position = $marker->position();

    $self->debug( "marker position is ". $marker->position());
#     print("add_element: \$o_position is $o_position\n");
#     print("add_element: \$marker is $marker\n");

    my $position;
    unless ( $o_position->isa('Bio::Map::LinkagePosition') ||
	     $o_position->isa('Bio::Map::OrderedPosition')
	     ) {
	$self->warn("You really should use a Linkage Position for this object. This insures that there is only one position. Trying anyway...");	
	my @p = ( $o_position->each_position_value($self));
	$position = shift @p;
	if( ! defined $position ) {
	    $self->throw("This marker ($marker) does not have a position in this map ($self)");
	}
    } else {
	$position = $o_position->order;
    }

    if ($self->{'_elements'}[$position]) {
	$self->warn("Replacing the marker in position $position because in a linkage map the position is a key.");
    }	
    $self->{'_elements'}[$position] = $marker;
}

1;
