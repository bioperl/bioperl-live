#
# BioPerl module for Bio::MapIO::mapmaker
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::MapIO::mapmaker - A Mapmaker Map reader

=head1 SYNOPSIS

# do not use this object directly it is accessed through the Bio::MapIO system 

    use Bio::MapIO;
    my $mapio = Bio::MapIO->new(-format => "mapmaker",
			                      -file   => "mapfile.map");
    while ( my $map = $mapio->next_map ) {  # get each map
	    foreach my $marker ( $map->each_element ) {
	       # loop through the markers associated with the map
	    }
    }

=head1 DESCRIPTION

This object contains code for parsing and processing Mapmaker output
and creating L<Bio::Map::MapI> objects from it.

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

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::MapIO::mapmaker;
use strict;

use Bio::Map::SimpleMap;
use Bio::Map::LinkagePosition;
use Bio::Map::Marker;

use base qw(Bio::MapIO);

=head2 next_map

 Title   : next_map
 Usage   : my $map = $factory->next_map;
 Function: Get one or more map objects from the Mapmaker input
 Returns : Bio::Map::MapI
 Args    : none

See L<Bio::Map::MapI>

=cut

sub next_map{
   my ($self) = @_;
   my $map = Bio::Map::SimpleMap->new(-name  => '',
												  -units => 'cM',
												  -type  => 'Genetic');

	# Mapmaker input can be free-form, like the result of a copy-paste
	# from a terminal, with no particular format before or after the 
	# map data. The $in_map variable is a flag that's set to 1 when 
	# we're reading map data lines and set back to 0 when we're finished.
   my ($in_map,$runningDistance);

   while ( defined ($_ = $self->_readline()) ) {
		if ( /^\s+Markers\s+Distance/ ) {
			$in_map = 1;
			next;
		} 
		next unless $in_map;
 
		s/ +/\t/;
		my ($number,$name,$distance) = split;
		$runningDistance += $distance unless ($distance =~ /-+/);
		$runningDistance = '0.0' if ($runningDistance == 0 || $distance =~ /-+/);

		my $pos = Bio::Map::LinkagePosition->new(-order => $number,
															 -map   => $map,
															 -value => $runningDistance );
		my $marker = Bio::Map::Marker->new(-name     => $name,
													 -position => $pos );
		
		if ($distance =~ /-+/) { # last marker
			$in_map = 0;
			return $map;
		}  
	}
}

=head2 write_map

 Title   : write_map
 Usage   : $factory->write_map($map);
 Function: Write a map out through the factory
 Returns : none
 Args    : Bio::Map::MapI

=cut

sub write_map{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

1;
