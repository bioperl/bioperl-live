# $Id$
#
# BioPerl module for Bio::MapIO::mapmaker
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
    my $mapio = new Bio::MapIO(-format => "mapmaker",
			       -file   => "mapfile.map");
    while( my $map = $mapio->next_map ) { 
	# get each map
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

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::MapIO::mapmaker;
use vars qw(@ISA);
use strict;

use Bio::MapIO;
use Bio::Map::SimpleMap;
use Bio::Map::LinkagePosition;
use Bio::Map::Marker;

@ISA = qw(Bio::MapIO );

=head2 next_map

 Title   : next_tree
 Usage   : my $map = $factory->next_map;
 Function: Get a map from the factory
 Returns : L<Bio::Map::MapI>
 Args    : none

=cut

sub next_map{
   my ($self) = @_;
   my ($ready,$map) = (0,new Bio::Map::SimpleMap('-name'  => '',
						 '-units' => 'cM',
						 '-type'  => 'Genetic'));
   my @markers;
   my $runningDistance = 0;
   while( defined($_ = $self->_readline()) ) {
       if ( $ready || /^\s+Markers\s+Distance/ ) { 
	   unless ( $ready ) { $ready = 1; next }
       } else { next }

       last if ( /-{5,}/); # map terminator is ------- 
       s/ +/\t/;
       my ($number,$name,$distance) = split;
       $runningDistance += $distance;
       $runningDistance = '0.0' if $runningDistance == 0;
#       print "$_|$number-$name-$distance---------";
       my $pos = new Bio::Map::LinkagePosition (-order => $number,
						-map => $map,
						-value => $runningDistance
						);
       my $marker = new Bio::Map::Marker(-name=> $name,
					 -position => $pos,
					 );
       $marker->position($pos);
#       use Data::Dumper; print Dumper($marker); exit;
#       print $marker->position->value, "\n";
#       use Data::Dumper; print Dumper($pos);
#	 $map->add_element(new Bio::Map::Marker('-name'=> $name,
#						'-position' => $pos,
#						));   
   }
#   return undef if( ! $ready );
   return $map;
}

=head2 write_map

 Title   : write_tree
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
