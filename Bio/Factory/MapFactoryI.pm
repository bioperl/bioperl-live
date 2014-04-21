#
# BioPerl module for Bio::Factory::MapFactoryI
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

Bio::Factory::MapFactoryI - A Factory for getting markers

=head1 SYNOPSIS

    # get a Map Factory somehow likely from Bio::MapIO system

    while( my $map = $mapin->next_map ) {
	print "map name is ", $map->name, " length is ", 
	    $map->length, " ", $map->units, "\n";
	$mapout->write_map($map);
    }

=head1 DESCRIPTION

This interface describes the necessary minimum methods for getting
Maps from a data stream.  It also supports writing Map data back to a
stream.

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


package Bio::Factory::MapFactoryI;
use strict;

use base qw(Bio::Root::RootI);

=head2 next_map

 Title   : next_map
 Usage   : my $map = $factory->next_map;
 Function: Get a map from the factory
 Returns : L<Bio::Map::MapI>
 Args    : none

=cut

sub next_map{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 write_map

 Title   : write_map
 Usage   : $factory->write_map($map);
 Function: Write a map out through the factory
 Returns : none
 Args    : L<Bio::Map::MapI>

=cut

sub write_map{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

1;
