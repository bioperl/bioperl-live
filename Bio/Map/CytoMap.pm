#
# BioPerl module for Bio::Map::CytoMap
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::CytoMap - A Bio::MapI compliant map implementation handling cytogenic bands 

=head1 SYNOPSIS

    use Bio::Map::CytoMap;
    my $map = Bio::Map::CytoMap->new(-name => 'human1',
				      -species => $human);
    foreach my $marker ( @markers ) { # get a list of markers somewhere
	$map->add_element($marker);
    }

=head1 DESCRIPTION

This is the simple implementation of cytogenetic maps based on
L<Bio::Map::MapI>.  It handles the essential storage of name, species,
type, and units as well as in memory representation of the elements of
a map.

For CytoMaps type is hard coded to be 'cytogeneticmap' and
units are set to '' but can be set to something else.

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

=head1 AUTHOR - Heikki Lehvaslaiho

Email heikki-at-bioperl-dot-org

=head1 CONTRIBUTORS

Jason Stajich      jason@bioperl.org
Lincoln Stein      lstein@cshl.org
Sendu Bala         bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Map::CytoMap;
use vars qw($MAPCOUNT);
use strict;


use base qw(Bio::Map::SimpleMap);
BEGIN { $MAPCOUNT = 1; }

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Map::CytoMap->new();
 Function: Builds a new Bio::Map::CytoMap object
 Returns : Bio::Map::CytoMap
 Args    : -name    => name of map (string)
           -species => species for this map (Bio::Species) [optional]
           -elements=> elements to initialize with
                       (arrayref of Bio::Map::MappableI objects) [optional]

           -uid     => Unique Id

=cut

sub new {
    my ($class, @args) = @_;
	
    my $self = $class->SUPER::new(@args);
	
    $self->{'_uid'} = $MAPCOUNT++;
    my ($uid) = $self->_rearrange([qw(UID)], @args);
    defined $uid && $self->unique_id($uid);
	
    return $self;
}

=head2 type

 Title   : type
 Usage   : my $type = $map->type
 Function: Get hard-coded Map type
 Returns : String coding Map type (always 'cyto')
 Args    : none

=cut

sub type {
   return 'cyto';
}

=head2 length

 Title   : length
 Usage   : my $length = $map->length();
 Function: Retrieves the length of the map,
 Returns : 0 since length is not calculatable for cytogenetic maps
 Args    : none

=cut

sub length {
   return 0;
}

1;
