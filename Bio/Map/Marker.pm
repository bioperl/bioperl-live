# BioPerl module for Bio::Map::Marker
#
# Cared for by Chad Matsalla <bioinformatics1@dieselwurks.com>
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::Marker - An object representing a generic marker.

=head1 SYNOPSIS

$o_usat = new Bio::Map::Marker(-name=>'Chad Super Marker 2',
	-position => Bio::Map::PositionI-child);

=head1 DESCRIPTION

This object handles the notion of a generic marker. This marker will have
a name and a position.

This object is intended to be used by a marker parser like Mapmaker.pm and
then blessed into the proper type of marker (ie Microsatellite) by the
calling script.

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

Heikki Lehvaslaiho heikki@ebi.ac.uk
Lincoln Stein      lstein@cshl.org
Jason Stajich      jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::Marker;
use vars qw(@ISA);
use strict;
use Bio::Root::Root;

@ISA = qw(Bio::Root::Root Bio::Map::MarkerI);

=head2 new

 Title   : new
 Usage   : $o_marker = new Bio::Map::Marker( -name => 'Whizzy marker',
	                                     -position => $position);
 Function: Builds a new Bio::Map::Marker object
 Returns : Bio::Map::Marker
 Args    :
	-name    => name of this microsatellite (optional, string,
		default 'Unknown microsatellite')
	-position => position of this marker (optional,
		Bio::Map::PositionI-inherited object, no default)

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($name, $position) = $self->_rearrange([qw(NAME POSITION)], @args);
    if ($name) { $self->name($name); } else {$self->name('Unnamed marker'); }
    $position && $self->position($position);
    return bless $self;
}


=head2 position

 Title   : position
 Usage   : my $position = $mappable->position(); 
 Function: Get/Set the Bio::Map::PositionI for a mappable element
 Returns : Bio::Map::PositionI
 Args    : Bio::Map::PositionI

=cut

sub position {
    my ($self,$o_position) = @_;
    if ($o_position) {
	$self->{'_position'} = $o_position;
    }
    return $self->{'_position'};	
}

=head2 name($new_name)

 Title   : name($new_name)
 Usage   : my $name = $o_usat->name($new_name) _or_
	my $name = $o_usat->name()
 Function: Get/Set the name for this Microsatellite
 Returns : A scalar representing the current name of this Microsatellite
 Args    : If provided, the current name of this Microsatellite
	will be set to $new_name.

=cut

sub name {
    my ($self,$name) = @_;
    if ($name) {
	$self->{'_name'} = $name;
    }
    return $self->{'_name'};	
}

1;


