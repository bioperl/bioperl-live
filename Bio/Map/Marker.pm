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
			       -position => [ [$map, $position] ] );

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
use Bio::Map::MarkerI;

@ISA = qw(Bio::Root::Root Bio::Map::MarkerI);

=head2 new

 Title   : new
 Usage   : $o_marker = new Bio::Map::Marker( -name => 'Whizzy marker',
	                                     -position => $position);
 Function: Builds a new Bio::Map::Marker object
 Returns : Bio::Map::Marker
 Args    :
           -name    => name of this microsatellite 
                       [optional], string,default 'Unknown'
          
           -positions => map position for this marker, [optional]
                Bio::Map::PositionI-inherited obj, no default)

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($name, $position) = $self->_rearrange([qw(NAME POSITION)], @args);
    if ($name) { $self->name($name); } else {$self->name('Unnamed marker'); }
    $position && $self->position($position); 
    return $self;
}

=head2 position

 Title   : position
 Usage   : my $position = $mappable->position($map); OR
           $mappable->position($map,$position); OR
 Function: Get/Set the Bio::Map::PositionI for a mappable element
           in a specific Map
 Returns : Bio::Map::PositionI
 Args    : $map  - Bio::Map::MapI # Map we are talking about
           $position - Bio::Map::PositionI # Position we want to set

=cut

sub position {
    my ($self,$o_position) = @_;
    
    if ($o_position) {
	if( !ref($o_position) ||
	    ! $o_position->isa('Bio::Map::PositionI') ) {
	    $self->warn("Must specify a Bio::Map::PositionI when trying to set the position not ". ref($o_position));
	    return;
	}
	$self->{'_position'} = $o_position;
    }
    return $self->{'_position'};
}

=head2 name($new_name)

 Title   : name($new_name)
 Usage   : $o_usat->name($new_name) _or_
	   my $name = $o_usat->name()
 Function: Get/Set the name for this Microsatellite
 Returns : A scalar representing the current name of this Microsatellite
 Args    : If provided, the current name of this Microsatellite
	   will be set to $new_name.

=cut

sub name {
    my ($self,$name) = @_;
    my $last = $self->{'_name'};
    if ($name) {
	$self->{'_name'} = $name;
    }
    return $last;	
}

=head2 equals

 Title   : equals
 Usage   : if( $mappable->equals($mapable2)) ...
 Function: Test if a position is equal to another position
 Returns : boolean
 Args    : Bio::Map::MappableI

=cut

sub equals {
    my ($self,$compare) = @_;
    return 0 unless defined $compare;
    if( $compare->isa('Bio::Map::MappableI') ){
	return ($self->position->equals($compare->position));
    } elsif( $compare->isa('Bio::Map::PositionI') ) {
	return ($self->position->equals($compare));
    } else { 
	$self->warn("Can only run equals with Bio::Map::MappableI or Bio::Map::PositionI"); 
    }
    return 0;
}

=head2 less_than

 Title   : less_than
 Usage   : if( $mappable->less_than($m2) ) ...
 Function: Tests if a position is less than another position
 Returns : boolean
 Args    : Bio::Map::MappableI

=cut

sub less_than {
    my ($self,$compare) = @_;
    return 0 unless defined $compare;
    if( $compare->isa('Bio::Map::MappableI') ){
	return ($self->position->less_than($compare->position));
    } elsif( $compare->isa('Bio::Map::PositionI') ) {
	return ($self->position->less_than($compare));
    } else { 
	$self->warn("Can only run less_than with Bio::Map::MappableI or Bio::Map::PositionI"); 
    }
    return 0;
}

=head2 greater_than

 Title   : greater_than
 Usage   : if( $mappable->greater_than($m2) ) ...
 Function: Tests if position is greater than another position
 Returns : boolean
 Args    : Bio::Map::MappableI

=cut

sub greater_than {
    my ($self,$compare) = @_;
    return 0 unless defined $compare;
    if( $compare->isa('Bio::Map::MappableI') ){
	return ($self->position->greater_than($compare->position));
    } elsif( $compare->isa('Bio::Map::PositionI') ) {
	return ($self->position->greater_than($compare));
    } else { 
	$self->warn("Can only run greater_than with Bio::Map::MappableI or Bio::Map::PositionI"); 
    }
    return 0;
}

1;


