# $Id$
#
# bioperl module for Bio::Structure::Atom
#
# Cared for by Kris Boulez <kris.boulez@algonomics.com>
#
# Copyright Kris Boulez
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Structure::Atom - Bioperl structure Object, describes an Atom

=head1 SYNOPSIS

=head1 DESCRIPTION

This object stores a Bio::Structure::Atom

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org             - General discussion
  http://bio.perl.org/MailList.html - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Kris Boulez

Email kris.boulez@algonomics.com

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Structure::Atom;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::Structure::Residue;
@ISA = qw(Bio::Root::Root);


=head2 new()

 Title   : new()
 Usage   : $struc = Bio::Structure::Atom->new( 
                                           -id  => 'human_id',
                                           );

 Function: Returns a new Bio::Structure::Atom object from basic 
	constructors. Probably most called from Bio::Structure::IO.
 Returns : a new Bio::Structure::Atom object

=cut


sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my($id, $x, $y, $z) =
        $self->_rearrange([qw(
			      ID
			      X
			      Y
			      Z
                              )],
                          @args);

    $id		&& $self->id($id);
    $x		&& $self->x($x);
    $y		&& $self->y($y);
    $z		&& $self->z($z);

    return $self;
}



=head2 x()

 Title   : x
 Usage   : $x = $atom->x($x);
 Function: Set/gets the X coordinate for an Atom
 Returns : The value for the X coordinate of the Atom (This is just a number,
 	it is expected to be in Angstrom, but no garantees)
 Args    : The X coordinate as a number

=cut

sub x {
	my ($self,$value) = @_;
	if( defined $value) {
		# do we want to check if $value contains really a number ?
		$self->{'x'} = $value;
   	}
	return $self->{'x'};
}


=head2 y()

 Title   : y
 Usage   : $y = $atom->y($y);
 Function: Set/gets the Y coordinate for an Atom
 Returns : The value for the Y coordinate of the Atom (This is just a number,
 	it is eypected to be in Angstrom, but no garantees)
 Args    : The Y coordinate as a number

=cut

sub y {
	my ($self,$value) = @_;
	if( defined $value) {
		# do we want to check if $value contains really a number ?
		$self->{'y'} = $value;
   	}
	return $self->{'y'};
}


=head2 z()

 Title   : z
 Usage   : $z = $atom->z($z);
 Function: Set/gets the Z coordinate for an Atom
 Returns : The value for the Z coordinate of the Atom (This is just a number,
 	it is ezpected to be in Angstrom, but no garantees)
 Args    : The Z coordinate as a number

=cut

sub z {
	my ($self,$value) = @_;
	if( defined $value) {
		# do we want to check if $value contains really a number ?
		$self->{'z'} = $value;
   	}
	return $self->{'z'};
}


=head2 xyz()

 Title   : xyz
 Usage   : ($x,$y,$z) = $atom->xyz;
 Function: Gets the XYZ coordinates for an Atom
 Returns : A list with the value for the XYZ coordinate of the Atom 
 Args    : 

=cut

sub xyz {
	my ($self) = @_;

	return ($self->x, $self->y, $self->z);
}


=head2 residue()

 Title   : residue
 Usage   : 
 Function:  No code here, all parent/child stuff via Entry
 Returns : 
 Args    : 

=cut

sub residue {
	my($self, $value) = @_;

	$self->throw("all parent/child stuff via Entry\n");
}


=head2 icode()

 Title   : icode
 Usage   : $icode = $atom->icode($icode)
 Function: Sets/gets the icode
 Returns : Returns the icode for this atom
 Args    : reference to an Atom

=cut

sub icode {
	my($self, $value) = @_;

	if (defined $value) {
		$self->{'icode'} = $value;
	}
	return $self->{'icode'};
}


=head2 occupancy()

 Title   : occupancy
 Usage   : $occupancy = $atom->occupancy($occupancy)
 Function: Sets/gets the occupancy
 Returns : Returns the occupancy for this atom
 Args    : reference to an Atom

=cut

sub occupancy {
	my($self, $value) = @_;

	if (defined $value) {
		$self->{'occupancy'} = $value;
	}
	return $self->{'occupancy'};
}


=head2 tempfactor()

 Title   : tempfactor
 Usage   : $tempfactor = $atom->tempfactor($tempfactor)
 Function: Sets/gets the tempfactor
 Returns : Returns the tempfactor for this atom
 Args    : reference to an Atom

=cut

sub tempfactor {
	my($self, $value) = @_;

	if (defined $value) {
		$self->{'tempfactor'} = $value;
	}
	return $self->{'tempfactor'};
}


=head2 element()

 Title   : element
 Usage   : $element = $atom->element($element)
 Function: Sets/gets the element
 Returns : Returns the element for this atom
 Args    : reference to an Atom

=cut

sub element {
	my($self, $value) = @_;

	if (defined $value) {
		$self->{'element'} = $value;
	}
	return $self->{'element'};
}


=head2 charge()

 Title   : charge
 Usage   : $charge = $atom->charge($charge)
 Function: Sets/gets the charge
 Returns : Returns the charge for this atom
 Args    : reference to an Atom

=cut

sub charge {
	my($self, $value) = @_;

	if (defined $value) {
		$self->{'charge'} = $value;
	}
	return $self->{'charge'};
}


=head2 id()

 Title   : id
 Usage   : $atom->id("CZ2")
 Function: Gets/sets the ID for this atom
 Returns : the ID
 Args    : the ID

=cut

sub id {
        my ($self, $value) = @_;;
        if (defined $value) {
	        $self->{'id'} = $value;
        }
        return $self->{'id'};
}

sub DESTROY {
	my $self =  shift;
	
	# dummy, nothing needs to be done here
}

#
# from here on only private methods
#

=head2 _remove_residue()

 Title   : _remove_residue
 Usage   : 
 Function: Removes the Residue this Atom is atttached to.
 Returns : 
 Args    : 

=cut

sub _remove_residue {
	my ($self) = shift;

	$self->throw("no code here at the moment\n");
}


=head2 _grandparent()

 Title   : _grandparent
 Usage   : 
 Function: get/set a symbolic reference to our grandparent
 Returns : 
 Args    : 

=cut

sub _grandparent {
	my($self,$symref) = @_;

	if (ref($symref)) {
		$self->throw("Thou shall only pass strings in here, no references $symref\n");
	}
	if (defined $symref) {
		$self->{'grandparent'} = $symref;
	}
	return $self->{'grandparent'};
}


1;
