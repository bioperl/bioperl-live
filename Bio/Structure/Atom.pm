#
# bioperl module for Bio::Structure::Atom
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

  #add synopsis here

=head1 DESCRIPTION

This object stores a Bio::Structure::Atom

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Kris Boulez

Email kris.boulez@algonomics.com

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Structure::Atom;
use strict;

use Bio::Structure::Residue;
use base qw(Bio::Root::Root);


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


=head2 serial()

 Title   : serial
 Usage   : $serial = $atom->serial($serial)
 Function: Sets/gets the serial number
 Returns : Returns the serial number for this atom
 Args    : reference to an Atom

=cut

sub serial {
	my($self, $value) = @_;

	if (defined $value) {
		$self->{'serial'} = $value;
	}
	return $self->{'serial'};
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


=head2 segID()

 Title   : segID
 Usage   : $segID = $atom->segID($segID)
 Function: Sets/gets the segID
 Returns : Returns the segID for this atom
 Args    : reference to an Atom

=cut

sub segID {
	my($self, $value) = @_;

	if (defined $value) {
		$self->{'segID'} = $value;
	}
	return $self->{'segID'};
}


=head2 pdb_atomname()

 Title   : pdb_atomname
 Usage   : $pdb_atomname = $atom->pdb_atomname($pdb_atomname)
 Function: Sets/gets the pdb_atomname (atomname used in the PDB file)
 Returns : Returns the pdb_atomname for this atom
 Args    : reference to an Atom

=cut

sub pdb_atomname {
	my($self, $value) = @_;

	if (defined $value) {
		$self->{'pdb_atomname'} = $value;
	}
	return $self->{'pdb_atomname'};
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


=head2 sigx()

 Title   : sigx
 Usage   : $sigx = $atom->sigx($sigx)
 Function: Sets/gets the sigx
 Returns : Returns the sigx for this atom
 Args    : reference to an Atom

=cut

sub sigx {
	my($self, $value) = @_;

	if (defined $value) {
		$self->{'sigx'} = $value;
	}
	return $self->{'sigx'};
}


=head2 sigy()

 Title   : sigy
 Usage   : $sigy = $atom->sigy($sigy)
 Function: Sets/gets the sigy
 Returns : Returns the sigy for this atom
 Args    : reference to an Atom

=cut

sub sigy {
	my($self, $value) = @_;

	if (defined $value) {
		$self->{'sigy'} = $value;
	}
	return $self->{'sigy'};
}


=head2 sigz()

 Title   : sigz
 Usage   : $sigz = $atom->sigz($sigz)
 Function: Sets/gets the sigz
 Returns : Returns the sigz for this atom
 Args    : reference to an Atom

=cut

sub sigz {
	my($self, $value) = @_;

	if (defined $value) {
		$self->{'sigz'} = $value;
	}
	return $self->{'sigz'};
}


=head2 sigocc()

 Title   : sigocc
 Usage   : $sigocc = $atom->sigocc($sigocc)
 Function: Sets/gets the sigocc
 Returns : Returns the sigocc for this atom
 Args    : reference to an Atom

=cut

sub sigocc {
	my($self, $value) = @_;

	if (defined $value) {
		$self->{'sigocc'} = $value;
	}
	return $self->{'sigocc'};
}


=head2 sigtemp()

 Title   : sigtemp
 Usage   : $sigtemp = $atom->sigtemp($sigtemp)
 Function: Sets/gets the sigtemp
 Returns : Returns the sigtemp for this atom
 Args    : reference to an Atom

=cut

sub sigtemp {
	my($self, $value) = @_;

	if (defined $value) {
		$self->{'sigtemp'} = $value;
	}
	return $self->{'sigtemp'};
}


=head2 aniso()

 Title   : aniso
 Usage   : $u12 = $atom->aniso("u12", $u12)
 Function: Sets/gets the anisotropic temperature factors
 Returns : Returns the requested factor for this atom
 Args    : reference to an Atom, name of the factor, value for the factor

=cut

sub aniso {
	my($self, $name, $value) = @_;

	if ( !defined $name) {
		$self->throw("You need to supply a name of the anisotropic temp factor you want to get");
	}
	if (defined $value) {
		$self->{$name} = $value;
	}
	return $self->{$name};
}

# placeholders 
sub u11 {
	my ($self, $name, $value) = @_;
	$self->aniso($name,$value);
}
sub u22 {
	my ($self, $name, $value) = @_;
	$self->aniso($name,$value);
}
sub u33 {
	my ($self, $name, $value) = @_;
	$self->aniso($name,$value);
}
sub u12 {
	my ($self, $name, $value) = @_;
	$self->aniso($name,$value);
}
sub u13 {
	my ($self, $name, $value) = @_;
	$self->aniso($name,$value);
}
sub u23 {
	my ($self, $name, $value) = @_;
	$self->aniso($name,$value);
}
sub sigu11 {
	my ($self, $name, $value) = @_;
	$self->aniso($name,$value);
}
sub sigu22 {
	my ($self, $name, $value) = @_;
	$self->aniso($name,$value);
}
sub sigu33 {
	my ($self, $name, $value) = @_;
	$self->aniso($name,$value);
}
sub sigu12 {
	my ($self, $name, $value) = @_;
	$self->aniso($name,$value);
}
sub sigu13 {
	my ($self, $name, $value) = @_;
	$self->aniso($name,$value);
}
sub sigu23 {
	my ($self, $name, $value) = @_;
	$self->aniso($name,$value);
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
