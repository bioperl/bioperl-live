# $Id $
#
# bioperl module for Bio::Structure::Residue
#
# Cared for by Kris Boulez <kris.boulez@algonomics.com>
#
# Copyright Kris Boulez
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Structure::Residue - Bioperl structure Object, describes a Residue

=head1 SYNOPSIS

=head1 DESCRIPTION

This object stores a Bio::Structure::Residue

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

package Bio::Structure::Residue;
use vars qw(@ISA);
use strict;

use Bio::Root::RootI;
use Bio::Structure::Chain;
use Bio::Structure::Atom;
@ISA = qw(Bio::Root::RootI);


=head2 new()

 Title   : new()
 Usage   : $residue = Bio::Structure::Residue->new( 
                                           -id  => 'human_id',
                                           );

 Function: Returns a new Bio::Structure::Residue object from basic 
	constructors. Probably most called from Bio::Structure::IO.
 Returns : a new Bio::Structure::Residue object

=cut


sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my($id, $atom ) =
        $self->_rearrange([qw(
			      ID
			      ATOM
                              )],
                          @args);

    $id      && $self->id($id);

    $self->{'atom'} = [];

    # the 'smallest' (and only) item that can be added to a residue is an atom 

    $atom  && $self->atom($atom);

    return $self;
}



=head2 atom()

 Title   : atom 
 Usage   : @atoms  = $residue->atom($atom);
 Function: Connects an (or a list of) Atom objects to a Bio::Structure::Residue.
 	To add an atom (and keep the existing ones) use add_atom()
 	It returns a list of Atom objects.
 Returns : list of Bio::Structure::Atom objects
 Args    : One Atom or a reference to an array of Atom objects

=cut

sub atom {
	my ($self,$value) = @_;
	if( defined $value) {
                if( (ref($value) eq "ARRAY") ||
                      ($value->isa('Bio::Structure::Atom')) ) {
                        # remove existing ones, tell they've become orphan
		        $self->_remove_atoms;
			# add the new ones
			$self->add_atom($value);
                }
		else {
			$self->throw("Supplied a $value to atom , we want a Bio::Structure::Atom or a list of these\n");
		}
   	}
	return @{ $self->{'atom'} };
}


=head2 add_atom()

 Title   : add_atom
 Usage   : $atom->add_atom($atom);
 Function: Adds a (or a list of) Atom objects to a Bio::Structure::Atom.
 Returns : 
 Args    : One Atom or a reference to an array of Atom objects

=cut

sub add_atom {
	my($self,$value) = @_;
	if (defined $value) {
		if (ref($value) eq "ARRAY") {
    			# if the user passed in a reference to an array
			for my $a ( @{$value} ) {
				if( ! $a->isa('Bio::Structure::Atom') ) {
					$self->throw("$a is not a Atom\n");
				}
				push @{$self->{'atom'}}, $a;
				# tell the child who his parent is
				$a->residue($self);
			}
		}
		elsif ( $value->isa('Bio::Structure::Atom') ) { 
			push @{$self->{'atom'}}, $value;
			$value->residue($self);
		}
		else {
			$self->throw("Supplied a $value to add_atom, we want a Atom or list of Atoms\n");
		}
	}
	else {
		$self->warn("You need to supply an argument of type Atom to add_atom\n");
	}
}


=head2 chain()

 Title   : chain
 Usage   : $chain = $residue->chain($chain)
 Function: Sets the Chain this Residue belongs to
 Returns : Returns the Chain this Residue belongs to
 Args    : reference to a Chain

=cut

sub chain {
	my($self, $value) = @_;

	if (defined $value) {
		if (! $value->isa('Bio::Structure::Chain') ) {
			$self->throw("Need to supply a Bio::Structure::Chain
				to chain(), not a $value\n");
		}
		$self->{'chain'} = $value;
	}

	return $self->{'chain'};
}


=head2 id()

 Title   : id
 Usage   : $residue->id("TRP-35")
 Function: Gets/sets the ID for this residue
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


=head2 DESTROY()

 Title   : DESTROY
 Usage   : 
 Function: destructor ( get rid of circular references )
 Returns : 
 Args    : 

=cut

sub DESTROY {
	my $self = shift;

	# no need to explicitely call Atom destructor as it's the lowest level
	$self->{'atom'} = [];

}


#
# from here on only private methods
#

=head2 _remove_atoms()

 Title   : _remove_atoms
 Usage   : 
 Function: Removes the atoms attached to a Residue. Tells the atoms they
 	don't belong to this Residue any more
 Returns : 
 Args    : 

=cut

sub _remove_atoms {
	my ($self) = shift;

	for my $atom ( $self->atom ) {
		$atom->_remove_residue;
	}
	$self->{'atom'} = [];
}


=head2 _remove_chain()

 Title   : _remove_chain
 Usage   : 
 Function: Removes the Chain this Residue is atttached to.
 Returns : 
 Args    : 

=cut

sub _remove_chain {
	my ($self) = shift;

	$self->{'chain'} = undef;
}


1;
