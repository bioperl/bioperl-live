#
# bioperl module for Bio::Structure::Residue
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

Bio::Structure::Residue - Bioperl structure Object, describes a Residue

=head1 SYNOPSIS

  #add synopsis here

=head1 DESCRIPTION

This object stores a Bio::Structure::Residue

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


package Bio::Structure::Residue;
use strict;

use Bio::Structure::Chain;
use Bio::Structure::Atom;
use base qw(Bio::Root::Root);


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

    $atom && $self->throw("add atoms via an Entry object entry->add_atom(residue,atom)\n");

    return $self;
}



=head2 atom()

 Title   : atom 
 Usage   : 
 Function:  nothing useful untill I get symbolic references to do what I want
 Returns : 
 Args    : 

=cut

sub atom {
	my ($self,$value) = @_;

	$self->throw("no code down here, go see an Entry object nearby\n");
}


=head2 add_atom()

 Title   : add_atom
 Usage   : 
 Function:  nothing useful untill I get symbolic references to do what I want
 Returns : 
 Args    : 

=cut

sub add_atom {
	my($self,$value) = @_;

	$self->throw("nothing here, use a method on an Entry object\n");
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

	$self->throw("use an Entry based method please\n");
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

	# no specific destruction for now
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

	$self->throw("no code here\n");
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
