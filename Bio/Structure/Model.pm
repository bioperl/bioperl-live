#
# bioperl module for Bio::Structure::Model
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

Bio::Structure::Model - Bioperl structure Object, describes a Model

=head1 SYNOPSIS

  #add synopsis here

=head1 DESCRIPTION

This object stores a Bio::Structure::Chain

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

package Bio::Structure::Model;
use strict;

use Bio::Structure::Entry;
use Bio::Structure::Chain;
use base qw(Bio::Root::Root);


=head2 new()

 Title   : new()
 Usage   : $struc = Bio::Structure::Model->new( 
                                           -id  => 'human_id',
                                           );

 Function: Returns a new Bio::Structure::Model object from basic 
	constructors. Probably most called from Bio::Structure::IO.
 Returns : a new Bio::Structure::Model object

=cut



sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my($id, $chain, $residue ) =
        $self->_rearrange([qw(
			      ID
			      CHAIN
			      RESIDUE
                              )],
                          @args);

    $id      && $self->id($id);

    $chain && $self->throw("you have to add chain via an Entry object\n");

    $residue && $self->throw("you have to add residues via an Entry object\n");

    return $self;
}



=head2 chain()

 Title   : chain
 Usage   : 
 Function: will eventually allow parent/child navigation not via an Entry object
 Returns : 
 Args    : 

=cut

sub chain {
	my ($self,$value) = @_;

	$self->throw("go via an Entry object\n");
}


=head2 add_chain()

 Title   : add_chain
 Usage   : 
 Function:  will eventually allow parent/child navigation not via an Entry object
 Returns : 
 Args    : 

=cut

sub add_chain {
	my ($self,$value) = @_;

	$self->throw("go via an Entry object for now\n");
}

=head2 entry()

 Title   : entry
 Usage   : 
 Function:  will eventually allow parent/child navigation not via an Entry object
 Returns : 
 Args    : 

=cut

sub entry {
	my($self) = @_;

	$self->throw("Model::entry go via an Entry object please\n");
}


=head2 id()

 Title   : id
 Usage   : $model->id("model 5")
 Function: Gets/sets the ID for this model
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

=head2 residue()

 Title   : residue
 Usage   : 
 Function:  will eventually allow parent/child navigation not via an Entry object
 Returns : 
 Args    : 

=cut

sub residue {
	my ($self, @args) = @_;

	$self->throw("need to go via Entry object or learn symbolic refs\n");
}


=head2 add_residue()

 Title   : add_residue
 Usage   : 
 Function:  will eventually allow parent/child navigation not via an Entry object
 Returns : 
 Args    : 

=cut

sub add_residue {
	my ($self, @args) = @_;

	$self->throw("go via entry->add_residue(chain, residue)\n");
}



sub DESTROY {
	my $self = shift;

	# no specific DESTROY for now
}

#
# from here on only private methods
#

=head2 _remove_chains()

 Title   : _remove_chains
 Usage   : 
 Function: Removes the chains attached to a Model. Tells the chains they
 	don't belong to this Model any more
 Returns : 
 Args    : 

=cut

sub _remove_chains {
	my ($self) = shift;

	$self->throw("use Entry methods pleae\n");
}


=head2 _remove_entry()

 Title   : _remove_entry
 Usage   : 
 Function: Removes the Entry this Model is atttached to.
 Returns : 
 Args    : 

=cut

sub _remove_entry {
	my ($self) = shift;

	$self->throw("use a method based on an Entry object\n");
}


=head2 _create_default_chain()

 Title   : _create_default_chain
 Usage   : 
 Function: Creates a default Chain for this Model. Typical situation
 	in an X-ray structure where there is only one chain
 Returns : 
 Args    : 

=cut

sub _create_default_chain {
	my ($self) = shift;

	my $chain = Bio::Structure::Chain->new(-id => "default");
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
