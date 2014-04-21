#
# bioperl module for Bio::Structure::Chain
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

Bio::Structure::Chain - Bioperl structure Object, describes a chain

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

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Structure::Chain;
use strict;

use Bio::Structure::Entry;
use Bio::Structure::Model;
use base qw(Bio::Root::Root);


=head2 new()

 Title   : new()
 Usage   : $struc = Bio::Structure::Chain->new( 
                           -id  => 'human_id',
                           -accession_number => 'AL000012',
                                           );

 Function: Returns a new Bio::Structure::Chain object from basic 
	        constructors. Usually called from Bio::Structure::IO.
 Returns : a new Bio::Structure::Chain object

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my($id, $residue ) =
        $self->_rearrange([qw(
			      ID
			      RESIDUE
                              )],
                          @args);

    $id      && $self->id($id);
    $self->{'residue'} = [];
    # the 'smallest' item that can be added to a chain is a residue. 
    $residue && $self->throw("use a method based on an Entry object for now");
    return $self;
}



=head2 residue()

 Title   : residue 
 Usage   : 
 Function:  nothing useful until I get symbolic references to do what I want
 Returns : 
 Args    : 

=cut

sub residue {
	my ($self,$value) = @_;

	$self->throw("use a method on an Entry object to do what you want");
}


=head2 add_residue()

 Title   : add_residue
 Usage   : 
 Function: nothing useful until I get symbolic references to do what I want
 Returns : 
 Args    : 

=cut

sub add_residue {
	my($self,$value) = @_;

	$self->throw("you want entry->add_residue(chain, residue)\n");
}

=head2 model()

 Title   : model
 Usage   : 
 Function: nothing useful until I get symbolic references to do what I want
 Returns : 
 Args    : 

=cut

sub model {
	my($self, $value) = @_;

	$self->throw("go via a Entry object please\n");
}


=head2 id()

 Title   : id
 Usage   : $chain->id("chain B")
 Function: Gets/sets the ID for this chain
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
	my $self = shift;
	
	# no specific destruction for now
}


#
# from here on only private methods
#

=head2 _remove_residues()

 Title   : _remove_residues
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub _remove_residues {
	my ($self) = shift;

	$self->throw("nothing usefull in here, go see Entry\n");
}


=head2 _remove_model()

 Title   : _remove_model
 Usage   : 
 Function: Removes the Model this Chain is atttached to.
 Returns : 
 Args    : 

=cut

sub _remove_model {
	my ($self) = shift;

	$self->throw("go see an Entry object, nothing here\n");
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
