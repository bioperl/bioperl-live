# $Id $
#
# bioperl module for Bio::Structure::Model
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

=head1 DESCRIPTION

This object stores a Bio::Structure::Chain

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

package Bio::Structure::Model;
use vars qw(@ISA);
use strict;

use Bio::Root::RootI;
use Bio::Structure::Entry;
use Bio::Structure::Chain;
@ISA = qw(Bio::Root::RootI);


=head2 new()

 Title   : new()
 Usage   : $struc = Bio::Structure::Model->new( 
                                           -id  => 'human_id',
                                           -accession_number => 'AL000012',
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

    $self->{'chain'} = [];

    # the 'smallest' item that can be added to a model is a chain. If
    # chain is not present, a default is taken

    if($chain) {
	    $self->chain($chain);
    } else {  # create first chain
	    $chain = Bio::Structure::Chain->new();
	    $self->chain($chain);
    }

    $residue  && $self->residue($residue);

    return $self;
}



=head2 chain()

 Title   : chain
 Usage   : @chains  = $model->chain($chain);
 Function: Connects a (or a list of) Chain objects to a Bio::Structure::Model.
 	To add a Chain (and keep the existing ones) use add_chain()
 	It returns a list of Chain objects.
 Returns : list of Bio::Structure::Chain objects
 Args    : One Chain or a reference to an array of Chain objects

=cut

sub chain {
	my ($self,$value) = @_;
	if( defined $value) {
                if( (ref($value) eq "ARRAY") ||
                      ($value->isa('Bio::Structure::Chain')) ) {
                        # remove existing ones, tell they've become orphan
                        $self->_remove_chains;
                        # add the new ones
                        $self->add_chain($value);
		}
		else {
			$self->throw("Supplied a $value to chain , we want a Bio::Structure::Chain or a list of these\n");
		}
   	}
	return @{ $self->{'chain'} };
}


=head2 add_chain()

 Title   : add_chain
 Usage   : $model->add_chain($chain);
 Function: Adds a (or a list of) Chain objects to a Bio::Structure::Model.
 Returns : 
 Args    : One Chain or a reference to an array of Chain objects

=cut

sub add_chain {
	my($self,$value) = @_;
	if (defined $value) {
		if (ref($value) eq "ARRAY") {
    			# if the user passed in a reference to an array
			for my $c ( @{$value} ) {
				if( ! $c->isa('Bio::Structure::Chain') ) {
					$self->throw("$c is not a Model\n");
				}
				push @{$self->{'chain'}}, $c;
				# tell the child who his parent is
				$c->model($self);
			}
		}
		elsif ( $value->isa('Bio::Structure::Chain') ) { 
			push @{$self->{'chain'}}, $value;
			$value->model($self);
		}
		else {
			$self->throw("Supplied a $value to add_chain, we want a Chain or list of Chains\n");
		}
	}
	else {
		$self->warn("You need to supply an argument of type Chain to add_chain\n");
	}
}

=head2 entry()

 Title   : entry
 Usage   : $entry = $model->entry($entry)
 Function: Sets the Entry this Model belongs to
 Returns : Returns the Entry this Model belongs to
 Args    : reference to an Entry

=cut

sub entry {
	my($self, $value) = @_;

	if (defined $value) {
		if (! $value->isa('Bio::Structure::Entry') ) {
			$self->throw("Need to supply a Bio::Structure::Entry
				to entry(), not a $value\n");
		}
		$self->{'entry'} = $value;
	}

	return $self->{'entry'};
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

	for my $chain ( $self->chain ) {
		$chain->_remove_model;
	}
	$self->{'chain'} = [];
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

	$self->{'entry'} = undef;
}

1;
