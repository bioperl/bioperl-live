# $Id $
#
# bioperl module for Bio::Structure::Entry
#
# Cared for by Kris Boulez <kris.boulez@algonomics.com>
#
# Copyright Kris Boulez
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Structure::Entry - Bioperl structure Object, describes the whole entry

=head1 SYNOPSIS

=head1 DESCRIPTION

This object stores a whole Bio::Structure entry. It can consist of one or
more models (Bio::Structure::Model), which in turn consist of one or more
chains (Bio::Structure::Chain). A chain is composed of residues 
(Bio::Structure::Residue) and a residue consists of atoms (Bio::Structure::Atom)
If no specific model or chain is chosen, the first one is choosen.

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

package Bio::Structure::Entry;
use vars qw(@ISA);
use strict;

use Bio::Root::RootI;
use Bio::StructureI;
use Bio::Structure::Model;
use Bio::Structure::Chain;
@ISA = qw(Bio::Root::RootI Bio::StructureI);


=head2 new()

 Title   : new()
 Usage   : $struc = Bio::Structure::Entry->new( 
                                           -id  => 'human_id',
                                           -accession_number => 'AL000012',
                                           );

 Function: Returns a new Bio::Structure::Entry object from basic 
	constructors. Probably most called from Bio::Structure::IO.
 Returns : a new Bio::Structure::Model object

=cut




sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my($id, $model, $chain, $residue ) =
        $self->_rearrange([qw(
			      ID
			      MODEL
			      CHAIN
			      RESIDUE
                              )],
                          @args);

    $id      && $self->id($id);

    $self->{'model'} = [];

    # the 'smallest' item that can be added to an Entry is a residue. If
    # neither model nor chain are present, a default is taken

    if($model) {
	    $self->model($model);
    } else {  # create first model
	    $model = Bio::Structure::Model->new();
	    $self->model($model);
    }

    if($chain) {
		for my $m ($self->model) { # add this chain on all models
		    $m->chain($chain);
    		}
    }

    $residue  && $self->residue($residue);

    return $self;
}



=head2 model()

 Title   : model
 Usage   : @models  = $structure->model($model);
 Function: Connects a (or a list of) Model objects to a Bio::Structure::Entry.
 	To add a Model (and keep the existing ones) use add_model()
 	It returns a list of Model objects.
 Returns : list of Bio::Structure::Model objects
 Args    : One Model or a reference to an array of Model objects

=cut

sub model {
	my ($self,$value) = @_;
	if( defined $value) {
		if( (ref($value) eq "ARRAY") ||
		      ($value->isa('Bio::Structure::Model')) ) {
			# remove existing ones, tell they've become orphan
			$self->_remove_models;
			# add the new ones
			$self->add_model($value);
		}
		else {
			$self->throw("Supplied a $value to model, we want a Bio::Structure::Model or a list of these\n");
		}
   	}
	return @{ $self->{'model'} };
}



=head2 add_model()

 Title   : add_model
 Usage   : $structure->add_model($model);
 Function: Adds a (or a list of) Model objects to a Bio::Structure::Entry.
 Returns : 
 Args    : One Model or a reference to an array of Model objects

=cut

sub add_model {
	my($self,$value) = @_;
	if (defined $value) {
		if (ref($value) eq "ARRAY") {
    			# if the user passed in a reference to an array
			for my $m ( @{$value} ) {
				if( ! $m->isa('Bio::Structure::Model') ) {
					$self->throw("$m is not a Model\n");
				}
				push @{$self->{'model'}}, $m;
				# tell the child who his parent is
				$m->entry($self);
			}
		}
		elsif ( $value->isa('Bio::Structure::Model') ) { 
			push @{$self->{'model'}}, $value;
			$value->entry($self);
		}
		else {
			$self->throw("Supplied a $value to add_model, we want a Model or list of Models\n");
		}
	}
	else {
		$self->warn("You need to supply an argument of type Model to add_model\n");
	}
}


=head2 id()

 Title   : id
 Usage   : $entry->id("identity");
 Function: Gets/sets the ID 
 Returns : 
 Args    : 

=cut

sub id {
	my ($self, $value) = @_;
	if (defined $value) {
		$self->{'id'} = $value;
	}
	return $self->{'id'};
}


sub DESTROY {
	my $self = shift;

	for my $model ($self->model) {
		next unless (defined $model);
		$model->DESTROY;
	}
	$self->{'model'} = [];
}

#
# from here on only private methods
#

=head2 _remove_models()

 Title   : _remove_models
 Usage   : 
 Function: Removes the models attached to an Entry. Tells the models they
 	don't belong to this Entry any more
 Returns : 
 Args    : 

=cut

sub _remove_models {
	my ($self) = shift;

	for my $model ( $self->model ) {
		$model->_remove_entry;
	}
	$self->{'model'} = [];
}



1;
