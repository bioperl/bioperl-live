# $Id $
#
# bioperl module for Bio::Structure::Chain
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

package Bio::Structure::Chain;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::Structure::Entry;
use Bio::Structure::Model;
@ISA = qw(Bio::Root::Root);


=head2 new()

 Title   : new()
 Usage   : $struc = Bio::Structure::Chain->new( 
                                           -id  => 'human_id',
                                           -accession_number => 'AL000012',
                                           );

 Function: Returns a new Bio::Structure::Chain object from basic 
	constructors. Probably most called from Bio::Structure::IO.
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

    $residue  && $self->residue($residue);

    return $self;
}



=head2 residue()

 Title   : residue 
 Usage   : @residues  = $chain->residue($residue);
 Function: Connects a (or a list of) Residue objects to a Bio::Structure::Chain.
 	To add a chain (and keep the existing ones) use add_residue
 	It returns a list of Residue objects.
 Returns : list of Bio::Structure::Residue objects
 Args    : One Residue or a reference to an array of Residue objects

=cut

sub residue {
	my ($self,$value) = @_;
	if( defined $value) {
                if( (ref($value) eq "ARRAY") ||
                      ($value->isa('Bio::Structure::Residue')) ) {
                        # remove existing ones, tell they've become orphan
		        $self->_remove_residues;
			# add the new ones
			$self->add_residue($value);
                }
		else {
			$self->throw("Supplied a $value to residue , we want a Bio::Structure::Residue or a list of these\n");
		}
   	}
	return @{ $self->{'residue'} };
}


=head2 add_residue()

 Title   : add_residue
 Usage   : $chain->add_residue($residue);
 Function: Adds a (or a list of) Residue objects to a Bio::Structure::Chain.
 Returns : 
 Args    : One Residue or a reference to an array of Residue objects

=cut

sub add_residue {
	my($self,$value) = @_;
	if (defined $value) {
		if (ref($value) eq "ARRAY") {
    			# if the user passed in a reference to an array
			for my $r ( @{$value} ) {
				if( ! $r->isa('Bio::Structure::Residue') ) {
					$self->throw("$r is not a Residue\n");
				}
				push @{$self->{'residue'}}, $r;
				# tell the child who his parent is
				$r->chain($self);
			}
		}
		elsif ( $value->isa('Bio::Structure::Residue') ) { 
			push @{$self->{'residue'}}, $value;
			$value->chain($self);
		}
		else {
			$self->throw("Supplied a $value to add_chain, we want a Residue or list of Residues\n");
		}
	}
	else {
		$self->warn("You need to supply an argument of type Residue to add_residue\n");
	}
}

=head2 model()

 Title   : model
 Usage   : $model = $chain->model($model)
 Function: Sets the Model this Chain belongs to
 Returns : Returns the Model this Chain belongs to
 Args    : reference to a Model

=cut

sub model {
	my($self, $value) = @_;

	if (defined $value) {
		if (! $value->isa('Bio::Structure::Model') ) {
			$self->throw("Need to supply a Bio::Structure::Model
				to model(), not a $value\n");
		}
		$self->{'model'} = $value;
	}

	return $self->{'model'};
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
	
	# get rid of circular references
	for my $res ($self->residue) {
		next unless (defined $res);
		$res->DESTROY;
	}
	$self->{'residue'} = [];
}


#
# from here on only private methods
#

=head2 _remove_residues()

 Title   : _remove_residues
 Usage   : 
 Function: Removes the residues attached to an Chain. Tells the residues they
 	don't belong to this Chain any more
 Returns : 
 Args    : 

=cut

sub _remove_residues {
	my ($self) = shift;

	for my $residue ( $self->residue ) {
		$residue->_remove_chain;
	}
	$self->{'residue'} = [];
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

	$self->{'model'} = undef;
}


1;
