# $Id$
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

  #add synopsis here

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

use Bio::Root::Root;
use Bio::StructureI;
use Bio::Structure::Model;
use Bio::Structure::Chain;
use Bio::Annotation::Collection;
use Tie::RefHash;

@ISA = qw(Bio::Root::Root Bio::StructureI);


=head2 new()

 Title   : new()
 Usage   : $struc = Bio::Structure::Entry->new( 
                                           -id  => 'structure_id',
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

    # where to store parent->child relations (1 -> 1..n)
    #  value to this hash will be an array ref
    #  by using Tie::RefHash we can store references in this hash
    $self->{'p_c'} = ();
    tie %{ $self->{'p_c'} } , "Tie::RefHash";
    
    # where to store child->parent relations (1 -> 1)
    $self->{'c_p'} = ();
    tie %{ $self->{'c_p'} } , "Tie::RefHash";

    $id      && $self->id($id);

    $self->{'model'} = [];
    $model   && $self->model($model);


    if($chain) {
	if ( ! defined($self->model) ) { # no model yet, create default one
		$self->_create_default_model;
	}
	for my $m ($self->model) { # add this chain on all models
		$m->chain($chain);
    	}
    }

    $residue  && $self->residue($residue);

    # taken from Bio::Seq (or should we just inherit Bio::Seq and override some methods)
    my $ann = Bio::Annotation::Collection->new;
    $self->annotation($ann);

    return $self;
}


=head2 model()

 Title   : model
 Function: Connects a (or a list of) Model objects to a Bio::Structure::Entry.
 	To add a Model (and keep the existing ones) use add_model()
 	It returns a list of Model objects.
 Returns : list of Bio::Structure::Model objects
 Args    : One Model or a reference to an array of Model objects

=cut

sub model {
	my ($self, $model) = @_;
	
	if( defined $model) {
		if( (ref($model) eq "ARRAY") ||
		      ($model->isa('Bio::Structure::Model')) ) {
			# remove existing ones, tell they've become orphan
			my @obj = $self->model;
			if (@obj) {
				for my $m (@obj) {
					$self->_remove_from_graph($m);
					$self->{'model'} = [];
				}
			}
			# add the new ones
			$self->add_model($self,$model);
		}
		else {
			$self->throw("Supplied a $model to model, we want a Bio::Structure::Model or a list of these\n");
		}
   	}
	# give back list of models via general get method
	$self->get_models($self);
}



=head2 add_model()

 Title   : add_model
 Usage   : $structure->add_model($model);
 Function: Adds a (or a list of) Model objects to a Bio::Structure::Entry.
 Returns : 
 Args    : One Model or a reference to an array of Model objects

=cut

sub add_model {
	my($self,$entry,$model) = @_;

	# if only one argument and it's a model, change evrything one place
	# this is for people calling $entry->add_model($model);
	if ( !defined $model && ref($entry) =~ /^Bio::Structure::Model/) {
		$model = $entry;
		$entry = $self;
	}
	# $self and $entry are the same here, but it's used for uniformicity
	if ( !defined($entry) || ref($entry) !~ /^Bio::Structure::Entry/) {
		$self->throw("first argument to add_model needs to be a Bio::Structure::Entry object\n");
	}
	if (defined $model) {
		if (ref($model) eq "ARRAY") {
    			# if the user passed in a reference to an array
			for my $m ( @{$model} ) {
				if( ! $m->isa('Bio::Structure::Model') ) {
					$self->throw("$m is not a Model\n");
				}
				if ( $self->_parent($m) ) {
					$self->throw("$m already assigned to a parent\n");
				}
				push @{$self->{'model'}}, $m;
				# create a stringified version of our ref
				# not used untill we get symbolic ref working
				#my $str_ref = "$self";
				#$m->_grandparent($str_ref);
			}
		}
		elsif ( $model->isa('Bio::Structure::Model') ) { 
			if ( $self->_parent($model) ) { # already assigned to a parent
				$self->throw("$model already assigned\n");
			}
			push @{$self->{'model'}}, $model;
			# create a stringified version of our ref
			#my $str_ref = "$self";
			#$model->_grandparent($str_ref);
		}
		else {
			$self->throw("Supplied a $model to add_model, we want a Model or list of Models\n");
		}
	}

	my $array_ref = $self->{'model'};
	return $array_ref ? @{$array_ref} : ();
}


=head2 get_models()

 Title   : get_models
 Usage   : $structure->get_models($structure);
 Function: general get method for models attached to an Entry
 Returns : a list of models attached to this entry
 Args    : an Entry

=cut

sub get_models {
	my ($self, $entry) = @_;

	# self and entry can be the same
	if ( !defined $entry) {
		$entry = $self;
	}
	# pass through to add_model
	$self->add_model($entry);
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


=head2 chain()

 Title   : chain
 Usage   : @chains  = $structure->chain($chain);
 Function: Connects a (or a list of) Chain objects to a Bio::Structure::Entry.
 Returns : list of Bio::Structure::Residue objects
 Args    : One Residue or a reference to an array of Residue objects

=cut

sub chain {
	my ($self, $chain) = @_;

	if ( ! $self->model ) {
		$self->_create_default_model;
	}
	my @models = $self->model;
	my $first_model = $models[0];

	if ( defined $chain) {
	
		if( (ref($chain) eq "ARRAY") ||
		      ($chain->isa('Bio::Structure::Chain')) ) {
			# remove existing ones, tell they've become orphan
			my @obj = $self->get_chains($first_model);
			if (@obj) {
				for my $c (@obj) {
					$self->_remove_from_graph($c);
				}
			}
			# add the new ones
			$self->add_chain($first_model,$chain);
		}
		else {
			$self->throw("Supplied a $chain to chain, we want a Bio::Structure::Chain or a list of these\n");
		}
	}
	$self->get_chains($first_model);
}


=head2 add_chain()

 Title   : add_chain
 Usage   : @chains  = $structure->add_chain($add_chain);
 Function: Adds a (or a list of) Chain objects to a Bio::Structure::Entry.
 Returns : 
 Args    : 

=cut

sub add_chain {
	my($self, $model, $chain) = @_;

	if (ref($model) !~ /^Bio::Structure::Model/) {
		$self->throw("add_chain: first argument needs to be a Model object ($model)\n");
	}
	if (defined $chain) {
		if (ref($chain) eq "ARRAY") {
    			# if the user passed in a reference to an array
			for my $c ( @{$chain} ) {
				if( ! $c->isa('Bio::Structure::Chain') ) {
					$self->throw("$c is not a Chain\n");
				}
				if ( $self->_parent($c) ) {
					$self->throw("$c already assigned to a parent\n");
				}
				$self->_parent($c, $model);
				$self->_child($model, $c);
				# stringify $self ref
				#my $str_ref = "$self";
				#$c->_grandparent($str_ref);
			}
		}
		elsif ( $chain->isa('Bio::Structure::Chain') ) { 
			if ( $self->_parent($chain) ) { # already assigned to parent
				$self->throw("$chain already assigned to a parent\n");
			}
			$self->_parent($chain,$model);
			$self->_child($model, $chain);
			# stringify $self ref
			#my $str_ref = "$self";
			#$chain->_grandparent($str_ref);
		}
		else {
			$self->throw("Supplied a $chain to add_chain, we want a Chain or list of Chains\n");
		}
	}
	my $array_ref = $self->_child($model);
	return $array_ref ? @{$array_ref} : ();
}


=head2 get_chains()

 Title   : get_chains
 Usage   : $entry->get_chains($model);
 Function: general get method for chains attached to a Model
 Returns : a list of chains attached to this model
 Args    : a Model

=cut

sub get_chains {
	my ($self, $model) = @_;

	if (! defined $model) {
		$model = ($self->get_models)[0];
	}
	# pass through to add_chain
	$self->add_chain($model);
}


=head2 residue()

 Title   : residue
 Usage   : @residues  = $structure->residue($residue);
 Function: Connects a (or a list of) Residue objects to a Bio::Structure::Entry.
 Returns : list of Bio::Structure::Residue objects
 Args    : One Residue or a reference to an array of Residue objects

=cut

sub residue {
	my ($self, $residue) = @_;

	if ( ! $self->model ) {
		my $m = $self->_create_default_model;
		$self->add_model($self,$m);
	}
	my @models = $self->model;
	my $first_model = $models[0];
	
	if ( ! $self->get_chains($first_model) ) {
		my $c = $self->_create_default_chain;
		$self->add_chain($first_model, $c);
	}
	my @chains = $self->get_chains($first_model);
	my $first_chain = $chains[0];

	if( defined $residue) {
		if( (ref($residue) eq "ARRAY") ||
		      ($residue->isa('Bio::Structure::Residue')) ) {
			# remove existing ones, tell they've become orphan
			my @obj = $self->get_residues($first_chain);
			if (@obj) {
				for my $r (@obj) {
					$self->_remove_from_graph($r);
				}
			}
			# add the new ones
			$self->add_residue($first_chain,$residue);
		}
		else {
			$self->throw("Supplied a $residue to residue, we want a Bio::Structure::Residue or a list of these\n");
		}
   	}
	$self->get_residues($first_chain);
}


=head2 add_residue()

 Title   : add_residue
 Usage   : @residues  = $structure->add_residue($residue);
 Function: Adds a (or a list of) Residue objects to a Bio::Structure::Entry.
 Returns : list of Bio::Structure::Residue objects
 Args    : One Residue or a reference to an array of Residue objects

=cut

sub add_residue {
	my($self,$chain,$residue) = @_;

	if (ref($chain) !~ /^Bio::Structure::Chain/) {
		$self->throw("add_residue: first argument needs to be a Chain object\n");
	}
	if (defined $residue) {
		if (ref($residue) eq "ARRAY") {
    			# if the user passed in a reference to an array
			for my $r ( @{$residue} ) {
				if( ! $r->isa('Bio::Structure::Residue') ) {
					$self->throw("$r is not a Residue\n");
				}
				if ( $self->_parent($r) ) {
					$self->throw("$r already belongs to a parent\n");
				}
				$self->_parent($r, $chain);
				$self->_child($chain, $r);
				# stringify
				my $str_ref = "$self";
				$r->_grandparent($str_ref);
			}
		}
		elsif ( $residue->isa('Bio::Structure::Residue') ) { 
			if ( $self->_parent($residue) ) {
				$self->throw("$residue already belongs to a parent\n");
			}
			$self->_parent($residue, $chain);
			$self->_child($chain, $residue);
			# stringify
			my $str_ref = "$self";
			$residue->_grandparent($str_ref);
		}
		else {
			$self->throw("Supplied a $residue to add_residue, we want a Residue or list of Residues\n");
		}
	}
	my $array_ref = $self->_child($chain);
	return $array_ref ? @{$array_ref} : ();
}


=head2 get_residues()

 Title   : get_residues
 Usage   : $structure->get_residues($chain);
 Function: general get method for residues attached to a Chain
 Returns : a list of residues attached to this chain
 Args    : a chain

=cut

sub get_residues {
	my ($self, $chain) = @_;

	if ( !defined $chain) {
		$self->throw("get_residues needs a Chain as argument");
	}
	# pass through to add_residue
	$self->add_residue($chain);
}


=head2 add_atom()

 Title   : add_atom
 Usage   : @atoms  = $structure->add_atom($residue,$atom);
 Function: Adds a (or a list of) Atom objects to a Bio::Structure::Residue.
 Returns : list of Bio::Structure::Atom objects
 Args    : a residue and an atom

=cut

sub add_atom {
	my($self,$residue,$atom) = @_;

	if (ref($residue) !~ /^Bio::Structure::Residue/) {
		$self->throw("add_atom: first argument needs to be a Residue object\n");
	}
	if (defined $atom) {
		if (ref($atom) eq "ARRAY") {
    			# if the user passed in a reference to an array
			for my $a ( @{$atom} ) {
				if( ! $a->isa('Bio::Structure::Atom') ) {
					$self->throw("$a is not an Atom\n");
				}
				if ( $self->_parent($a) ) {
					$self->throw("$a already belongs to a parent\n");
				}
				$self->_parent($a, $residue);
				$self->_child($residue, $a);
				# stringify
				#my $str_ref = "$self";
				#$r->_grandparent($str_ref);
			}
		}
		#elsif ( $atom->isa('Bio::Structure::Atom') ) { 
		elsif ( ref($atom) =~ /^Bio::Structure::Atom/ ) { 
			if ( $self->_parent($atom) ) {
				$self->throw("$atom already belongs to a parent\n");
			}
			$self->_parent($atom, $residue);
			$self->_child($residue, $atom);
			# stringify
			#my $str_ref = "$self";
			#$atom->_grandparent($str_ref);
		}
	}
	my $array_ref = $self->_child($residue);
	return $array_ref ? @{$array_ref} : ();
}


=head2 get_atoms()

 Title   : get_atoms
 Usage   : $structure->get_atoms($residue);
 Function: general get method for atoms attached to a Residue
 Returns : a list of atoms attached to this residue
 Args    : a residue

=cut

sub get_atoms {
	my ($self, $residue) = @_;

	if ( !defined $residue) {
		$self->throw("get_atoms needs a Residue as argument");
	}
	# pass through to add_atom
	$self->add_atom($residue);
}


=head2 parent()

 Title   : parent
 Usage   : $structure->parent($residue);
 Function: returns the parent of the argument
 Returns : the parent of the argument
 Args    : a Bio::Structure object

=cut

=head2 conect()

 Title   : conect
 Usage   : $structure->conect($source);
 Function: get/set method for conect
 Returns : a list of serial numbers for atoms connected to source
 	(together with $entry->get_atom_by_serial($model, $serial) this should be OK for now)
 Args    : the serial number for the source atom

=cut

sub conect {
	my ($self, $source, $serial) = @_;
	
	if ( !defined $source ) {
		$self->throw("You need to supply at least a source to conect");
	}
	if ( defined $serial ) {
		if ( !exists(${$self->{'conect'}}{$source}) || ref(${$self->{'conect'}}{$source} !~ /^ARRAY/ ) ) {
			${$self->{'conect'}}{$source} = [];
		}
		push @{ ${$self->{'conect'}}{$source} }, $serial;
	}
	return @{ ${$self->{'conect'}}{$source} };
}

=head2 master()

 Title   : master
 Usage   : $structure->master($source);
 Function: get/set method for master
 Returns : the master line
 Args    : the master line for this entry

=cut

sub master {
	my ($self, $value) = @_;
	if (defined $value) {
		$self->{'master'} = $value;
	}
	return $self->{'master'};
}


=head2 seqres()

 Title   : seqres
 Usage   : $seqobj = $structure->seqres("A");
 Function: gets a sequence object containing the sequence from the SEQRES record.
 	    if a chain-ID is given , the sequence for this chain is given, if none
	    is provided the first chain is choosen
 Returns : a Bio::PrimarySeq
 Args    : the chain-ID of the chain you want the sequence from

=cut

sub seqres {
	my ($self, $chainid) = @_;
	my $s_u = "x4 A1 x7 A3 x1 A3 x1 A3 x1 A3 x1 A3 x1 A3 x1 A3 x1 A3 x1 A3 x1 A3 x1 A3 x1 A3 x1 A3";
	my $seq;
	if ( !defined $chainid) {
		my $m = ($self->get_models($self))[0];
		my $c = ($self->get_chains($m))[0];
		$chainid = $c->id;
	}
	my $seqres = ($self->annotation->get_Annotations("seqres"))[0];
	my $seqres_string = $seqres->as_text;
$self->debug("seqres : $seqres_string\n");
	$seqres_string =~ s/^Value: //;
	my $pos = 0;
	while ($pos < length $seqres_string) {
		my $st = substr($seqres_string, $pos, 63);
		$pos += 63;
$self->debug("seqres: $st\n");
		my ($chain,@res) = unpack $s_u, $st;
		$chain = 'default' if ($chain =~ /^\s*$/);
		next if ($chain ne $chainid);
		for my $res (@res) {
			next if $res =~ /^\s*$/;
			$seq .= $self->_three_to_one($res);
		}
	}
	my $pseq = Bio::PrimarySeq->new;
	$pseq->seq($seq);
	my $id = $self->id . "_" . $chainid;
	$pseq->id($id);
	return $pseq;
}

=head2 get_atom_by_serial()

 Title   : get_atom_by_serial
 Usage   : $structure->get_atom_by_serial($module, $serial);
 Function: get the Atom for a  for get_atom_by_serial
 Returns : the Atom object with this serial number in the model
 Args    : Model on which to work, serial number for atom
 	(if only a number is supplied, the first model is chosen)

=cut

sub get_atom_by_serial {
	my ($self, $model, $serial) = @_;

	if ($model =~ /^\d+$/ && !defined $serial) { # only serial given
		$serial = $model;
		my @m = $self->get_models($self);
		$model = $m[0];
	}
	if ( !defined $model || ref($model) !~ /^Bio::Structure::Model/ ) {
		$self->throw("Could not find (first) model\n");
	}
	if ( !defined $serial || ($serial !~ /^\d+$/) ) {
		$self->throw("The serial number you provided looks fishy ($serial)\n");
	}
	for my $chain ($self->get_chains($model) ) {
		for my $residue ($self->get_residues($chain) ) {
			for my $atom ($self->get_atoms($residue) ) {
				# this could get expensive, do we cache ???
				next unless ($atom->serial == $serial);
				return $atom;
			}
		}
	}
} 

sub parent {
	my ($self, $obj) = @_;
	
	if ( !defined $obj) {
		$self->throw("parent: you need to supply an argument to get the parent from\n");
	}

	# for now we pass on to _parent, untill we get the symbolic ref thing working.
	$self->_parent($obj);
}

sub DESTROY {
	my $self = shift;

	#print STDERR "DESTROY on $self being called\n";

##	for my $pc (keys %{ $self->{'p_c'} } ) {
##		next unless ( defined ${ $self->{'p_c'} }{$pc} );
##		delete ${$self->{'p_c'}}{$pc};
##	}
##	for my $cp (keys %{ $self->{'c_p'} } ) {
##		next unless ( defined ${ $self->{'c_p'} }{$cp} );
##		delete ${$self->{'c_p'}}{$cp};
##	}
	%{ $self->{'p_c'} } = ();
	%{ $self->{'c_p'} } = ();
}

# copied from Bio::Seq.pm
#
=head2 annotation

 Title   : annotation
 Usage   : $obj->annotation($seq_obj)
 Function:
 Example :
 Returns : value of annotation
 Args    : newvalue (optional)


=cut

sub annotation {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'annotation'} = $value;
    }
    return $obj->{'annotation'};

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

	;
}


=head2 _create_default_model()

 Title   : _create_default_model
 Usage   : 
 Function: Creates a default Model for this Entry. Typical situation
 	in an X-ray structure where there is only one model
 Returns : 
 Args    : 

=cut

sub _create_default_model {
	my ($self) = shift;

	my $model = Bio::Structure::Model->new(-id => "default");
	return $model;
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
	return $chain;
}



=head2 _parent()

 Title   : _parent
 Usage   : This is an internal function only. It is used to have one 
 	place that keeps track of which object has which other object 
	as parent. Thus allowing the underlying modules (Atom, Residue,...)
	to have no knowledge about all this (and thus removing the possibility
	of reference cycles).
	This method hides the details of manipulating references to an anonymous
	hash.
 Function: To get/set an objects parent 
 Returns : a reference to the parent if it exist, undef otherwise. In the 
 	current implementation each node should have a parent (except Entry).
 Args    : 

=cut

# manipulating the c_p hash

sub _parent {
	no strict "refs";
	my ($self, $key, $value) = @_;
	
	if ( (!defined $key) || (ref($key) !~ /^Bio::/) ) {
		$self->throw("First argument to _parent needs to be a reference to a Bio:: object ($key)\n");
	}
	if ( (defined $value) && (ref($value) !~ /^Bio::/) ) {
		$self->throw("Second argument to _parent needs to be a reference to a Bio:: object\n");
	}
	# no checking here for consistency of key and value, needs to happen in caller
	
	if (defined $value) {
		# is this value already in, shout
		if (exists( ${ $self->{'c_p'} }{$key} ) ) {
			$self->throw("_parent: $key already has a parent ${$self->{'c_p'}}{$key}\n");
		}
		${$self->{'c_p'}}{$key} = $value;
	}
	return ${$self->{'c_p'}}{$key}; 
}


=head2 _child()

 Title   : _child
 Usage   : This is an internal function only. It is used to have one 
 	place that keeps track of which object has which other object 
	as child. Thus allowing the underlying modules (Atom, Residue,...)
	to have no knowledge about all this (and thus removing the possibility
	to have no knowledge about all this (and thus removing the possibility
	of reference cycles).
	This method hides the details of manipulating references to an anonymous
	hash.
 Function: To get/set an object's child(ren) 
 Returns : a reference to an array of child(ren) if it exist, undef otherwise. 
 Args    : 

=cut

# manipulating the p_c hash
sub _child {
	my ($self, $key, $value) = @_;
	
	if ( (!defined $key) || (ref($key) !~ /^Bio::/) ) {
		$self->throw("First argument to _child needs to be a reference to a Bio:: object\n");
	}
	if ( (defined $value) && (ref($value) !~ /^Bio::/) ) {
		$self->throw("Second argument to _child needs to be a reference to a Bio:: object\n");
	}
	# no checking here for consistency of key and value, needs to happen in caller
	
	if (defined $value) {
		if ( !exists(${$self->{'p_c'}}{$key}) || ref(${$self->{'p_c'}}{$key}) !~ /^ARRAY/ ) {
			${$self->{'p_c'}}{$key} = [];
		}
		push @{ ${$self->{'p_c'}}{$key} }, $value;
	}
	return  ${$self->{'p_c'}}{$key}; 
}



=head2 _remove_from_graph()

 Title   : _remove_from_graph
 Usage   : This is an internal function only. It is used to remove from
 	the parent/child graph. We only remove the links from object to
	his parent. Not the ones from object to its children.
 Function: To remove an object from the parent/child graph
 Returns : 
 Args    : the object to be orphaned

=cut

sub _remove_from_graph {
	my ($self, $object) = @_;
	
	if ( !defined($object) && ref($object) !~ /^Bio::/) {
		$self->throw("_remove_from_graph needs a Bio object as argument");
	}
	if ( $self->_parent($object) ) {
		my $dad = $self->_parent($object);
		# if we have a parent, remove me as being a child
		for my $k (0 .. $#{$self->_child($dad)}) {
			if ($object eq ${$self->{'p_c'}{$dad}}[$k]) {
				splice(@{$self->{'p_c'}{$dad}}, $k,1);
			}
		}
		delete( $self->{'c_p'}{$object});
	}
}

			
sub _print_stats_pc {
	# print stats about the parent/child hashes
	my ($self) =@_;
	my $pc = scalar keys %{$self->{'p_c'}};
	my $cp = scalar keys %{$self->{'c_p'}};
	my $now_time = Time::HiRes::time();
	$self->debug("pc stats: P_C $pc C_P $cp $now_time\n");
}

sub _three_to_one {
	# convert three letter AA codes to one letter, return X if unknown
	my ($self, $three) = @_;
	my $one;
	my %three_to_one = (
		"ALA" => "A",
		"ARG" => "R",
		"ASN" => "N",
		"ASP" => "D",
		"CYS" => "C",
		"GLN" => "Q",
		"GLU" => "E",
		"GLY" => "G",
		"HIS" => "H",
		"ILE" => "I",
		"LEU" => "L",
		"LYS" => "K",
		"MET" => "M",
		"PHE" => "F",
		"PRO" => "P",
		"SER" => "S",
		"THR" => "T",
		"TRP" => "W",
		"TYR" => "Y",
		"VAL" => "V"
	);
	if(exists $three_to_one{$three}) {
		$one = $three_to_one{$three};
	} 
	else {
		$one = "X";
	}
	return $one;
}


	


1;
