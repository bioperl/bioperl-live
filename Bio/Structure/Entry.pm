#
# bioperl module for Bio::Structure::Entry
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

Bio::Structure::Entry - Bioperl structure Object, describes the whole entry

=head1 SYNOPSIS

  #add synopsis here

=head1 DESCRIPTION

This object stores a whole Bio::Structure entry. It can consist of one
or more models (L<Bio::Structure::Model>), which in turn consist of one 
or more chains (L<Bio::Structure::Chain>). A chain is composed of residues 
(L<Bio::Structure::Residue>) and a residue consists of atoms 
(L<Bio::Structure::Atom>). If no specific model or chain is chosen, the 
first one is chosen.

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

The rest of the documentation details each of the object methods. Internal 
methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Structure::Entry;
use strict;

use Bio::Structure::Model;
use Bio::Structure::Chain;
use Bio::Annotation::Collection;
use Tie::RefHash;

use base qw(Bio::Root::Root Bio::Structure::StructureI);

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
			      RESIDUE )], @args);

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
 Returns : List of Bio::Structure::Model objects
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
 Returns : The ID
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
 Function: Connects a Chain or a list of Chain objects to a Bio::Structure::Entry.
 Returns : List of Bio::Structure::Chain objects
 Args    : A Chain or a reference to an array of Chain objects

=cut

sub chain {
	my ($self, $chain) = @_;

	if ( ! $self->model ) {
		$self->_create_default_model;
	}
	my @models = $self->model;
	my $first_model = $models[0];

	if ( defined $chain) {
		
		if( (ref($chain) eq "ARRAY") || ($chain->isa('Bio::Structure::Chain')) ) {
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
 Usage   : @chains  = $structure->add_chain($model,$chain);
 Function: Adds one or more Chain objects to a Bio::Structure::Entry.
 Returns : List of Chain objects associated with the Model
 Args    : A Model object and a Chain object or a reference to an array of 
           of Chain objects

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
 Function: General get method for Chains attached to a Model
 Returns : A list of Chains attached to this model
 Args    : A Model

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
 Returns : List of Bio::Structure::Residue objects
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
 Usage   : @residues  = $structure->add_residue($chain,$residue);
 Function: Adds one or more Residue objects to a Bio::Structure::Entry.
 Returns : List of Bio::Structure::Residue objects
 Args    : A Chain object and a Residue object or a reference to an array of 
           Residue objects

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
 Function: General get method for Residues attached to a Chain
 Returns : A list of residues attached to this Chain
 Args    : A Chain

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
 Returns : List of Bio::Structure::Atom objects
 Args    : A Residue and an Atom

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
 Function: General get method for Atoms attached to a Residue
 Returns : A list of Atoms attached to this Residue
 Args    : A Residue

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
 Function: Returns the parent of the argument
 Returns : The parent of the argument
 Args    : A Bio::Structure object

=cut

=head2 connect

 Title   : connect
 Usage   : 
 Function: Alias to conect()
 Returns : 
 Args    : 

=cut

sub connect {
	my $self = shift;
	return $self->conect(@_);
}

=head2 conect()

 Title   : conect
 Usage   : $structure->conect($source);
 Function: Get/set method for conect
 Returns : A list of serial numbers for Atoms connected to source
 	        (together with $entry->get_atom_by_serial($model, $serial),
           this should be OK for now)
 Args    : The source, the serial number for the source Atom, and the type

=cut

sub conect {
	my ($self, $source, $serial, $type) = @_;
	
	if ( !defined $source ) {
		$self->throw("You need to supply at least a source to connect");
	}
	if ( defined $serial && defined $type ) {
		if ( !exists(${$self->{'conect'}}{$source}) || 
			  ref(${$self->{'conect'}}{$source} !~ /^ARRAY/ ) ) {
			${$self->{'conect'}}{$source} = [];
		}
		# we also need to store type, a conect object might be better 
		my $c = $serial . "_" . $type;
		push @{ ${$self->{'conect'}}{$source} }, $c;
	}
	# Bug 1894
	return () if ( !exists $self->{'conect'}{$source} || 
					  !defined $self->{'conect'}{$source} );
	return @{ ${$self->{'conect'}}{$source} };
}

=head2 get_all_connect_source

 Title   : get_all_connect_source
 Usage   : 
 Function: Alias to get_all_conect_source()
 Returns : 
 Args    : 

=cut

sub get_all_connect_source {
	my $self = shift;
	return get_all_conect_source(@_);
}

=head2 get_all_conect_source()

 Title   : get_all_conect_source
 Usage   : @sources = $structure->get_all_conect_source;
 Function: Get all the sources for the conect records
 Returns : A list of serial numbers for atoms connected to source
 	        (together with $entry->get_atom_by_serial($model, $serial), 
           this should be OK for now)
 Args    : 
 Notes   : This is a bit of a kludge, but it is the best for now. Conect info might need
 	        to go in a separate object

=cut

sub get_all_conect_source {
	my ($self) = shift;
	my (@sources);

	for my $source (sort {$a<=>$b} keys %{$self->{'conect'}}) {
		push @sources, $source;
	}
	return @sources;
}


=head2 master()

 Title   : master
 Usage   : $structure->master($source);
 Function: Get/set method for master
 Returns : The master line
 Args    : The master line for this entry

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
 Function: Gets a sequence object containing the sequence from the SEQRES record.
 	        if a chain-ID is given, the sequence for this chain is given, if none
	        is provided the first chain is chosen
 Returns : A Bio::PrimarySeq
 Args    : The chain-ID of the chain you want the sequence from

=cut

sub seqres {
	my ($self, $chainid) = @_;
	my $s_u = "x3 A1 x7 A3 x1 A3 x1 A3 x1 A3 x1 A3 x1 A3 x1 A3 x1 A3 x1 A3 x1 A3 x1 A3 x1 A3 x1 A3";
	my (%seq_ch);
	if ( !defined $chainid) {
		my $m = ($self->get_models($self))[0];
		my $c = ($self->get_chains($m))[0];
		$chainid = $c->id;
	}
	my $seqres = ($self->annotation->get_Annotations("seqres"))[0];
	my $seqres_string = $seqres->as_text;
	$self->debug("seqres : $seqres_string\n");
	$seqres_string =~ s/^Value: //;
	# split into lines of 62 long
	my @l = unpack("A62" x (length($seqres_string)/62), $seqres_string);
	for my $line (@l) {
		# get out chain_id and sequence
		# we use a1, as A1 strips all spaces :(
		my ($chid, $seq) = unpack("x3 a1 x7 A51", $line);
		if ($chid eq " ") {
			$chid = "default";
		}
		$seq =~ s/(\w+)/\u\L$1/g;	# ALA -> Ala  (for SeqUtils)
		$seq =~ s/\s//g; 		# strip all spaces
		$seq_ch{$chid} .= $seq;
		$self->debug("seqres : $chid $seq_ch{$chid}\n");
	}
	# do we have a seqres for this chainid
	if(! exists $seq_ch{$chainid} ) {
		$self->warn("There is no SEQRES known for chainid \"$chainid\"");
		return;
	}

	# this will break for non-protein structures (about 10% for now) XXX KB
	my $pseq = Bio::PrimarySeq->new(-alphabet => 'protein');
	$pseq = Bio::SeqUtils->seq3in($pseq,$seq_ch{$chainid});
	my $id = $self->id . "_" . $chainid;
	$pseq->id($id);
	return $pseq;
}


=head2 get_atom_by_serial()

 Title   : get_atom_by_serial
 Usage   : $structure->get_atom_by_serial($model,$serial);
 Function: Get the Atom by serial
 Returns : The Atom object with this serial number in the model
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

	%{ $self->{'p_c'} } = ();
	%{ $self->{'c_p'} } = ();
}

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
 	        do not belong to this Entry any more
 Returns : 
 Args    : 

=cut

#

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
 Returns : A reference to the parent if it exist, undef otherwise. In the 
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
	if (defined ( $self->{'c_p'}->{$key}) && 
	    exists ( $self->{'c_p'}->{$key})
	    ) {
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
 Function: To get/set an the children of an object 
 Returns : A reference to an array of child(ren) if they exist, undef otherwise. 
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
 Args    : The object to be orphaned

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


1;
