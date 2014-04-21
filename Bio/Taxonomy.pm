#
# BioPerl module for Bio::Taxonomy
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Juguang Xiao
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Taxonomy - representing Taxonomy.

=head1 SYNOPSIS

  # NB: This module is deprecated. Use Bio::Taxon in combination with
  # Bio::Tree::Tree methods instead.

  use Bio::Taxonomy;

  # CREATION: You can either create an instance by assigning it,
  # or fetch it through factory.

  # Create the nodes first. See Bio::Taxonomy::Node for details.
  my $node_species_sapiens = Bio::Taxonomy::Node->new(
      -object_id => 9606, # or -ncbi_taxid. Requird tag
      -names => {
          'scientific' => ['sapiens'],
          'common_name' => ['human']
      },
      -rank => 'species'  # Required tag
  );
  my $node_genus_Homo = Bio::Taxonomy::Node->new(
      -object_id => 9605,
      -names => { 'scientific' => ['Homo'] },
      -rank => 'genus'
  );
  my $node_class_Mammalia = Bio::Taxonomy::Node->new(
      -object_id => 40674,
      -names => {
          'scientific' => ['Mammalia'],
          'common' => ['mammals']
      },
      -rank => 'class'
  );
  my $taxonomy = Bio::Taxonomy->new;
  $taxonomy->add_node($node_class_Mammalia);
  $taxonomy->add_node($node_species_sapiens);
  $taxonomy->add_node($node_genus_Homo);

  # OR you can fetch it through a factory implementing
  # Bio::Taxonomy::FactoryI
  my $factory;

  my $taxonomy = $factory->fetch_by_ncbi_taxid(40674);

  # USAGE

  # In this case, binomial returns a defined value.
  my $binomial = $taxonomy->binomial;

  # 'common_names' refers to the lowest-rank node's common names, in
  # array.
  my @common_names = $taxonomy->common_names;

  # 'get_node', will return undef if the rank is no defined in
  # taxonomy object.  It will throw error if the rank string is not
  # defined, say 'species lah'.
  my $node = $taxonomy->get_node('class');
  my @nodes = $taxonomy->get_all_nodes;

  # Also, you can search for parent and children nodes, if taxonomy
  # comes with factory.

  my $parent_taxonomy = $taxonomy->get_parent

=head1 DESCRIPTION

Bio::Taxonomy object represents any rank-level in taxonomy system,
rather than Bio::Species which is able to represent only
species-level.

There are two ways to create Taxonomy object, e.g.
1) instantiate an object and assign all nodes on your own code; and
2) fetch an object by factory.

=head2 Creation by instantiation

The abstraction of Taxonomy is actually a hash in data structure
term. The keys of the hash are the rank names, such as 'genus' and
'species', and the values are the instances of Bio::Taxonomy::Node.

=head2 Creation by Factory fetching

NCBI Taxonomy system is well accepted as the standard. The Taxonomy
Factories in bioperl access this system, through HTTP to NCBI Entrez,
dump file, and advanced biosql database.

Bio::Taxonomy::FactoryI defines all methods that all implementations
must obey.

$factory-E<gt>fetch is a general method to fetch Taxonomy by either
NCBI taxid or any types of names.

$factory-E<gt>fetch_parent($taxonomy), returns a Taxonomy that is
one-step higher rank of the taxonomy specified as argument.

$factory-E<gt>fetch_children($taxonomy), reports an array of Taxonomy
those are one-step lower rank of the taxonomy specified as the
argument.

=head2 Usage of Taxonomy object

##

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 CONTACT

Juguang Xiao, juguang@tll.org.sg

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# code begins...


package Bio::Taxonomy;
use strict;


use base qw(Bio::Root::Root);


=head2 new

 Title   : new
 Usage   : my $obj = Bio::Taxonomy->new();
 Function: Builds a new Bio::Taxonomy object
 Returns : Bio::Taxonomy
 Args    : -method  -> method used to decide classification
                       (none|trust|lookup)
           -ranks   -> what ranks are there

=cut


sub new {
   my ($class,@args) = @_;

   my $self = $class->SUPER::new(@args);
   $self->warn("Bio::Taxonomy is deprecated. Use Bio::Taxon in combination with Bio::Tree::Tree instead.");

   $self->{'_method'}='none';
   $self->{'_ranks'}=[];
   $self->{'_rank_hash'}={};
    $self->{_hierarchy} = {}; # used to store the nodes, with ranks as keys.
   my ($method,$ranks,$order) = $self->_rearrange([qw(METHOD RANKS ORDER)], @args);

   if ($method) {
      $self->method($method);
   }

   if (defined $ranks &&
      (ref($ranks) eq "ARRAY") ) {
      $self->ranks(@$ranks);
   } else {
      # default ranks
      # I think these are in the right order, but not sure:
      # some parvorder|suborder and varietas|subspecies seem
      # to be at the same level - any taxonomists?
      # I don't expect that these will actually be used except as a way
      # to find what ranks there are in taxonomic use
      $self->ranks(('root',
        'superkingdom', 'kingdom',
        'superphylum', 'phylum', 'subphylum',
        'superclass', 'class', 'subclass', 'infraclass',
        'superorder', 'order', 'suborder', 'parvorder', 'infraorder',
        'superfamily', 'family', 'subfamily',
        'tribe', 'subtribe',
        'genus', 'subgenus',
        'species group', 'species subgroup', 'species', 'subspecies',
        'varietas', 'forma', 'no rank'));
   }

   return $self;
}


=head2 method

 Title   : method
 Usage   : $obj = taxonomy->method($method);
 Function: set or return the method used to decide classification
 Returns : $obj
 Args    : $obj

=cut


sub method {
   my ($self,$value) = @_;
   if (defined $value && $value=~/none|trust|lookup/) {
       $self->{'_method'} = $value;
   }
   return $self->{'_method'};
}


=head2 classify

 Title   : classify
 Usage   : @obj[][0-1] = taxonomy->classify($species);
 Function: return a ranked classification
 Returns : @obj of taxa and ranks as word pairs separated by "@"
 Args    : Bio::Species object

=cut


sub classify {
   my ($self,$value) = @_;
   my @ranks;

   if (! $value->isa('Bio::Species') ) {
      $self->throw("Trying to classify $value which is not a Bio::Species object");
   }

   my @classes=reverse($value->classification);

   if ($self->method eq 'none') {
      for (my $i=0; $i < @classes-2; $i++) {
         ($ranks[$i][0],$ranks[$i][1])=($classes[$i],'no rank');
      }
      push @ranks,[$classes[-2],'genus'];
      push @ranks,[$value->binomial,'species'];
   } elsif ($self->method eq 'trust') {
      if (scalar(@classes)==scalar($self->ranks)) {
         for (my $i=0; $i < @classes; $i++) {
            if ($self->rank_of_number($i) eq 'species') {
               push @ranks,[$value->binomial,$self->rank_of_number($i)];
            } else {
               push @ranks,[$classes[$i],$self->rank_of_number($i)];
            }
         }
      } else {
         $self->throw("Species object and taxonomy object cannot be reconciled");
      }
   } elsif ($self->method eq 'lookup') {
      # this will lookup a DB for the rank of a taxon name
      # I imagine that some kind of Bio::DB class will be need to
      # be given to the taxonomy object to act as an DB interface
      # (I'm not sure how useful this is though - if you have a DB of
      # taxonomy - why would you be doing things this way?)
      $self->throw_not_implemented();
   }

   return @ranks;
}


=head2 level_of_rank

 Title   : level_of_rank
 Usage   : $obj = taxonomy->level_of_rank($obj);
 Function: returns the level of a rank name
 Returns : $obj
 Args    : $obj

=cut


sub level_of {
   my ($self,$value) = @_;

   return $self->{'_rank_hash'}{$value};
}


=head2 rank_of_number

 Title   : rank_of_number
 Usage   : $obj = taxonomy->rank_of_number($obj);
 Function: returns the rank name of a rank level
 Returns : $obj
 Args    : $obj

=cut


sub rank_of_number {
   my ($self,$value) = @_;

   return ${$self->{'_ranks'}}[$value];
}


=head2 ranks

 Title   : ranks
 Usage   : @obj = taxonomy->ranks(@obj);
 Function: set or return all ranks
 Returns : @obj
 Args    : @obj

=cut


sub ranks {
   my ($self,@value) = @_;

   # currently this makes no uniqueness sanity check (this should be done)
   # I am think that adding a way of converting multiple 'no rank' ranks
   # to unique 'no rank #' ranks so that the level of a 'no rank' is
   # abstracted way from the user - I'm not sure of the value of this

   if (@value) {
      $self->{'_ranks'}=\@value;
   }

   for (my $i=0; $i <= @{$self->{'_ranks'}}-1; $i++) {
      $self->{'_rank_hash'}{$self->{'_ranks'}[$i]}=$i unless $self->{'_ranks'}[$i] eq 'no rank';
   }

   return @{$self->{'_ranks'}};
}

=head2 add_node

  Title:    add_node
  Usage:    $obj->add_node($node[, $node2, ...]);
  Function: add one or more Bio::Taxonomy::Node objects
  Returns:  None
  Args:     any number of Bio::Taxonomy::Node(s)

=cut

sub add_node {
    my ($self, @nodes) = @_;
    foreach(@nodes){
        $self->throw("A Bio::Taxonomy::Node object needed")
            unless($_->isa('Bio::Taxonomy::Node'));
        my ($node, $rank) = ($_, $_->rank);
        if(exists $self->{_hierarchy}->{$rank}){
#            $self->throw("$rank has been defined");
#            print STDERR "RANK:$rank\n";
#            return;
        }
        $self->{_hierarchy}->{$rank} = $node;
    }
}

=head2 binomial

  Title   : binomial
  Usage   : my $val = $obj->binomial;
  Function: returns the binomial name if this taxonomy reachs species level
  Returns : the binomial name
            OR undef if taxonmy does not reach species level
  Args    : [No arguments]

=cut

sub binomial {
    my $self = shift;
    return $self->get_node('species')->scientific_name;
    my $genus = $self->get_node('genus');
    my $species = $self->get_node('species');
    return ($species && $genus) ? "$species $genus" : undef;
}

=head2 get_node

  Title   : get_node
  Usage   : $node = $taxonomy->get_node('species');
  Function: get a Bio::Taxonomy::Node object according to rank name
  Returns : a Bio::Taxonomy::Node object or undef if null
  Args    : a vaild rank name

=cut

sub get_node {
    my ($self, $rank) = @_;
    unless(grep /$rank/, keys %{$self->{_hierarchy}}){
        $self->throw("'$rank' is not in the rank list");
    }
    return (exists $self->{_hierarchy}->{$rank})?
        $self->{_hierarchy}->{$rank} : undef;
}

=head2 classification

  Title   : classification
  Usage   : @names = $taxonomy->classification;
  Function: get the classification names of one taxonomy
  Returns : array of names
  Args    : [No arguments]

=cut

sub classification {
    my $self = shift;
    my %rank_hash = %{$self->{_rank_hash}};
    my %hierarchy = %{$self->{_hierarchy}};
    my @ordered_nodes = sort {
        ($rank_hash{$a} <=> $rank_hash{$b})
    } keys %hierarchy;
    return map {$hierarchy{$_}->scientific_name} @ordered_nodes;
}

1;
