# $Id$
#
# BioPerl module for Bio::Taxonomy::Node
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Taxonomy::Node - A node in a represented taxonomy

=head1 SYNOPSIS

  use Bio::Taxonomy::Node;
  # typically you will get a Node from a Bio::DB::Taxonomy object
  # but here is how you initialize one
  my $node = new Bio::Taxonomy::Node(-name      => $name,
                                     -object_id => $oid,
                                     -parent_id => $pid,
                                     -rank   => $rank,
                                     -division  => $div,
                                     -dbh       => $dbh);

  my $dbh = new Bio::DB::Taxonomy(-source   => 'flatfile',
                                  -directory=> '/tmp',
                                  -nodesfile=> '/path/to/nodes.dmp',
                                  -namesfile=> '/path/to/names.dmp');
  my $hum_node = $dbh->get_Taxonomy_Node(-name => 'Homo sapiens');
  my $hum_node2= $dbh->get_Taxonomy_Node(-taxonid => '9606');

  print "rank is ", $hum_node->rank, "\n";
  print "classification is ", join(" ", $hum_node->classification),"\n"; 
  print "division is ", $node->division, "\n";
  my $mmu_node = $dbh->get_Taxonomy_Node(-name => 'Mus musculus');
  my @mmu_lineage = $mmu->get_Lineage_Nodes;

=head1 DESCRIPTION

This is the next generation (for Bioperl) of representing Taxonomy
information.  Previously all information was managed by a single
object called Bio::Species.  This new implementation allows
representation of the intermediate nodes not just the species nodes
and can relate their connections.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 CONTRIBUTORS

Juguang Xiao,     juguang@tll.org.sg
Gabriel Valiente, valiente@lsi.upc.edu
Sendu Bala,       bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Taxonomy::Node;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::IdentifiableI;
use Bio::DB::Taxonomy;

@ISA = qw(Bio::Root::Root Bio::IdentifiableI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Taxonomy::Node();
 Function: Builds a new Bio::Taxonomy::Node object 
 Returns : an instance of Bio::Taxonomy::Node
 Args    : -dbh        => a reference to a Bio::DB::Taxonomy object [defaults to
                          one using -source => 'entrez']
           -name       => a string representing the node name
           -object_id  => unique identifier - typically NCBI Taxid
           -ncbi_taxid => alias for -object_id, only use one of these
           -parent_id  => parent id (unique identifier for parent)
           -rank       => node rank (one of 'species', 'genus', etc)
           -classification => Arrayref of full lineage classification
           -common_names   => array ref of all common names
           -genetic_code   => genetic code table number
           -mito_genetic_code => mitochondrial genetic code table number
           -create_date       => date created (where available)
           -update_date       => date last updated
           -pub_date          => date published
           
           #*** deprecated?
           -organelle      => organelle type where appropriate
           -division       => 'primates', 'rodents', or 'inv', etc
           -sub_species    => if this is a subspecies, provide that as well,
                              most relavent for virii and others
           -variant        => provide variant info

=cut

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($name, $uniqueid, $parentid,
	$rank, $div, $dbh, $classification, $ncbitaxid, $commonname, $commonnames,
	$organelle, $sub_species, $variant, $gcode, $mitocode,
	$createdate, $updatedate, $pubdate) = 
	    $self->_rearrange([ qw ( NAME OBJECT_ID PARENT_ID RANK DIVISION
				     DBH
				     CLASSIFICATION
				     NCBI_TAXID
				     COMMON_NAME
                     COMMON_NAMES
				     ORGANELLE
				     SUB_SPECIES
				     VARIANT
				     GENETIC_CODE
				     MITO_GENETIC_CODE
				     CREATE_DATE
				     UPDATE_DATE
				     PUB_DATE
				     )],
			      @args);
    
    if( defined $ncbitaxid && defined $uniqueid && $ncbitaxid ne $uniqueid  ) {
	$self->warn("Only provide one of -object_id or -ncbi_taxid, using $uniqueid\n");
    } elsif( ! defined $uniqueid && defined $ncbitaxid ) { 
	$uniqueid = $ncbitaxid;
    }
    $uniqueid && $self->object_id($uniqueid);
    defined $parentid && $self->parent_id($parentid);
    
    defined $rank && $self->rank($rank);
    defined $name && $self->node_name($name);
    
    my @common_names;
    if ($commonnames) {
        $self->throw("-common_names takes only an array reference") unless ref($commonnames) eq 'ARRAY';
        @common_names = @{$commonnames};
        if ($commonname) {
            my %c_names = map { $_ => 1 } @common_names;
            unless (exists $c_names{$commonname}) {
                unshift(@common_names, $commonname);
            }
        }
    }
    @common_names > 0 && $self->common_names(@common_names);
    
    defined $gcode      && $self->genetic_code($gcode);
    defined $mitocode   && $self->mitochondrial_genetic_code($mitocode);
    defined $createdate && $self->create_date($createdate);
    defined $updatedate && $self->update_date($updatedate);
    defined $pubdate    && $self->pub_date($pubdate);
    
    $self->db_handle($dbh || Bio::DB::Taxonomy->new(-source => 'entrez'));
    
    if( defined $classification ) {
    if( ref($classification) !~ /ARRAY/ ) {
        $self->warn("Classification can only be initialize with an arrayref\n");
    }
    $self->classification(@$classification);
    }
    
    #*** deprecated?
    defined $organelle && $self->name('organelle', $organelle);
    defined $div  && $self->division($div);
    defined $sub_species && $self->name('sub_species', $sub_species);
    defined $variant && $self->name('variant',$variant);
    
    return $self;
}

=head1 Bio::IdentifiableI interface 

Also see L<Bio::IdentifiableI>

=head2 version

 Title   : version
 Usage   : $obj->version($newval)
 Function: 
 Example : 
 Returns : value of version (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub version {
    my $self = shift;
    return $self->{'version'} = shift if @_;
    return $self->{'version'};
}

=head2 authority

 Title   : authority
 Usage   : $obj->authority($newval)
 Function: 
 Example : 
 Returns : value of authority (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub authority {
    my $self = shift;
    return $self->{'authority'} = shift if @_;
    return $self->{'authority'};
}

=head2 namespace

 Title   : namespace
 Usage   : $obj->namespace($newval)
 Function: 
 Example : 
 Returns : value of namespace (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub namespace {
    my $self = shift;
    return $self->{'namespace'} = shift if @_;
    return $self->{'namespace'};
}

=head1 Bio::Taxonomy::Node implementation

=head2 db_handle

 Title   : db_handle
 Usage   : $obj->db_handle($newval)
 Function: Get/Set Bio::DB::Taxonomy Handle
 Returns : value of db_handle (a scalar) (Bio::DB::Taxonomy object)
 Args    : on set, new value (a scalar or undef, optional) Bio::DB::Taxonomy object

Also see L<Bio::DB::Taxonomy>

=cut

sub db_handle {
    my $self = shift;
    if( @_ ) {
	my $v = shift;
	# until we establish some other higher level TaxonomyDB interface
	if( ! ref($v) || ! $v->isa('Bio::DB::Taxonomy') ) {
	    $self->throw("Must have provided a valid Bio::DB::Taxonomy object");
	}
	$self->{'db_handle'} = $v;
    }
    return $self->{'db_handle'};
}

=head2 rank

 Title   : rank
 Usage   : $obj->rank($newval)
 Function: Get/set rank of this Node, 'species', 'genus', 'order', etc...
 Returns : value of rank (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub rank {
    my $self = shift;
    return $self->{'rank'} = shift if @_;
    return $self->{'rank'};
}

=head2 object_id

 Title   : object_id
 Usage   : $obj->object_id($newval)
 Function: Get/Set object id (NCBI Taxonomy ID in most cases); ncbi_taxid() is
           a synonym of this method.
 Returns : value of object_id (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub object_id {
    my $self = shift;
    return $self->{'_object_id'} = shift if @_;
    return $self->{'_object_id'};
}

*ncbi_taxid = \&object_id;

=head2 parent_id

 Title   : parent_id
 Usage   : $obj->parent_id($newval)
 Function: Get/Set parent ID, (NCBI Taxonomy ID in most cases);
           parent_taxon_id() is a synonym of this method.
 Returns : value of parent_id (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub parent_id {
    my $self = shift;
    return $self->{'parent_id'} = shift if @_;
    return $self->{'parent_id'};
}

*parent_taxon_id = \&parent_id;

=head2 genetic_code

 Title   : genetic_code
 Usage   : $obj->genetic_code($newval)
 Function: Get/set genetic code table
 Returns : value of genetic_code (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub genetic_code {
    my $self = shift;
    return $self->{'genetic_code'} = shift if @_;
    return $self->{'genetic_code'};
}

=head2 mitochondrial_genetic_code

 Title   : mitochondrial_genetic_code
 Usage   : $obj->mitochondrial_genetic_code($newval)
 Function: Get/set mitochondrial genetic code table
 Returns : value of mitochondrial_genetic_code (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub mitochondrial_genetic_code {
    my $self = shift;
    return $self->{'mitochondrial_genetic_code'} = shift if @_;
    return $self->{'mitochondrial_genetic_code'};
}

=head2 create_date

 Title   : create_date
 Usage   : $obj->create_date($newval)
 Function: Get/Set Date this node was created (in the database)
 Returns : value of create_date (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub create_date {
    my $self = shift;
    return $self->{'create_date'} = shift if @_;
    return $self->{'create_date'};
}

=head2 update_date

 Title   : update_date
 Usage   : $obj->update_date($newval)
 Function: Get/Set Date this node was updated (in the database)
 Returns : value of update_date (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub update_date {
    my $self = shift;
    return $self->{'update_date'} = shift if @_;
    return $self->{'update_date'};
}

=head2 pub_date

 Title   : pub_date
 Usage   : $obj->pub_date($newval)
 Function: Get/Set Date this node was published (in the database)
 Returns : value of pub_date (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub pub_date {
    my $self = shift;
    return $self->{'pub_date'} = shift if @_;
    return $self->{'pub_date'};
}

=head2 get_Parent_Node

 Title   : get_Parent_Node
 Usage   : my $parentnode = $node->get_Parent_Node()
 Function: Retrieve the full Parent node from the database
 Returns : Bio::Taxonomy::Node
 Args    : none

=cut

sub get_Parent_Node {
   my ($self) = @_;

   if( ! $self->db_handle ||
       ! defined $self->parent_id ) {
       $self->warn("Cannot get the parent node for ".$self->node_name.
		   " because parent_id or db handle is not defined\n");
       return;
   }

   my $node = $self->db_handle->get_Taxonomy_Node(-taxonid => $self->parent_id);
   unless ( defined $node ) {
       $self->warn("Could not find node for parent id ". $self->parent_id);
       return;
   }
   return $node;
}

=head2 get_Children_Nodes

 Title   : get_Children_Nodes
 Usage   : my @nodes = $node->get_Children_Nodes();
 Function: Get the children of a node as L<Bio::Taxonomy::Node> objects
 Returns : Array of L<Bio::Taxonomy::Node> objects
 Args    : none

=cut

sub get_Children_Nodes {
   my ($self) = @_;
   if( ! $self->db_handle ||
       ! defined $self->object_id ) {
       $self->warn("Cannot get the children nodes for ".$self->name('common').
		   " because object_id or db handle is not defined\n");
       return;
   }
   my @nodes = $self->db_handle->get_Children_Taxids($self->object_id);
   my @children;
   for my $n ( @nodes ) {
       my $node = $self->db_handle->get_Taxonomy_Node(-taxonid => $n);
       unless( defined $node ) {
	   $self->warn("Could not find node for id $n child of ".$self->object_id);
       } else {
	   push @children, $node;
       }
   }
   return @children;
}

=head2 get_Lineage_Nodes

 Title   : get_Lineage_Nodes
 Usage   : my @nodes = $node->get_Lineage_Nodes();
 Function: Get the full lineage of a node as L<Bio::Taxonomy::Node> objects
 Returns : Array of L<Bio::Taxonomy::Node> objects
 Args    : none

=cut

sub get_Lineage_Nodes {
   my ($self) = @_;
   if( ! $self->db_handle ) {
       $self->warn("Cannot get the lineage nodes for ".$self->node_name.
		   " because db handle is not defined\n");
       return;
   }
   my $node = $self;
   my @lineage;
   while ($node->node_name ? $node->node_name ne "root" : 1) {
      $node = $node->get_Parent_Node;
      $node || last;
      unshift @lineage, $node;
   }
   return @lineage;
}

=head2 get_LCA_Node

 Title   : get_LCA_Node
 Usage   : my $lca = $node1->get_LCA_Node($node2);
 Function: Get the most recent common ancestor of two nodes as a L<Bio::Taxonomy::Node> object
 Returns : L<Bio::Taxonomy::Node> object
 Args    : Bio::Taxonomy::Node object

=cut

sub get_LCA_Node {
   my ($self, $node) = @_;
   if ($self eq $node || ($self->object_id && $node->object_id && $self->object_id == $node->object_id)) {
      return $self;
   } else {
      my @PATH1 = $self->get_Lineage_Nodes;
      my @PATH2 = $node->get_Lineage_Nodes;
      my ($root1, $root2, $lca);
      do {
         $root1 = shift @PATH1;
         $root2 = shift @PATH2;
         
         # node_name isn't necessarily unique in the taxonomy, so must compare
         # on object_id
         my ($r1_oid, $r2_oid) = ($root1->object_id, $root2->object_id);
         unless (defined $r1_oid && defined $r2_oid) {
            $self->warn("One of the lineages had a node with no object_id, can't calculate the common ancestor");
            last;
         }
         $lca = $root1 if $root1->object_id eq $root2->object_id;
      } while ($root1->object_id eq $root2->object_id);
      return $lca;
   }
}

=head2 name

  Title:    name
  Usage:    $obj->name('scientific', 'Homo sapiens');
            $obj->name('common', 'human', 'man');
            my @names = @{$obj->name('common')};
  Function: Get/set the names. node_name(), scientific_name() and common_names()
            are shorthands to name('scientific'), name('scientific') and
            name('common') respectively.
  Returns:  names (a array reference)
  Args:     Arg1 => the name_class. You can assign any text, but the words
                'scientific' and 'common' have the special meaning, as
                scientific name and common name, respectively. 'scientific' is
                treated specially, allowing only the first value in the Arg2
                list to be set.
            Arg2 .. => list of names

=cut

sub name {
    my ($self, $name_class, @names) = @_;
    $self->throw('No name class specified') unless defined $name_class;
    
    if (@names) {
        if ($name_class =~ /scientific/i) {
            delete $self->{'_names_hash'}->{$name_class};
            @names = (shift(@names));
        }
        push @{$self->{'_names_hash'}->{$name_class}}, @names;
    }
    return $self->{'_names_hash'}->{$name_class} || return;
}

=head2 node_name

 Title   : node_name
 Usage   : $obj->node_name($newval)
 Function: Get/set the name of this node, typically the scientific name of the
           node, eg. 'Primate' or 'Homo'; scientific_name() is a synonym of this
           method.
 Returns : value of node_name (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub node_name {
    my $self = shift;
    my @v = @{$self->name('scientific', @_) || []};
    return shift @v;
}

*scientific_name = \&node_name;

=head2 binomial

 Title   : binomial
 Usage   : $binomial = $self->binomial();
           $binomial = $self->binomial('FULL');
 Function: Returns a string "Genus species", or "Genus species subspecies",
           if the first argument is 'FULL' (and the species has a subspecies).
 Args    : Optionally the string 'FULL' to get the full name including
           the subspecies.

=cut

sub binomial {
    my( $self, $full ) = @_;
    my $rank = $self->rank || '';
    if( $rank ne 'species' ) {
        $self->warn("Asking for binomial for a $rank node which is not a species, you may not get what you are expecting (genus, species is not available)\n");
    }
    
    # do we already have the binomial?
    my $sci_name = $self->scientific_name || '';
    if ($rank eq 'species' && $sci_name =~ /\w+\s+\w+/) {
        return $sci_name;
    }
    
    my( $species, $genus ) = $self->classification();
    unless( defined $species) {
        $species = 'sp.';
        $self->warn("requested binomial but classification was not set");
    }
    $genus = ''   unless( defined $genus);
    my $bi = "$genus $species";
    if (defined($full) && $full =~ /full/i) { 
        my $ssp = $self->sub_species;
        $bi .= " $ssp" if $ssp;
    }
    return $bi;
}

=head2 common_names

 Title   : common_names
 Usage   : $obj->common_names($newval)
 Function: Get/set the other names of this node, typically the genbank common
           name and others, eg. 'Human' and 'man'. common_name() is a synonym
           of this method.
 Returns : array of names
 Args    : on set, new list of names (scalars, optional)

=cut

sub common_names {
    my $self = shift;
    my @v = @{$self->name('common', @_) || []};
    return ( wantarray ) ? @v : shift @v;
}

*common_name = \&common_names;

=head2 classification

 Title   : classification
 Usage   : $self->classification(@class_array);
           @classification = $self->classification();
 Function: Fills/returns the classification list in the object.
 Returns : Classification array (scalars)
 Args    : [optional] array of vals to set classification to

=cut

sub classification {
   my ($self,@vals) = @_;
   my $p;
   if( @vals ) {
       $self->{'_classification'} = [@vals];
       return @vals;
   } elsif ( defined $self->{'_classification'} ) {
       return @{$self->{'_classification'}};
   } else {
       foreach my $node ($self->get_Lineage_Nodes, $self) {
            my $name = $node->node_name || next;
            if (($self->rank && $self->rank ne 'no rank') || $self->show_all) {
                unshift(@vals, $name);
            }
       }
       $self->{'_classification'} = \@vals;
   }
   return @vals;
}

=head2 show_all

 Title   : show_all
 Usage   : $obj->show_all($newval)
 Function: Boolean flag whether or not we should show all intermediete
           nodes that do not have actual ranks.
 Returns : boolean
 Args    : boolean

=cut

sub show_all{
    my $self = shift;
    return $self->{'show_all'} = shift if @_;
    return $self->{'show_all'};
}

=head2 division

 Title   : division
 Usage   : $obj->division($newval)
 Function: Get/set the division this node belongs to, eg. 'Primates' or
           'Bacteria'.
 Returns : value of division (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub division {
    my $self = shift;
    my @v = @{$self->name('division',@_) || []};
    return wantarray ? @v : shift @v;
}

=head2 species

 Title   : species
 Usage   : $self->species($species);
           $species = $self->species();
 Function: Get or set the scientific species name.
 Example : $self->species('Homo sapiens');
 Returns : Scientific species name as string
 Args    : Scientific species name as string

=cut

sub species {
    my($self, $species) = @_;
    
    my $rank = $self->rank || '';
    if ($rank ne 'species') {
        $self->warn("Getting or setting species for a '$rank' node which is not a species will probably result in a messed up classification()");
    }
    
    if (defined $species) {
        $self->{'_classification'}[0] = $species;
    }
    return $self->{'_classification'}[0];
}

=head2 genus

 Title   : genus
 Usage   : $self->genus( $genus );
           $genus = $self->genus();
 Function: Get or set the scientific genus name.
 Example : $self->genus( 'Homo' );
 Returns : Scientific genus name as string
 Args    : Scientific genus name as string

=cut

sub genus {
    my($self, $genus) = @_;
    
    my $rank = $self->rank || '';
    if ($rank ne 'species') {
        $self->warn("Getting or setting genus for a '$rank' node which is not a species will probably result in a messed up classification()");
    }
    
    if (defined $genus) {
        $self->{'_classification'}[1] = $genus;
    }
    return $self->{'_classification'}[1];
}

=head2 sub_species

 Title   : sub_species
 Usage   : $obj->sub_species($newval)
 Function:
 Returns : value of sub_species
 Args    : newvalue (optional)

=cut

sub sub_species {
    my $self = shift;
    my @v = @{$self->name('sub_species', @_) || []};
    return wantarray ? @v : shift @v;
}

=head2 organelle

 Title   : organelle
 Usage   : $self->organelle( $organelle );
           $organelle = $self->organelle();
 Function: Get or set the organelle name
 Example : $self->organelle('Chloroplast')
 Returns : The organelle name in a string
 Args    : String, which is the organelle name

=cut

sub organelle {
    my $self = shift;
    my @v = @{$self->name('organelle', @_) || []};
    return wantarray ? @v : shift @v;
}

=head2 variant

 Title   : variant
 Usage   : $obj->variant($newval)
 Function: Get/set variant information for this species object (strain,
           isolate, etc).
 Example : 
 Returns : value of variant (a scalar)
 Args    : new value (a scalar or undef, optional)


=cut

sub variant{
    my $self = shift;
    my @v = @{$self->name('variant',@_) || []};
    return wantarray ? @v : shift @v;
}

1;