#
# BioPerl module for Bio::Taxon
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Sendu Bala, based heavily on a module by Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Taxon - A node in a represented taxonomy

=head1 SYNOPSIS

  use Bio::Taxon;

  # Typically you will get a Taxon from a Bio::DB::Taxonomy object
  # but here is how you initialize one
  my $taxon = Bio::Taxon->new(-name      => $name,
                              -id        => $id,
                              -rank      => $rank,
                              -division  => $div);

  # Get one from a database
  my $dbh = Bio::DB::Taxonomy->new(-source   => 'flatfile',
                                   -directory=> '/tmp',
                                   -nodesfile=> '/path/to/nodes.dmp',
                                   -namesfile=> '/path/to/names.dmp');
  my $human = $dbh->get_taxon(-name => 'Homo sapiens');
  $human = $dbh->get_taxon(-taxonid => '9606');

  print "id is ", $human->id, "\n"; # 9606
  print "rank is ", $human->rank, "\n"; # species
  print "scientific name is ", $human->scientific_name, "\n"; # Homo sapiens
  print "division is ", $human->division, "\n"; # Primates

  my $mouse = $dbh->get_taxon(-name => 'Mus musculus');

  # You can quickly make your own lineages with the list database
  my @ranks = qw(superkingdom class genus species);
  my @h_lineage = ('Eukaryota', 'Mammalia', 'Homo', 'Homo sapiens');
  my $list_dbh = Bio::DB::Taxonomy->new(-source => 'list', -names => \@h_lineage,
                                                           -ranks => \@ranks);
  $human = $list_dbh->get_taxon(-name => 'Homo sapiens');
  my @names = $human->common_names; # @names is empty
  $human->common_names('woman');
  @names = $human->common_names; # @names contains woman

  # You can switch to another database when you need more information
  my $entrez_dbh = Bio::DB::Taxonomy->new(-source => 'entrez');
  $human->db_handle($entrez_dbh);
  @names = $human->common_names; # @names contains woman, human, man

  # Since Bio::Taxon implements Bio::Tree::NodeI, we have access to those
  # methods (and can manually create our own taxa and taxonomy without the use
  # of any database)
  my $homo = $human->ancestor;

  # Though be careful with each_Descendent - unless you add_Descendent()
  # yourself, you won't get an answer because unlike for ancestor(), Bio::Taxon
  # does not ask the database for the answer. You can ask the database yourself
  # using the same method:
  ($human) = $homo->db_handle->each_Descendent($homo);

  # We can also take advantage of Bio::Tree::Tree* methods:
  # a) some methods are available with just an empty tree object
  use Bio::Tree::Tree;
  my $tree_functions = Bio::Tree::Tree->new();
  my @lineage = $tree_functions->get_lineage_nodes($human);
  my $lineage = $tree_functions->get_lineage_string($human);
  my $lca = $tree_functions->get_lca($human, $mouse);

  # b) for other methods, create a tree using your Taxon object
  my $tree = Bio::Tree::Tree->new(-node => $human);
  my @taxa = $tree->get_nodes;
  $homo = $tree->find_node(-rank => 'genus');

  # Normally you can't get the lca of a list-database derived Taxon and an
  # entrez or flatfile-derived one because the two different databases might
  # have different roots and different numbers of ranks between the root and the
  # taxa of interest. To solve this, make a tree of the Taxon with the more
  # detailed lineage and splice out all the taxa that won't be in the lineage of
  # your other Taxon:
  my $entrez_mouse = $entrez_dbh->get_taxon(-name => 'Mus musculus');
  my $list_human = $list_dbh->get_taxon(-name => 'Homo sapiens');
  my $mouse_tree = Bio::Tree::Tree->new(-node => $entrez_mouse);
  $mouse_tree->splice(-keep_rank => \@ranks);
  $lca = $mouse_tree->get_lca($entrez_mouse, $list_human);

=head1 DESCRIPTION

This is the next generation (for Bioperl) of representing Taxonomy
information. Previously all information was managed by a single
object called Bio::Species. This new implementation allows
representation of the intermediate nodes not just the species nodes
and can relate their connections.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 CONTRIBUTORS

Jason Stajich,    jason-at-bioperl-dot-org (original Bio::Taxonomy::Node)
Juguang Xiao,     juguang@tll.org.sg
Gabriel Valiente, valiente@lsi.upc.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


package Bio::Taxon;
use strict;
use Scalar::Util qw(blessed);

use Bio::DB::Taxonomy;

use base qw(Bio::Tree::Node Bio::IdentifiableI);


=head2 new

 Title   : new
 Usage   : my $obj = Bio::Taxonomy::Node->new();
 Function: Builds a new Bio::Taxonomy::Node object 
 Returns : an instance of Bio::Taxonomy::Node
 Args    : -dbh               => a reference to a Bio::DB::Taxonomy object
                                 [no default]
           -name              => a string representing the taxon name
                                 (scientific name)
           -id                => human readable id - typically NCBI taxid
           -ncbi_taxid        => same as -id, but explicitly say that it is an
                                 NCBI taxid
           -rank              => node rank (one of 'species', 'genus', etc)
           -common_names      => array ref of all common names
           -division          => 'Primates', 'Rodents', etc
           -genetic_code      => genetic code table number
           -mito_genetic_code => mitochondrial genetic code table number
           -create_date       => date created in database
           -update_date       => date last updated in database
           -pub_date          => date published in database

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($name, $id, $objid, $rank, $div, $dbh, $ncbitaxid, $commonname,
        $commonnames, $gcode, $mitocode, $createdate, $updatedate, $pubdate,
        $parent_id) = $self->_rearrange([qw(NAME ID OBJECT_ID RANK DIVISION DBH
                                            NCBI_TAXID COMMON_NAME COMMON_NAMES
                                            GENETIC_CODE MITO_GENETIC_CODE
                                            CREATE_DATE UPDATE_DATE PUB_DATE
                                            PARENT_ID)], @args);
    
    if (defined $id && (defined $ncbitaxid && $ncbitaxid ne $id || defined $objid && $objid ne $id)) {
        $self->warn("Only provide one of -id, -object_id or -ncbi_taxid, using $id\n");
    }
    elsif(!defined $id) { 
        $id = $objid || $ncbitaxid;
    }
    defined $id && $self->id($id);
    $self->{_ncbi_tax_id_provided} = 1 if $ncbitaxid;
    
    defined $rank && $self->rank($rank);
    defined $name && $self->node_name($name);
    
    my @common_names;
    if ($commonnames) {
        $self->throw("-common_names takes only an array reference") unless $commonnames
            && ref($commonnames) eq 'ARRAY';
        @common_names = @{$commonnames};
    }
    if ($commonname) {
        my %c_names = map { $_ => 1 } @common_names;
        unless (exists $c_names{$commonname}) {
            unshift(@common_names, $commonname);
        }
    }
    @common_names > 0 && $self->common_names(@common_names);
    
    defined $gcode      && $self->genetic_code($gcode);
    defined $mitocode   && $self->mitochondrial_genetic_code($mitocode);
    defined $createdate && $self->create_date($createdate);
    defined $updatedate && $self->update_date($updatedate);
    defined $pubdate    && $self->pub_date($pubdate);
    defined $div        && $self->division($div);
    defined $dbh        && $self->db_handle($dbh);
    
    # deprecated and will issue a warning when method called,
    # eventually to be removed completely as option
    defined $parent_id  && $self->parent_id($parent_id);
    
    # some things want to freeze/thaw Bio::Species objects, but
    # _root_cleanup_methods contains a CODE ref, delete it.
    delete $self->{_root_cleanup_methods};
    
    return $self;
}


=head1 Bio::IdentifiableI interface 

Also see L<Bio::IdentifiableI>

=head2 version

 Title   : version
 Usage   : $taxon->version($newval)
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
 Usage   : $taxon->authority($newval)
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
 Usage   : $taxon->namespace($newval)
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
 Usage   : $taxon->db_handle($newval)
 Function: Get/Set Bio::DB::Taxonomy Handle
 Returns : value of db_handle (a scalar) (Bio::DB::Taxonomy object)
 Args    : on set, new value (a scalar, optional) Bio::DB::Taxonomy object

Also see L<Bio::DB::Taxonomy>

=cut

sub db_handle {
    my $self = shift;
    if (@_) {
        my $db = shift;
        
        if (! ref($db) || ! $db->isa('Bio::DB::Taxonomy')) {
            $self->throw("Must provide a valid Bio::DB::Taxonomy object to db_handle()");
        }
        if (!$self->{'db_handle'} || ($self->{'db_handle'} && $self->{'db_handle'} ne $db)) {
            my $new_self = $self->_get_similar_taxon_from_db($self, $db);
            $self->_merge_taxa($new_self) if $new_self;
        }
        
        # NB: The Bio::DB::Taxonomy modules access this data member directly
        # to avoid calling this method and going infinite
        $self->{'db_handle'} = $db;
    }
    return $self->{'db_handle'};
}


=head2 rank

 Title   : rank
 Usage   : $taxon->rank($newval)
 Function: Get/set rank of this Taxon, 'species', 'genus', 'order', etc...
 Returns : value of rank (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub rank {
    my $self = shift;
    return $self->{'rank'} = shift if @_;
    return $self->{'rank'};
}


=head2 id

 Title   : id
 Usage   : $taxon->id($newval)
 Function: Get/Set id (NCBI Taxonomy ID in most cases); object_id() and
           ncbi_taxid() are synonyms of this method.
 Returns : id (a scalar)
 Args    : none to get, OR scalar to set

=cut

sub id {
    my $self = shift;
    return $self->SUPER::id(@_);
}

*object_id = \&id;


=head2 ncbi_taxid

 Title   : ncbi_taxid
 Usage   : $taxon->ncbi_taxid($newval)
 Function: Get/Set the NCBI Taxonomy ID; This actually sets the id() but only
           returns an id when ncbi_taxid has been explictely set with this
           method.
 Returns : id (a scalar)
 Args    : none to get, OR scalar to set

=cut

sub ncbi_taxid {
    my ($self, $id) = @_;
    
    if ($id) {
        $self->{_ncbi_tax_id_provided} = 1;
        return $self->SUPER::id($id);
    }
    
    if ($self->{_ncbi_tax_id_provided}) {
        return $self->SUPER::id;
    }
    return;
}


=head2 parent_id

 Title   : parent_id
 Usage   : $taxon->parent_id()
 Function: Get parent ID, (NCBI Taxonomy ID in most cases);
           parent_taxon_id() is a synonym of this method.
 Returns : value of parent_id (a scalar)
 Args    : none
 Status  : deprecated

=cut

sub parent_id {
    my $self = shift;
    if (@_) {
        $self->warn("You can no longer set the parent_id - use ancestor() instead");
    }
    my $ancestor = $self->ancestor() || return;
    return $ancestor->id;
}

*parent_taxon_id = \&parent_id;


=head2 genetic_code

 Title   : genetic_code
 Usage   : $taxon->genetic_code($newval)
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
 Usage   : $taxon->mitochondrial_genetic_code($newval)
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
 Usage   : $taxon->create_date($newval)
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
 Usage   : $taxon->update_date($newval)
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
 Usage   : $taxon->pub_date($newval)
 Function: Get/Set Date this node was published (in the database)
 Returns : value of pub_date (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub pub_date {
    my $self = shift;
    return $self->{'pub_date'} = shift if @_;
    return $self->{'pub_date'};
}


=head2 ancestor

 Title   : ancestor
 Usage   : my $ancestor_taxon = $taxon->ancestor()
 Function: Retrieve the ancestor taxon. Normally the database is asked what the
           ancestor is.

           If you manually set the ancestor (or you make a Bio::Tree::Tree with
           this object as an argument to new()), the database (if any) will not
           be used for the purposes of this method.

           To restore normal database behaviour, call ancestor(undef) (which
           would remove this object from the tree), or request this taxon again
           as a new Taxon object from the database.

 Returns : Bio::Taxon
 Args    : none

=cut

sub ancestor {
    my $self = shift;
    my $ancestor = $self->SUPER::ancestor(@_);
    if ($ancestor) {
        return $ancestor;
    }
    my $dbh = $self->db_handle;
    #*** could avoid the db lookup if we knew our current id was definitely
    #    information from the db...
    my $definitely_from_dbh = $self->_get_similar_taxon_from_db($self);
    return $dbh->ancestor($definitely_from_dbh);
}


=head2 get_Parent_Node

 Title   : get_Parent_Node
 Function: Synonym of ancestor()
 Status  : deprecated

=cut

sub get_Parent_Node {
    my $self = shift;
    $self->warn("get_Parent_Node is deprecated, use ancestor() instead");
    return $self->ancestor(@_);
}


=head2 each_Descendent

 Title   : each_Descendent
 Usage   : my @taxa = $taxon->each_Descendent();
 Function: Get all the descendents for this Taxon (but not their descendents,
           ie. not a recursive fetchall). get_Children_Nodes() is a synonym of
           this method.

           Note that this method never asks the database for the descendents;
           it will only return objects you have manually set with
           add_Descendent(), or where this was done for you by making a
           Bio::Tree::Tree with this object as an argument to new().

           To get the database descendents use
           $taxon->db_handle->each_Descendent($taxon).

 Returns : Array of Bio::Taxon objects
 Args    : optionally, when you have set your own descendents, the string
           "height", "creation", "alpha", "revalpha", or coderef to be used to
           sort the order of children nodes.

=cut


# implemented by Bio::Tree::Node

=head2 get_Children_Nodes

 Title   : get_Children_Nodes
 Function: Synonym of each_Descendent()
 Status  : deprecated

=cut

sub get_Children_Nodes {
    my $self = shift;
    $self->warn("get_Children_Nodes is deprecated, use each_Descendent() instead");
    return $self->each_Descendent(@_);
}


=head2 name

  Title:    name
  Usage:    $taxon->name('scientific', 'Homo sapiens');
            $taxon->name('common', 'human', 'man');
            my @names = @{$taxon->name('common')};
  Function: Get/set the names. node_name(), scientific_name() and common_names()
            are shorthands to name('scientific'), name('scientific') and
            name('common') respectively.
  Returns:  names (a array reference)
  Args:     Arg1 => the name_class. You can assign any text, but the words
                'scientific' and 'common' have the special meaning, as
                scientific name and common name, respectively. 'scientific' and
                'division' are treated specially, allowing only the first value
                in the Arg2 list to be set.
            Arg2 ... => list of names

=cut

sub name {
    my ($self, $name_class, @names) = @_;
    $self->throw('No name class specified') unless defined $name_class;
    
    if (@names) {
        if ($name_class =~ /scientific|division/i) {
            delete $self->{'_names_hash'}->{$name_class};
            @names = (shift(@names));
        }
        push @{$self->{'_names_hash'}->{$name_class}}, @names;
    }
    return $self->{'_names_hash'}->{$name_class} || return;
}


=head2 node_name

 Title   : node_name
 Usage   : $taxon->node_name($newval)
 Function: Get/set the name of this taxon (node), typically the scientific name
           of the taxon, eg. 'Primate' or 'Homo'; scientific_name() is a synonym
           of this method.
 Returns : value of node_name (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub node_name {
    my $self = shift;
    my @v = @{$self->name('scientific', @_) || []};
    return pop @v;
}

*scientific_name = \&node_name;


=head2 common_names

 Title   : common_names
 Usage   : $taxon->common_names($newval)
 Function: Get/add the other names of this taxon, typically the genbank common
           name and others, eg. 'Human' and 'man'. common_name() is a synonym
           of this method.
 Returns : array of names in list context, one of those names in scalar context
 Args    : on add, new list of names (scalars, optional)

=cut

sub common_names {
    my $self = shift;
    my @v = @{$self->name('common', @_) || []};
    return ( wantarray ) ? @v : pop @v;
}

*common_name = \&common_names;


=head2 division

 Title   : division
 Usage   : $taxon->division($newval)
 Function: Get/set the division this taxon belongs to, eg. 'Primates' or
           'Bacteria'.
 Returns : value of division (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub division {
    my $self = shift;
    my @v = @{$self->name('division',@_) || []};
    return pop @v;
}


# get a node from the database that is like the supplied node
sub _get_similar_taxon_from_db {
    #*** not really happy with this having to be called so much; there must be
    #    a better way...
    my ($self, $taxon, $db) = @_;
    $self->throw("Must supply a Bio::Taxon") unless ref($taxon) && $taxon->isa("Bio::Taxon");
    ($self->id || $self->node_name) || return;
    $db ||= $self->db_handle || return;
    if (!blessed($db) || !$db->isa('Bio::DB::Taxonomy')) {
        $self->throw("DB handle is not a Bio::DB::Taxonomy: got $db in node ".$self->node_name)
    }
    my $db_taxon = $db->get_taxon(-taxonid => $taxon->id) if $taxon->id;
    unless ($db_taxon) {
        my @try_ids = $db->get_taxonids($taxon->node_name) if $taxon->node_name;
        
        my $own_rank = $taxon->rank || 'no rank';
        foreach my $try_id (@try_ids) {
            my $try = $db->get_taxon(-taxonid => $try_id);
            my $try_rank = $try->rank || 'no rank';
            if ($own_rank eq 'no rank' || $try_rank eq 'no rank' || $own_rank eq $try_rank) {
                $db_taxon = $try;
                last;
            }
        }
    }
    
    return $db_taxon;
}


# merge data from supplied Taxon into self
sub _merge_taxa {
    my ($self, $taxon) = @_;
    $self->throw("Must supply a Bio::Taxon object") unless ref($taxon) && $taxon->isa('Bio::Taxon');
    return if ($taxon eq $self);
    
    foreach my $attrib (qw(scientific_name version authority namespace genetic_code mitochondrial_genetic_code create_date update_date pub_date division id)) {
        my $own = $self->$attrib();
        my $his = $taxon->$attrib();
        if (!$own && $his) {
            $self->$attrib($his);
        }
    }
    
    my $own = $self->rank || 'no rank';
    my $his = $taxon->rank || 'no rank';
    if ($own eq 'no rank' && $his ne 'no rank') {
        $self->rank($his);
    }
    
    my %own_cnames = map { $_ => 1 } $self->common_names;
    my %his_cnames = map { $_ => 1 } $taxon->common_names;
    foreach (keys %his_cnames) {
        unless (exists $own_cnames{$_}) {
            $self->common_names($_);
        }
    }
    
    #*** haven't merged the other things in names() hash, could do above much easier with direct access to object data
}


=head2 remove_Descendent

 Title   : remove_Descendent
 Usage   : $node->remove_Descedent($node_foo);
 Function: Removes a specific node from being a Descendent of this node
 Returns : nothing
 Args    : An array of Bio::Node::NodeI objects which have been previously
           passed to the add_Descendent call of this object.

=cut

sub remove_Descendent {
    # need to override this method from Bio::Tree::Node since it casually
    # throws away nodes if they don't branch
    my ($self,@nodes) = @_;
    my $c= 0;
    foreach my $n ( @nodes ) {
        if ($self->{'_desc'}->{$n->internal_id}) {
            $self->{_removing_descendent} = 1;
            $n->ancestor(undef);
            $self->{_removing_descendent} = 0;
            $self->{'_desc'}->{$n->internal_id}->ancestor(undef);
            delete $self->{'_desc'}->{$n->internal_id};
            $c++;
        }
    }
    return $c;
}


1;
