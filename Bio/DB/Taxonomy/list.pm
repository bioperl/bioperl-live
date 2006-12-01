# $Id$
#
# BioPerl module for Bio::DB::Taxonomy::list
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Sendu Bala
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Taxonomy::list - An implementation of Bio::DB::Taxonomy
that accepts lists of words to build a database

=head1 SYNOPSIS

  use Bio::DB::Taxonomy;

  my @names = ('Eukaryota', 'Mammalia', 'Homo', 'Homo sapiens');
  my @ranks = qw(superkingdom class genus species);
  my $db = new Bio::DB::Taxonomy(-source => 'list', -names => \@names,
                                                    -ranks => \@ranks);

  @names = ('Eukaryota', 'Mammalia', 'Mus', 'Mus musculus');
  $db->add_lineage(-names => \@names, -ranks => \@ranks);

=head1 DESCRIPTION

This is an implementation which uses supplied lists of words to create a
database from which you can extract Bio::Taxon objects.

=head1 TODO

It is possible this module could do something like store the data it builds
up to disc. Would that be useful?
At any rate, this is why the module is called 'list' and not 'in_memory' or
similar.

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

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::Taxonomy::list;
use strict;
use Bio::Taxon;

use base qw(Bio::DB::Taxonomy);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::DB::Taxonomy::list();
 Function: Builds a new Bio::DB::Taxonomy::list object 
 Returns : an instance of Bio::DB::Taxonomy::list
 Args    : optional, as per the add_lineage() method.

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    $self->{db} = {};
    $self->add_lineage(@args) if @args;
    
    return $self;
}

=head2 add_lineage

 Title   : add_lineage
 Usage   : $db->add_lineage(-names => \@names)
 Function: Add a lineage to the database, where the lineage is described by
           a list of scientific names in the order root->leaf. The rank of each
           name can optionally be described by supplying an additional list
           of rank names in the same order (eg. superkingdom->species).
 Returns : n/a
 Args    : -names => [] : array ref of scientific names, REQUIRED
           -ranks => [] : array ref of rank names, same order as above, OPTIONAL

=cut

sub add_lineage {
    my $self = shift;
    my ($names, $ranks) = $self->_rearrange([qw (NAMES RANKS)], @_);
    $self->throw("-names must be supplied and its value must be an array reference") unless $names && ref($names) eq 'ARRAY';
    my @names = @{$names};
    
    my @ranks;
    if ($ranks) {
        $self->throw("-ranks must be an array reference") unless ref($ranks) eq 'ARRAY';
        $self->throw("The -names and -ranks lists must be of equal length") unless @{$names} == @{$ranks};
        @ranks = @{$ranks};
    }
    else {
        for (0..$#names) {
            push(@ranks, 'no rank');
        }
    }
    
    # This is non-trivial because names are not guaranteed unique in a taxonomy,
    # and neither are name&rank combinations. Furthermore, different name&rank
    # combinations can actually refer to the same taxon, eg. when one time
    # 'Homo'&'genus' is supplied, while another time 'Homo'&'no rank'.
    #
    # name&rank&ancestor could well be unique (or good enough 99.9999% of the
    # time), but we have the added complication that lineages could sometimes be
    # supplied with differing numbers of taxa. Ideally we want realise that
    # the first of these two lineages shares all its nodes with the second:
    # ('Mammalia', 'Homo', 'Homo sapiens')
    # ('Mammalia', 'Hominidae', 'Homo', 'Homo sapiens')
    #
    # Clearly with limited information we can't do a perfect job, but we can try
    # and do a reasonable one.
    
    
    #...
    
    
    # All that said, let's just do the trivial implementation now and see how
    # bad it is! (assume names are unique, always have the same ancestor)
    
    my %names;
    foreach my $i (0..$#names) {
        my $name = $names[$i];
        $names{$name}++;
        if ($names{$name} > 1 && $name ne $names[$i - 1]) {
            $self->throw("The lineage '".join(', ', @names)."' had two non-consecutive nodes with the same name. Can't cope!");
        }
    }
    
    my $ancestor_node_id;
    my @node_ids;
    for my $i (0..$#names) {
        my $name = $names[$i];
        my $rank = $ranks[$i];
        
        # this is a new node with a new id if we haven't seen this name before,
        # or if the ancestor of this node in this supplied lineage has the
        # same name as this node (like '... Pinus, Pinus, Pinus densiflora').
        my $db_name = $name eq $names[$i - 1] ? $name.'_'.$rank : $name;
        if (! exists $self->{db}->{name_to_id}->{$db_name} || $name eq $names[$i - 1]) {
            my $next_num = ++$self->{db}->{node_ids};
            $self->{db}->{name_to_id}->{$db_name} = 'list'.$next_num; # so definitely not confused with ncbi taxonomy ids
        }
        my $node_id = $self->{db}->{name_to_id}->{$db_name};
        
        unless (exists $self->{db}->{node_data}->{$node_id}) {
            $self->{db}->{node_data}->{$node_id} = [($name, '')];
        }
        my $node_data = $self->{db}->{node_data}->{$node_id};
        
        if (!$node_data->[1] || ($node_data->[1] eq 'no rank' && $rank ne 'no rank')) {
            $node_data->[1] = $rank;
        }
        
        if ($ancestor_node_id) {
            if ($self->{db}->{ancestors}->{$node_id} && $self->{db}->{ancestors}->{$node_id} ne $ancestor_node_id) {
                $self->throw("This lineage (".join(', ', @names).") and a previously computed lineage share a node name but have different ancestries for that node. Can't cope!");
            }
            $self->{db}->{ancestors}->{$node_id} = $ancestor_node_id;
        }
        
        $ancestor_node_id = $node_id;
        push(@node_ids, $node_id);
    }
    
    # go through the lineage in reverse so we can remember the children
    my $child_id;
    foreach my $node_id (reverse @node_ids) {
        unless ($child_id) {
            $child_id = $node_id;
            next;
        }
        
        $self->{db}->{children}->{$node_id}->{$child_id} = 1;
    }
    
    #*** would prefer to use Digest::MD5 or similar for the hash keys, but this
    #    needs to work for everyone without hassle
    #
    #my $rank_list_id;
    #if (exists $DATABASE->{rank_lists}->{"@ranks"}) {
    #    $rank_list_id = ${$DATABASE->{rank_lists}->{"@ranks"}}[1];
    #}
    #else {
    #    $DATABASE->{rank_lists}->{"@ranks"} = [\@ranks, ++$DATABASE->{rank_id}];
    #    $DATABASE->{rank_id_to_list}->{$DATABASE->{rank_id}} = "@ranks";
    #}
    #
    ## have we already added this lineage?
    #if (exists $DATABASE->{name_lists}->{"@names"}) {
    #    foreach my $this_rank_id (@{${$DATABASE->{name_lists}->{"@names"}}[1]}) {
    #        return if $this_rank_id == $rank_list_id;
    #    }
    #}
    #else {
    #    $DATABASE->{name_lists}->{"@names"} = [\@names, []];
    #}
    #
    #push(@{${$DATABASE->{name_lists}->{"@names"}}[1]}, $rank_list_id);
    #
    #*** ideally we would also avoid the next step if new lineage is a branch
    #    of a longer existing lineage in the database
    #
    # compute the whole taxonomic tree from scratch, so that we aren't dependant
    # on the order lineages are added
    #$self->_compute_tree;
}

=head2 Bio::DB::Taxonomy Interface implementation

=cut

=head2 get_taxon

 Title   : get_taxon
 Usage   : my $taxon = $db->get_taxon(-taxonid => $taxonid)
 Function: Get a Bio::Taxon object from the database.
 Returns : Bio::Taxon object
 Args    : just a single value which is the database id, OR named args:
           -taxonid => taxonomy id (to query by taxonid; NB: these are not
                       NCBI taxonomy ids but 'list' pre-fixed ids unique to the
                       list database)
            OR
           -name    => string (to query by a taxonomy name)

=cut

sub get_taxon {
    my $self = shift;
    my ($taxonid, $name);
    
    if (@_ > 1) {
        ($taxonid, $name) = $self->_rearrange([qw(TAXONID NAME)],@_);
        if ($name) {
            ($taxonid, my @others) = $self->get_taxonids($name);
            $self->warn("There were multiple ids ($taxonid @others) matching '$name', using '$taxonid'") if @others > 0;
        }
    }
    else {
        $taxonid = shift;
    }
    
    my $node = $self->{db}->{node_data}->{$taxonid} || return;
    my ($sci_name, $rank) = @{$node};
    
    my $taxon = new Bio::Taxon(
                        -name         => $sci_name,
                        -object_id    => $taxonid, # since this is NOT a real ncbi taxid, set it as simply the object id
                        -rank         => $rank );
    # we can't use -dbh or the db_handle() method ourselves or we'll go
    # infinite on the merge attempt
    $taxon->{'db_handle'} = $self;
    
    $self->_handle_internal_id($taxon, 1);
    
    return $taxon;
}

*get_Taxonomy_Node = \&get_taxon;

=head2 get_taxonids

 Title   : get_taxonids
 Usage   : my @taxonids = $db->get_taxonids('Homo sapiens');
 Function: Searches for a taxonid (generated by the list module) based on a
           query string. Note that multiple taxonids can match to the same
           supplied name.
 Returns : array of integer ids in list context, one of these in scalar context
 Args    : string representing taxon's name

=cut

sub get_taxonids {
    my ($self, $query) = @_;
    my $id = $self->{db}->{name_to_id}->{$query} || return;
    return $id;
}

*get_taxonid = \&get_taxonids;

=head2 ancestor

 Title   : ancestor
 Usage   : my $ancestor_taxon = $db->ancestor($taxon)
 Function: Retrieve the full ancestor taxon of a supplied Taxon from the
           database. 
 Returns : Bio::Taxon
 Args    : Bio::Taxon (that was retrieved from this database)

=cut

sub ancestor {
    my ($self, $taxon) = @_;
    $taxon || return; # for bug 2092, or something similar to it at least: shouldn't need this!
    $self->throw("Must supply a Bio::Taxon") unless ref($taxon) && $taxon->isa('Bio::Taxon');
    $self->throw("The supplied Taxon must belong to this database") unless $taxon->db_handle && $taxon->db_handle eq $self;
    my $id = $taxon->id || $self->throw("The supplied Taxon is missing its id!");
    
    my $ancestor_id = $self->{db}->{ancestors}->{$id} || return;
    return $self->get_taxon($ancestor_id);
}

=head2 each_Descendent

 Title   : each_Descendent
 Usage   : my @taxa = $db->each_Descendent($taxon);
 Function: Get all the descendents of the supplied Taxon (but not their
           descendents, ie. not a recursive fetchall).
 Returns : Array of Bio::Taxon objects
 Args    : Bio::Taxon (that was retrieved from this database)

=cut

sub each_Descendent {
    my ($self, $taxon) = @_;
    $self->throw("Must supply a Bio::Taxon") unless ref($taxon) && $taxon->isa('Bio::Taxon');
    $self->throw("The supplied Taxon must belong to this database") unless $taxon->db_handle && $taxon->db_handle eq $self;
    my $id = $taxon->id || $self->throw("The supplied Taxon is missing its id!");
    
    my @children_ids = keys %{$self->{db}->{children}->{$id} || {}};
    my @children;
    foreach my $child_id (@children_ids) {
        push(@children, $self->get_taxon($child_id) || next);
    }
    
    return @children;
}

=head2 Helper methods 

=cut

# look at all the lineages we have and work out the overall tree
#sub _compute_tree {
#    my $self = shift;
#    #tba
#}

1;
