#
# BioPerl module for Bio::DB::Taxonomy::list
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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
  my $db = Bio::DB::Taxonomy->new(-source => 'list', -names => \@names,
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

  https://redmine.open-bio.org/projects/bioperl/

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
 Usage   : my $obj = Bio::DB::Taxonomy::list->new();
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
    # supplied with differing numbers of taxa. Ideally we want to realise that
    # the first of these two lineages shares all its nodes with the second:
    # ('Mammalia', 'Homo', 'Homo sapiens')
    # ('Mammalia', 'Hominidae', 'Homo', 'Homo sapiens')
    #
    # Clearly with limited information we can't do a perfect job, but we can try
    # and do a reasonable one.
    
    
    #...
    
    
    # All that said, let's just do the trivial implementation now and see how
    # bad it is! (assumes ranks are unique except for 'no rank')
    
    
    my $first_lineage = $self->{db}->{node_ids} ? 0 : 1;
    
    my $ancestor_node_id;
    my @node_ids;
    for my $i (0..$#names) {
        my $name = $names[$i];
        my $rank = $ranks[$i];
        
        # This is a new node with a new id if we haven't seen this name before.
        # It's also always a new node if this is the first lineage going into
        # the db.
        #
        # We need to handle, however, situations in the future where we try to
        # merge in a new lineage but we have non-unique names in the lineage
        # and possible missing classes in some lineages
        # (eg.
        # '... Anophelinae, Anopheles, Anopheles, Angusticorn, Anopheles...'
        # merged with
        # '... Anophelinae, Anopheles, Angusticorn, Anopheles...'),
        # but still need the 'tree' to be correct
        
        my $is_new = 0;
        if ($first_lineage || ! exists $self->{db}->{name_to_id}->{$name}) {
            $is_new = 1;
        }
        
        my $node_id;
        unless ($is_new) {
            my @same_named = @{$self->{db}->{name_to_id}->{$name}};
            
            # look for the node that is consistent with this lineage
            SAME_NAMED: foreach my $s_id (@same_named) {
                my $this_ancestor_id;
                if ($ancestor_node_id) {
                    $this_ancestor_id = $self->{db}->{ancestors}->{$s_id};
                    if ($ancestor_node_id eq $this_ancestor_id) {
                        $node_id = $s_id;
                        last SAME_NAMED;
                    }
                }
                
                if ($names[$i + 1]) {
                    my $my_child_name = $names[$i + 1];
                    my @children_ids = keys %{$self->{db}->{children}->{$s_id} || {}};
                    foreach my $c_id (@children_ids) {
                        my $this_child_name = $self->{db}->{node_data}->{$c_id}->[0];
                        if ($my_child_name eq $this_child_name) {
                            
                            if ($ancestor_node_id) {
                                my @s_ancestors;
                                while ($this_ancestor_id = $self->{db}->{ancestors}->{$this_ancestor_id}) {
                                    if ($ancestor_node_id eq $this_ancestor_id) {
                                        $node_id = $s_id;
                                        $ancestor_node_id = $self->{db}->{ancestors}->{$s_id};
                                        push(@node_ids, @s_ancestors, $ancestor_node_id);
                                        last SAME_NAMED;
                                    }
                                    unshift(@s_ancestors, $this_ancestor_id);
                                }
                            }
                            else {
                                #$self->warn("This new lineage (@names) doesn't start at the same root as the existing lineages.".
                                #            "\nI'm assuming '$name' corresponds to node $s_id");
                                $node_id = $s_id;
                                last SAME_NAMED;
                            }
                        }
                    }
                }
            }
            
            $node_id || $is_new++;
        }
        
        if ($is_new) {
            my $next_num = ++$self->{db}->{node_ids};
            # 'list' so definitely not confused with ncbi taxonomy ids
            $node_id = 'list'.$next_num;
            push(@{$self->{db}->{name_to_id}->{$name}}, $node_id);
        }
        
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
        $child_id = $node_id;
    }
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
    
    my $taxon = Bio::Taxon->new(
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
    return @{$self->{db}->{name_to_id}->{$query} || []};
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

1;
