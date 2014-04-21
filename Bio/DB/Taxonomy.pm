#
# BioPerl module for Bio::DB::Taxonomy
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Taxonomy - Access to a taxonomy database

=head1 SYNOPSIS

  use Bio::DB::Taxonomy;
  my $db = Bio::DB::Taxonomy->new(-source => 'entrez');
  # use NCBI Entrez over HTTP
  my $taxonid = $db->get_taxonid('Homo sapiens');

  # get a taxon
  my $taxon = $db->get_taxon(-taxonid => $taxonid);

=head1 DESCRIPTION

This is a front end module for access to a taxonomy database.

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

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 CONTRIBUTORS

Sendu Bala: bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::Taxonomy;
use vars qw($DefaultSource $TAXON_IIDS);
use strict;
use Bio::Tree::Tree;

use base qw(Bio::Root::Root);

$DefaultSource = 'entrez';
$TAXON_IIDS = {};


=head2 new

 Title   : new
 Usage   : my $obj = Bio::DB::Taxonomy->new(-source => 'entrez');
 Function: Builds a new Bio::DB::Taxonomy object.
 Returns : an instance of Bio::DB::Taxonomy
 Args    : -source => which database source 'entrez' (NCBI taxonomy online),
                      'flatfile' (local NCBI taxonomy), 'greengenes' (local
                      GreenGenes taxonomy), 'silva' (local Silva taxonomy), or
                      'list' (Do-It-Yourself taxonomy)

=cut

sub new {
    my($class,@args) = @_;

    if( $class =~ /Bio::DB::Taxonomy::(\S+)/ ) {
        my ($self) = $class->SUPER::new(@args);
        $self->_initialize(@args);
        return $self;
    } else { 
        my %param = @args;
        @param{ map { lc $_ } keys %param } = values %param; # lowercase keys
        my $source = $param{'-source'} || $DefaultSource;

        $source = "\L$source"; # normalize capitalization to lower case

        # normalize capitalization
        return unless( $class->_load_tax_module($source) );
        return "Bio::DB::Taxonomy::$source"->new(@args);
    }
}


# empty for now
sub _initialize { }


=head2 get_num_taxa

 Title   : get_num_taxa
 Usage   : my $num = $db->get_num_taxa();
 Function: Get the number of taxa stored in the database.
 Returns : A number
 Args    : None

=cut

sub get_num_taxa {
    shift->throw_not_implemented();
}


=head2 get_taxon

 Title   : get_taxon
 Usage   : my $taxon = $db->get_taxon(-taxonid => $taxonid);
 Function: Get a Bio::Taxon object from the database.
 Returns : Bio::Taxon object
 Args    : just a single value which is the database id, OR named args:
           -taxonid => taxonomy id (to query by taxonid)
             OR
           -name    => string (to query by a taxonomy name: common name, 
                       scientific name, etc)

=cut

sub get_taxon {
    shift->throw_not_implemented();
}

*get_Taxonomy_Node = \&get_taxon;


=head2 get_taxonids

 Title   : get_taxonids
 Usage   : my @taxonids = $db->get_taxonids('Homo sapiens');
 Function: Searches for a taxonid (typically ncbi_taxon_id) based on a query
           string. Note that multiple taxonids can match to the same supplied
           name.
 Returns : array of integer ids in list context, one of these in scalar context
 Args    : string representing the taxon's name

=cut

sub get_taxonids {
    shift->throw_not_implemented();
}

*get_taxonid = \&get_taxonids;
*get_taxaid  = \&get_taxonids;


=head2 get_tree

 Title   : get_tree
 Usage   : my $tree = $db->get_tree(@species_names);
 Function: Generate a tree comprised of the full lineages of all the supplied
           species names. The nodes for the requested species are given
           name('supplied') values corresponding to the supplied name, such that
           they can be identified if the real species name in the database
           (stored under node_name()) is different. The nodes are also given an
           arbitrary branch length of 1.
 Returns : Bio::Tree::Tree
 Args    : A list of species names (strings) to include in the tree.

=cut

sub get_tree {
    my ($self, @species_names) = @_;
    
    # the full lineages of the species are merged into a single tree
    my $tree;
    for my $name (@species_names) {
        my @ids = $self->get_taxonids($name);
        if (not scalar @ids) {
            $self->throw("Could not find species $name in the taxonomy");
        }
        for my $id (@ids) {
            my $node = $self->get_taxon(-taxonid => $id);
            $node->name('supplied', $name);
            if ($tree) {
                $tree->merge_lineage($node);
            } else {
                $tree = Bio::Tree::Tree->new(-verbose => $self->verbose, -node => $node);
            }
        }
    }

    # add arbitrary branch length
    for my $node ($tree->get_nodes) {
        $node->branch_length(1);
    }
    
    return $tree;
}


=head2 ancestor

 Title   : ancestor
 Usage   : my $ancestor_taxon = $db->ancestor($taxon);
 Function: Retrieve the full ancestor taxon of a supplied Taxon from the
           database. 
 Returns : Bio::Taxon
 Args    : Bio::Taxon (that was retrieved from this database)

=cut

sub ancestor {
    shift->throw_not_implemented();
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
    shift->throw_not_implemented();
}


=head2 get_all_Descendents

 Title   : get_all_Descendents
 Usage   : my @taxa = $db->get_all_Descendents($taxon);
 Function: Like each_Descendent(), but do a recursive fetchall
 Returns : Array of Bio::Taxon objects
 Args    : Bio::Taxon (that was retrieved from this database)

=cut

sub get_all_Descendents {
    my ($self, $taxon) = @_;
    my @taxa;
    foreach my $desc_taxon ($self->each_Descendent($taxon)) {
        push @taxa, ($desc_taxon, $self->get_all_Descendents($desc_taxon));
    }
    return @taxa;
}


=head2 _load_tax_module

 Title   : _load_tax_module
 Usage   : *INTERNAL Bio::DB::Taxonomy stuff*
 Function: Loads up (like use) a module at run time on demand

=cut

sub _load_tax_module {
    my ($self, $source) = @_;
    my $module = "Bio::DB::Taxonomy::" . $source;
    my $ok;

    eval { $ok = $self->_load_module($module) };
    if ( $@ ) {
        print STDERR $@;
        print STDERR <<END;
$self: $source cannot be found
Exception $@
For more information about the Bio::DB::Taxonomy system please see
the Bio::DB::Taxonomy docs.  This includes ways of checking for 
formats at compile time, not run time.
END
  ;
    }
    return $ok;
}


=head2 _handle_internal_id

 Title   : _handle_internal_id
 Usage   : *INTERNAL Bio::DB::Taxonomy stuff*
 Function: Add an internal ID to a taxon object, ensuring that the taxon gets
           the same internal ID, regardless of which database it is retrieved
           from.
 Returns : The assigned internal ID
 Args    : * A Bio::Taxon
           * An optional boolean to decide whether or not to try and do the job
             using scientific name & rank in addition to taxon ID. This is
             useful if your IDs are not comparable to that of other databases,
             e.g. if they are arbitrary, as in the case of Bio::DB::Taxonomy::list.
             CAVEAT: will handle ambiguous names within a database fine, but not
             across multiple databases.

=cut

sub _handle_internal_id {
    my ($self, $taxon, $try_name) = @_;
    $self->throw("Must supply a Bio::Taxon") unless ref($taxon) && $taxon->isa('Bio::Taxon');

    my $taxid = $taxon->id              || return;
    my $name  = $taxon->scientific_name || '';
    my $rank  = $taxon->rank            || 'no rank';
    my $dbh   = $try_name ? $taxon->db_handle : 'any';

    my $iid = $TAXON_IIDS->{taxids}->{$dbh}->{$taxid};
    if ( (not defined $iid) && $try_name && $name && exists $TAXON_IIDS->{names}->{$name}) {
        # Search for a suitable IID based on species name and ranks
        my %test_ranks = map {$_ => undef} ($rank, 'no rank');
        SEARCH: while (my ($test_rank, undef) = each %test_ranks) {
            # Search at the specified rank first, then with 'no rank'
            while ( my ($test_iid, $test_info) = each %{$TAXON_IIDS->{names}->{$name}->{$rank}} ) {
                while (my ($test_db, $test_taxid) = each %$test_info) {
                    if ( ($test_db eq $dbh) && not($test_taxid eq $taxid) ) {
                        # Taxa are different (same database, different taxid)
                        next;
                    }
                    # IID is acceptable since taxa are from different databases,
                    # or from the same database but have the same taxid
                    $iid = $test_iid;
                    $TAXON_IIDS->{taxids}->{$dbh}->{$taxid} = $iid;
                    last SEARCH;
                }
            }
        }
    }

    if (defined $iid) {
        # Assign Bio::DB::Taxonomy IID with risky Bio::Tree::Node internal method
        $taxon->_creation_id($iid);
    } else {
        # Register new IID in Bio::DB::Taxonomy
        $iid = $taxon->internal_id;
        $TAXON_IIDS->{taxids}->{$dbh}->{$taxid} = $iid;
        if ($name) {
            $TAXON_IIDS->{names}->{$name}->{$rank}->{$iid}->{$taxon->db_handle} = $taxid
        }
    }

    return $iid;

}


1;
