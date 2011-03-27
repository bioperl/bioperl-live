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

  https://redmine.open-bio.org/projects/bioperl/

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
 Args    : -source => which database source 'entrez' or 'flatfile' or 'list'

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

      $source = "\L$source";	# normalize capitalization to lower case

      # normalize capitalization
      return unless( $class->_load_tax_module($source) );
      return "Bio::DB::Taxonomy::$source"->new(@args);
  }
}

# empty for now
sub _initialize { }

=head2 get_taxon

 Title   : get_taxon
 Usage   : my $taxon = $db->get_taxon(-taxonid => $taxonid)
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
 Args    : string representing taxon's name

=cut

sub get_taxonids {
    shift->throw_not_implemented();
}

*get_taxonid = \&get_taxonids;
*get_taxaid = \&get_taxonids;

=head2 get_tree

 Title   : get_tree
 Usage   : my $tree = $db->get_tree(@species_names)
 Function: Generate a tree comprised of the full lineages of all the supplied
           species names. The nodes for the requested species are given
           name('supplied') values corresponding to the supplied name, such that
           they can be identified if the real species name in the database
           (stored under node_name()) is different.
 Returns : Bio::Tree::Tree
 Args    : a list of species names (strings)

=cut

sub get_tree {
    my ($self, @species_names) = @_;
    
    # the full lineages of the species are merged into a single tree
    my $tree;
    foreach my $name (@species_names) {
        my $ncbi_id = $self->get_taxonid($name);
        if ($ncbi_id) {
            my $node = $self->get_taxon(-taxonid => $ncbi_id);
            $node->name('supplied', $name);
            
            if ($tree) {
                $tree->merge_lineage($node);
            }
            else {
                $tree = Bio::Tree::Tree->new(-verbose => $self->verbose, -node => $node);
            }
        }
        else {
            $self->throw("No taxonomy database node for species ".$name);
        }
    }
    
    return $tree;
}

=head2 ancestor

 Title   : ancestor
 Usage   : my $ancestor_taxon = $db->ancestor($taxon)
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
 Function: Tries to ensure that when a taxon is requested from any database,
           the Taxon object returned will have the same internal id regardless
           of database.
 Args    : Bio::Taxon, and optionally true value to try and do the job using
           scientific name & rank if your ids aren't comparable to other dbs.

=cut

sub _handle_internal_id {
    my ($self, $taxon, $try_name) = @_;
    $self->throw("Must supply a Bio::Taxon") unless ref($taxon) && $taxon->isa('Bio::Taxon');
    my $taxid = $taxon->id || return;
    my $sci_name = $taxon->scientific_name || '';
    my $rank = $taxon->rank || 'no rank';
    
    if ($try_name && $sci_name && defined $TAXON_IIDS->{names}->{$sci_name}) {
        if (defined $TAXON_IIDS->{names}->{$sci_name}->{$rank}) {
            $TAXON_IIDS->{taxids}->{$taxid} = $TAXON_IIDS->{names}->{$sci_name}->{$rank};
        }
        elsif ($rank eq 'no rank') {
            # pick the internal id of one named rank taxa at random
            my ($iid) = values %{$TAXON_IIDS->{names}->{$sci_name}};
            $TAXON_IIDS->{taxids}->{$taxid} = $iid;
        }
    }
    
    if (defined $TAXON_IIDS->{taxids}->{$taxid}) {
        # a little dangerous to use this internal method of Bio::Tree::Node;
        # but it is how internal_id() is set
        $taxon->_creation_id($TAXON_IIDS->{taxids}->{$taxid});
    }
    else {
        $TAXON_IIDS->{taxids}->{$taxid} = $taxon->internal_id;
        $TAXON_IIDS->{names}->{$sci_name}->{$rank} = $taxon->internal_id if $sci_name;
    }
}

1;
