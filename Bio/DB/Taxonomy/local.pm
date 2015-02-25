#
# BioPerl module for Bio::DB::Taxonomy::local
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Chris Fields <cjfields-at-cpan-dot-org>
#
# Copyright Chris Fields, munged heavily from Bio::DB::Taxonomy::flatfile by
#           Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Taxonomy::local - Create and use a taxonomy from local files

=head1 SYNOPSIS

  use Bio::DB::Taxonomy;

  my $db = Bio::DB::Taxonomy->new(-source    => 'local' ,  
                                  -adaptor   => 'DB_File',
                                  -taxonomy  => 'NCBI',      # NCBI taxonomy (default)
                                  -create    => 1,           # create database
                                  -adaptor_args => {         # adaptor-spec. args
                                    -directory  => 'NCBI-Taxonomy',
                                    -nodesfile  => 'nodes.dmp',
                                    -namesfile  => 'names.dmp',
                                  },
                                  );

=head1 DESCRIPTION

This is an implementation of Bio::DB::Taxonomy which stores and accesses 
taxonomy information stored locally.  This could be via a DBI, memory, or
locally indexed flat files.

The intent is a new implementation that abstracts the type of data being loaded
and the API for loading and retrieval of information into separate classes. In
order to implement support for a specific taxonomy format under any backend, one
should implement a format-specific Loader (which parses the raw data into a
usable common format), and a database interface. It is mainly meant to act as a
replacement for Bio::DB::Taxonomy::flatfile that's flexible (allowing for
alternative backends such as SQLite or Neo4J, and possibly allowing different
taxonomies, such as Silva or Greengenes).

=head1 TODO

* Initial implementation focused on NCBI and DB_File (rewrite of flatfile)

* Reimplementation focused on NCBI and SQLite or similar

* Make a short-hand DSN-like call available

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

=head1 AUTHOR - Chris Fields

Email cjfields-at-cpan-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::Taxonomy::local;

use strict;

use base qw(Bio::DB::Taxonomy);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::DB::Taxonomy::flatfile->new();
 Function: Builds a new Bio::DB::Taxonomy::flatfile object 
 Returns : an instance of Bio::DB::Taxonomy::flatfile
 Args    : -directory => name of directory where index files should be created
           -nodesfile => name of file containing nodes (nodes.dmp from NCBI)
           -namesfile => name of the file containing names(names.dmp from NCBI)
           -force     => 1 to replace current indexes even if they exist

=cut

sub new {
    my ( $class, @args ) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    my ($adaptor, $tax, $create, $data_files) =
      $self->_rearrange( [qw(ADAPTOR TAXONOMY CREATE DATA_FILES)], @args );
    
    return $self;
}

=head2 Bio::DB::Taxonomy interface implementation

=head2 get_num_taxa

 Title   : get_num_taxa
 Usage   : my $num = $db->get_num_taxa();
 Function: Get the number of taxa stored in the database.
 Returns : A number
 Args    : None

=cut

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

=head2 get_taxonids

 Title   : get_taxonids
 Usage   : my @taxonids = $db->get_taxonids('Homo sapiens');
 Function: Searches for a taxonid (typically ncbi_taxon_id) based on a query
           string. Note that multiple taxonids can match to the same supplied
           name.
 Returns : array of integer ids in list context, one of these in scalar context
 Args    : string representing taxon's name

=cut

=head2 get_Children_Taxids

 Title   : get_Children_Taxids
 Usage   : my @childrenids = $db->get_Children_Taxids 
 Function: Get the ids of the children of a node in the taxonomy
 Returns : Array of Ids
 Args    : Bio::Taxon or a taxon_id
 Status  : deprecated (use each_Descendent())

=cut

=head2 ancestor

 Title   : ancestor
 Usage   : my $ancestor_taxon = $db->ancestor($taxon)
 Function: Retrieve the full ancestor taxon of a supplied Taxon from the
           database. 
 Returns : Bio::Taxon
 Args    : Bio::Taxon (that was retrieved from this database)

=cut

=head2 each_Descendent

 Title   : each_Descendent
 Usage   : my @taxa = $db->each_Descendent($taxon);
 Function: Get all the descendents of the supplied Taxon (but not their
           descendents, ie. not a recursive fetchall).
 Returns : Array of Bio::Taxon objects
 Args    : Bio::Taxon (that was retrieved from this database)

=cut

=head2 Helper methods 

=cut

# internal method which does the indexing
# connect the internal db handle

=head2 index_directory

 Title   : index_directory
 Funtion : Get/set the location that index files are stored. (this module
           will index the supplied database)
 Usage   : $obj->index_directory($newval)
 Returns : value of index_directory (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

1;
