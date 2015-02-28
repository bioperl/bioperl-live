#
# BioPerl module for Bio::DB::Taxonomy::flatfile
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Chris Fields <cjfields-at-cpan-dot-org>
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Taxonomy::sqlite - SQLite-based implementation of Bio::DB::Taxonomy::flatfile

=head1 SYNOPSIS

  use Bio::DB::Taxonomy;

  my $db = Bio::DB::Taxonomy->new(-source    => 'sqlite' ,
                                  -nodesfile => 'nodes.dmp',
                                  -namesfile => 'names.dmp');

=head1 DESCRIPTION

This is an implementation of Bio::DB::Taxonomy which stores and accesses the
NCBI taxonomy using a simple SQLite3 database stored locally on disk.

The required database files, nodes.dmp and names.dmp can be obtained from
ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

=head1 TODO

Beyond completing the implementation and optimization, this will
likely be rolled into a more flexible backend at some future point.

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

package Bio::DB::Taxonomy::sqlite;

use strict;
use DB_File;
use Bio::Taxon;
use File::Spec::Functions;
use DBI;

use constant SEPARATOR => ':';

our $DEFAULT_INDEX_DIR     = $Bio::Root::IO::TEMPDIR;    # /tmp
our $DEFAULT_DB_NAME       = 'taxonomy.sqlite';
our $DEFAULT_NODE_INDEX    = 'nodes';
our $DEFAULT_NAME2ID_INDEX = 'names2id';
our $DEFAULT_ID2NAME_INDEX = 'id2names';
our $DEFAULT_PARENT_INDEX  = 'parents';

our @DIVISIONS = (
    [qw(BCT Bacteria)],
    [qw(INV Invertebrates)],
    [qw(MAM Mammals)],
    [qw(PHG Phages)],
    [qw(PLN Plants)],                                    # (and fungi)
    [qw(PRI Primates)],
    [qw(ROD Rodents)],
    [qw(SYN Synthetic)],
    [qw(UNA Unassigned)],
    [qw(VRL Viruses)],
    [qw(VRT Vertebrates)],
    [qw(ENV 'Environmental samples')]
);

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
    my ( $dir, $nodesfile, $namesfile, $db, $force ) =
      $self->_rearrange( [qw(DIRECTORY NODESFILE NAMESFILE DB FORCE)], @args );

    $self->index_directory( $dir || $DEFAULT_INDEX_DIR );
    
    $self->db_name( $db || $DEFAULT_DB_NAME );
    
    if ($nodesfile) {
        $self->_build_index( $nodesfile, $namesfile, $force );
    }

    $self->_db_connect;
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

sub get_num_taxa {
    my ($self) = @_;
    
    my $ct = $self->_dbh_do(<<SQL);
    SELECT COUNT(*) FROM taxon
SQL
    
    return shift @{$ct};
}

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
    my ($self) = shift;
    my ( $taxonid, $name );

    if ( @_ > 1 ) {
        ( $taxonid, $name ) = $self->_rearrange( [qw(TAXONID NAME)], @_ );
        if ($name) {
            ( $taxonid, my @others ) = $self->get_taxonids($name);
            $self->warn(
"There were multiple ids ($taxonid @others) matching '$name', using '$taxonid'"
            ) if @others > 0;
        }
    }
    else {
        $taxonid = shift;
    }

    return unless $taxonid;

    $taxonid =~ /^\d+$/ || return;
    my $node = $self->{'_nodes'}->[$taxonid] || return;
    length($node) || return;
    my ( $taxid, undef, $rank, $code, $divid, $gen_code, $mito ) =
      split( SEPARATOR, $node );
    last unless defined $taxid;
    my ($taxon_names) = $self->{'_id2name'}->[$taxid];
    my ( $sci_name, @common_names ) = split( SEPARATOR, $taxon_names );

    my $taxon = Bio::Taxon->new(
        -name         => $sci_name,
        -common_names => [@common_names],
        -ncbi_taxid =>
          $taxid,    # since this is a real ncbi taxid, explicitly set it as one
        -rank              => $rank,
        -division          => $DIVISIONS[$divid]->[1],
        -genetic_code      => $gen_code,
        -mito_genetic_code => $mito
    );

    # we can't use -dbh or the db_handle() method ourselves or we'll go
    # infinite on the merge attempt
    $taxon->{'db_handle'} = $self;

    $self->_handle_internal_id($taxon);

    return $taxon;
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
    my ( $self, $query ) = @_;
    my $ids = $self->{'_name2id'}->{ lc($query) };
    unless ($ids) {
        if ( $query =~ /_/ ) {

            # try again converting underscores to spaces
            $query =~ s/_/ /g;
            $ids = $self->{'_name2id'}->{ lc($query) };
        }
        $ids || return;
    }
    my @ids = split( SEPARATOR, $ids );
    return wantarray() ? @ids : shift @ids;
}

*get_taxonid = \&get_taxonids;

=head2 get_Children_Taxids

 Title   : get_Children_Taxids
 Usage   : my @childrenids = $db->get_Children_Taxids 
 Function: Get the ids of the children of a node in the taxonomy
 Returns : Array of Ids
 Args    : Bio::Taxon or a taxon_id
 Status  : deprecated (use each_Descendent())

=cut

sub get_Children_Taxids {
    my ( $self, $node ) = @_;
    $self->warn(
        "get_Children_Taxids is deprecated, use each_Descendent instead");
    my $id;
    if ( ref($node) ) {
        if ( $node->can('object_id') ) {
            $id = $node->object_id;
        }
        elsif ( $node->can('ncbi_taxid') ) {
            $id = $node->ncbi_taxid;
        }
        else {
            $self->warn(
                "Don't know how to extract a taxon id from the object of type "
                  . ref($node)
                  . "\n" );
            return;
        }
    }
    else { $id = $node }
    my @vals = $self->{'_parentbtree'}->get_dup($id);
    return @vals;
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
    my ( $self, $taxon ) = @_;
    $self->throw("Must supply a Bio::Taxon")
      unless ref($taxon) && $taxon->isa('Bio::Taxon');
    $self->throw("The supplied Taxon must belong to this database")
      unless $taxon->db_handle && $taxon->db_handle eq $self;
    my $id =
      $taxon->id || $self->throw("The supplied Taxon is missing its id!");

    my $node = $self->{'_nodes'}->[$id];
    if ( length($node) ) {
        my ( undef, $parent_id ) = split( SEPARATOR, $node );
        $parent_id || return;
        $parent_id eq $id && return;    # one of the roots
        return $self->get_taxon($parent_id);
    }
    return;
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
    my ( $self, $taxon ) = @_;
    $self->throw("Must supply a Bio::Taxon")
      unless ref($taxon) && $taxon->isa('Bio::Taxon');
    $self->throw("The supplied Taxon must belong to this database")
      unless $taxon->db_handle && $taxon->db_handle eq $self;
    my $id =
      $taxon->id || $self->throw("The supplied Taxon is missing its id!");

    my @desc_ids = $self->{'_parentbtree'}->get_dup($id);
    my @descs;
    foreach my $desc_id (@desc_ids) {
        push( @descs, $self->get_taxon($desc_id) || next );
    }
    return @descs;
}

=head2 Helper methods 

=cut

=head2 index_directory

 Title   : index_directory
 Funtion : Get/set the location that index files are stored. (this module
           will index the supplied database)
 Usage   : $obj->index_directory($newval)
 Returns : value of index_directory (a scalar)
 Args    : on set, new value (a scalar or undef, optional)
 Note    : kept for backwards compatibility with older DB_File implementation

=cut

sub index_directory {
    my $self = shift;
    return $self->{'index_directory'} = shift if @_;
    return $self->{'index_directory'};
}

=head2 db_name

 Title   : db_name
 Funtion : Get/set the name of the SQLite3 database where data is stored
 Usage   : $obj->index_directory($newval)
 Returns : value of index_directory (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

# TODO: this may need some disambiguation w/ index_directory above; for now we
# assume this doesn't have a full path name (though I see no reason why this
# shouldn't allow that)

sub db_name {
    my $self = shift;
    return $self->{'db_name'} = shift if @_;
    return $self->{'db_name'};
}

# internal method which does the indexing
sub _build_index {
    my ( $self, $nodesfile, $namesfile, $force ) = @_;

    # TODO: need to disambiguate using index_directory here since we only have
    # one file.  Mayeb ignore it in favor of having full path for db_name?
    my ($dir, $db_name) = ($self->index_directory, $self->db_name);
    
    # TODO: we're ignoring index_directory for now, may add support for this
    # down the way
    my $dbh = DBI->connect("dbi:SQLite:dbname=$db_name","","") or die $!;

    $dbh->do('PRAGMA synchronous = 0');      # Non transaction safe!!!
    $dbh->do('PRAGMA cache_size = 100000');  # TODO: make this so user can define the cache size?

    if (! -e $db_name || $force) {
        
        $self->debug("Loading taxon table data\n");
        $self->_init_db($dbh);
        open my $NODES, '<', $nodesfile
            or $self->throw("Could not read node file '$nodesfile': $!");
    
        my $sth = $dbh->prepare_cached(<<SQL);
    INSERT INTO taxon (taxon_id, parent_id, rank, code, division_id, gencode_id, mito_id) VALUES (?,?,?,?,?,?,?)
SQL
        $dbh->do("BEGIN");
        while (<$NODES>) {
            next if /^\s*$/;
            chomp;
            my ($taxid,$parent,$rank,$code,$divid,undef,$gen_code,undef,$mito) = split(/\t\|\t/,$_);
            next if $taxid == 1;
            if ($parent == 1) {
                $parent = $taxid;
            }
            
            $sth->execute($taxid, $parent, $rank, $code, $divid, $gen_code, $mito) or die $sth->errstr.": TaxID $taxid";
        }
        $dbh->do("COMMIT");
        
        # TODO:index parent_id
        close $NODES;
        
        $self->debug("Loading name table data\n");
        open my $NAMES, '<', $namesfile
            or $self->throw("Could not read names file '$namesfile': $!");
    
        my $sth = $dbh->prepare_cached(<<SQL) or $self->throw($dbh->errstr);
    INSERT INTO names (taxon_id, name, uniq_name, class) VALUES (?,?,?,?)
SQL
        $dbh->do("BEGIN");
        while (<$NAMES>) {
            next if /^$/;
            chomp;
            my ($taxid, $name, $unique_name, $class) = split(/\t\|\t/,$_);
            # don't include the fake root node 'root' or 'all' with id 1
            next if $taxid == 1;
    
            $class =~ s/\s+\|\s*$//;
            my $lc_name = lc($name);
            my $orig_name = $name;
            
            # NOTE: From original implementation
            # ----
            # unique names aren't always in the correct column, sometimes they
            # are uniqued by adding bracketed rank names to the normal name;
            # store the uniqued version then fix the name for normal use
            # ----
            # TODO: the above no longer seems to be the case (e.g. all 'unique
            # names' in that column appear to be truly unique as of 2015/2/28).
            # Can add a column constraint to catch this, but I think (beyond
            # much data munging) we can now skip this step and simply load the
            # table in.

            # TODO: maybe FTS4 table for full-text search?
            
            # As of 2015/2/28 it appears this conditional is never true
            if ($lc_name =~ /\(class\)$/) { # it seems that only rank of class is ever used in this situation
                $name =~ s/\s+\(class\)$//;
                $lc_name = lc($name);
            }
            
            $sth->execute($taxid, $lc_name, lc($unique_name), $class) or $self->throw($sth->errstr);
        }
        $dbh->do("COMMIT");
        $dbh->do("PRAGMA foreign_keys = ON");
        close $NAMES;
        $self->{dbh} = $dbh;
        $self->{'_initialized'} = 1;
    }
    1;
}

# connect the internal db handle
sub _db_connect {
    my $self = shift;
    return if $self->{'_initialized'};

    my ($dir, $db_name) = ($self->index_directory, $self->db_name);
    
    # TODO: we're ignoring index_directory for now, may add support for this
    # down the way
    $self->{dbh} = DBI->connect("dbi:SQLite:dbname=$db_name","","") or die $!;
    $self->{'_initialized'} = 1;
}

sub _init_db {
    my ($self, $dbh) = @_;
    my $schema = $self->taxon_schema();
    # TODO: set up handler parameters here
    for my $table (sort keys %$schema) {
        $dbh->do("DROP TABLE IF EXISTS $table") or $self->throw($dbh->errstr);
        $dbh->do("CREATE TABLE $table ".$schema->{$table}) or $self->throw($dbh->errstr);
    }
    1;
}

sub _dbh_do {
    my ($self, $sql) = @_;
    # TODO: more sanity checks
    my $rows = $self->{dbh}->do($sql) or $self->throw( $self->{dbh}->errstr );
    return $rows;
}

# TODO: check data size, this is a ballpark estimate (could be reduced)
sub taxon_schema {
    my $self   = shift;
    return {
    taxon   => <<SCHEMA,
    (
        taxon_id        INTEGER PRIMARY KEY NOT NULL,
        parent_id       INTEGER,
        left_id         INTEGER,
        right_id        INTEGER,
        name            VARCHAR(100),
        rank            VARCHAR(25),
        code            VARCHAR(5),
        division_id     INTEGER,
        gencode_id      INTEGER,
        mito_id         INTEGER,
        FOREIGN KEY(parent_id) REFERENCES taxon(taxon_id)
    )
SCHEMA
    
    names   => <<SCHEMA,
    (
        name_id         INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
        taxon_id        INTEGER,
        name            VARCHAR(50),
        uniq_name       VARCHAR(50),
        class           VARCHAR(25),
        FOREIGN KEY(taxon_id) REFERENCES taxon(taxon_id)
    )
SCHEMA
    };
}

sub DESTROY {
    my $self = shift;
    
    undef $self->{_dbh};
    
    my $default_temp = quotemeta $DEFAULT_INDEX_DIR;
    if ($self->{index_directory} =~ m/^$default_temp/) {
        unlink catfile($self->index_directory(),$self->db_name());
    }
}

1;

