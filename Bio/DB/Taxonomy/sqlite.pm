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

  my $db = Bio::DB::Taxonomy->new(-source    => 'sqlite',
                                  -db        => 'mytax.db'  # default 'taxonomy.sqlite'
                                  -nodesfile => 'nodes.dmp',
                                  -namesfile => 'names.dmp');

=head1 DESCRIPTION

This is an implementation of Bio::DB::Taxonomy which stores and accesses the
NCBI taxonomy using a simple SQLite3 database stored locally on disk.

With this implementation, one can do the same basic searches as with the 'flatfile'
database.  A test lookup of 1000 NCBI TaxIDs with full lineage information took
about 2 seconds on my older MacBook Pro laptop with an on-disk implementation.  

A few key differences:

=over 4

=item * You can use typical SQL syntax to run a query search; for instance, if you want you can run:

   @ids = sort $db->get_taxonids('Chloroflexi%');

=item * In-memory database is allowed

  my $db = Bio::DB::Taxonomy->new(-source    => 'sqlite',
                                  -db        => ':memory:',
                                  -nodesfile => 'nodes.dmp',
                                  -namesfile => 'names.dmp');

=back

The required database files, nodes.dmp and names.dmp can be obtained from
ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

=head1 TODO

=over 4

=item * Small optimizations, such as optimizing name lookups

=item * Possibly use L<recursive CTE|http://www.sqlite.org/lang_with.html> to do lineage lookups 

=item * Clean up SQL (still kind of a mess right now)

=item * Check compat. with other NCBI-specific L<Bio::DB::Taxonomy> implementations

=item * Plan out feasibility of allowing other backends (Neo4J, other DBI, etc)

=item * Optionally calculate left/right ID values for TaxID nodes

=back

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

use 5.010;
use strict;
use DB_File;
use Bio::Taxon;
use File::Spec::Functions;
use Data::Dumper;
use DBI;

use constant SEPARATOR => ':';

our $DEFAULT_INDEX_DIR     = $Bio::Root::IO::TEMPDIR;    # /tmp
our $DEFAULT_CACHE_SIZE    = 0;    # /tmp
our $DEFAULT_DB_NAME       = 'taxonomy.sqlite';

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

# TODO: get rid of globals!
sub new {
    my ( $class, @args ) = @_;

    my $self = $class->SUPER::new(@args);
    
    my ( $dir, $nodesfile, $namesfile, $db, $force, $cs ) =
      $self->_rearrange( [qw(DIRECTORY NODESFILE NAMESFILE DB FORCE CACHE_SIZE)], @args );

    $self->index_directory( $dir || $DEFAULT_INDEX_DIR );
    
    $self->db_name( $db || $DEFAULT_DB_NAME );
    
    $self->cache_size($cs // $DEFAULT_CACHE_SIZE);
    
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
    
    my $ct = $self->_dbh_fetch(<<SQL);
    SELECT COUNT(*) FROM taxon
SQL
    
    return @{$ct}[0];
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

    $taxonid =~ /^\d+$/ || $self->throw("TaxID must be integer, got [$taxonid]");
    
    my ( $parent_id, $rank, $code, $divid, $gen_code, $mito, $nm, $uniq, $class );
    # single join or two calls?
    my $sth = $self->_prepare_cached(<<SQL);
    SELECT tax.parent_id, tax.rank, tax.code, tax.division_id, tax.gencode_id, tax.mito_id, names.name, names.uniq_name, names.class
    FROM taxon as tax, names
    WHERE
        tax.taxon_id = ?
    AND
        names.taxon_id = tax.taxon_id
SQL
    
    $sth->bind_columns(\$parent_id, \$rank, \$code, \$divid, \$gen_code, \$mito, \$nm, \$uniq, \$class);
    
    $sth->execute($taxonid) or $self->throw($sth->errstr);
    
    my ($sci_name, @common_names);
    
    while ($sth->fetch) {
        if ($class eq 'scientific name') {
            $sci_name = $nm;
        } else {
            push @common_names, $nm;
        }
    }
        
    my $taxon = Bio::Taxon->new(
        -name         => $sci_name,
        -common_names => [@common_names],
        -ncbi_taxid   => $taxonid,
        -parent_id    => $parent_id,   
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
    
    # TODO: note we're not cleaning the query here, so you could technically
    # have a fuzzy match (or Bobby Tables someone)
    
    # TODO: OR'd match seems poor optimally
    my $taxids = $self->{dbh}->selectcol_arrayref(<<SQL);
    SELECT DISTINCT taxon_id FROM names
    WHERE
        name LIKE "$query"
    OR
        uniq_name LIKE "$query"
SQL

    return wantarray() ? @{$taxids} : @{$taxids}[0];
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
    $self->deprecated(); # ?
    #$self->warn(
    #    "get_Children_Taxids is deprecated, use each_Descendent instead");
    #my $id;
    #if ( ref($node) ) {
    #    if ( $node->can('object_id') ) {
    #        $id = $node->object_id;
    #    }
    #    elsif ( $node->can('ncbi_taxid') ) {
    #        $id = $node->ncbi_taxid;
    #    }
    #    else {
    #        $self->warn(
    #            "Don't know how to extract a taxon id from the object of type "
    #              . ref($node)
    #              . "\n" );
    #        return;
    #    }
    #}
    #else { $id = $node }
    #my @vals = $self->{'_parentbtree'}->get_dup($id);
    #return @vals;
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
    
    # TODO:
    # Note here we explicitly set the parent ID, but use a separate method to
    # check whether it is defined. Mixing back-end databases, even if from the
    # same source, should still work (since a different backend wouldn't
    # explicitly set the parent_id)
    
    if ($taxon->trusted_parent_id) {
        # this is the failsafe when we hit the root node
        if ($taxon->parent_id eq $id) {
            return;
        }
        return $self->get_taxon(-taxonid => $taxon->parent_id);
    } else {
        # TODO: would there be any other option?
        return;
    }
}

# TODO: this may act as a drop-in for a recursive CTE lookup

#=head2 ancestors
#
# Title   : ancestors
# Usage   : my @ancestor_taxa = $db->ancestors($taxon)
# Function: Retrieve the full ancestor taxon of a supplied Taxon from the
#           database. 
# Returns : List of Bio::Taxon
# Args    : Bio::Taxon (that was retrieved from this database)
#
#=cut

#sub ancestors { ... }

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
      unless $taxon->db_handle && $taxon->db_handle eq $self;  # yikes
    
    my $id =
      $taxon->id || $self->throw("The supplied Taxon is missing its id!");
    
    #my ( $parent_id, $rank, $code, $divid, $gen_code, $mito, $nm, $uniq, $class );
    # single join or two calls?
    
    # probably not optimal, maybe set up as a cached statement with bindings?
    my $desc_ids = $self->{dbh}->selectcol_arrayref(<<SQL) or $self->throw($self->{dbh}->errstr);
    SELECT tax.taxon_id
    FROM taxon as tax
    WHERE
        tax.parent_id = $id
SQL
    
    return unless ref $desc_ids eq 'ARRAY';
    
    my @descs;
    foreach my $desc_id (@$desc_ids) {
        push( @descs, $self->get_taxon($desc_id) || next );
    }
    return @descs;
}

=head2 Helper methods 

=cut

=head2 index_directory

 Title   : index_directory
 Function : Get/set the location that index files are stored. (this module
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
 Function : Get/set the name of the SQLite3 database where data is stored
 Usage   : $obj->db_name($newval)
 Returns : value of db_name (a scalar)
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

=head2 cache_size

 Title   : cache_size
 Function : Get/set the cachesize used for loading the SQLite3 database
 Usage   : $obj->cache_size($newval)
 Returns : value of cache_size (a scalar)
 Args    : on set, new value (a scalar or undef, optional)
 Note    : we do no checking on whether this value is an integer (SQLite does this for use)

=cut

sub cache_size {
    my $self = shift;
    return $self->{'cache_size'} = shift if defined($_[0]);
    return $self->{'cache_size'};
}

# internal method which does the indexing
sub _build_index {
    my ( $self, $nodesfile, $namesfile, $force ) = @_;

    # TODO: need to disambiguate using index_directory here since we only have
    # one file.  Mayeb ignore it in favor of having full path for db_name?
    my ($dir, $db_name) = ($self->index_directory, $self->db_name);
    if (! -e $db_name || $force) {
        
        # TODO: we're ignoring index_directory for now, may add support for this
        # down the way
        my $dbh = DBI->connect("dbi:SQLite:dbname=$db_name","","") or die $!;
        
        $self->debug("Running SQLite version:".$dbh->{sqlite_version}."\n");
    
        #$dbh->do('PRAGMA synchronous = 0');      # Non transaction safe!!!
        
        if ($self->cache_size) {
            my $cs = $self->cache_size;
            $self->debug("Setting cache size $cs\n");
            $dbh->do("PRAGMA cache_size = $cs") 
        }

        $self->debug("Loading taxon table data\n");
        $self->_init_db($dbh);
        open my $NODES, '<', $nodesfile
            or $self->throw("Could not read node file '$nodesfile': $!");
    
        # TODO: this has the really unnecessary 'OR IGNORE' option added,
        # apparently b.c the test data expects to handle cases where the TaxID
        # is repeated in this table (which should never happen in this table). I
        # will likely change this to throw under those circumstances
        
        my $sth = $dbh->prepare_cached(<<SQL);
    INSERT OR IGNORE INTO taxon (taxon_id, parent_id, rank, code, division_id, gencode_id, mito_id) VALUES (?,?,?,?,?,?,?)
SQL
        $dbh->do("BEGIN");
        while (<$NODES>) {
            next if /^\s*$/;
            chomp;
            my ($taxid,$parent,$rank,$code,$divid,undef,$gen_code,undef,$mito) = split(/\t\|\t/,$_);
            next if $taxid == 1;
            if ($parent == 1) {
                $parent = undef;
            }
            
            $sth->execute($taxid, $parent, $rank, $code, $divid, $gen_code, $mito) or die $sth->errstr.": TaxID $taxid";
        }
        $dbh->do("COMMIT") or $self->throw($dbh->errstr);
        
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
            
            #if ($name =~ /\(class\)$/) { # it seems that only rank of class is ever used in this situation
            #    $name =~ s/\s+\(class\)$//;
            #}
            
            $sth->execute($taxid, $name, $unique_name, $class) or $self->throw($sth->errstr);
        }
        close $NAMES;

        $dbh->do("COMMIT");
        
        $self->debug("Creating taxon index\n");
        $dbh->do("CREATE INDEX parent_idx ON taxon (parent_id)") or $self->throw($dbh->errstr);
        $self->debug("Creating name index\n");
        $dbh->do("CREATE INDEX name_idx ON names (name)") or $self->throw($dbh->errstr);
        $self->debug("Creating taxon name table index\n");
        $dbh->do("CREATE INDEX taxon_name_idx ON names (taxon_id)") or $self->throw($dbh->errstr);

        $dbh->do("PRAGMA foreign_keys = ON");
        
        #$dbh->do('PRAGMA synchronous = 1');
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
    my $dbh = DBI->connect("dbi:SQLite:dbname=$db_name","","") or die $!;
    $dbh->do("PRAGMA foreign_keys = ON");
    if ($self->cache_size) {
        my $cs = $self->cache_size;
        $self->debug("Setting cache size $cs\n");
        $dbh->do("PRAGMA cache_size = $cs") 
    }
    $self->{dbh} = $dbh;

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

sub _dbh_fetch {
    my ($self, $sql) = @_;
    # TODO: more sanity checks
    my $rows = $self->{dbh}->selectrow_arrayref($sql) or $self->throw( $self->{dbh}->errstr );
    return $rows;
}

sub _prepare_cached {
    my ($self, $sql) = @_;
    # TODO: more sanity checks
    my $sth = $self->{dbh}->prepare_cached($sql) or $self->throw( $self->{dbh}->errstr );
    $sth;
}


# TODO: check data size, this is a ballpark estimate (could be reduced)
sub taxon_schema {
    my $self   = shift;
    return {
    taxon   => <<SCHEMA,
    (
        taxon_id        INTEGER UNIQUE PRIMARY KEY NOT NULL,
        parent_id       INTEGER,
        left_id         INTEGER,
        right_id        INTEGER,
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
    undef $self->{dbh};
}

1;

