#
# BioPerl module for Bio::DB::Taxonomy::flatfile
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Taxonomy::flatfile - Use the NCBI taxonomy from local indexed flat files

=head1 SYNOPSIS

  use Bio::DB::Taxonomy;

  my $db = Bio::DB::Taxonomy->new(-source    => 'flatfile' ,
                                  -nodesfile => 'nodes.dmp',
                                  -namesfile => 'names.dmp');

=head1 DESCRIPTION

This is an implementation of Bio::DB::Taxonomy which stores and accesses the
NCBI taxonomy using flat files stored locally on disk and indexed using the
DB_File module RECNO data structure for fast retrieval.

The required database files, nodes.dmp and names.dmp can be obtained from
ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

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

Email jason-at-bioperl-dot-org

=head1 CONTRIBUTORS

Sendu Bala: bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::Taxonomy::flatfile;

use vars qw($DEFAULT_INDEX_DIR $DEFAULT_NODE_INDEX $DEFAULT_NAME2ID_INDEX
            $DEFAULT_ID2NAME_INDEX $DEFAULT_PARENT_INDEX @DIVISIONS);

use strict;
use DB_File;
use Bio::Taxon;
use File::Spec::Functions;

use constant SEPARATOR => ':';

$DEFAULT_INDEX_DIR     = $Bio::Root::IO::TEMPDIR; # /tmp
$DEFAULT_NODE_INDEX    = 'nodes';
$DEFAULT_NAME2ID_INDEX = 'names2id';
$DEFAULT_ID2NAME_INDEX = 'id2names';
$DEFAULT_PARENT_INDEX  = 'parents';

$DB_BTREE->{'flags'} = R_DUP; # allow duplicate values in DB_File BTREEs

@DIVISIONS =   ([qw(BCT Bacteria)],
                [qw(INV Invertebrates)],
                [qw(MAM Mammals)],
                [qw(PHG Phages)],
                [qw(PLN Plants)], # (and fungi)
                [qw(PRI Primates)],
                [qw(ROD Rodents)],
                [qw(SYN Synthetic)],
                [qw(UNA Unassigned)],
                [qw(VRL Viruses)],
                [qw(VRT Vertebrates)],
                [qw(ENV 'Environmental samples')]);

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
  my($class, @args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($dir,$nodesfile,$namesfile,$force) =
      $self->_rearrange([qw(DIRECTORY NODESFILE NAMESFILE FORCE)], @args);
  
  $self->index_directory($dir || $DEFAULT_INDEX_DIR);
  if ( $nodesfile ) {
          $self->_build_index($nodesfile,$namesfile,$force);
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
    if (not exists $self->{_num_taxa}) {
        my $num = 0;
        while ( my ($parent, undef) = each %{$self->{_parent2children}} ) {
           $num++;
        }
        $self->{_num_taxa} = $num;
    }
    return $self->{_num_taxa};
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
    
    return unless $taxonid;
    
    $taxonid =~ /^\d+$/ || return;
    my $node = $self->{'_nodes'}->[$taxonid] || return;
    length($node) || return;
    my ($taxid, undef, $rank, $code, $divid, $gen_code, $mito) = split(SEPARATOR,$node);
    last unless defined $taxid;
    my ($taxon_names) = $self->{'_id2name'}->[$taxid];
    my ($sci_name, @common_names) = split(SEPARATOR, $taxon_names);
    
    my $taxon = Bio::Taxon->new(
                        -name         => $sci_name,
                        -common_names => [@common_names],
                        -ncbi_taxid   => $taxid, # since this is a real ncbi taxid, explicitly set it as one
                        -rank         => $rank,
                        -division     => $DIVISIONS[$divid]->[1],
                        -genetic_code => $gen_code,
                        -mito_genetic_code => $mito );
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
    my ($self, $query) = @_;
    my $ids = $self->{'_name2id'}->{lc($query)};
    unless ($ids) {
        if ($query =~ /_/) {
            # try again converting underscores to spaces
            $query =~ s/_/ /g;
            $ids = $self->{'_name2id'}->{lc($query)};
        }
        $ids || return;
    }
    my @ids = split(SEPARATOR, $ids);
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
   my ($self, $node) = @_;
   $self->warn("get_Children_Taxids is deprecated, use each_Descendent instead");
   my $id;
   if( ref($node) ) {
       if( $node->can('object_id') ) {
           $id = $node->object_id;
       } elsif( $node->can('ncbi_taxid') ) {
           $id = $node->ncbi_taxid;
       } else { 
           $self->warn("Don't know how to extract a taxon id from the object of type ".ref($node)."\n");
           return;
       }
   } else { $id = $node }
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
    my ($self, $taxon) = @_;
    $self->throw("Must supply a Bio::Taxon") unless ref($taxon) && $taxon->isa('Bio::Taxon');
    $self->throw("The supplied Taxon must belong to this database") unless $taxon->db_handle && $taxon->db_handle eq $self;
    my $id = $taxon->id || $self->throw("The supplied Taxon is missing its id!");
    
    my $node = $self->{'_nodes'}->[$id];
    if (length($node)) {
        my (undef, $parent_id) = split(SEPARATOR,$node);
        $parent_id || return;
        $parent_id eq $id && return; # one of the roots
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
    my ($self, $taxon) = @_;
    $self->throw("Must supply a Bio::Taxon") unless ref($taxon) && $taxon->isa('Bio::Taxon');
    $self->throw("The supplied Taxon must belong to this database") unless $taxon->db_handle && $taxon->db_handle eq $self;
    my $id = $taxon->id || $self->throw("The supplied Taxon is missing its id!");

    my @desc_ids = $self->{'_parentbtree'}->get_dup($id);
    my @descs;
    foreach my $desc_id (@desc_ids) {
        push(@descs, $self->get_taxon($desc_id) || next);
    }
    return @descs;
}


=head2 Helper methods 

=cut

# internal method which does the indexing
sub _build_index {
    my ($self, $nodesfile, $namesfile, $force) = @_;
    
    my $dir = $self->index_directory;
    my $nodeindex         = catfile($dir, $DEFAULT_NODE_INDEX);
    my $name2idindex      = catfile($dir, $DEFAULT_NAME2ID_INDEX);
    my $id2nameindex      = catfile($dir, $DEFAULT_ID2NAME_INDEX);
    my $parent2childindex = catfile($dir, $DEFAULT_PARENT_INDEX);
    $self->{'_nodes'}           = [];
    $self->{'_id2name'}         = [];
    $self->{'_name2id'}         = {};
    $self->{'_parent2children'} = {};
    
    if (! -e $nodeindex || $force) {
        my (%parent2children,@nodes);
        open my $NODES, '<', $nodesfile
            or $self->throw("Could not read node file '$nodesfile': $!");
        
        unlink $nodeindex;
        unlink $parent2childindex;
        my $nh = tie ( @nodes, 'DB_File', $nodeindex, O_RDWR|O_CREAT, 0644, $DB_RECNO) || 
            $self->throw("Cannot open file '$nodeindex': $!");
        my $btree = tie( %parent2children, 'DB_File', $parent2childindex, O_RDWR|O_CREAT, 0644, $DB_BTREE) || 
            $self->throw("Cannot tie to file '$parent2childindex': $!");
        
        while (<$NODES>) {
            next if /^$/;
            chomp;
            my ($taxid,$parent,$rank,$code,$divid,undef,$gen_code,undef,$mito) = split(/\t\|\t/,$_);
            # don't include the fake root node 'root' with id 1; we essentially have multiple roots here
            next if $taxid == 1;
            if ($parent == 1) {
                $parent = $taxid;
            }

            # keep this stringified
            $nodes[$taxid] = join(SEPARATOR, ($taxid,$parent,$rank,$code,$divid,$gen_code,$mito));
            $btree->put($parent,$taxid);
        }
        close $NODES;
        
        $nh = $btree = undef;
        untie @nodes ;
        untie %parent2children;
    }
    
    if ((! -e $name2idindex || -z $name2idindex) || (! -e $id2nameindex || -z $id2nameindex) || $force) { 
        open my $NAMES, '<', $namesfile
            or $self->throw("Could not read names file '$namesfile': $!");
        
        unlink $name2idindex;
        unlink $id2nameindex;
        my (@id2name,%name2id);
        my $idh = tie (@id2name, 'DB_File', $id2nameindex, O_RDWR|O_CREAT, 0644, $DB_RECNO) || 
            $self->throw("Cannot tie to file '$id2nameindex': $!");
        my $nameh = tie ( %name2id, 'DB_File', $name2idindex, O_RDWR|O_CREAT, 0644, $DB_HASH) || 
            $self->throw("Cannot tie to file '$name2idindex': $!");
        
        while (<$NAMES>) {
            next if /^$/;
            chomp; 
            my ($taxid, $name, $unique_name, $class) = split(/\t\|\t/,$_);
            # don't include the fake root node 'root' or 'all' with id 1
            next if $taxid == 1;

            $class =~ s/\s+\|\s*$//;
            my $lc_name = lc($name);
            my $orig_name = $name;
            
            # unique names aren't always in the correct column, sometimes they
            # are uniqued by adding bracketed rank names to the normal name;
            # store the uniqued version then fix the name for normal use
            if ($lc_name =~ /\(class\)$/) { # it seems that only rank of class is ever used in this situation
                $name2id{$lc_name} = $taxid;
                $name =~ s/\s+\(class\)$//;
                $lc_name = lc($name);
            }
            
            # handle normal names which aren't necessarily unique
            my $taxids = $name2id{$lc_name} || '';
            my %taxids = map { $_ => 1 } split(SEPARATOR, $taxids);
            unless (exists $taxids{$taxid}) {
                $taxids{$taxid} = 1;
                $name2id{$lc_name} = join(SEPARATOR, keys %taxids);
            }
            
            # store unique names in name2id
            if ($unique_name) {
                $name2id{lc($unique_name)} = $taxid;
            }
            
            # store all names in id2name array
            my $names = $id2name[$taxid] || '';
            my @names = split(SEPARATOR, $names);
            if ($class && $class eq 'scientific name') {
                # the scientific name should be the first name stored
                unshift(@names, $name);
                push(@names, $orig_name) if ($orig_name ne $name);
                push(@names, $unique_name) if $unique_name;
            }
            else {
                # all other ('common' in this simplification) names get added after
                push(@names, $name);
                push(@names, $orig_name) if ($orig_name ne $name);
                push(@names, $unique_name) if $unique_name;
            }
            $id2name[$taxid] = join(SEPARATOR, @names);
        }
        close $NAMES;
        
        $idh = $nameh = undef;
        untie( %name2id);
        untie( @id2name);
    }
}


# connect the internal db handle
sub _db_connect {
    my $self = shift;
    return if $self->{'_initialized'};

    my $dir = $self->index_directory;
    my $nodeindex         = catfile($dir, $DEFAULT_NODE_INDEX);
    my $name2idindex      = catfile($dir, $DEFAULT_NAME2ID_INDEX);
    my $id2nameindex      = catfile($dir, $DEFAULT_ID2NAME_INDEX);
    my $parent2childindex = catfile($dir, $DEFAULT_PARENT_INDEX);
    $self->{'_nodes'}           = [];
    $self->{'_id2name'}         = [];
    $self->{'_name2id'}         = {};
    $self->{'_parent2children'} = {};
    
    if( ! -e $nodeindex ||
        ! -e $name2idindex || 
        ! -e $id2nameindex ) {
        $self->warn("Index files have not been created");
        return 0;
    }
    tie ( @{$self->{'_nodes'}}, 'DB_File', $nodeindex, O_RDWR,undef, $DB_RECNO) 
        || $self->throw("$! $nodeindex");
    tie (@{$self->{'_id2name'}}, 'DB_File', $id2nameindex,O_RDWR, undef, 
        $DB_RECNO) || $self->throw("$! $id2nameindex");
    
    tie ( %{$self->{'_name2id'}}, 'DB_File', $name2idindex, O_RDWR,undef, 
        $DB_HASH) || $self->throw("$! $name2idindex");
    $self->{'_parentbtree'} = tie( %{$self->{'_parent2children'}},
                                   'DB_File', $parent2childindex, 
                                   O_RDWR, 0644, $DB_BTREE);

    $self->{'_initialized'} = 1;
}


=head2 index_directory

 Title   : index_directory
 Funtion : Get/set the location that index files are stored. (this module
           will index the supplied database)
 Usage   : $obj->index_directory($newval)
 Returns : value of index_directory (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut


sub index_directory {
    my $self = shift;
    return $self->{'index_directory'} = shift if @_;
    return $self->{'index_directory'};
}


sub DESTROY {
    my $self = shift;
    # Destroy all filehandle references
    # to be able to remove temporary files
    undef $self->{_id2name};
    undef $self->{_name2id};
    undef $self->{_nodes};
    undef $self->{_parent2children};
    undef $self->{_parentbtree};
    unlink catfile($self->{index_directory},'id2names');
    unlink catfile($self->{index_directory},'names2id');
    unlink catfile($self->{index_directory},'nodes');
    unlink catfile($self->{index_directory},'parents');
}

1;
