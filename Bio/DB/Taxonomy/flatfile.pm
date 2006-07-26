# $Id$
#
# BioPerl module for Bio::DB::Taxonomy::flatfile
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Taxonomy::flatfile - An implementation of Bio::DB::Taxonomy
which uses local flat files

=head1 SYNOPSIS

  use Bio::DB::Taxonomy;

  my $db = new Bio::DB::Taxonomy(-source => 'flatfile'
                                 -nodesfile => $nodesfile,
                                 -namesfile => $namefile);

=head1 DESCRIPTION

This is an implementation which uses local flat files and the DB_File
module RECNO data structures to manage a local copy of the NCBI
Taxonomy database.

Required database files can be obtained from
ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

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

Describe contact details here

=head1 CONTRIBUTORS

Sendu Bala: bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::Taxonomy::flatfile;
use vars qw(@ISA $DEFAULT_INDEX_DIR $DEFAULT_NODE_INDEX 
	    $DEFAULT_NAME2ID_INDEX $DEFAULT_ID2NAME_INDEX
	    $NCBI_TAXONOMY_HOSTNAME $DEFAULT_PARENT_INDEX
	    $NCBI_TAXONOMY_FILE @DIVISIONS);
use strict;
use Bio::DB::Taxonomy;
use Bio::Taxonomy::Node;
use Bio::Species;
use DB_File;

use constant SEPARATOR => ':';

$DEFAULT_INDEX_DIR = '/tmp';
$DEFAULT_NODE_INDEX = 'nodes';
$DEFAULT_NAME2ID_INDEX = 'names2id';
$DEFAULT_ID2NAME_INDEX = 'id2names';
$DEFAULT_PARENT_INDEX = 'parents';
$NCBI_TAXONOMY_HOSTNAME = 'ftp.ncbi.nih.gov';
$NCBI_TAXONOMY_FILE = '/pub/taxonomy/taxdump.tar.gz';

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

@ISA = qw( Bio::DB::Taxonomy );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::DB::Taxonomy::flatfile();
 Function: Builds a new Bio::DB::Taxonomy::flatfile object 
 Returns : an instance of Bio::DB::Taxonomy::flatfile
 Args    : -directory => name of directory where index files should be created
           -nodesfile => name of file containing nodes (nodes.dmp from NCBI)
           -namesfile => name of the file containing names(names.dmp from NCBI)
           -force     => 1 replace current indexes even if they exist

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($dir,$nodesfile,$namesfile,$force) = $self->_rearrange([qw
	  (DIRECTORY NODESFILE NAMESFILE FORCE)], @args);
  
  $self->index_directory($dir || $DEFAULT_INDEX_DIR);
  if ( $nodesfile ) {
	  $self->_build_index($nodesfile,$namesfile,$force);
  }

  $self->_db_connect;
  return $self;
}

=head2 Bio::DB::Taxonomy Interface implementation

=cut

=head2 get_Taxonomy_Node

 Title   : get_Taxonomy_Node
 Usage   : my $species = $db->get_Taxonomy_Node(-taxonid => $taxaid)
 Function: Get a Bio::Taxonomy::Node object for a taxonid
 Returns : Bio::Taxonomy::Node object
 Args    : -taxonid => taxonomy id (to query by taxonid)
            OR
           -name   => string (to query by a taxonomy name: common name, 
                              species, genus, etc)

See L<Bio::Taxonomy::Taxon>

=cut

sub get_Taxonomy_Node {
   my ($self) = shift;
   my (%item,$taxonid,$name);

   if( @_ > 1 ) {
       ($taxonid,$name) = $self->_rearrange([qw(TAXONID NAME)],@_);
       if ($name) {
            ($taxonid, my @others) = $self->get_taxonid($name);
            $self->warn("There were multiple ids ($taxonid @others) matching '$name', using '$taxonid'") if @others > 0;
       }
   } else {  
       $taxonid = shift;
   }
   
   my (@fields,$node,$taxonnode);
   my $first = 1;
   my @classification;
   while( defined ($node = $self->{'_nodes'}->[$taxonid]) && length($node) ) {
       my ($taxid,$parent,$rank,$code,$divid,$gen_code,$mito) = split(SEPARATOR,$node);
       my ($taxon_names) = $self->{'_id2name'}->[$taxid];
       my ($sci_name, @common_names) = split(SEPARATOR, $taxon_names);
       
       if ($first) {
	   $taxonnode = new Bio::Taxonomy::Node(-dbh       => $self,
						-name      => $sci_name,
                        -common_names => [@common_names],
						-object_id => $taxid,
						-parent_id => $parent,
						-rank      => $rank,
						-division  => $DIVISIONS[$divid]->[1],
                        -genetic_code => $gen_code,
                        -mito_genetic_code => $mito );
	   $first = 0;
       }
       
       push @fields, $sci_name if ($rank && $rank ne 'no rank');
       last if ! defined $parent || $parent == 1 || ! $taxid;
       $taxonid = $parent;
   }
   $taxonnode->classification(@fields) if defined $taxonnode;
   
   return $taxonnode;
}

=head2 get_taxonids

 Title   : get_taxonids
 Usage   : my @taxonids = $db->get_taxonids('Homo sapiens');
 Function: Searches for a taxonid (typically ncbi_taxon_id) based on a query
           string. Note that multiple taxonids can match to the same supplied
           name.
 Returns : array of integer ids in list context, one of these in scalar context
 Args    : string representing taxanomic (node) name

=cut

sub get_taxonids {
    my ($self, $query) = @_;
    my $ids = $self->{'_name2id'}->{lc($query)} || return;
    my @ids = split(SEPARATOR, $ids);
    return wantarray() ? @ids : shift @ids;
}

*get_taxonid = \&get_taxonids;

=head2 get_Children_Taxids

 Title   : get_Children_Taxids
 Usage   : my @childrenids = $db->get_Children_Taxids 
 Function: Get the children of a node in the taxonomy
 Returns : Array of Ids
 Args    : Bio::Taxonomy::Node or a taxon_id

See L<Bio::Taxonomy::Node>

=cut

sub get_Children_Taxids {
   my ($self,$node) = @_;
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

=head2 Helper methods 

=cut

# internal method which does the indexing
sub _build_index {
    my ($self,$nodesfile,$namesfile,$force) = @_;
    
    my ($dir) = ($self->index_directory);
    my $nodeindex = "$dir/$DEFAULT_NODE_INDEX";
    my $name2idindex = "$dir/$DEFAULT_NAME2ID_INDEX";
    my $id2nameindex = "$dir/$DEFAULT_ID2NAME_INDEX";
    my $parent2childindex = "$dir/$DEFAULT_PARENT_INDEX";
    $self->{'_nodes'}    = [];
    $self->{'_id2name'} = [];
    $self->{'_name2id'} = {};
    $self->{'_parent2children'} = {};
    
    if (! -e $nodeindex || $force) {
        my (%parent2children,@nodes);
        open(NODES,$nodesfile) || 
            $self->throw("Cannot open node file '$nodesfile' for reading");
        
        unlink $nodeindex;
        unlink $parent2childindex;
        my $nh = tie ( @nodes, 'DB_File', $nodeindex, O_RDWR|O_CREAT, 0644, $DB_RECNO) || 
            $self->throw("Cannot open file '$nodeindex': $!");	
        my $btree = tie( %parent2children, 'DB_File', $parent2childindex, O_RDWR|O_CREAT, 0644, $DB_BTREE) || 
            $self->throw("Cannot open file '$parent2childindex': $!");	
        
        while (<NODES>) {
            chomp;
            my ($taxid,$parent,$rank,$code,$divid,undef,$gen_code,undef,$mito) = split(/\t\|\t/,$_);
            # keep this stringified
            $nodes[$taxid] = join(SEPARATOR, ($taxid,$parent,$rank,$code,$divid,$gen_code,$mito));
            $btree->put($parent,$taxid);
        }
        close(NODES);
        
        $nh = $btree = undef;
        untie @nodes ;
        untie %parent2children;
    }
    
    if ((! -e $name2idindex || -z $name2idindex) || (! -e $id2nameindex || -z $id2nameindex) || $force) { 
        open(NAMES,$namesfile) || 
            $self->throw("Cannot open names file '$namesfile' for reading");
        
        unlink $name2idindex;
        unlink $id2nameindex;
        my (@id2name,%name2id);
        my $idh = tie (@id2name, 'DB_File', $id2nameindex, O_RDWR|O_CREAT, 0644, $DB_RECNO) || 
            $self->throw("Cannot open file '$id2nameindex': $!");
        my $nameh = tie ( %name2id, 'DB_File', $name2idindex, O_RDWR|O_CREAT, 0644, $DB_HASH) || 
            $self->throw("Cannot open file '$name2idindex': $!");
        
        while (<NAMES>) {
            chomp;	    
            my ($taxid, $name, $unique_name, $class) = split(/\t\|\t/,$_);
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
        close(NAMES);
        
        $idh = $nameh = undef;
        untie( %name2id);
        untie( @id2name);
    }
}

# connect the internal db handle
sub _db_connect {
    my $self = shift;
    return if $self->{'_initialized'};
    
    $self->{'_nodes'}   = [];
    $self->{'_id2name'} = [];
    $self->{'_name2id'} = {};
    
    my ($dir) = ($self->index_directory);
    my $nodeindex = "$dir/$DEFAULT_NODE_INDEX";
    my $name2idindex = "$dir/$DEFAULT_NAME2ID_INDEX";
    my $id2nameindex = "$dir/$DEFAULT_ID2NAME_INDEX";
    my $parent2childindex = "$dir/$DEFAULT_PARENT_INDEX";
    
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
    $self->{'_initialized'}  = 1;
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

1;