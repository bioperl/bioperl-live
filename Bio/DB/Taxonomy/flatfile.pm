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

Bio::DB::Taxonomy::flatfile - An implementation of Bio::DB::Taxonomy which uses local flat files

=head1 SYNOPSIS

  use Bio::DB::Taxonomy;

  my $db = new Bio::DB::Taxonomy(-source => 'flatfile'
                                 -nodesfile => $nodesfile,
                                 -namesfile => $namefile);

=head1 DESCRIPTION

This is an implementation which uses local flat files and the
DB_File module RECNO data structures to manage a local
copy of the NCBI Taxonomy database.

File can be obtained from ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::DB::Taxonomy::flatfile;
use vars qw(@ISA $DEFAULT_INDEX_DIR $DEFAULT_NODE_INDEX 
	    $DEFAULT_NAME2ID_INDEX $DEFAULT_ID2NAME_INDEX
	    $NCBI_TAXONOMY_HOSTNAME
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
$NCBI_TAXONOMY_HOSTNAME = 'ftp.ncbi.nih.gov';
$NCBI_TAXONOMY_FILE = '/pub/taxonomy/taxdump.tar.gz';

@DIVISIONS = ([qw(BCT Bacteria)],
	      [qw(INV Invertebrates)],
	      [qw(MAM Mammals)],
	      [qw(PHG Phages)],
	      [qw(PLN Plants)], # (and fungi)
	      [qw(PRI Primates)],
	      [qw(ROD Rodents)],
	      [qw(SYN Synthetic)],
	      [qw(UNA Unassigned)],
	      [qw(VRL Viruses)],
	      [qw(VRT Vertebrates)]
	      );

@ISA = qw(Bio::DB::Taxonomy );

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
  my ($dir,$nodesfile,$namesfile,$force) = $self->_rearrange([qw(DIRECTORY
								 NODESFILE
								 NAMESFILE
								 FORCE)],
							     @args);
  
  $self->index_directory($dir || $DEFAULT_INDEX_DIR);
  if ( $nodesfile ) {
      $self->_build_index($nodesfile,$namesfile,$force);
  }

  $self->_db_connect;
  return $self;
}

=head2 Bio::DB::Taxonomy Interface implementation

=head2 get_Taxonomy_Node

 Title   : get_Taxonomy_Node
 Usage   : my $species = $db->get_Taxonomy_Node(-taxonid => $taxaid)
 Function: Get a Bio::Taxonomy::Taxon object for a taxonid
 Returns : Bio::Taxonomy::Taxon object
 Args    : -taxonid => taxonomy id (to query by taxonid)
            OR
           -name   => string (to query by a taxonomy name: common name, 
                              species, genus, etc)


=cut

sub get_Taxonomy_Node{
   my ($self) = shift;
   my (%item,$taxonid,$name);

   if( @_ > 1 ) {
       ($taxonid,$name) = $self->_rearrange([qw(TAXONID
						NAME)],@_);
       if( $name ) {
	   ($taxonid) = $self->get_taxonid($name);
       }
   } else {  
       $taxonid = shift;
   }
   my $orig_taxonid = $taxonid;
   my (@fields,$node,$taxonnode);
   my $first = 1;
   while( defined ($node = $self->{'_nodes'}->[$taxonid]) ) {
       my ($taxid,$parent,$rank,$code,$divid) = split(SEPARATOR,$node);
       my ($taxon_name) = $self->{'_id2name'}->[$taxid];
       push @fields, $taxon_name if ($rank && $rank ne 'no rank') ;
       if( $first ) {	   
	   $taxonnode = new Bio::Taxonomy::Node(-dbh       => $self,
						-name      => $taxon_name,
						-object_id => $taxid,
						-parent_id => $parent,
						-rank      => $rank,
						-division  => $DIVISIONS[$divid]->[0]);
	   $first = 0;
       }

       last if $parent == 1 || ! $parent || ! $taxid;
       $taxonid = $parent;
   }

   my $speciesnode = new Bio::Species(-ncbi_taxid     => $orig_taxonid,
#				      -common_name    => $item{'CommonName'},
#				      -division       => $item{'Division'});
				      -classification => [@fields],
				      );
   return $taxonnode;
}

=head2 get_taxonid

 Title   : get_taxonid
 Usage   : my $taxonid = $db->get_taxonid('Homo sapiens');
 Function: Searches for a taxonid (typically ncbi_taxon_id) 
           based on a query string 
 Returns : Integer ID
 Args    : String representing species/node name 


=cut

sub get_taxonid {
    my ($self,$query) = @_;
    my $id = $self->{'_name2id'}->{lc($query)};
    if( $id ) { 
	my ($taxid,$name,$fullname) = split(/:/,$id);
	return $taxid;
    }
    return 0;
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
    $self->{'_nodes'}    = [];
    $self->{'_id2name'} = [];
    $self->{'_name2id'} = {};
        
    if( ! -e $nodeindex || $force ) {
	open(NODES,$nodesfile) || 
	    $self->throw("Cannot open node file '$nodesfile' for reading");
	
	unlink $nodeindex;
	tie ( @{$self->{'_nodes'}}, 'DB_File', $nodeindex, O_RDWR|O_CREAT, 
	      0644, $DB_RECNO) || 
		  $self->throw("Cannot open file '$nodeindex': $!");	
	while(<NODES>) {
	    chomp;
	    my ($taxid,$parent,$rank,$code,$divid) = split(/\t\|\t/,$_);
	    # keep this stringified
	    $self->{'_nodes'}->[$taxid] = join(SEPARATOR, 
					       ($taxid,$parent,$rank,
						$code,$divid));
	}
	close(NODES);
	undef $self->{'_nodes'};
	untie( @{$self->{'_nodes'}} );
    }
    if( ! -e $name2idindex || ! -e $id2nameindex || $force ) { 
	open(NAMES,$namesfile) || 
	    $self->throw("Cannot open names file '$namesfile' for reading");

	unlink $name2idindex;
	unlink $id2nameindex;

	tie (@{$self->{'_id2name'}}, 'DB_File', $id2nameindex, 
	     O_RDWR|O_CREAT, 0644, $DB_RECNO) || 
		 $self->throw("Cannot open file '$id2nameindex': $!");
	
	tie ( %{$self->{'_name2id'}}, 'DB_File', $name2idindex, O_RDWR|O_CREAT, 
	      0644, $DB_HASH) || 
		  $self->throw("Cannot open file '$name2idindex': $!");
	
	while(<NAMES>) {
	    chomp;	    
	    my ($taxid,$name,$uniquename,$class) = split(/\t\|\t/,$_);
	    $class =~ s/\s+\|\s*$//;
	    $uniquename = $name unless $uniquename;
	    my $idx = lc($name);
	    $self->{'_name2id'}->{$idx} = join(SEPARATOR,
					       ($taxid, $name,$uniquename,
						$class));
	    if( $class && $class eq 'scientific name' ) {
		# only store the id2name lookup when it is the "proper" name
		$self->{'_id2name'}->[$taxid] = $uniquename;
	    }
	}
	close(NAMES);
	undef $self->{'_id2name'};
	undef $self->{'_name2id'};
	untie( %{$self->{'_name2id'}} );
	untie( @{$self->{'_id2name'}} );
    }
}

# connect the internal db handle and 

sub _db_connect {
    my $self = shift;
    return if $self->{'_initialized'};

#    undef $self->{'_nodes'};
#    undef $self->{'_id2name'};
#    undef $self->{'_name2id'};
#    untie( %{$self->{'_name2id'}} );
#    untie( @{$self->{'_id2name'}} );
#    untie( @{$self->{'_nodes'}} );
    
    $self->{'_nodes'}   = [];
    $self->{'_id2name'} = [];
    $self->{'_name2id'} = {};

    my ($dir) = ($self->index_directory);
    my $nodeindex = "$dir/$DEFAULT_NODE_INDEX";
    my $name2idindex = "$dir/$DEFAULT_NAME2ID_INDEX";
    my $id2nameindex = "$dir/$DEFAULT_ID2NAME_INDEX";
    if( ! -e $nodeindex ||
	! -e $name2idindex || 
	! -e $id2nameindex ) {
	$self->warn("Index files have not been created");
	return 0;
    }
    tie ( @{$self->{'_nodes'}}, 'DB_File', 
	  $nodeindex, O_RDONLY,0644, $DB_RECNO) 
	|| $self->throw("$! $nodeindex");
    tie (@{$self->{'_id2name'}}, 'DB_File', $id2nameindex,O_RDONLY, 0644, 
	 $DB_RECNO) || $self->throw("$! $id2nameindex");
    
    tie ( %{$self->{'_name2id'}}, 'DB_File', $name2idindex, O_RDONLY,0644, 
	  $DB_HASH) || $self->throw("$! $name2idindex");
    $self->{'_initialized'}  = 1;
}


=head2 index_directory

 Title   : index_directory
 Usage   : $obj->index_directory($newval)
 Function: 
 Example : 
 Returns : value of index_directory (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub index_directory {
    my $self = shift;

    return $self->{'index_directory'} = shift if @_;
    return $self->{'index_directory'};
}
1;
