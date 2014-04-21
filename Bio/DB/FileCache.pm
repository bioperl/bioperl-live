#
# POD documentation - main docs before the code
#
#

=head1 NAME

Bio::DB::FileCache - In file cache for BioSeq objects

=head1 SYNOPSIS



  $cachedb = Bio::DB::FileCache->new($real_db);

  #
  # $real_db is a Bio::DB::RandomAccessI database
  #

  $seq = $cachedb->get_Seq_by_id('ROA1_HUMAN');

  #
  # $seq is a Bio::Seq object
  #

  # more control provided with named-parameter form

  $cachedb = Bio::DB::FileCache->new( -seqdb => $real_db,
				      -file  => $path,
				      -keep  => $flag,
				    );
=head1 DESCRIPTION

This is a disk cache system which saves the objects returned by
Bio::DB::RandomAccessI on disk.  The disk cache grows without limit,
while the process is running, but is automatically unlinked at process
termination unless the -keep flag is set.

This module requires DB_File and Storable.

=head1 CONTACT

Lincoln Stein E<lt>lstein@cshl.orgE<gt>

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::FileCache;

use DB_File;
use Storable qw(freeze thaw);
use Fcntl qw(O_CREAT O_RDWR O_RDONLY);
use File::Temp 'tmpnam';

use strict;


use base qw(Bio::Root::Root Bio::DB::SeqI);

use Bio::Seq::RichSeq;
use Bio::Location::Split;
use Bio::Location::Fuzzy;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::Species;
use Bio::Annotation::Collection;

=head2 new

 Title   : new
 Usage   : $db = Bio::DB::FileCache->new(
                 -seqdb => $db,   # Bio::DB::RandomAccessI database
                 -file  => $path, # path to index file
                 -keep  => $flag, # don't unlink index file
          )
 Function: creates a new on-disk cache
 Returns : a Bio::DB::RandomAccessI database
 Args    : as above
 Throws  : "Must be a randomaccess database" exception
           "Could not open primary index file" exception

If no index file is specified, will create a temporary file in your
system's temporary file directory.  The name of this temporary file
can be retrieved using file_name().

=cut

#'
sub new {
    my ($class,@args) = @_;

    my $self = Bio::Root::Root->new();
    bless $self,$class;

    my ($seqdb,$file_name,$keep) = $self->_rearrange([qw(SEQDB FILE
							 KEEP)],@args);

    if( !defined $seqdb || !ref $seqdb ||
	! $seqdb->isa('Bio::DB::RandomAccessI') ) {
       $self->throw("Must be a randomaccess database not a [$seqdb]");
    }

    $self->seqdb($seqdb);
    $file_name ||= tmpnam();
    $self->file_name($file_name);
    $self->keep($keep);

    $self->_open_database($file_name);
    return $self;
}

=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : $seq = $db->get_Seq_by_id('ROA1_HUMAN')
 Function: Gets a Bio::Seq object by its name
 Returns : a Bio::Seq object
 Args    : the id (as a string) of a sequence
 Throws  : "id does not exist" exception


=cut

sub get_Seq_by_id{
   my ($self,$id) = @_;

   # look in the cache first
   my $obj = $self->_get('id' => $id);
   return $obj if defined $obj;

   # get object from seqdb
   $obj = $self->seqdb->get_Seq_by_id($id);
   $self->_store('id' => $id, $obj);

   return $obj;
}

=head2 get_Seq_by_acc

 Title   : get_Seq_by_acc
 Usage   : $seq = $db->get_Seq_by_acc('X77802');
 Function: Gets a Bio::Seq object by accession number
 Returns : A Bio::Seq object
 Args    : accession number (as a string)
 Throws  : "acc does not exist" exception


=cut

sub get_Seq_by_acc{
   my ($self,$acc) = @_;

   # look in the cache first
   my $obj = $self->_get('acc' => $acc);
   return $obj if defined $obj;

   # get object from seqdb
   $obj = $self->seqdb->get_Seq_by_acc($acc);
   $self->_store('acc' => $acc, $obj);

   return $obj;
}

=head2 seqdb

 Title   : seqdb
 Usage   : $seqdb = $db->seqdb([$seqdb])
 Function: gets/sets the Bio::DB::RandomAccessI database
 Returns : a Bio::DB::RandomAccessI database
 Args    : new sequence database (optional)
 Throws  : nothing

=cut

sub seqdb {
    my ($self, $seqdb) = @_;
    if ($seqdb) {
        $self->{'seqdb'} = $seqdb;
    } else {
        return $self->{'seqdb'};
    }
}

=head2 file_name

 Title   : file_name
 Usage   : $path = $db->file_name([$file_name])
 Function: gets/sets the name of the cache file
 Returns : a path
 Args    : new cache file name (optional)
 Throws  : nothing

It probably isn't useful to set the cache file name after you've
opened it.

=cut

#'

sub file_name {
  my $self = shift;
  my $d = $self->{file_name};
  $self->{file_name} = shift if @_;
  $d;
}

=head2 keep

 Title   : keep
 Usage   : $keep = $db->keep([$flag])
 Function: gets/sets the value of the "keep" flag
 Returns : current value
 Args    : new value (optional)
 Throws  : nothing

The keep flag will cause the index file to be unlinked when the
process exits.  Since on some operating systems (Unix, OS/2) the
unlinking occurs during the new() call immediately after opening the
file, it probably isn't safe to change this value.

=cut

#'
sub keep {
  my $self = shift;
  my $d = $self->{keep};
  $self->{keep} = shift if @_;
  $d;
}

=head2 db

 Title   : db
 Usage   : $db->db
 Function: returns tied hash to index database
 Returns : a Berkeley DB tied hashref
 Args    : none
 Throws  : nothing

=cut

sub db { shift->{db} }

=head2 flush

 Title   : flush
 Usage   : $db->flush
 Function: flushes the cache
 Returns : nothing
 Args    : none
 Throws  : nothing

=cut

sub flush {
  my $db = shift->db or return;
  %{$db} = ();
}

sub _get {
  my $self = shift;
  my ($type,$id) = @_;
  my $serialized = $self->db->{"${type}_${id}"};
  my $obj = thaw($serialized);
  $obj;
}

sub _store {
  my $self = shift;
  my ($type,$id,$obj) = @_;
  if( ! defined $obj ) {
      # bug #1628
      $self->debug("tried to store an undefined value for $id, skipping");
      return;
  }
  my $serialized = freeze($obj);
  $self->db->{"${type}_${id}"} = $serialized;
}

=head2 get_Seq_by_version

 Title   : get_Seq_by_version
 Usage   : $seq = $db->get_Seq_by_version('X77802.1');
 Function: Gets a Bio::Seq object by sequence version
 Returns : A Bio::Seq object
 Args    : accession.version (as a string)
 Throws  : "acc.version does not exist" exception

=cut

sub get_Seq_by_version{
   my ($self,@args) = @_;
   $self->throw("Not implemented it");
}

sub DESTROY {
  my $self = shift;
  unlink $self->file_name unless $self->keep;
}


sub _open_database {
  my $self = shift;
  my $file = shift;
  my $flags = O_CREAT|O_RDWR;
  my %db;
  tie(%db,'DB_File',$file,$flags,0666,$DB_BTREE)
    or $self->throw("Could not open primary index file");
  $self->{db} = \%db;
  unlink $file unless $self->keep;
}

## End of Package

1;
