#
#
# BioPerl module for Bio::DB::Flat::BDB
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Lincoln Stein <lstein@cshl.org>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Flat::BDB - Interface for BioHackathon standard BDB-indexed flat file

=head1 SYNOPSIS

  #You should not be using this module directly.

See L<Bio::DB::Flat>.

=head1 DESCRIPTION

This object provides the basic mechanism to associate positions in
files with primary and secondary name spaces. Unlike
Bio::Index::Abstract (see L<Bio::Index::Abstract>), this is specialized
to work with the BerkeleyDB-indexed "common" flat file format worked
out at the 2002 BioHackathon.

This object is the guts to the mechanism, which will be used by the
specific objects inheriting from it.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via
email or the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Lincoln Stein

Email - lstein@cshl.org

=head1 SEE ALSO

L<Bio::DB::Flat>,

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with an "_" (underscore).

=cut


# Let the code begin...

package Bio::DB::Flat::BDB;

use strict;
use DB_File;
use IO::File;
use Fcntl qw(O_CREAT O_RDWR O_RDONLY);
use File::Spec;
use Bio::SeqIO;
use Bio::DB::RandomAccessI;
use Bio::Root::Root;
use Bio::Root::IO;

use base qw(Bio::DB::Flat);

sub _initialize {
  my $self = shift;
  my ($max_open) = $self->_rearrange(['MAXOPEN'],@_);
  $self->{bdb_maxopen} = $max_open || 32;
}

# return a filehandle seeked to the appropriate place
# this only works with the primary namespace
sub _get_stream {
  my ($self,$id) = @_;
  my ($filepath,$offset,$length) = $self->_lookup_primary($id)
    or $self->throw("Unable to find a record for $id in the flat file index");
  my $fh = $self->_fhcache($filepath)
    or $self->throw("couldn't open $filepath: $!");
  seek($fh,$offset,0) or $self->throw("can't seek on $filepath: $!");
  $fh;
}

# return records corresponding to the indicated index
# if there are multiple hits will return a list in list context,
# otherwise will throw an exception
sub fetch_raw {
  my ($self,$id,$namespace) = @_;

  # secondary lookup
  if (defined $namespace && $namespace ne $self->primary_namespace) {
    my @hits = $self->_lookup_secondary($namespace,$id);
    $self->throw("Multiple records correspond to $namespace=>$id but function called in a scalar context")
      unless wantarray;
    return map {$self->_read_record(@$_)} @hits;
  }

  # primary lookup
  my @args = $self->_lookup_primary($id)
    or $self->throw("Unable to find a record for $id in the flat file index");
  return $self->_read_record(@args);
}

# create real live Bio::Seq object
sub get_Seq_by_id {
  my $self = shift;
  my $id   = shift;
  my $fh   = eval {$self->_get_stream($id)} or return;
  my $seqio =
    $self->{bdb_cached_parsers}{fileno $fh} ||= Bio::SeqIO->new( -Format => $self->file_format,
								 -fh     => $fh);
  return $seqio->next_seq;
}

# fetch array of Bio::Seq objects
sub get_Seq_by_acc {
  my $self = shift;
  unshift @_,'ACC' if @_==1;
  my ($ns,$key) = @_;
  my @primary_ids = $self->expand_ids($ns => $key);
  $self->throw("more than one sequences correspond to this accession")
      if @primary_ids > 1 && ! wantarray;
  my @rc = map {$self->get_Seq_by_id($_)} @primary_ids;
  return wantarray ? @rc : $rc[0];
}

# fetch array of Bio::Seq objects
sub get_Seq_by_version {
  my $self = shift;
  unshift @_,'VERSION' if @_==1;
  my ($ns,$key) = @_;
  my @primary_ids = $self->expand_ids($ns => $key);
  $self->throw("more than one sequences correspond to this accession")
    if @primary_ids > 1 && !wantarray;
  my @rc = map {$self->get_Seq_by_id($_)} @primary_ids;
  return wantarray ? @rc : $rc[0];
}

=head2 get_PrimarySeq_stream

 Title   : get_PrimarySeq_stream
 Usage   : $stream = get_PrimarySeq_stream
 Function: Makes a Bio::DB::SeqStreamI compliant object
           which provides a single method, next_primary_seq
 Returns : Bio::DB::SeqStreamI
 Args    : none


=cut

sub get_PrimarySeq_stream {
  my $self = shift;
  my @files  = $self->files || 0;
  my $out = Bio::SeqIO::MultiFile->new( -format => $self->file_format ,
					-files  => \@files);
  return $out;
}

sub get_all_primary_ids {
  my $self = shift;
  my $db   = $self->primary_db;
  return keys %$db;
}

=head2 get_all_primary_ids

 Title   : get_all_primary_ids
 Usage   : @ids = $seqdb->get_all_primary_ids()
 Function: gives an array of all the primary_ids of the
           sequence objects in the database.
 Example :
 Returns : an array of strings
 Args    : none

=cut

# this will perform an ID lookup on a (possibly secondary)
# id, returning all the corresponding ids
sub expand_ids {
  my $self = shift;
  my ($ns,$key) = @_;
  return $key unless defined $ns;
  return $key if $ns eq $self->primary_namespace;
  my $db   = $self->secondary_db($ns)
    or $self->throw("invalid secondary namespace $ns");
  my $record = $db->{$key} or return;  # nothing there
  return $self->unpack_secondary($record);
}

# build index from files listed
sub build_index {
  my $self  = shift;
  my @files = @_;
  my $count = 0;
  for my $file (@files) {
    $file = File::Spec->rel2abs($file)
      unless File::Spec->file_name_is_absolute($file);
    $count += $self->_index_file($file);
  }
  $self->write_config;
  $count;
}

sub _index_file {
  my $self = shift;
  my $file = shift;

  my $fileno = $self->_path2fileno($file);
  defined $fileno or $self->throw("could not create a file number for $file");

  my $fh     = $self->_fhcache($file) or $self->throw("could not open $file for indexing: $!");
  my $offset = 0;
  my $count  = 0;

  while (!eof($fh)) {
    my ($ids,$adjustment)  = $self->parse_one_record($fh) or next;
    $adjustment ||= 0;  # prevent uninit variable warning
    my $pos = tell($fh) + $adjustment;
    $self->_store_index($ids,$file,$offset,$pos-$offset);
    $offset = $pos;
    $count++;
  }
  $count;
}

=head2 To Be Implemented in Subclasses

The following methods MUST be implemented by subclasses.

=cut

=head2 May Be Overridden in Subclasses

The following methods MAY be overridden by subclasses.

=cut

sub default_primary_namespace {
  return "ACC";
}

sub default_secondary_namespaces {
  return;
}

sub _read_record {
  my $self = shift;
  my ($filepath,$offset,$length) = @_;
  my $fh = $self->_fhcache($filepath)
    or $self->throw("couldn't open $filepath: $!");
  seek($fh,$offset,0) or $self->throw("can't seek on $filepath: $!");
  my $record;
  read($fh,$record,$length) or $self->throw("can't read $filepath: $!");
  $record
}

# return a list in the form ($filepath,$offset,$length)
sub _lookup_primary {
  my $self    = shift;
  my $primary = shift;
  my $db     = $self->primary_db
    or $self->throw("no primary namespace database is open");

  my $record = $db->{$primary} or return;  # nothing here

  my($fileid,$offset,$length) = $self->unpack_primary($record);
  my $filepath = $self->_fileno2path($fileid)
    or $self->throw("no file path entry for fileid $fileid");
  return ($filepath,$offset,$length);
}

# return a list of array refs in the form [$filepath,$offset,$length]
sub _lookup_secondary {
  my $self = shift;
  my ($namespace,$secondary) = @_;
  my @primary = $self->expand_ids($namespace=>$secondary);
  return map {[$self->_lookup_primary($_)]} @primary;
}

# store indexing information into a primary & secondary record
# $namespaces is one of:
#     1. a scalar corresponding to the primary name
#     2. a hashref corresponding to namespace=>id identifiers
#              it is valid for secondary id to be an arrayref
sub _store_index {
  my $self = shift;
  my ($keys,$filepath,$offset,$length) = @_;
  my ($primary,%secondary);

  if (ref $keys eq 'HASH') {
    my %valid_secondary = map {$_=>1} $self->secondary_namespaces;
    while (my($ns,$value) = each %$keys) {
      if ($ns eq $self->primary_namespace) {
	$primary = $value;
      } else {
	$valid_secondary{$ns} or $self->throw("invalid secondary namespace $ns");
	push @{$secondary{$ns}},$value;
      }
    }
    $primary or $self->throw("no primary namespace ID provided");
  } else {
    $primary = $keys;
  }

  $self->throw("invalid primary ID; must be a scalar") 
    if ref($primary) =~ /^(ARRAY|HASH)$/;  # but allow stringified objects

  $self->_store_primary($primary,$filepath,$offset,$length);
  for my $ns (keys %secondary) {
    my @ids = ref $secondary{$ns} ? @{$secondary{$ns}} : $secondary{$ns};
    $self->_store_secondary($ns,$_,$primary) foreach @ids;
  }

  1;
}

# store primary index
sub _store_primary {
  my $self = shift;
  my ($id,$filepath,$offset,$length) = @_;

  my $db = $self->primary_db
    or $self->throw("no primary namespace database is open");
  my $fileno = $self->_path2fileno($filepath);
  defined $fileno or $self->throw("could not create a file number for $filepath");

  my $record = $self->pack_primary($fileno,$offset,$length);
  $db->{$id} = $record or return;  # nothing here
  1;
}

# store a primary index name under a secondary index
sub _store_secondary {
  my $self = shift;
  my ($secondary_ns,$secondary_id,$primary_id) = @_;

  my $db   = $self->secondary_db($secondary_ns)
    or $self->throw("invalid secondary namespace $secondary_ns");

  # first get whatever secondary ids are already stored there
  my @primary = $self->unpack_secondary($db->{$secondary_id});
  # uniqueify
  my %unique  = map {$_=>undef} @primary,$primary_id;

  my $record = $self->pack_secondary(keys %unique);
  $db->{$secondary_id} = $record;
}

# get output file handle
sub _outfh {
  my $self = shift;
#### XXXXX FINISH #####
#  my $
}

# unpack a primary record into fileid,offset,length
sub unpack_primary {
  my $self = shift;
  my $index_record = shift;
  return split "\t",$index_record;
}

# unpack a secondary record into a list of primary ids
sub unpack_secondary {
  my $self = shift;
  my $index_record = shift or return;
  return split "\t",$index_record;
}

# pack a list of fileid,offset,length into a primary id record
sub pack_primary {
  my $self = shift;
  my ($fileid,$offset,$length) = @_;
  return join "\t",($fileid,$offset,$length);
}

# pack a list of primary ids into a secondary id record
sub pack_secondary {
  my $self = shift;
  my @secondaries = @_;
  return join "\t",@secondaries;
}

sub primary_db {
  my $self = shift;
  # lazy opening
  $self->_open_bdb unless exists $self->{bdb_primary_db};
  return $self->{bdb_primary_db};
}

sub secondary_db {
  my $self = shift;
  my $secondary_namespace = shift
    or $self->throw("usage: secondary_db(\$secondary_namespace)");
  $self->_open_bdb unless exists $self->{bdb_primary_db};
  return $self->{bdb_secondary_db}{$secondary_namespace};
}

sub _open_bdb {
  my $self = shift;

  my $flags = $self->write_flag ? O_CREAT|O_RDWR : O_RDONLY;

  my $primary_db = {};
  tie(%$primary_db,'DB_File',$self->_catfile($self->_primary_db_name),$flags,0666,$DB_BTREE)
    or $self->throw("Could not open primary index file: $! (did you remember to use -write_flag=>1?)");
  $self->{bdb_primary_db} = $primary_db;

  for my $secondary ($self->secondary_namespaces) {
    my $secondary_db = {};
    tie(%$secondary_db,'DB_File',$self->_catfile($self->_secondary_db_name($secondary)),$flags,0666,$DB_BTREE)
      or $self->throw("Could not open primary index file");
    $self->{bdb_secondary_db}{$secondary} = $secondary_db;
  }

  1;
}

sub _primary_db_name {
  my $self = shift;
  my $pns  = $self->primary_namespace or $self->throw('no primary namespace defined');
  return "key_$pns";
}

sub _secondary_db_name {
  my $self  = shift;
  my $sns   = shift;
  return "id_$sns";
}

sub _fhcache {
  my $self  = shift;
  my $path  = shift;
  my $write = shift;

  if (!$self->{bdb_fhcache}{$path}) {
    $self->{bdb_curopen} ||= 0;
    if ($self->{bdb_curopen} >= $self->{bdb_maxopen}) {
      my @lru = sort {$self->{bdb_cacheseq}{$a} <=> $self->{bdb_cacheseq}{$b};} keys %{$self->{bdb_fhcache}};
      splice(@lru, $self->{bdb_maxopen} / 3);
      $self->{bdb_curopen} -= @lru;
      for (@lru) { delete $self->{bdb_fhcache}{$_} }
    }
    if ($write) {
      my $modifier = $self->{bdb_fhcache_seenit}{$path}++ ? '>' : '>>';
      $self->{bdb_fhcache}{$path} = IO::File->new("${modifier}${path}") or return;
    } else {
      $self->{bdb_fhcache}{$path} = IO::File->new($path) or return;
    }
    $self->{bdb_curopen}++;
  }
  $self->{bdb_cacheseq}{$path}++;
  $self->{bdb_fhcache}{$path}
}

1;
