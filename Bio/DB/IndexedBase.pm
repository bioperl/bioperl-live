#
# BioPerl module for Bio::DB::IndexedBase
#
# You may distribute this module under the same terms as perl itself
#


=head1 NAME

Bio::DB::IndexedBase - Base class for modules using indexed sequence files

=head1 SYNOPSIS

  use Bio::DB::Fasta;

  # create database from directory of fasta files
  my $db       = Bio::DB::Fasta->new('/path/to/fasta/files');
  my @ids     = $db->ids;

  # simple access (for those without Bioperl)
  my $seq      = $db->seq('CHROMOSOME_I',4_000_000 => 4_100_000);
  my $revseq   = $db->seq('CHROMOSOME_I',4_100_000 => 4_000_000);
  my $length   = $db->length('CHROMOSOME_I');
  my $header   = $db->header('CHROMOSOME_I');
  my $alphabet = $db->alphabet('CHROMOSOME_I');

  # Bioperl-style access
  my $obj     = $db->get_Seq_by_id('CHROMOSOME_I');
  my $seq     = $obj->seq;                            # string
  my $subseq  = $obj->subseq(4_000_000 => 4_100_000); # string
  my $trunc   = $obj->trunc(4_000_000 => 4_100_000);  # object
  my $length  = $obj->length;
  # etc

  # Bio::SeqIO-style access
  my $stream  = $db->get_PrimarySeq_stream;
  while (my $seq = $stream->next_seq) {
    # Bio::PrimarySeqI stuff
  }

  my $fh = Bio::DB::Fasta->newFh('/path/to/fasta/files');
  while (my $seq = <$fh>) {
    # Bio::PrimarySeqI stuff
  }

  # tied hash access
  tie %sequences,'Bio::DB::Fasta','/path/to/fasta/files';
  print $sequences{'CHROMOSOME_I:1,20000'};

=head1 DESCRIPTION

Bio::DB::Fasta provides indexed access to one or more Fasta files.  It
provides random access to each sequence entry, and to subsequences
within each entry, allowing you to retrieve portions of very large
sequences without bringing the entire sequence into memory.

When you initialize the module, you point it at a single fasta file or
a directory of multiple such files.  The first time it is run, the
module generates an index of the contents of the file or directory
using the AnyDBM_File module (BerkeleyDB preferred, followed by GDBM_File,
NDBM_File, and SDBM_File).  Thereafter it uses the index file to find
the file and offset for any requested sequence.  If one of the source
fasta files is updated, the module reindexes just that one file.  (You
can also force reindexing manually).  For improved performance, the
module keeps a cache of open filehandles, closing less-recently used
ones when the cache is full.

The fasta files may contain any combination of nucleotide and protein
sequences; during indexing the module guesses the molecular type.
Entries may have any line length up to 65,536 characters, and
different line lengths are allowed in the same file.  However, within
a sequence entry, all lines must be the same length except for the
last.

An error will be thrown if this is not the case.

The module uses /^E<gt>(\S+)/ to extract the primary ID of each sequence
from the Fasta header.  During indexing, you may pass a callback routine to
modify this primary ID.  For example, you may wish to extract a
portion of the gi|gb|abc|xyz nonsense that GenBank Fasta files use.
The original header line can be recovered later.

This module was developed for use with the C. elegans and human
genomes, and has been tested with sequence segments as large as 20
megabases.  Indexing the C. elegans genome (100 megabases of genomic
sequence plus 100,000 ESTs) takes ~5 minutes on my 300 MHz pentium
laptop. On the same system, average access time for any 200-mer within
the C. elegans genome was E<lt>0.02s.

*Berkeley DB can be obtained free from www.sleepycat.com. After it is
installed you will need to install the BerkeleyDB Perl module.

=head1 DATABASE CREATION AND INDEXING

The two constructors for this class are new() and newFh().  The former
creates a Bio::DB::Fasta object which is accessed via method calls.
The latter creates a tied filehandle which can be used Bio::SeqIO
style to fetch sequence objects in a stream fashion.  There is also a
tied hash interface.

=over 2

=item $db = Bio::DB::Fasta-E<gt>new($path [,%options])

Create a new Bio::DB::IndexedBase object from the file or files
indicated by $fasta_path.  Indexing will be performed automatically if
needed.  If successful, new() will return the database accessor
object.  Otherwise it will return undef.

$fasta_path may be an individual Fasta file, or may refer to a
directory containing one or more of such files.  Following the path,
you may pass a series of name=E<gt>value options or a hash with these
same name=E<gt>value pairs.  Valid options are:

 Option Name   Description               Default
 -----------   -----------               -------

 -glob         Glob expression to        *
               search for files in
               in directories.

 -makeid       A code subroutine for     None
               transforming Fasta IDs.

 -maxopen      Maximum size of           32
               filehandle cache.

 -debug        Turn on status            0
               messages.

 -reindex      Force the index to be     0
               rebuilt.

 -dbmargs      Additional arguments      None
               to pass to the DBM
               routines when tied
               (scalar or array ref).

-dbmargs can be used to control the format of the index.  For example,
you can pass $DB_BTREE to this argument so as to force the IDs to be
sorted and retrieved alphabetically.  Note that you must use the same
arguments every time you open the index!

-reindex can be used to force the index to be recreated from scratch.

The -makeid option gives you a chance to modify sequence IDs during
indexing.  The option value should be a code reference that will
take a scalar argument and return a scalar result, like this:

  $db = Bio::DB::Fasta->new("file.fa",-makeid=>\&make_my_id);

  sub make_my_id {
    my $description_line = shift;
    # get a different id from the header, e.g.
    $description_line =~ /(\S+)$/;
    return $1;
  }

make_my_id() will be called with the full fasta id line (including the
"E<gt>" symbol!).  For example:

 >A12345.3 Predicted C. elegans protein egl-2

By default, this module will use the regular expression /^E<gt>(\S+)/
to extract "A12345.3" for use as the ID.  If you pass a -makeid
callback, you can extract any portion of this, such as the "egl-2"
symbol.

The -makeid option is ignored after the index is constructed.

=item $fh = Bio::DB::Fasta-E<gt>newFh($fasta_path [,%options])

Create a tied filehandle opened on a Bio::DB::Fasta object.  Reading
from this filehandle with E<lt>E<gt> will return a stream of sequence objects,
Bio::SeqIO style.

=item $db->index_dir($dir)

Set the index directory and index the files within

=item $index_name = $db-E<gt>index_name

Return the path to the index file.

=item $path = $db-E<gt>path

Return the path to the Fasta file(s).

=back

=head1 OBJECT METHODS

The following object methods are provided.

=over 10

=item $raw_seq = $db-E<gt>seq($id [,$start, $stop])

Return the raw sequence (a string) given an ID and optionally a start
and stop position in the sequence.  In the case of DNA sequence, if
$stop is less than $start, then the reverse complement of the sequence
is returned (this violates Bio::Seq conventions).

For your convenience, subsequences can be indicated with any of the
following compound IDs:

   $db->seq("$id:$start,$stop")

   $db->seq("$id:$start..$stop")

   $db->seq("$id:$start-$stop")

=item $length = $db-E<gt>length($id)

Return the length of the indicated sequence.

=item $header = $db-E<gt>header($id)

Return the header line for the ID, including the initial "E<gt>".

=item $type = $db-E<gt>alphabet($id)

Return the molecular type of the indicated sequence.  One of "dna",
"rna" or "protein".

=item $filename = $db-E<gt>file($id)

Return the name of the file in which the indicated sequence can be
found.

=item $offset = $db-E<gt>offset($id)

Return the offset of the indicated sequence from the beginning of the
file in which it is located.  The offset points to the beginning of
the sequence, not the beginning of the header line.

=item $header_length = $db-E<gt>headerlen($id)

Return the length of the header line for the indicated sequence.

=item $header_offset = $db-E<gt>header_offset($id)

Return the offset of the header line for the indicated sequence from
the beginning of the file in which it is located.

=back

For BioPerl-style access, the following methods are provided:

=over 4

=item $seq = $db-E<gt>get_Seq_by_id($id)

Return a Bio::PrimarySeq::Fasta object, which obeys the
Bio::PrimarySeqI conventions.  For example, to recover the raw DNA or
protein sequence, call $seq-E<gt>seq().

Note that get_Seq_by_id() does not bring the entire sequence into
memory until requested.  Internally, the returned object uses the
accessor to generate subsequences as needed.

=item $seq = $db-E<gt>get_Seq_by_acc($id)

=item $seq = $db-E<gt>get_Seq_by_primary_id($id)

These methods all do the same thing as get_Seq_by_id().

=item $stream = $db-E<gt>get_PrimarySeq_stream()

Return a Bio::DB::Fasta::Stream object, which supports a single method
next_seq(). Each call to next_seq() returns a new
Bio::PrimarySeq::Fasta object, until no more sequences remain.

=back

See L<Bio::PrimarySeqI> for methods provided by the sequence objects
returned from get_Seq_by_id() and get_PrimarySeq_stream().

=head1 TIED INTERFACES

This module provides two tied interfaces, one which allows you to
treat the sequence database as a hash, and the other which allows you
to treat the database as an I/O stream.

=head2 Creating a Tied Hash

The tied hash interface is very straightforward

=over 1

=item $obj = tie %db,'Bio::DB::Fasta','/path/to/fasta/files' [,@args]

Tie %db to Bio::DB::Fasta using the indicated path to the Fasta files.
The optional @args list is the same set of named argument/value pairs
used by Bio::DB::Fasta-E<gt>new().

If successful, tie() will return the tied object.  Otherwise it will
return undef.

=back

Once tied, you can use the hash to retrieve an individual sequence by
its ID, like this:

  my $seq = $db{CHROMOSOME_I};

You may select a subsequence by appending the comma-separated range to
the sequence ID in the format "$id:$start,$stop".  For example, here
is the first 1000 bp of the sequence with the ID "CHROMOSOME_I":

  my $seq = $db{'CHROMOSOME_I:1,1000'};

(The regular expression used to parse this format allows sequence IDs
to contain colons.)

When selecting subsequences, if $start E<gt> stop, then the reverse
complement will be returned for DNA sequences.

The keys() and values() functions will return the sequence IDs and
their sequences, respectively.  In addition, each() can be used to
iterate over the entire data set:

 while (my ($id,$sequence) = each %db) {
    print "$id => $sequence\n";
 }

When dealing with very large sequences, you can avoid bringing them
into memory by calling each() in a scalar context.  This returns the
key only.  You can then use tied(%db) to recover the Bio::DB::Fasta
object and call its methods.

 while (my $id = each %db) {
    print "$id => $db{$sequence:1,100}\n";
    print "$id => ",tied(%db)->length($id),"\n";
 }

You may, in addition invoke Bio::DB::Fasta the FIRSTKEY and NEXTKEY tied
hash methods directly.

=over 2

=item $id = $db-E<gt>FIRSTKEY

Return the first ID in the database.

=item $id = $db-E<gt>NEXTKEY($id)

Given an ID, return the next ID in sequence.

=back

This allows you to write the following iterative loop using just the
object-oriented interface:

 my $db = Bio::DB::Fasta->new('/path/to/fasta/files');
 for (my $id=$db->FIRSTKEY; $id; $id=$db->NEXTKEY($id)) {
    # do something with sequence
 }

=head2 Creating a Tied Filehandle

The Bio::DB::Fasta-E<gt>newFh() method creates a tied filehandle from
which you can read Bio::PrimarySeq::Fasta sequence objects
sequentially.  The following bit of code will iterate sequentially
over all sequences in the database:

 my $fh = Bio::DB::Fasta->newFh('/path/to/fasta/files');
 while (my $seq = <$fh>) {
   print $seq->id,' => ',$seq->length,"\n";
 }

When no more sequences remain to be retrieved, the stream will return
undef.

=head1 BUGS

When a sequence is deleted from one of the Fasta files, this deletion
is not detected by the module and removed from the index.  As a
result, a "ghost" entry will remain in the index and will return
garbage results if accessed.

Currently, the only way to accommodate deletions is to rebuild the
entire index, either by deleting it manually, or by passing
-reindex=E<gt>1 to new() when initializing the module.

=head1 SEE ALSO

L<bioperl>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut


package Bio::DB::IndexedBase;

BEGIN {
  @AnyDBM_File::ISA = qw(DB_File GDBM_File NDBM_File SDBM_File)
}

use strict;
use IO::File;
use AnyDBM_File;
use Fcntl;
use File::Spec;
use File::Basename qw(dirname);

use base qw(Bio::Root::Root);

*seq = *sequence = \&subseq;
*ids = \&get_all_ids;
*get_seq_by_primary_id = *get_Seq_by_acc = \&get_Seq_by_id;

use constant STRUCT    =>'NNnnCa*';
use constant STRUCTBIG =>'QQnnCa*'; # 64-bit file offset and seq length

use constant NA        => 0;
use constant DNA       => 1;
use constant RNA       => 2;
use constant PROTEIN   => 3;

use constant DIE_ON_MISSMATCHED_LINES => 1; # if you want


=head2 new

 Title   : new
 Usage   : my $db = Bio::DB::Fasta->new( $path, @options);
 Function: initialize a new Bio::DB::Fasta object
 Returns : new Bio::DB::Fasta object
 Args    : a single file, or path to dir, or arrayref of files

These are optional arguments to pass in as well.

 -glob         Glob expression to use    *.{fa,fasta,fast,FA,FASTA,FAST}
               for searching for Fasta
               files in directories.

 -makeid       A code subroutine for     None
               transforming Fasta IDs.

 -maxopen      Maximum size of           32
               filehandle cache.

 -debug        Turn on status            0
               messages.

 -reindex      Force the index to be     0
               rebuilt.

 -dbmargs      Additional arguments      None
               to pass to the DBM
               routines when tied
               (scalar or array ref).

=cut

##sub new {
##
##}


=head2 newFh

 Title   : newFh
 Usage   : my $fh = Bio::DB::Qual->newFh('/path/to/qual/files', %options);
 Function: Get a new Fh for a file or directory containing several files
 Returns : filehandle object
 Args    : Fasta filename and options

=cut

sub newFh {
    my $class = shift;
    my $self  = $class->new(@_);
    require Symbol;
    my $fh = Symbol::gensym or return;
    tie $$fh,'Bio::DB::Indexed::Stream',$self or return;
    return $fh;
}


=head2 dbmargs

 Title   : dbmargs
 Usage   : my @args = $db->dbmargs;
 Function: gets stored dbm arguments
 Returns : array
 Args    : none

=cut

sub dbmargs {
    my $self = shift;
    my $args = $self->{dbmargs} or return;
    return ref($args) eq 'ARRAY' ? @$args : $args;
}


sub _open_index {
    my ($self, $index_file, $write) = @_;
    my %offsets;
    my $flags = $write ? O_CREAT|O_RDWR : O_RDONLY;
    my @dbmargs = $self->dbmargs;
    tie %offsets, 'AnyDBM_File', $index_file, $flags, 0644, @dbmargs 
        or $self->throw( "Could not open index file $index_file: $!");
    return \%offsets;
}


sub _close_index {
    my ($self, $index) = @_;
    untie %$index;
    return 1;
}


=head2 set_pack_method

 Title   : set_pack_method
 Usage   : $db->set_pack_method( @files )
 Function: Determines whether data packing uses 32 or 64 bit integers
 Returns : 1 for success
 Args    : one or more file paths

=cut

sub set_pack_method {
    my $self = shift;
    # Find the maximum file size:
    my ($maxsize) = sort { $b <=> $a } map { -s $_ } @_;
    my $fourGB    = (2 ** 32) - 1;

    if ($maxsize > $fourGB) {
        # At least one file exceeds 4Gb - we will need to use 64 bit ints
        $self->{packmeth}   = \&_packBig;
        $self->{unpackmeth} = \&_unpackBig;
    } else {
        $self->{packmeth}   = \&_pack;
        $self->{unpackmeth} = \&_unpack;
    }
    return 1;
}


sub _pack {
    return pack STRUCT, @_;
}


sub _packBig {
    return pack STRUCTBIG, @_;
}


sub _unpack {
    return unpack STRUCT, shift;
}


sub _unpackBig {
    return unpack STRUCTBIG, shift;
}


sub _type {
  my ($self, $string) = @_;
  return $string =~ m/^[gatcnGATCN*-]+$/   ? DNA
         : $string =~ m/^[gaucnGAUCN*-]+$/ ? RNA
         : PROTEIN;
}


=head2 index_dir

 Title   : index_dir
 Usage   : $db->index_dir($dir)
 Function: Index the files in the given directory that match -glob
 Returns : hashref of offsets
 Args    : dirname
           boolean to force a reload of all files

=cut

sub index_dir {
    my ($self, $dir, $force_reindex) = @_;
    my @files = glob( File::Spec->catfile($dir, $self->{glob}) );
    $self->throw("No suitable files found in $dir") if scalar @files == 0;
    my $index = $self->index_name( File::Spec->catfile($dir, 'directory.index') );
    my $offsets = $self->index_files(\@files, $force_reindex);

    return $offsets;
}


=head2 index_file

 Title   : index_file
 Usage   : $db->index_file($filename)
 Function: Index the given file
 Returns : hashref of offsets
 Args    : filename
           boolean to force reloading a file

=cut

sub index_file {
    my ($self, $file, $force_reindex) = @_;
    my $index = $self->index_name( "$file.index" );
    my $offsets = $self->index_files([$file], $force_reindex);
    return $self->{offsets};
}


=head2 index_files

 Title   : index_files
 Usage   : $db->index_files(\@files)
 Function: Index the given files
 Returns : hashref of offsets
 Args    : arrayref of filnames
           boolean to force a reload of all files

=cut

sub index_files {
    my ($self, $files, $force_reindex) = @_;

    $self->set_pack_method( @$files );

    # get name of index file, set it if it does not exist
    my $index = $self->index_name;
    if (not defined $index) {
        $self->index_name( 'files.index' );
    }

    # if caller has requested reindexing, unlink the index file.
    unlink $index if $force_reindex;

    # get the modification time of the index
    my $indextime = (stat $index)[9] || 0;

    # list recently updated files
    my $modtime = 0;
    my @updated;
    for my $file (@$files) {
        my $m = (stat $file)[9] || 0;
        if ($m > $modtime) {
           $modtime = $m;
        }
        if ($m > $indextime) {
           push @updated, $file;
        }
    }

    my $reindex      = $force_reindex || (scalar @updated > 1);
    $self->{offsets} = $self->_open_index($index, $reindex) or return;

    if ($reindex) {
        # reindex contents of changed files
        $self->{indexing} = $index;
        for my $file (@updated) {
            $self->calculate_offsets($file, $self->{offsets});
        }
        delete $self->{indexing};
    }

    # closing and reopening might help corrupted index file problem on Windows
    $self->_close_index($self->{offsets});

    return $self->{offsets} = $self->_open_index($index);
}


=head2 index_name

 Title   : index_name
 Usage   : my $indexname = $db->index_name($path);
 Function: Get or set the path of the index file
 Returns : string
 Args    : none, or path to set

=cut

sub index_name {
    my ($self, $path) = @_;
    if (defined $path) {
        $self->{index_name} = $path;
    }
    return $self->{index_name};
}


sub path {
    return shift->{dirname};
}


sub _check_linelength {
    my ($self, $linelength) = @_;
    return unless defined $linelength;
    $self->throw(
        "Each line of the qual file must be less than 65,536 characters.Line ".
        "$. is $linelength chars."
    ) if $linelength > 65535;
}


sub _caloffset {
    my ($self, $id, $a) = @_;
    $a--;
    my ($offset,$seqlength,$linelength,$firstline,$file)
        = &{$self->{unpackmeth}}($self->{offsets}{$id});
    $a = 0            if $a < 0;
    $a = $seqlength-1 if $a >= $seqlength;
    my $tl = $self->{offsets}{__termination_length};
    return $offset + $linelength * int($a/($linelength-$tl)) + $a % ($linelength-$tl);
}


sub _fhcache {
    my ($self, $path) = @_;
    if (!$self->{fhcache}{$path}) {
        if ($self->{curopen} >= $self->{maxopen}) {
            my @lru = sort {$self->{cacheseq}{$a} <=> $self->{cacheseq}{$b};}
                keys %{$self->{fhcache}};
            splice(@lru, $self->{maxopen} / 3);
            $self->{curopen} -= @lru;
            for (@lru) {
                delete $self->{fhcache}{$_};
            }
        }
        $self->{fhcache}{$path} = IO::File->new($path) || return;
        binmode $self->{fhcache}{$path};
        $self->{curopen}++;
    }
    $self->{cacheseq}{$path}++;
    return $self->{fhcache}{$path};
}


sub _fileno2path {
    my ($self, $no) = @_;
    return $self->{offsets}{"__file_$no"};
}


sub _path2fileno {
    my ($self, $path) = @_;
    if ( !defined $self->{offsets}{"__path_$path"} ) {
        my $fileno = ($self->{offsets}{"__path_$path"} = 0+ $self->{fileno}++);
        $self->{offsets}{"__file_$fileno"} = $path;
    }
    return $self->{offsets}{"__path_$path"}
}


#-------------------------------------------------------------
# Tied hash logic
#

sub TIEHASH {
  my $self = shift;
  return $self->new(@_);
}

sub FETCH {
  return shift->subseq(@_);
}

sub STORE {
  shift->throw("Read-only database");
}

sub DELETE {
  shift->throw("Read-only database");
}

sub CLEAR {
  shift->throw("Read-only database");
}

sub EXISTS {
  return defined shift->offset(@_);
}

sub FIRSTKEY {
  return tied(%{shift->{offsets}})->FIRSTKEY(@_);
}

sub NEXTKEY {
  return tied(%{shift->{offsets}})->NEXTKEY(@_);
}

sub DESTROY {
  my $self = shift;
  if ($self->{indexing}) {  # killed prematurely, so index file is no good!
    warn "indexing was interrupted, so deleting $self->{indexing}";
    unlink $self->{indexing};
  }
  return 1;
}


#-------------------------------------------------------------
# stream-based access to the database
#

package Bio::DB::Indexed::Stream;
use base qw(Tie::Handle Bio::DB::SeqI);


sub new {
    my ($class, $db) = @_;
    my $key = $db->FIRSTKEY;
    return bless {
        db  => $db,
        key => $key
    }, $class;
}

sub next_seq {
    my $self = shift;
    my ($key, $db) = @{$self}{'key', 'db'};
    return if not defined $key;
    while ($key =~ /^__/) {
        $key = $db->NEXTKEY($key);
        return if not defined $key;
    }
    my $value = $db->get_Seq_by_id($key);
    $self->{key} = $db->NEXTKEY($key);
    return $value;
}

sub TIEHANDLE {
    my ($class, $db) = @_;
    return $class->new($db);
}

sub READLINE {
    my $self = shift;
    return $self->next_seq;
}


1;

__END__
