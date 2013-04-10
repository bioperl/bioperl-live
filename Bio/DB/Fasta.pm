#
# BioPerl module for Bio::DB::Fasta
#
# You may distribute this module under the same terms as perl itself
#

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Fasta -- Fast indexed access to a directory of fasta files

=head1 SYNOPSIS

  use Bio::DB::Fasta;

  # create database from directory of fasta files
  my $db      = Bio::DB::Fasta->new('/path/to/fasta/files');

  # simple access (for those without Bioperl)
  my $seq      = $db->seq('CHROMOSOME_I',4_000_000 => 4_100_000);
  my $revseq   = $db->seq('CHROMOSOME_I',4_100_000 => 4_000_000);
  my @ids     = $db->ids;
  my $length   = $db->length('CHROMOSOME_I');
  my $alphabet = $db->alphabet('CHROMOSOME_I');
  my $header   = $db->header('CHROMOSOME_I');

  # Bioperl-style access
  my $db      = Bio::DB::Fasta->new('/path/to/fasta/files');

  my $obj     = $db->get_Seq_by_id('CHROMOSOME_I');
  my $seq     = $obj->seq; # sequence string
  my $subseq  = $obj->subseq(4_000_000 => 4_100_000); # string
  my $trunc   = $obj->trunc(4_000_000 => 4_100_000); # seq object
  my $length  = $obj->length;
  # (etc)

  # Bio::SeqIO-style access
  my $stream  = Bio::DB::Fasta->new('/path/to/files')->get_PrimarySeq_stream;
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
using the AnyDBM module (Berkeley DB* preferred, followed by GDBM_File,
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

=item $db = Bio::DB::Fasta-E<gt>new($fasta_path [,%options])

Create a new Bio::DB::Fasta object from the Fasta file or files
indicated by $fasta_path.  Indexing will be performed automatically if
needed.  If successful, new() will return the database accessor
object.  Otherwise it will return undef.

$fasta_path may be an individual Fasta file, or may refer to a
directory containing one or more of such files.  Following the path,
you may pass a series of name=E<gt>value options or a hash with these
same name=E<gt>value pairs.  Valid options are:

 Option Name   Description               Default
 -----------   -----------               -------

 -glob         Glob expression to use    *.{fa,fasta,fast,FA,FASTA,FAST,dna}
               for searching for Fasta
	            files in directories. 

 -makeid       A code subroutine for     None
	            transforming Fasta IDs.

 -maxopen      Maximum size of		     32
	            filehandle cache.

 -debug        Turn on status		        0
	            messages.

 -reindex      Force the index to be     0
               rebuilt.

 -dbmargs      Additional arguments      none
               to pass to the DBM
               routines when tied
               (scalar or array ref).

-dbmargs can be used to control the format of the index.  For example,
you can pass $DB_BTREE to this argument so as to force the IDs to be
sorted and retrieved alphabetically.  Note that you must use the same
arguments every time you open the index!

-reindex can be used to force the index to be recreated from scratch.

=item $fh = Bio::DB::Fasta-E<gt>newFh($fasta_path [,%options])

Create a tied filehandle opened on a Bio::DB::Fasta object.  Reading
from this filehandle with E<lt>E<gt> will return a stream of sequence objects,
Bio::SeqIO style.

=back

The -makeid option gives you a chance to modify sequence IDs during
indexing.  The option value should be a code reference that will
take a scalar argument and return a scalar result, like this:

  $db = Bio::DB::Fasta->new("file.fa",-makeid=>\&make_my_id);

  sub make_my_id {
    my $description_line = shift;
    # get a different id from the fasta header, e.g.
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

=item $type  = $db-E<gt>alphabet($id)

Return the molecular type of the indicated sequence.  One of "dna",
"rna" or "protein".

=item $filename  = $db-E<gt>file($id)

Return the name of the file in which the indicated sequence can be
found.

=item $offset    = $db-E<gt>offset($id)

Return the offset of the indicated sequence from the beginning of the
file in which it is located.  The offset points to the beginning of
the sequence, not the beginning of the header line.

=item $header_length = $db-E<gt>headerlen($id)

Return the length of the header line for the indicated sequence.

=item $header_offset = $db-E<gt>header_offset($id)

Return the offset of the header line for the indicated sequence from
the beginning of the file in which it is located.

=item $index_name  = $db-E<gt>index_name

Return the path to the index file.

=item $path = $db-E<gt>path

Return the path to the Fasta file(s).

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

#'
package Bio::DB::Fasta;

BEGIN {
  @AnyDBM_File::ISA = qw(DB_File GDBM_File NDBM_File SDBM_File)
}

use strict;
use IO::File;
use AnyDBM_File;
use Fcntl;
use File::Glob ':glob';
use File::Basename qw(basename dirname);

use base qw(Bio::DB::SeqI Bio::Root::Root);

*seq = *sequence = \&subseq;
*ids = \&get_all_ids;
*get_seq_by_primary_id = *get_Seq_by_acc  = \&get_Seq_by_id;

use constant STRUCT =>'NNnnCa*';
use constant STRUCTBIG =>'QQnnCa*'; # 64-bit file offset and seq length
use constant DNA     => 1;
use constant RNA     => 2;
use constant PROTEIN => 3;
use constant DIE_ON_MISSMATCHED_LINES => 1; # if you want 

# Bio::DB-like object
# providing fast random access to a directory of FASTA files

=head2 new

 Title   : new
 Usage   : my $db = Bio::DB::Fasta->new( $path, @options);
 Function: initialize a new Bio::DB::Fasta object
 Returns : new Bio::DB::Fasta object
 Args    : path to dir of fasta files or a single filename

These are optional arguments to pass in as well.

 -glob         Glob expression to use    *.{fa,fasta,fast,FA,FASTA,FAST}
               for searching for Fasta
	             files in directories. 

 -makeid       A code subroutine for     none
	             transforming Fasta IDs.

 -maxopen      Maximum size of		       32
	             filehandle cache.

 -debug        Turn on status		         0
	             messages.

 -reindex      Force the index to be     0
               rebuilt.

 -dbmargs      Additional arguments      none
               to pass to the DBM
               routines when tied
               (scalar or array ref).

=cut

sub new {
  my $class = shift;
  my $path  = shift;
  my %opts  = @_;

  my $self = bless { debug      => $opts{-debug},
	  makeid     => $opts{-makeid},
	  glob       => $opts{-glob}    || '*.{fa,fasta,FA,FASTA,fast,FAST,dna,FNA,fna,FAA,faa,FSA,fsa}',
	  maxopen    => $opts{-maxopen} || 32,
	  dbmargs    => $opts{-dbmargs} || undef,
	  fhcache    => {},
	  cacheseq   => {},
	  curopen    => 0,
	  openseq    => 1,
	  dirname    => undef,
	  offsets    => undef,
		   }, $class;
  my ($offsets,$dirname);

  if (-d $path) {
    # because Win32 glob() is broken with respect to long file names
    # that contain whitespace.
    $path = Win32::GetShortPathName($path)
      if $^O =~ /^MSWin/i && eval 'use Win32; 1';
    $offsets = $self->index_dir($path,$opts{-reindex}) or return;
    $dirname = $path;
  } elsif (-f _) {
    $offsets = $self->index_file($path,$opts{-reindex});
    $dirname = dirname($path);
  } else {
    $self->throw( "$path: Invalid file or dirname");
  }
  @{$self}{qw(dirname offsets)} = ($dirname,$offsets);

  $self;
}

=head2 newFh

 Title   : newFh
 Function: gets a new Fh for a file
 Example : internal method
 Returns : GLOB 
 Args    :

=cut

sub newFh {
  my $class = shift;
  my $self  = $class->new(@_);
  require Symbol;
  my $fh = Symbol::gensym or return;
  tie $$fh,'Bio::DB::Fasta::Stream',$self or return;
  $fh;
}

sub _open_index {
  my $self = shift;
  my ($index,$write) = @_;
  my %offsets;
  my $flags = $write ? O_CREAT|O_RDWR : O_RDONLY;
  my @dbmargs = $self->dbmargs;
  eval {
      tie %offsets,'AnyDBM_File',$index,$flags,0644,@dbmargs 
	  or die "Can't open sequence index file $index: $!";
  };
  warn $@ if $@;
  return \%offsets;
}

sub _close_index {
  my $self = shift;
  my $index = shift;
  untie %$index;
}

=head2 index_dir

 Title   : index_dir
 Usage   : $db->index_dir($dir)
 Function: set the index dir and load all files in the dir
 Returns : hashref of seq offsets in each file
 Args    : dirname, boolean to force a reload of all files

=cut

sub index_dir {
  my $self = shift;
  my $dir  = shift;
  my $force_reindex = shift;

  # find all fasta files
  my @files = glob("$dir/$self->{glob}");
  return unless @files;

  # get name of index
  my $index = $self->index_name($dir,1);

  # if caller has requested reindexing, then unlink
  # the index file.
  unlink $index if $force_reindex;

  # get the modification time of the index
  my $indextime   = 0;
  for my $suffix('','.pag','.dir') {
    $indextime ||= (stat("${index}${suffix}"))[9];
  }
  $indextime ||= 0;  # prevent some uninit variable warnings

  # get the most recent modification time of any of the contents
  my $modtime = 0;
  my %modtime;
  $self->set_pack_method( @files );
  foreach (@files) {
    my $m = (stat($_))[9];
    $modtime{$_} = $m;
    $modtime = $m if defined $m && $modtime < $m;
  }

  my $reindex      = $force_reindex || $indextime < $modtime;
  $self->{offsets} = $self->_open_index($index,$reindex) or return;

  # no indexing needed
  return $self->{offsets} unless $reindex;

  # otherwise reindex contents of changed files
  $self->{indexing} = $index;
  foreach (@files) {
    next if( defined $indextime && $modtime{$_} <= $indextime);
    $self->calculate_offsets($_,$self->{offsets});
  }
  delete $self->{indexing};

  # we've been having troubles with corrupted index files on Windows systems,
  # so possibly closing and reopening will help
  $self->_close_index($self->{offsets});

  return $self->{offsets}  = $self->_open_index($index);
}

=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : my $seq = $db->get_Seq_by_id($id)
 Function: Bio::DB::RandomAccessI method implemented
 Returns : Bio::PrimarySeqI object
 Args    : id

=cut

sub get_Seq_by_id {
  my $self = shift;
  my $id   = shift;
  return unless exists $self->{offsets}{$id};
  return Bio::PrimarySeq::Fasta->new($self,$id);
}

=head2 set_pack_method

 Title   : set_pack_method
 Usage   : $db->set_pack_method( @files )
 Function: Determines whether data packing uses 32 or 64 bit integers
 Returns :
 Args    : one or more file paths

=cut

sub set_pack_method {
  my $self = shift;
  # Find the maximum file size:eq
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
}

=head2 index_file

 Title   : index_file
 Usage   : $db->index_file($filename)
 Function: (re)loads a sequence file and indexes sequences offsets in the file
 Returns : seq offsets in the file
 Args    : filename, 
           boolean to force reloading a file

=cut


sub index_file {
  my $self = shift;
  my $file = shift;
  my $force_reindex = shift;

  $self->set_pack_method( $file );
  my $index = $self->index_name($file);
  # if caller has requested reindexing, then unlink the index
  unlink $index if $force_reindex;

  # get the modification time of the index
  my $indextime = (stat($index))[9] || 0;
  my $modtime   = (stat($file))[9]  || 0;

  my $reindex = $force_reindex || $indextime < $modtime;
  my $offsets = $self->_open_index($index,$reindex) or return;
  $self->{offsets} = $offsets;

  return $self->{offsets} unless $reindex;

  $self->{indexing} = $index;
  $self->calculate_offsets($file,$offsets);
  delete $self->{indexing};
  return $self->{offsets};
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

=head2 index_name

 Title   : index_name
 Usage   : my $indexname = $db->index_name($path,$isdir);
 Function: returns the name of the index for a specific path 
 Returns : string
 Args    : path to check, 
           boolean if it is a dir

=cut

sub index_name {
  my $self  = shift;
  my ($path,$isdir) = @_;
  unless ($path) {
    my $dir = $self->{dirname} or return;
    return $self->index_name($dir,-d $dir);
  } 
  return "$path/directory.index" if $isdir;
  return "$path.index";
}

=head2 calculate_offsets

 Title   : calculate_offsets
 Usage   : $db->calculate_offsets($filename,$offsets);
 Function: calculates the sequence offsets in a file based on id
 Returns : offset hash for each file
 Args    : file to process
           $offsets - hashref of id to offset storage

=cut

sub calculate_offsets {
  my $self = shift;
  my ($file,$offsets) = @_;
  my $base = $self->path2fileno(basename($file));

  my $fh = IO::File->new($file) or $self->throw( "Can't open $file: $!");
  binmode $fh;
  warn "indexing $file\n" if $self->{debug};
  my ($offset,@id,$linelength,$type,$firstline,$count,
      $termination_length,$seq_lines,$last_line,%offsets);
  my ($l3_len,$l2_len,$l_len)=(0,0,0);

  while (<$fh>) {		# don't try this at home
    $termination_length ||= /\r\n$/ ? 2 : 1; # account for crlf-terminated Windows files
    next unless /\S/;
    if (index($_, ">") == 0) {
        if (/^>(\S+)/) {
          print STDERR "indexed $count sequences...\n" 
        if $self->{debug} && (++$count%1000) == 0;
          
        
          my $pos = tell($fh);
          if (@id) {
        my $seqlength    = $pos - $offset - length($_);
        $seqlength      -= $termination_length * $seq_lines;
        my $ppos = &{$self->{packmeth}}($offset,$seqlength,
                           $linelength,$firstline,
                           $type,$base);
        for my $id (@id) { $offsets->{$id}  = $ppos }
          }
          @id = ref($self->{makeid}) eq 'CODE' ? $self->{makeid}->($_) : $1;
          ($offset,$firstline,$linelength) = ($pos,length($_),0);
          $self->_check_linelength($linelength);
          ($l3_len,$l2_len,$l_len)=(0,0,0);
          $seq_lines = 0;
        } else {
          # catch bad header lines, bug 3172
          $self->throw("FASTA header doesn't match '>(\\S+)': $_")
        }
    } else {
      $l3_len= $l2_len; $l2_len= $l_len; $l_len= length($_); # need to check every line :(
      if (DIE_ON_MISSMATCHED_LINES &&
	  $l3_len>0 && $l2_len>0 && $l3_len!=$l2_len) {
	  my $fap= substr($_,0,20)."..";
	  $self->throw("Each line of the fasta entry must be the same length except the last.
   Line above #$. '$fap' is $l2_len != $l3_len chars.");
      }
      $linelength ||= length($_);
      $type       ||= $self->_type($_);
      $seq_lines++;
    }
    $last_line = $_;
  }

  $self->_check_linelength($linelength);
  # deal with last entry
  if (@id) {
    my $pos = tell($fh);
    my $seqlength   = $pos - $offset;
    if ($linelength == 0) { # yet another pesky empty chr_random.fa file
      $seqlength = 0;
    } else {
      if ($last_line !~ /\s$/) {
	$seq_lines--;
      }
      $seqlength -= $termination_length * $seq_lines;
    };
    my $ppos = &{$self->{packmeth}}($offset,$seqlength,
					   $linelength,$firstline,
					   $type,$base);
    for my $id (@id) { $offsets->{$id}  = $ppos }
  }
  $offsets->{__termination_length} = $termination_length;
  return \%offsets;
}

=head2 get_all_ids

 Title   : get_all_ids
 Usage   : my @ids = $db->get_all_ids
 Function: gets all the stored ids in all indexes
 Returns : list of ids
 Args    : none

=cut

sub get_all_ids  { grep {!/^__/} keys %{shift->{offsets}} }

sub offset {
  my $self = shift;
  my $id   = shift;
  my $offset = $self->{offsets}{$id} or return;
  (&{$self->{unpackmeth}}($offset))[0];
}

sub length {
  my $self = shift;
  my $id   = shift;
  my $offset = $self->{offsets}{$id} or return;
  (&{$self->{unpackmeth}}($offset))[1];
}

sub linelen {
  my $self = shift;
  my $id   = shift;
  my $offset = $self->{offsets}{$id} or return;
  (&{$self->{unpackmeth}}($offset))[2];
}

sub headerlen {
  my $self = shift;
  my $id   = shift;
  my $offset = $self->{offsets}{$id} or return;
  (&{$self->{unpackmeth}}($offset))[3];
}

sub alphabet {
  my $self = shift;
  my $id   = shift;
  my $offset = $self->{offsets}{$id} or return;
  my $type = (&{$self->{unpackmeth}}($offset))[4];
  return $type == DNA ? 'dna'
         : $type == RNA ? 'rna'
         : 'protein';

}

sub path { shift->{dirname} } 

sub header_offset {
    my $self = shift;
    my $id   = shift;
    return unless $self->{offsets}{$id};
    return $self->offset($id) - $self->headerlen($id);
}

sub file {
  my $self = shift;
  my $id   = shift;
  my $offset = $self->{offsets}{$id} or return;
  $self->fileno2path((&{$self->{unpackmeth}}($offset))[5]);
}

sub fileno2path {
  my $self = shift;
  my $no   = shift;
  return $self->{offsets}{"__file_$no"};
}

sub path2fileno {
  my $self = shift;
  my $path = shift;
  if ( !defined $self->{offsets}{"__path_$path"} ) {
    my $fileno  = ($self->{offsets}{"__path_$path"} = 0+ $self->{fileno}++);
    $self->{offsets}{"__file_$fileno"} = $path;
  }
  return $self->{offsets}{"__path_$path"}
}

sub _check_linelength {
  my $self       = shift;
  my $linelength = shift;
  return unless defined $linelength;
  $self->throw("Each line of the fasta file must be less than 65,536 characters.  Line $. is $linelength chars.")	if $linelength > 65535.

}

=head2 subseq

 Title   : subseq
 Usage   : $seqdb->subseq($id,$start,$stop);
 Function: returns a subseq of a sequence in the db
 Returns : subsequence data
 Args    : id of sequence, starting point, ending point

=cut

sub subseq {
  my ($self,$id,$start,$stop) = @_;
  if ($id =~ /^(.+):([\d_]+)(?:,|-|\.\.)([\d_]+)$/) {
    ($id,$start,$stop) = ($1,$2,$3);
    $start =~ s/_//g;
    $stop =~ s/_//g;
  }
  $start ||= 1;
  $stop  ||= $self->length($id);

  my $reversed;
  if (defined $stop && $start > $stop) {
    ($start,$stop) = ($stop,$start);
    $reversed++;
  }

  my $data;

  my $fh = $self->fh($id) or return;
  my $filestart = $self->caloffset($id,$start);
  my $filestop  = $self->caloffset($id,$stop);

  seek($fh,$filestart,0);
  read($fh,$data,$filestop-$filestart+1);
  $data =~ s/\n//g;
  $data =~ s/\r//g;
  if ($reversed) {
    $data = reverse $data;
    $data =~ tr/gatcGATC/ctagCTAG/;
  }
  $data;
}

sub fh {
  my $self = shift;
  my $id   = shift;
  my $file = $self->file($id) or return;
  $self->fhcache("$self->{dirname}/$file") or $self->throw( "Can't open file $file");
}

sub header {
  my $self = shift;
  my $id   = shift;
  my ($offset,$seqlength,$linelength,$firstline,$type,$file) 
    = &{$self->{unpackmeth}}($self->{offsets}{$id}) or return;
  $offset -= $firstline;
  my $data;
  my $fh = $self->fh($id) or return;
  seek($fh,$offset,0);
  read($fh,$data,$firstline);
  chomp $data;
  substr($data,0,1) = '';
  $data;
}

sub caloffset {
  my $self = shift;
  my $id   = shift;
  my $a    = shift()-1;
  my ($offset,$seqlength,$linelength,$firstline,$type,$file) = &{$self->{unpackmeth}}($self->{offsets}{$id});
  $a = 0            if $a < 0;
  $a = $seqlength-1 if $a >= $seqlength;
  my $tl = $self->{offsets}{__termination_length};
  $offset + $linelength * int($a/($linelength-$tl)) + $a % ($linelength-$tl);
}

sub fhcache {
  my $self = shift;
  my $path = shift;
  if (!$self->{fhcache}{$path}) {
    if ($self->{curopen} >= $self->{maxopen}) {
      my @lru = sort {$self->{cacheseq}{$a} <=> $self->{cacheseq}{$b};} keys %{$self->{fhcache}};
      splice(@lru, $self->{maxopen} / 3);
      $self->{curopen} -= @lru;
      for (@lru) { delete $self->{fhcache}{$_} }
    }
    $self->{fhcache}{$path} = IO::File->new($path) or return;
    binmode $self->{fhcache}{$path};
    $self->{curopen}++;
  }
  $self->{cacheseq}{$path}++;
  $self->{fhcache}{$path}
}

sub _pack {
  pack STRUCT,@_;
}

sub _packBig {
  pack STRUCTBIG,@_;
}

sub _unpack {
  unpack STRUCT,shift;
}

sub _unpackBig {
  unpack STRUCTBIG,shift;
}

sub _type {
  shift;
  local $_ = shift;
  return /^[gatcnGATCN*-]+$/   ? DNA
         : /^[gaucnGAUCN*-]+$/ ? RNA
	 : PROTEIN;
}

=head2 get_PrimarySeq_stream

 Title   : get_PrimarySeq_stream
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub get_PrimarySeq_stream {
  my $self = shift;
  return Bio::DB::Fasta::Stream->new($self);
}

sub TIEHASH {
  my $self = shift;
  return $self->new(@_);
}

sub FETCH {
  shift->subseq(@_);
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
  defined shift->offset(@_);
}
sub FIRSTKEY { tied(%{shift->{offsets}})->FIRSTKEY(@_); }
sub NEXTKEY  { tied(%{shift->{offsets}})->NEXTKEY(@_);  }

sub DESTROY {
  my $self = shift;
  if ($self->{indexing}) {  # killed prematurely, so index file is no good!
    warn "indexing was interrupted, so unlinking $self->{indexing}";
    unlink $self->{indexing};
  }
}

#-------------------------------------------------------------
# Bio::PrimarySeqI compatibility
#
package Bio::PrimarySeq::Fasta;
use overload '""' => 'display_id';

use base qw(Bio::Root::Root Bio::PrimarySeqI);

sub new {
  my $class = shift;
  $class = ref($class) if ref $class;
  my ($db,$id,$start,$stop) = @_;
  return bless { db    => $db,
		 id    => $id,
		 start => $start || 1,
		 stop  => $stop  || $db->length($id)
	       },$class;
}

sub fetch_sequence { shift->seq(@_) }

sub seq {
  my $self = shift;
  return $self->{db}->seq($self->{id},$self->{start},$self->{stop});
}

sub subseq {
  my $self = shift;
  $self->trunc(@_)->seq();	
}

sub trunc {
  my $self = shift;
  my ($start,$stop) = @_;
  $self->throw("Stop cannot be smaller than start")  unless $start <= $stop;
  return $self->{start} <= $self->{stop} ?  $self->new($self->{db},
						       $self->{id},
						       $self->{start}+$start-1,
						       $self->{start}+$stop-1)
                                         :  $self->new($self->{db},
						       $self->{id},
						       $self->{start}-($start-1),
						       $self->{start}-($stop-1)
						      );  
	
}

sub is_circular {
  my $self = shift;
  return $self->{is_circular};
}

sub display_id {
  my $self = shift;
  return $self->{id};
}

sub accession_number {
  my $self = shift;
  return "unknown";
}

sub primary_id {
  my $self = shift;
  return overload::StrVal($self);
}

sub can_call_new { return 0 }

sub alphabet {
  my $self = shift;
  return $self->{db}->alphabet($self->{id});
}

sub revcom {
  my $self = shift;
  return $self->new(@{$self}{'db','id','stop','start'});
}

sub length {
  my $self = shift;
  #return $self->{db}->length($self->{id}); # wrong because ignores sequence start and stop values
  return length($self->seq);

}

sub description  { 
    my $self = shift;
    my $header = $self->{'db'}->header($self->{id});
    # remove the id from the header
    return (split(/\s+/,$header,2))[1];
}

*desc = \&description;

#-------------------------------------------------------------
# stream-based access to the database
#
package Bio::DB::Fasta::Stream;
use base qw(Tie::Handle Bio::DB::SeqI);


sub new {
  my $class = shift;
  my $db    = shift;
  my $key = $db->FIRSTKEY;
  return bless { db=>$db,key=>$key },$class;
}

sub next_seq {
  my $self = shift;
  my ($key,$db) = @{$self}{'key','db'};
  while ($key =~ /^__/) {
    $key = $db->NEXTKEY($key);
    return unless defined $key;
  }
  my $value = $db->get_Seq_by_id($key);
  $self->{key} = $db->NEXTKEY($key);
  $value;
}

sub TIEHANDLE {
  my $class = shift;
  my $db    = shift;
  return $class->new($db);
}
sub READLINE {
  my $self = shift;
  $self->next_seq;
}

1;

__END__

