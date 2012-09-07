#
# BioPerl module for Bio::DB::Qual
#
# You may distribute this module under the same terms as perl itself
#

=head1 NAME

Bio::DB::Qual - Fast indexed access to quality files

=head1 SYNOPSIS

  use Bio::DB::Qual;

  # create database from directory of qual files
  my $db      = Bio::DB::Qual->new('/path/to/qual/files');
  my @ids     = $db->ids;

  # simple access (for those without Bioperl)
  my @qual    = @{$db->qual('CHROMOSOME_I',4_000_000 => 4_100_000)};
  my @revqual = @{$db->qual('CHROMOSOME_I',4_100_000 => 4_000_000)};
  my $length  = $db->length('CHROMOSOME_I');
  my $header  = $db->header('CHROMOSOME_I');

  # Bioperl-style access
  my $obj     = $db->get_Qual_by_id('CHROMOSOME_I');
  my @qual    = @{$obj->qual};
  my @subqual = @{$obj->subqual(4_000_000 => 4_100_000)};
  my $length  = $obj->length;
  # etc

  # Bio::SeqIO-style access
  my $stream  = $db->get_PrimaryQual_stream;
  while (my $qual = $stream->next_seq) {
    # Bio::Seq::PrimaryQual operations
  }

  my $fh = Bio::DB::Qual->newFh('/path/to/qual/files');
  while (my $qual = <$fh>) {
    # Bio::Seq::PrimaryQual operations
  }

  # tied hash access
  tie %qualities,'Bio::DB::Qual','/path/to/qual/files';
  print $qualities{'CHROMOSOME_I:1,20000'};

=head1 DESCRIPTION

Bio::DB::Qual provides indexed access to one or more qual files. It provides
random access to each quality score entry without having to read the file from
the beginning. Access to subqualities (portions of a quality score) is provided,
although contrary to Bio::DB::Fasta, the full quality score has to be brought in
memory.

When you initialize the module, you point it at a single qual file or a
directory of multiple such files. The first time it is run, the module generates
an index of the contents of the file or directory using the AnyDBM module
(Berkeley DB* preferred, followed by GDBM_File, NDBM_File, and SDBM_File).
Thereafter it uses the index file to find the file and offset for any requested
quality score. If one of the source qual files is updated, the module reindexes
just that one file. (You can also force reindexing manually). For improved
performance, the module keeps a cache of open filehandles, closing less-recently
used ones when the cache is full.

The qual files may contain decimal quality scores. Entries may have any line
length up to 65,536 characters, and different line lengths are allowed in the
same file. However, within a quality score entry, all lines must be the same
length except for the last. An error will be thrown if this is not the case.

The module uses /^E<gt>(\S+)/ to extract the primary ID of each quality score
from the qual header. During indexing, you may pass a callback routine to modify
this primary ID.  For example, you may wish to extract a portion of the
gi|gb|abc|xyz prefixes that are commonly used. The original header line can be
recovered later.

*Berkeley DB can be obtained free from www.sleepycat.com. After it is installed
you will need to install the BerkeleyDB Perl module.

=head1 DATABASE CREATION AND INDEXING

The two constructors for this class are new() and newFh(). The former creates a
Bio::DB::Qual object which is accessed via method calls. The latter creates a
tied filehandle which can be used Bio::SeqIO-style to fetch quality score
objects in a data stream. There is also a tied hash interface.

=over 2

=item $db = Bio::DB::Qual-E<gt>new($qual_path [,%options])

Create a new Bio::DB::Qual object from the Qual file or files indicated by
$qual_path. Indexing will be performed automatically if needed. If successful,
new() will return the database accessor object. Otherwise it will return undef.

$qual_path may be an individual qual file, or may refer to a directory
containing one or more of such files. Following the path, you may pass a series
of name=E<gt>value options or a hash with these same name=E<gt>value pairs. 
Valid options are:

 Option Name   Description               Default
 -----------   -----------               -------

 -glob         Glob expression to use    *.{qa,QA,qual,QUAL}
               for searching for qual
               files in directories. 

 -makeid       A code subroutine for     None
               transforming qual IDs.

 -maxopen      Maximum size of           32
               filehandle cache.

 -debug        Turn on status            0
               messages.

 -reindex      Force the index to be     0
               rebuilt.

 -dbmargs      Additional arguments      none
               to pass to the DBM
               routines when tied
               (scalar or array ref).

 -index_name   Name of the file that     Auto
               holds the indexing
               information.

-dbmargs can be used to control the format of the index. For example, you can
pass $DB_BTREE to this argument so as to force the IDs to be sorted and
retrieved alphabetically. Note that you must use the same arguments every time
you open the index!

-reindex can be used to force the index to be recreated from scratch.

The -makeid option gives you a chance to modify quality score IDs during
indexing. The option value should be a code reference that will take a scalar
argument and return a scalar result, like this:

  $db = Bio::DB::Qual->new("file.qual",-makeid=>\&make_my_id);

  sub make_my_id {
    my $description_line = shift;
    # get a different id from the header, e.g.
    $description_line =~ /(\S+)$/;
    return $1;
  }

make_my_id() will be called with the full qual id line (including the "E<gt>"
symbol!). For example:

  >A12345.3 Predicted C. elegans protein egl-2

By default, this module will use the regular expression /^E<gt>(\S+)/ to extract
"A12345.3" for use as the ID.If you pass a -makeid callback, you can extract any
portion of this, such as the "egl-2" symbol.

The -makeid option is ignored after the index is constructed.

=item $fh = Bio::DB::Qual-E<gt>newFh($qual_path [,%options])

Create a tied filehandle opened on a Bio::DB::Qual object. Reading from this
filehandle with E<lt>E<gt> will return a stream of quality objects,
Bio::SeqIO-style.

=back

=item $index_name  = $db-E<gt>index_name

Return the path to the index file.

=item $path = $db-E<gt>path

Return the path to the Qual file(s).

=head1 OBJECT METHODS

The following object methods are provided.

=over 10

=item $raw_qual = $db-E<gt>qual($id [,$start, $stop])

Return a quality score array reference given an ID and optionally a start and
stop position (the quality value number) in the quality score. If $stop is less
than $start, then the reverse complement of the quality score is returned (this
violates Bio::Seq conventions).

For your convenience, subqualities can be indicated with any of the following
compound IDs:

   $db->qual("$id:$start,$stop")

   $db->qual("$id:$start..$stop")

   $db->qual("$id:$start-$stop")

=item $length = $db-E<gt>length($id)

Return the length of the indicated quality score, i.e. the number of quality
values.

=item $header = $db-E<gt>header($id)

Return the header line for the ID, including the initial "E<gt>".

=item $filename  = $db-E<gt>file($id)

Return the name of the file in which the indicated quality score can be found.

=item $offset    = $db-E<gt>offset($id)

Return the offset of the indicated quality score from the beginning of the file
in which it is located.  The offset points to the beginning of the quality
score, not the beginning of the header line.

=item $header_length = $db-E<gt>headerlen($id)

Return the length of the header line for the indicated quality score.

=item $header_offset = $db-E<gt>header_offset($id)

Return the offset of the header line for the indicated quality score from the
beginning of the file in which it is located.

=back

For BioPerl-style access, the following methods are provided:

=over 4

=item $qual = $db-E<gt>get_Qual_by_id($id)

Return a Bio::Seq::PrimaryQual object, which obeys the Bio::PrimarySeqI 
conventions. To recover the quality score, call $qual-E<gt>qual().

Note that get_Qual_by_id() does not bring the entire quality score into memory
until requested. Internally, the returned object uses the accessor to generate
subqualities as needed.

=item $qual = $db-E<gt>get_Qual_by_acc($id)

=item $qual = $db-E<gt>get_Qual_by_primary_id($id)

These methods all do the same thing as get_Qual_by_id().

=item $stream = $db-E<gt>get_PrimaryQual_stream()

Return a Bio::DB::Indexed::Stream object, which supports a single method 
next_seq(). Each call to next_seq() returns a new Bio::Seq::PrimaryQual object,
until no more quality scores remain.

=back

See L<Bio::Seq::PrimaryQual> and L<Bio::PrimarySeqI> for methods provided by the
quality objects returned from get_Qual_by_id() and get_PrimaryQual_stream().

=head1 TIED INTERFACES

This module provides two tied interfaces, one which allows you to treat the
quality score database as a hash, and the other which allows you to treat the
database as an I/O stream.

=head2 Creating a Tied Hash

The tied hash interface is very straightforward.

=over 1

=item $obj = tie %db,'Bio::DB::Qual','/path/to/qual/files' [,@args]

Tie %db to Bio::DB::Qual using the indicated path to the Qual files. The
optional @args list is the same set of named argument/value pairs used by
Bio::DB::Qual-E<gt>new().

If successful, tie() will return the tied object.  Otherwise it will return
undef.

=back

Once tied, you can use the hash to retrieve an individual quality score by its
ID, like this:

  my $qual = $db{CHROMOSOME_I};

You may select a subquality by appending the comma-separated range to the
quality score ID in the format "$id:$start,$stop". For example, here is the
first 1000 quality values of the quality score with ID "CHROMOSOME_I":

  my $qual = $db{'CHROMOSOME_I:1,1000'};

(The regular expression used to parse this format allows quality score IDs to
contain colons.)

When selecting subqualities, if $start E<gt> stop, then the reverse complement
will be returned.

The keys() and values() functions will return the IDs and their quality scores,
respectively. In addition, each() can be used to iterate over the entire data
set:

 while (my ($id,$quality) = each %db) {
    print "$id => $quality\n";
 }

When dealing with very large quality scores, you can avoid bringing them into
memory by calling each() in a scalar context. This returns the key only. You can
then use tied(%db) to recover the Bio::DB::Qual object and call its methods.

 while (my $id = each %db) {
    print "$id => $db{$quality:1,100}\n";
    print "$id => ",tied(%db)->length($id),"\n";
 }

You may, in addition invoke Bio::DB::Qual the FIRSTKEY and NEXTKEY tied hash
methods directly.

=over 2

=item $id = $db-E<gt>FIRSTKEY

Return the first ID in the database.

=item $id = $db-E<gt>NEXTKEY($id)

Given an ID, return the next quality score ID.

=back

This allows you to write the following iterative loop using just the object-
oriented interface:

 my $db = Bio::DB::Qual->new('/path/to/qual/files');
 for (my $id=$db->FIRSTKEY; $id; $id=$db->NEXTKEY($id)) {
    # do something with quality
 }

=head2 Creating a Tied Filehandle

The Bio::DB::Qual-E<gt>newFh() method creates a tied filehandle from which you
can read Bio::Seq::PrimaryQual quality score objects sequentially. The following
bit of code will iterate sequentially over all quality scores in the database:

 my $fh = Bio::DB::Qual->newFh('/path/to/qual/files');
 while (my $qual = <$fh>) {
   print $qual->id,' => ',$qual->length,"\n";
 }

When no more quality scores remain to be retrieved, the stream will return
undef.

=head1 LIMITATIONS

When a quality score is deleted from one of the qual files, this deletion is not
detected by the module and removed from the index. As a result, a "ghost" entry
will remain in the index and will return garbage results if accessed. Currently,
the only way to accommodate deletions is to rebuild the entire index, either by
deleting it manually, or by passing -reindex=E<gt>1 to new() when
initializing the module.

All quality score lines for a given quality score must have the same length
except for the last (not sure why there is this limitation). This is not
problematic for sequences but could be annoying for quality scores. A workaround
is to make sure that your quality scores fit on no more than 2 lines. Another
solution could be to padd them with blank spaces so that each line has the same
number of characters (maybe this padding should be implemented in
Bio::SeqIO::qual?).

=head1 AUTHOR

Florent E Angly E<lt>florent . angly @ gmail-dot-comE<gt>.  

Module largely based on and adapted from Bio::DB::Fasta by Lincoln Stein.

Copyright (c) 2007 Florent E Angly.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut


package Bio::DB::Qual;

use strict;
use IO::File;
use File::Spec;
use File::Basename qw(basename dirname);

use base qw(Bio::DB::IndexedBase Bio::DB::SeqI);

*qual = *quality = \&subqual;
*get_seq_by_primary_id = *get_Seq_by_acc = *get_Seq_by_id = *get_qual_by_primary_id = *get_qual_by_acc  = \&get_Qual_by_id;

use constant NA        => 0;
use constant DNA       => 1;
use constant RNA       => 2;
use constant PROTEIN   => 3;

use constant DIE_ON_MISSMATCHED_LINES => 1; # you can avoid dying if you want
                                            # but you're likely to get bad
                                            # results

my $termination_length;


=head2 new

 Title   : new
 Usage   : my $db = Bio::DB::Qual->new( $path, @options);
 Function: initialize a new Bio::DB::Qual object
 Returns : new Bio::DB::Qual object
 Args    : a single file, or path to dir, or arrayref of files

These are optional arguments to pass in as well (and their defaults).

 -glob         Glob expression to use    *.{qual,QUAL,qa,QA}
               for searching for qual
               files in directories. 

 -makeid       A code subroutine for     none
               transforming qual IDs.

 -maxopen      Maximum size of           32
               filehandle cache.

 -debug        Turn on status            0
               messages.

 -reindex      Force the index to be     0
               rebuilt.

 -dbmargs      Additional arguments      none
               to pass to the DBM
               routines when tied
               (scalar or array ref).

 -index_name   Name of the file that     Auto
               holds the indexing
               information.

=cut

sub new {
  my ($class, $path, %opts) = @_;
  $opts{-glob} ||= '*.{qual,QUAL,qa,QA}',
  my $self = Bio::DB::IndexedBase->new( $path, %opts );
  bless $self, __PACKAGE__;
  return $self;
}


=head2 _calculate_offsets

 Title   : _calculate_offsets
 Usage   : $db->_calculate_offsets($filename,$offsets);
 Function: calculates the quality score offsets in a file based on ID
 Returns : offset hash for this file
 Args    : file to process, $offsets - hashref of id to offset storage

=cut

sub _calculate_offsets {
    my ($self, $file, $offsets) = @_;
    my $base = $self->_path2fileno(basename($file));

    my $fh = IO::File->new($file) or $self->throw("Can't open $file: $!");
    binmode $fh;
    warn "Indexing $file\n" if $self->{debug};
    my ( $offset, $id, $linelength, $firstline, $count, $qual_lines, $last_line,
         %offsets );
    my ( $l3_len, $l2_len, $l_len ) = ( 0, 0, 0 );

    while (<$fh>) { # don't try this at home
        # account for crlf-terminated Windows files      
        $termination_length ||= /\r\n$/ ? 2 : 1;
        if (/^>(\S+)/) {
            print STDERR "indexed $count quality scores...\n" 
            if $self->{debug} && (++$count%1000) == 0;
            my $pos = tell($fh);
            if ($id) {
                my $qualstrlength = $pos - $offset - length($_);
                $qualstrlength -= $termination_length * $qual_lines;
                $offsets->{$id} = &{$self->{packmeth}}(
                    $offset,
                    $qualstrlength,
                    $linelength,
                    $firstline,
                    NA,
                    $base,
                );
            }
            $id = ref($self->{makeid}) eq 'CODE' ? $self->{makeid}->($_) : $1;
            ($offset, $firstline, $linelength) = ($pos, length($_), 0);
            $self->_check_linelength($linelength);
            ($l3_len, $l2_len, $l_len) = (0, 0, 0);
            $qual_lines = 0;
        } else {
            $l3_len = $l2_len;
            $l2_len = $l_len;
            $l_len = length($_);
            # need to check every line :(
            if (DIE_ON_MISSMATCHED_LINES &&
                $l3_len > 0 &&
                $l2_len > 0 &&
                $l3_len != $l2_len
            ) {
                my $fap = substr($_, 0, 20)."..";
                $self->throw("Each line of the qual entry must be the same ".
                "length except the last. Line above #$. '$fap' is $l2_len != ".
                "$l3_len chars.");
            }
            $linelength ||= length($_);
            $qual_lines++;
        }
        $last_line = $_;
    }

    $self->_check_linelength($linelength);
    # deal with last entry
    if ($id) {
        my $pos = tell($fh);
        my $qualstrlength = $pos - $offset;
      
        if ($linelength == 0) {
            $qualstrlength = 0;
        } else {
            if ($last_line !~ /\s$/) {
                $qual_lines--;
            }
            $qualstrlength -= $termination_length * $qual_lines;
        }
        $offsets->{$id} = &{$self->{packmeth}}(
            $offset,
            $qualstrlength,
            $linelength,
            $firstline,
            NA,
            $base,
        );
    }
    return \%offsets;
}


=head2 get_Qual_by_id

 Title   : get_Qual_by_id
 Usage   : my $qual = $db->get_Qual_by_id($id)
 Function: Bio::DB::RandomAccessI method implemented
 Returns : Bio::PrimarySeqI object
 Args    : id

=cut

sub get_Qual_by_id {
    my ($self, $id) = @_;
    $self->throw('Need to provide a sequence ID') if not defined $id;
    return if not exists $self->{offsets}{$id};
    return Bio::Seq::PrimaryQual::Qual->new($self,$id);
}


sub offset {
    my ($self, $id) = @_;
    $self->throw('Need to provide a sequence ID') if not defined $id;
    my $offset = $self->{offsets}{$id} or return;
    return (&{$self->{unpackmeth}}($offset))[0];
}


=head2 length

 Title   : length
 Usage   : $qualdb->length($seqid);
 Function: gets the number of quality values in a quality score
 Returns : scalar
 Args    : ID of a quality score

=cut

sub length {
    # the NUMBER of quality values
    my ($self, $id) = @_;
    $self->throw('Need to provide a sequence ID') if not defined $id;
    my $len = scalar(@{$self->subqual($id)});
    return $len;
}


sub lengthstr {
    # the length of the quality STRING
    my ($self, $id) = @_;
    $self->throw('Need to provide a sequence ID') if not defined $id;
    my $offset = $self->{offsets}{$id} or return;
    return (&{$self->{unpackmeth}}($offset))[1];
}


sub linelen {
    my ($self, $id) = @_;
    $self->throw('Need to provide a sequence ID') if not defined $id;
    my $offset = $self->{offsets}{$id} or return;
    return (&{$self->{unpackmeth}}($offset))[2];
}


sub headerlen {
    my ($self, $id) = @_;
    $self->throw('Need to provide a sequence ID') if not defined $id;
    my $offset = $self->{offsets}{$id} or return;
    return (&{$self->{unpackmeth}}($offset))[3];
}


sub header_offset {
    my ($self, $id) = @_;
    $self->throw('Need to provide a sequence ID') if not defined $id;
    return if not $self->{offsets}{$id};
    return $self->offset($id) - $self->headerlen($id);
}


sub file {
    my ($self, $id) = @_;
    $self->throw('Need to provide a sequence ID') if not defined $id;
    my $offset = $self->{offsets}{$id} or return;
    return $self->_fileno2path((&{$self->{unpackmeth}}($offset))[5]);
}


=head2 subqual

 Title   : subqual
 Usage   : my @qualarr = @{$qualdb->subqual($id,$start,$stop)};
 Function: returns a subqual of a quality score in the database
 Returns : subquality array reference
 Args    : id of quality score, starting quality value number, ending quality
           value number

=cut

sub subqual {
    my ($self, $id, $start, $stop) = @_;

    # Quality values in a quality score can have 1 or 2 digits and are separated
    # by one (or several?) spaces. Thus contrary to Bio::DB::Fasta, here there
    # is no easy way match the position of a quality value to its position in
    # the quality string.
    # As a consequence, if a subqual of the quality is requested, we still need
    # to grab the full quality string first - performance penalty for big
    # quality scores :(
    # I think there is no way around starting at the begining of the quality
    # score but maybe there is a resource-efficient way of starting at the
    # begining of the quality score and stopping when the the position of the
    # last quality value requested is reached??

    $self->throw('Need to provide a sequence ID') if not defined $id;

    # position of the quality values
    if ($id =~ /^(.+):([\d_]+)(?:,|-|\.\.)([\d_]+)$/) {
        ($id, $start, $stop) = ($1,$2,$3);
        $start =~ s/_//g;
        $stop  =~ s/_//g;
    }

    # position in quality string
    my $string_start = 1;
    my $string_stop = $self->lengthstr($id);

    # fetch full quality string
    my $fh = $self->fh($id) or return;
    my $filestart = $self->_caloffset($id, $string_start, $termination_length);
    my $filestop  = $self->_caloffset($id, $string_stop , $termination_length);
    seek($fh,$filestart,0);
    my $data;
    read($fh, $data, $filestop-$filestart+1);

    # process quality score
    $data =~ s/\n//g;
    $data =~ s/\r//g;
    my $reverse = 0;
    if ($stop && $start && $stop < $start) {
        $reverse = 1;
        my $tmp = $start;
        $start = $stop;
        $stop = $tmp; 
    }
    my $subqual = 0;
    $subqual = 1 if ( $start || $stop );
    my @data;
    if ( $subqual || $reverse ) {
        @data = split / /, $data, $stop+1;
        my $length = scalar(@data);
        $start = 1       if $start < 1;
        $stop  = $length if $stop  > $length;
        pop @data if ($stop != $length);
        splice @data, 0, $start-1;
        @data = reverse(@data) if $reverse;
        $data = join ' ', @data;
    } else {
        @data = split / /, $data;
    }

    return \@data;
}

sub fh {
    my ($self, $id) = @_;
    $self->throw('Need to provide a sequence ID') if not defined $id;
    my $file = $self->file($id) or return;
    return $self->_fhcache( File::Spec->catfile($self->{dirname},$file) ) or
        $self->throw( "Can't open file $file");
}

=head2 header

 Title   : header
 Usage   : $qualdb->header($id);
 Function: returns the header of a quality score in the database
 Returns : header string
 Args    : id of quality score

=cut

sub header {
    my ($self, $id) = @_;
    $self->throw('Need to provide a sequence ID') if not defined $id;
    my ($offset,$seqlength,$linelength,$firstline,$file) 
        = &{$self->{unpackmeth}}($self->{offsets}{$id}) or return;
    $offset -= $firstline;
    my $data;
    my $fh = $self->fh($id) or return;
    seek($fh,$offset,0);
    read($fh,$data,$firstline);
    chomp $data;
    substr($data,0,1) = '';
    return $data;
}


=head2 get_PrimaryQual_stream

 Title   : get_PrimaryQual_stream
 Usage   : $qualdb->get_PrimaryQual_stream
 Function: get a SeqIO-like stream of quality scores 
 Returns : stream object
 Args    : none

=cut

sub get_PrimaryQual_stream {
    my $self = shift;
    return Bio::DB::Indexed::Stream->new($self);
}


#-------------------------------------------------------------
# Tied hash overrides
#

sub FETCH {
    return shift->subqual(@_);
}


#-------------------------------------------------------------
# Bio::Seq::PrimaryQual compatibility
#
# Usage is the same as in Bio::Seq::PrimaryQual

package Bio::Seq::PrimaryQual::Qual;
use overload '""' => 'display_id';

use base qw(Bio::Root::Root Bio::Seq::PrimaryQual);

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($db, $id, $start, $stop) = $self->_rearrange(
                                 [qw(DATABASE ID START STOP)],
                                 @args);
    $self->{db}    = $db;
    $self->{id}    = $id;
    $self->{start} = $start || 1;
    $self->{stop}  = $stop  || $db->length($id);
    return $self;
}

sub qual {
    my $self = shift;
    my $qual = $self->{db}->qual($self->{id},$self->{start},$self->{stop});
    return $qual;
}

sub subqual {
    my ($self, $start, $stop) = @_;
    return $self->trunc($start, $stop)->qual;
}

sub trunc {
    my ($self, $start, $stop) = @_;
    $self->throw(
        "$stop is smaller than $stop. If you want to truncate and reverse ".
        "complement, you must call trunc followed by revcom."
    ) if $start > $stop;
    my ($left, $right);
    if ($self->{start} <= $self->{stop}) {
        $left = $self->{start}+$start-1;
        $right  = $self->{start}+$stop-1;
    } else {
        $left = $self->{start}-($start-1);
        $right  = $self->{start}-($stop-1);
    }
    my $obj = $self->new( -database => $self->{db},
                          -id       => $self->{id},
                          -start    => $left,
                          -stop     => $right
                        );
    return $obj;  
}

sub display_id {
    my $self = shift;
    return $self->{id};
}

sub primary_id {
    my $self = shift;
    return overload::StrVal($self);
}

sub length {
    # number of quality scores
    my $self = shift;
    return scalar(@{$self->qual});
}

sub description  { 
    my $self = shift;
    my $header = $self->{'db'}->header($self->{id});
    # remove the id from the header
    $header = (split(/\s+/,$header,2))[2];
    return $header;
}
*desc = \&description;


1;
