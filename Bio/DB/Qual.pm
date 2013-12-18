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
  my $db      = Bio::DB::Qual->new('/path/to/qual/files/');
  my @ids     = $db->get_all_primary_ids;

  # Simple access
  my @qualarr = @{$db->qual('CHROMOSOME_I',4_000_000 => 4_100_000)};
  my @revqual = @{$db->qual('CHROMOSOME_I',4_100_000 => 4_000_000)};
  my $length  = $db->length('CHROMOSOME_I');
  my $header  = $db->header('CHROMOSOME_I');

  # Access to sequence objects. See Bio::PrimarySeqI.
  my $obj     = $db->get_Qual_by_id('CHROMOSOME_I');
  my @qual    = @{$obj->qual};
  my @subqual = @{$obj->subqual(4_000_000 => 4_100_000)};
  my $length  = $obj->length;

  # Loop through sequence objects
  my $stream  = $db->get_PrimarySeq_stream;
  while (my $qual = $stream->next_seq) {
    # Bio::Seq::PrimaryQual operations
  }

  # Filehandle access
  my $fh = Bio::DB::Qual->newFh('/path/to/qual/files/');
  while (my $qual = <$fh>) {
    # Bio::Seq::PrimaryQual operations
  }

  # Tied hash access
  tie %qualities,'Bio::DB::Qual','/path/to/qual/files/';
  print $qualities{'CHROMOSOME_I:1,20000'};

=head1 DESCRIPTION

Bio::DB::Qual provides indexed access to a single Fasta file, several files,
or a directory of files. It provides random access to each quality score entry
without having to read the file from the beginning. Access to subqualities
(portions of a quality score) is provided, although contrary to Bio::DB::Fasta,
the full quality score has to be brought in memory. Bio::DB::Qual is based on
Bio::DB::IndexedBase. See this module's documentation for details.

The qual files should contain decimal quality scores. Entries may have any line
length up to 65,536 characters, and different line lengths are allowed in the
same file. However, within a quality score entry, all lines must be the same
length except for the last. An error will be thrown if this is not the case.

The module uses /^E<gt>(\S+)/ to extract the primary ID of each quality score
from the qual header. See -makeid in Bio::DB::IndexedBase to pass a callback
routine to reversibly modify this primary ID, e.g. if you wish to extract a
specific portion of the gi|gb|abc|xyz GenBank IDs.

=head1 DATABASE CREATION AND INDEXING

The object-oriented constructor is new(), the filehandle constructor is newFh()
and the tied hash constructor is tie(). They all allow to index a single Fasta
file, several files, or a directory of files. See Bio::DB::IndexedBase.

=head1 SEE ALSO

L<Bio::DB::IndexedBase>

L<Bio::DB::Fasta>

L<Bio::Seq::PrimaryQual>

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

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

For BioPerl-style access, the following methods are provided:

=head2 get_Seq_by_id

 Title   : get_Seq_by_id,  get_Seq_by_acc, get_Seq_by_version, get_Seq_by_primary_id,
           get_Qual_by_id, get_qual_by_acc, get_qual_by_version, get_qual_by_primary_id,
 Usage   : my $seq = $db->get_Seq_by_id($id);
 Function: Given an ID, fetch the corresponding sequence from the database.
 Returns : A Bio::PrimarySeq::Fasta object (Bio::PrimarySeqI compliant)
           Note that to save resource, Bio::PrimarySeq::Fasta sequence objects
           only load the sequence string into memory when requested using seq().
           See L<Bio::PrimarySeqI> for methods provided by the sequence objects
           returned from get_Seq_by_id() and get_PrimarySeq_stream().
 Args    : ID

=head2 get_PrimarySeq_stream

 Title   : get_Seq_stream, get_PrimarySeq_stream
 Usage   : my $stream = $db->get_Seq_stream();
 Function: Get a stream of Bio::PrimarySeq::Fasta objects. The stream supports a
           single method, next_seq(). Each call to next_seq() returns a new
           Bio::PrimarySeq::Fasta sequence object, until no more sequences remain.
 Returns : A Bio::DB::Indexed::Stream object
 Args    : None

=head1

For simple access, the following methods are provided:

=cut


package Bio::DB::Qual;

use strict;
use IO::File;
use File::Spec;

use base qw(Bio::DB::IndexedBase);

our $obj_class = 'Bio::Seq::PrimaryQual::Qual';
our $file_glob = '*.{qual,QUAL,qa,QA}';


=head2 new

 Title   : new
 Usage   : my $db = Bio::DB::Qual->new( $path, %options);
 Function: Initialize a new database object. When indexing a directory, files
           ending in .qual,qa are indexed by default.
 Returns : A new Bio::DB::Qual object
 Args    : A single file, or path to dir, or arrayref of files
           Optional arguments: see Bio::DB::IndexedBase

=cut


sub _calculate_offsets {
    # Bio::DB::IndexedBase calls this to calculate offsets
    my ($self, $fileno, $file, $offsets) = @_;

    my $fh = IO::File->new($file) or $self->throw("Could not open $file: $!");
    binmode $fh;
    warn "Indexing $file\n" if $self->{debug};
    my ($offset, @ids, $linelen, $headerlen, $count, $qual_lines, $last_line,
        $numres, %offsets);
    my ($l3_len, $l2_len, $l_len, $blank_lines) = (0, 0, 0, 0);

    my $termination_length = $self->{termination_length};
    while (my $line = <$fh>) {
        # Account for crlf-terminated Windows files      
        if (index($line, '>') == 0) {
            if ($line =~ /^>(\S+)/) {
                print STDERR "Indexed $count quality scores...\n" 
                    if $self->{debug} && (++$count%1000) == 0;
                $self->_check_linelength($linelen);
                my $pos = tell($fh);
                if (@ids) {
                    my $strlen = $pos - $offset - length($line);
                    $strlen   -= $termination_length * $qual_lines;
                    my $ppos = &{$self->{packmeth}}($offset, $strlen, $numres,
                        $linelen, $headerlen, Bio::DB::IndexedBase::NA, $fileno);
                    for my $id (@ids) {
                        $offsets->{$id} = $ppos;
                    }
                    $numres = 0;
                }
                @ids = $self->_makeid($line);
                ($offset, $headerlen, $linelen, $qual_lines) = ($pos, length $line, 0, 0);
                ($l3_len, $l2_len, $l_len, $blank_lines) = (0, 0, 0, 0);
            } else {
                # Catch bad header lines, bug 3172
                $self->throw("FASTA header doesn't match '>(\\S+)': $line");
            }
        } elsif ($line !~ /\S/) {
            # Skip blank line
            $blank_lines++;
            next;
        } else {
            # Need to check every line :(
            $l3_len = $l2_len;
            $l2_len = $l_len;
            $l_len  = length $line;
            if (Bio::DB::IndexedBase::DIE_ON_MISSMATCHED_LINES) {
                if ( ($l3_len > 0) && ($l2_len > 0) && ($l3_len != $l2_len) ) {
                    my $fap = substr($line, 0, 20)."..";
                    $self->throw("Each line of the qual entry must be the same ".
                        "length except the last. Line above #$. '$fap' is $l2_len".
                        " != $l3_len chars.");
                }
                if ($blank_lines) {
                    # Blank lines not allowed in entry
                    $self->throw("Blank lines can only precede header lines, ".
                        "found preceding line #$.");
                }
            }
            $linelen ||= length $line;
            $qual_lines++;
            $numres += scalar(split /\s+/, $line);
        }
        $last_line = $line;
    }

    # Process last entry
    $self->_check_linelength($linelen);
    my $pos = tell($fh);
    if (@ids) {
        my $strlen = $pos - $offset;
        if ($linelen == 0) {
            $strlen = 0;
        } else {
            if ($last_line !~ /\s$/) {
                $qual_lines--;
            }
            $strlen -= $termination_length * $qual_lines;
        }
        my $ppos = &{$self->{packmeth}}($offset, $strlen, $numres, $linelen,
            $headerlen, Bio::DB::IndexedBase::NA, $fileno);
        for my $id (@ids) {
            $offsets->{$id} = $ppos;
        }
    }

    return \%offsets;
}


# for backward compatibility
sub get_PrimaryQual_stream {
   my $self = shift;
   return $self->get_PrimarySeq_stream;
}


# for backward compatibility
sub get_Qual_by_id {
    my ($self, $id) = @_;
    return $self->get_Seq_by_id($id);
}

*get_qual_by_version = *get_qual_by_primary_id = *get_qual_by_acc = \&get_Qual_by_id;


=head2 qual

 Title   : qual, quality, subqual
 Usage   : # All quality scores
           my @qualarr = @{$qualdb->subqual($id)};
           # Subset of the quality scores
           my @subqualarr = @{$qualdb->subqual($id, $start, $stop, $strand)};
           # or...
           my @subqualarr = @{$qualdb->subqual($compound_id)};
 Function: Get a subqual of an entry in the database. For your convenience,
           the sequence to extract can be specified with any of the following
           compound IDs:
              $db->qual("$id:$start,$stop")
              $db->qual("$id:$start..$stop")
              $db->qual("$id:$start-$stop")
              $db->qual("$id:$start,$stop/$strand")
              $db->qual("$id:$start..$stop/$strand")
              $db->qual("$id:$start-$stop/$strand")
              $db->qual("$id/$strand")
           If $stop is less than $start, then the reverse complement of the
           sequence is returned. Avoid using it if possible since this goes
           against Bio::Seq conventions.
 Returns : Reference to an array of quality scores
 Args    : Compound ID of entry to retrieve
             or
           ID, optional start (defaults to 1), optional end (defaults to the
           number of quality scores for this sequence), and strand (defaults to
           1).

=cut

sub subqual {
    my ($self, $id, $start, $stop, $strand) = @_;

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
    ($id, $start, $stop, $strand) = $self->_parse_compound_id($id, $start, $stop, $strand);

    # Position in quality string
    my $string_start = 1;
    my $string_stop = $self->strlen($id);

    # Fetch full quality string
    my $fh = $self->_fh($id) or return;
    my $filestart = $self->_calc_offset($id, $string_start);
    my $filestop  = $self->_calc_offset($id, $string_stop );
    seek($fh, $filestart,0);
    my $data;
    read($fh, $data, $filestop-$filestart+1);

    # Process quality score
    $data =~ s/\n//g;
    $data =~ s/\r//g;
    my $subqual = 0;
    $subqual = 1 if ( $start || $stop );
    my @data;
    if ( $subqual || ($strand == -1) ) {
        @data = split / /, $data, $stop+1;
        my $length = scalar(@data);
        $start = 1       if $start < 1;
        $stop  = $length if $stop  > $length;
        pop @data if ($stop != $length);
        splice @data, 0, $start-1;
        @data = reverse(@data) if $strand == -1;
        $data = join ' ', @data;
    } else {
        @data = split / /, $data;
    }

    return \@data;
}

*qual = *quality = \&subqual;


=head2 header

 Title   : header
 Usage   : my $header = $db->header($id);
 Function: Get the header line (ID and description fields) of the specified entry.
 Returns : String
 Args    : ID of entry

=cut

sub header {
    my ($self, $id) = @_;
    $self->throw('Need to provide a sequence ID') if not defined $id;
    my ($offset, $headerlen) = (&{$self->{unpackmeth}}($self->{offsets}{$id}))[0,4];
    $offset -= $headerlen;
    my $data;
    my $fh = $self->_fh($id) or return;
    seek($fh, $offset, 0);
    read($fh, $data, $headerlen);
    # On Windows chomp remove '\n' but leaves '\r'
    # when reading '\r\n' in binary mode
    $data =~ s/\n//g;
    $data =~ s/\r//g;
    substr($data, 0, 1) = '';
    return $data;
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
    $self->{stop}  = $stop  || $db->length($id);
    $self->{start} = $start || ($self->{stop} > 0 ? 1 : 0); # handle 0-length seqs
    return $self;
}


sub qual {
    my $self = shift;
    my $qual = $self->{db}->qual($self->{id}, $self->{start}, $self->{stop});
    return $qual;
}


sub subqual {
    my ($self, $start, $stop) = @_;
    return $self->trunc($start, $stop)->qual;
}


sub trunc {
    # Override Bio::Seq::QualI trunc() method. This way, we create an object
    # that does not store the quality array in memory.
    my ($self, $start, $stop) = @_;
    $self->throw(
        "$stop is smaller than $stop. If you want to truncate and reverse ".
        "complement, you must call trunc followed by revcom."
    ) if $start > $stop;
    if ($self->{start} <= $self->{stop}) {
        $start = $self->{start}+$start-1;
        $stop  = $self->{start}+$stop-1;
    } else {
        $start = $self->{start}-($start-1);
        $stop  = $self->{start}-($stop-1);
    }
    my $obj = $self->new( -database => $self->{db},
                          -id       => $self->{id},
                          -start    => $start,
                          -stop     => $stop
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

sub revcom {
    # Override Bio::QualI revcom() with optimized method.
    my $self = shift;
    return $self->new(@{$self}{'db', 'id', 'stop', 'start'});
}

sub length {
    # Get length from quality location, not the quality array (too expensive)
    my $self = shift;
    return $self->{start} < $self->{stop}   ?
        $self->{stop}  - $self->{start} + 1 :
        $self->{start} - $self->{stop}  + 1 ;
}


sub description  { 
    my $self = shift;
    my $header = $self->{'db'}->header($self->{id});
    # remove the id from the header
    $header = (split(/\s+/, $header, 2))[2];
    return $header;
}
*desc = \&description;


1;
