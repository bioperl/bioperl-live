#
# BioPerl module for Bio::AlignIO::fasta
#
# Copyright Peter Schattner
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::fasta - fasta MSA Sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::AlignIO>
class.

=head1 DESCRIPTION

This object can transform L<Bio::SimpleAlign> objects to and from
fasta flat files.  This is for the fasta alignment format, not
for the FastA sequence analysis program.  To process the alignments from
FastA (FastX, FastN, FastP, tFastA, etc) use the Bio::SearchIO module.

=head1 FEEDBACK

=head2 Support

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and
reponsive experts will be able look at the problem and quickly
address it. Please include a thorough description of the problem
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHORS

Peter Schattner

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::fasta;
use strict;

use base qw(Bio::AlignIO);
our $WIDTH = 60;
use Bio::LocatableSeq;

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln
 Function: returns the next alignment in the stream.
 Returns : Bio::Align::AlignI object - returns 0 on end of file
            or on error
 Args    : -width => optional argument to specify the width sequence
           will be written (60 chars by default)

See L<Bio::Align::AlignI>

=cut

sub next_aln {
    my $self = shift;
    my ($width) = $self->_rearrange( [qw(WIDTH)], @_ );
    $self->width( $width || $WIDTH );

    my ($start, $end,      $name,     $seqname, $seq,  $seqchar,
        $entry, $tempname, $tempdesc, %align,   $desc, $maxlen
    );
    my $aln = Bio::SimpleAlign->new();

    while ( defined( $entry = $self->_readline ) ) {
        chomp $entry;
        if ( $entry =~ s/^>\s*(\S+)\s*// ) {
            $tempname = $1;
            chomp($entry);
            $tempdesc = $entry;
            if ( defined $name ) {
                $seqchar =~ s/\s//g;
                $seqname = $name;
                $start   = 1;
                $end     = $self->_get_len($seqchar);
                $seq     = Bio::LocatableSeq->new(
                    -seq         => $seqchar,
                    -display_id  => $seqname,
                    -description => $desc,
                    -start       => $start,
                    -end         => $end,
                    -alphabet    => $self->alphabet,
                );
                $aln->add_seq($seq);
                $self->debug("Reading $seqname\n");
            }
            $desc    = $tempdesc;
            $name    = $tempname;
            $desc    = $entry;
            $seqchar = "";
            next;
        }

        # removed redundant symbol validation
        # this is already done in Bio::PrimarySeq
        $seqchar .= $entry;
    }

    #  Next two lines are to silence warnings that
    #  otherwise occur at EOF when using <$fh>
    $name    = "" if ( !defined $name );
    $seqchar = "" if ( !defined $seqchar );
    $seqchar =~ s/\s//g;

    #  Put away last name and sequence
    if ( $name =~ /(\S+)\/(\d+)-(\d+)$/ ) {
        $seqname = $1;
        $start   = $2;
        $end     = $3;
    }
    else {
        $seqname = $name;
        $start   = 1;
        $end     = $self->_get_len($seqchar);
    }

    # This logic now also reads empty lines at the
    # end of the file. Skip this is seqchar and seqname is null
    unless ( length($seqchar) == 0 && length($seqname) == 0 ) {
        $seq = Bio::LocatableSeq->new(
            -seq         => $seqchar,
            -display_id  => $seqname,
            -description => $desc,
            -start       => $start,
            -end         => $end,
            -alphabet    => $self->alphabet,
        );
        $aln->add_seq($seq);
        $self->debug("Reading $seqname\n");
    }
    my $alnlen = $aln->length;
    foreach my $seq ( $aln->each_seq ) {
        if ( $seq->length < $alnlen ) {
            my ($diff) = ( $alnlen - $seq->length );
            $seq->seq( $seq->seq() . "-" x $diff );
        }
    }

    # no sequences means empty alignment (possible EOF)
    return $aln if $aln->num_sequences;
    return;
}


=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in fasta format
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object

See L<Bio::Align::AlignI>

=cut

sub write_aln {
    my ($self,@aln) = @_;
    my $width = $self->width;
    my ($seq,$desc,$rseq,$name,$count,$length,$seqsub);

    foreach my $aln (@aln) {
        if( ! $aln || ! $aln->isa('Bio::Align::AlignI')  ) {
            $self->warn("Must provide a Bio::Align::AlignI object when calling write_aln");
            next;
        }
        if( $self->force_displayname_flat ) {
            $aln->set_displayname_flat(1);
        }
        foreach $rseq ( $aln->each_seq() ) {
            $name = $aln->displayname($rseq->get_nse());
            $seq  = $rseq->seq();
            $desc = $rseq->description || '';
            $desc = ' '.$desc if $desc;
            $self->_print (">$name$desc\n") or return;
            $count = 0;
            $length = length($seq);
            if(defined $seq && $length > 0) {
                $seq =~ s/(.{1,$width})/$1\n/g;
            } else {
                $seq = "\n";
            }
            $self->_print($seq);
        }
    }
    $self->flush if $self->_flush_on_write && defined $self->_fh;
    return 1;
}

=head2 _get_len

 Title   : _get_len
 Usage   :
 Function: determine number of alphabetic chars
 Returns : integer
 Args    : sequence string

=cut

sub _get_len {
    my ($self,$seq) = @_;
    my $chars = $Bio::LocatableSeq::GAP_SYMBOLS.$Bio::LocatableSeq::FRAMESHIFT_SYMBOLS;
    $seq =~ s{[$chars]+}{}gi;
    return CORE::length($seq);
}

=head2 width

 Title   : width
 Usage   : $obj->width($newwidth)
           $width = $obj->width;
 Function: Get/set width of alignment
 Returns : integer value of width
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub width{
    my $self = shift;

    return $self->{'_width'} = shift if @_;
    return $self->{'_width'} || $WIDTH;
}

1;
