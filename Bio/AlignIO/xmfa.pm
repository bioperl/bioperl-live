#
# BioPerl module for Bio::AlignIO::xmfa
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::xmfa - XMFA MSA Sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::AlignIO> 
class.

=head1 DESCRIPTION

This object can transform L<Bio::SimpleAlign> objects from
XMFA flat file databases.  For more information, see:

  http://asap.ahabs.wisc.edu/mauve-aligner/mauve-user-guide/mauve-output-file-formats.html

This module is based on the AlignIO::fasta parser written by
Peter Schattner

=head1 TODO

Finish write_aln(), clean up code, allow LargeLocatableSeq (ie for
very large sequences a'la Mauve)

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

Chris Fields

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::xmfa;
use strict;

use base qw(Bio::AlignIO);
our $WIDTH = 60;

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
    my ($width) = $self->_rearrange([qw(WIDTH)],@_);
    $self->width($width || $WIDTH);

    my ($name, $tempname, $seqchar);
    my $aln = Bio::SimpleAlign->new();
    my $seqs = 0;
    # alignments
    while (defined (my $entry = $self->_readline) ) {
        chomp $entry;
        if ( index($entry, '=') == 0 ) {
            if (defined $name && $seqchar) {
                my $seq = $self->_process_seq($name, $seqchar);
                $aln->add_seq($seq);
            }
            if ($aln && $entry =~ m{score\s*=\s*(\d+)}) {
                $aln->score($1);
            }
            $seqchar = '';
            undef $name;
            last;
        } elsif ( $entry =~ m{^>.+$}xms) {
            if ( defined $name ) {
                my $seq = $self->_process_seq($name, $seqchar);
                $aln->add_seq($seq);
            }
            $seqchar = '';
            $name = $entry;
        } else {
            $seqchar .= $entry;
        }
    }
    
    # this catches last sequence if '=' is not present (Mauve)
    if ( defined $name ) {
        my $seq = $self->_process_seq($name, $seqchar);
        $aln->add_seq($seq);
    }
    $aln->num_sequences ? return $aln : return;
}

=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in xmfa format
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object

See L<Bio::Align::AlignI>

=cut

sub write_aln {
    my ($self,@aln) = @_;
    my $width = $self->width;
    my ($seq,$desc,$rseq,$name,$count,$length,$seqsub,$start,$end,$strand,$id);

    foreach my $aln (@aln) {
        if( ! $aln || ! $aln->isa('Bio::Align::AlignI')  ) { 
            $self->warn("Must provide a Bio::Align::AlignI object when calling write_aln");
            next;
        }
        #if( $self->force_displayname_flat ) {
        #    $aln->set_displayname_flat(1);
        #}
        my $seqct = 1;
        foreach $rseq ( $aln->each_seq() ) {
            ($start, $end, $strand, $id) = ($rseq->start, $rseq->end, $rseq->strand || 0,
                                            $rseq->display_id);
            $strand = ($strand == 1)  ? '+' :
                      ($strand == -1) ? '-' :
                      '';
            $name = sprintf("%d:%d-%d %s %s",$seqct,$start,$end,$strand,$id);
            $seq  = $rseq->seq();
            $desc = $rseq->description || '';
            $self->_print (">$name $desc\n") or return ;	
            $count = 0;
            $length = length($seq);
            if(defined $seq && $length > 0) {
            $seq =~ s/(.{1,$width})/$1\n/g;
            } else {
            $seq = "\n";
            }
            $self->_print($seq) || return 0;
            $seqct++;
        }
        my $alndesc = '';
        $alndesc = "score = ".$aln->score if ($aln->score);
        $self->_print("= $alndesc\n") || return 0;
        
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
    $seq =~ s/[^A-Z]//gi;
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

####### PRIVATE #######

sub _process_seq {
    my ($self, $entry, $seq) = @_;
    my ($start, $end, $strand, $seqname, $desc, $all);
    # put away last name and sequence
    if ( $entry =~ m{^>\s*\d+:(\d+)-(\d+)\s([+-]{1})(?:\s+(\S+)\s*(\S\.*)?)?} ) {
        ($start, $end, $seqname, $desc) = ($1, $2, $4, $5);
        $strand = ($3 eq '+')  ?  1  : -1;
    } else {
        $self->throw("Line does not comform to XMFA format:\n$entry");
    }
    my $seqobj = Bio::LocatableSeq->new(
					-nowarnonempty => 1,
					-strand      => $strand,
					-seq         => $seq,
					-display_id  => $seqname,
					-description => $desc || $all,
					-start       => $start,
					-end         => $end,
					-alphabet    => $self->alphabet,
					);
    $self->debug("Reading $seqname\n");
    return $seqobj;
}

1;
