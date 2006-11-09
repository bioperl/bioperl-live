# $Id$
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

  http://gel.ahabs.wisc.edu/docserver/mauve/files.stx

This module is based on the AlignIO::fasta parser written by
Peter Schattner

=head1 TODO

Finish write_aln(), clean up code

=head1 FEEDBACK

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

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

    my ($start, $end, $strand, $name, $seqname, $seq, $seqchar, $entry, 
         $tempname, $tempdesc, %align, $desc, $maxlen, $extra);
    my $aln = Bio::SimpleAlign->new();

    # alignments
    while (defined ($entry = $self->_readline) ) {
        chomp $entry;
        if ( index($entry, '=') == 0 ) {
            if ($entry =~ m{score\s*=\s*(\d+)}) {
                $aln->score($1);
            }
            last;
        }
        if ( $entry =~ s{^>(.*)$}{} ) {
            $tempname  = $1;
            chomp($entry);
            $tempdesc = $entry;
            if ( defined $name ) {
                # put away last name and sequence
                if ( $name =~ m{\d+:(\d+)-(\d+)\s([+-]{1})\s+(\S+)\s*(.*)} ) {
                    ($start, $end, $seqname, $extra) = ($1, $2, $4, $5);
                    $strand = ($3 eq '+')  ?  1  : 
                              ($3 eq '-')  ?  -1 :
                              0; 
                } else {
                    $self->throw("Does not comform to XMFA format");
                }
                $seq = new Bio::LocatableSeq(
                         -strand      => $strand,
                         -seq         => $seqchar,
                         -display_id  => $seqname,
                         -description => $extra,
                         -start       => $start,
                         -end         => $end,
                         );
                $aln->add_seq($seq);
                $self->debug("Reading $seqname\n");
            }
            $desc = $tempdesc;  
            $name = $tempname;
            $desc = $entry;
            $seqchar  = "";
            next;
        }
        $seqchar .= $entry;
    }

    #  Next two lines are to silence warnings that
    #  otherwise occur at EOF when using <$fh>
    $name = "" if (!defined $name);
    $seqchar="" if (!defined $seqchar);

    #  Put away last name and sequence
    if ( $name =~ m{\d+:(\d+)-(\d+)\s([+-]{1})\s+(\S+)\s*(.*)} ) {
        ($start, $end, $seqname, $extra) = ($1, $2, $4, $5);
        $strand = ($3 eq '+')  ?  1  : 
                  ($3 eq '-')  ?  -1 :
                   0;
    }

    #  If $end <= 0, we have either reached the end of
    #  file in <> or we have encountered some other error
    if ( $end <= 0 ) { 
        undef $aln; 
        return $aln;
    }

    # This logic now also reads empty lines at the 
    # end of the file. Skip this is seqchar and seqname is null
    unless ( length($seqchar) == 0 && length($seqname) == 0 ) {
        $seq = new Bio::LocatableSeq(-seq         => $seqchar,
                                     -strand      => $strand,
                                      -display_id  => $seqname,
                                      -description => $extra,
                                      -start       => $start,
                                      -end         => $end,
                                             );
        $aln->add_seq($seq);
        $self->debug("Reading $seqname\n");
    }
    my $alnlen = $aln->length;
    foreach my $seq ( $aln->each_seq ) {
        if ( $seq->length < $alnlen ) {
            my ($diff) = ($alnlen - $seq->length);
            $seq->seq( $seq->seq() . "-" x $diff);
        }
    }
    return $aln;
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
            ($start, $end, $strand, $id) = ($rseq->start, $rseq->end, $rseq->strand,
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

1;
