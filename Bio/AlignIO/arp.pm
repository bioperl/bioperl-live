# $Id$
#
# BioPerl module for Bio::AlignIO::arp
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::arp - ARP MSA Sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::AlignIO> 
class.

=head1 DESCRIPTION

This object can create L<Bio::SimpleAlign> objects from
ARP flat files.  These are typically configuration-like data files
for the program Arlequin.  For more information, see:

  http://lgb.unige.ch/arlequin/

Note that, at the moment, this only scans the allele sequence data
in the DATA section and inserts them into SimpleAlign objects.  ARP
files that contain other data (RFLP, etc.) will not be parsed; if the
DNA data is actually SNP data, then the LocatableSeq object instantiation
will throw an error.  

=head1 FEEDBACK

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHORS

Chris Fields (cjfields)

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::arp;
use strict;
use Data::Dumper;

use base qw(Bio::AlignIO);

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
    my $aln = Bio::SimpleAlign->new();
    
    my $data;
    my $sdflag = 0;
    my %arp;
    
    SCAN:
    while (defined ($data = $self->_readline) ) {
        if ($data =~ m{DataType=(\S+)}) {
            $self->datatype($1);
        } elsif ($data =~ m{SampleData=}) {
            $sdflag = 1;
            if (lc($self->datatype) ne 'dna') {
                $self->warn("ARP data does not contain sequence!");
                return;
            }
            next;
        } elsif ($data =~ m{SampleName="(.*)"}) {
            $aln->description($1);
        }
        DATA:
        while ($sdflag) {
            $data =~ s{(?:^\s+|\s+$)}{};
            my ($id, $score, $seq) = split m{\s+}, $data,2;
            # what to do with the score???
            my $temp;
            ($temp = $data) =~ s{[^A-Z]}{}gi;
            $self->debug($temp,"\n");
            my $newseq = new Bio::LocatableSeq(
                         -start       => 1,
                         -end         => CORE::length($temp),
                         -seq         => $seq,
                         -id          => $id,
                         );
            $aln->add_seq($newseq);
            $self->debug("Reading $id\n");
            $data = $self->_readline;
            if ($data =~ /^\s*}\s*$/) {
                last SCAN;
            }
        }
    }
    # alignments only returned if they contain sequences
    return $aln if $aln->no_sequences;
    return 0;
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
    $self->throw_not_implemented;
}

sub datatype {
    my $self = shift;
    return $self->{'_datatype'} = shift if @_;
    return $self->{'_datatype'};
}

1;
