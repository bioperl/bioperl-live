#
# BioPerl module for Bio::AlignIO::psi
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::psi - Read/Write PSI-BLAST profile alignment files

=head1 SYNOPSIS

This module will parse PSI-BLAST output of the format seqid XXXX  

=head1 DESCRIPTION

This is a parser for psi-blast blocks.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::AlignIO::psi;
use vars qw($BlockLen $IdLength);
use strict;

$BlockLen = 100; 
$IdLength = 13;

# Object preamble - inherits from Bio::Root::Root

use Bio::SimpleAlign;
use Bio::LocatableSeq;

use base qw(Bio::AlignIO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::AlignIO::psi->new();
 Function: Builds a new Bio::AlignIO::psi object 
 Returns : Bio::AlignIO::psi
 Args    :

=cut

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream
 Returns : Bio::Align::AlignI object
 Args    : NONE

See L<Bio::Align::AlignI>

=cut

sub next_aln {
    my ($self) = @_;
    my $aln;
    my %seqs;
    my @order;
    while( defined ($_ = $self->_readline ) ) {
	next if( /^\s+$/);
	if( !defined $aln ) {
	    $aln = Bio::SimpleAlign->new();
	}
	my ($id,$s) = split;
	push @order, $id if( ! defined $seqs{$id});
	$seqs{$id} .= $s;
    }
    foreach my $id ( @order) {
    my $gaps = $seqs{$id} =~ tr/-/-/;
	my $seq = Bio::LocatableSeq->new(-seq => $seqs{$id},
					-id  => $id,
					-start => 1,
					-end   => length($seqs{$id}) - $gaps,
					 -alphabet => $self->alphabet,
                    );
	$aln->add_seq($seq);
    }
    return $aln if defined $aln && $aln->num_sequences;
	return;
}

=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the NCBI psi-format object (.aln) into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Align::AlignI object

L<Bio::Align::AlignI>

=cut

sub write_aln {
	my ($self,$aln) = @_;
	unless( defined $aln && ref($aln) && 
			  $aln->isa('Bio::Align::AlignI') ) {
		$self->warn("Must provide a valid Bio::Align::AlignI to write_aln");
		return 0;
	}
	my $ct = 0;
	my @seqs = $aln->each_seq;
	my $len = 1;
	my $alnlen = $aln->length;
	my $idlen = $IdLength;
	my @ids = map { substr($_->display_id,0,$idlen) } @seqs;
	while( $len < ($alnlen + 1) ) {
		my $start = $len;
		my $end   = $len + $BlockLen;
		$end = $alnlen if ( $end > $alnlen ); 
		my $c = 0;
		foreach my $seq ( @seqs ) {
			$self->_print(sprintf("%-".$idlen."s %s\n",
										 $ids[$c++],
										 $seq->subseq($start,$end)));
		}
		$self->_print("\n");
		$len += $BlockLen+1;
	}
	$self->flush if $self->_flush_on_write && defined $self->_fh;
	return 1;
}

1;
