# $Id$
#
# BioPerl module for Bio::Align::Utilities
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Align::Utilities - A collection of utilities regarding converting
and manipulating alignment objects

=head1 SYNOPSIS

  use Bio::Align::Utilities qw(aa_to_dna_aln);
  # %dnaseqs is a hash of CDS sequences (spliced)


  # Even if the protein alignments are local make sure the start/end
  # stored in the LocatableSeq objects are to the full length protein.
  # The CoDing Sequence that is passed in should still be the full 
  # length CDS as the nt alignment will be generated.
  #
  my $dna_aln = aa_to_dna_aln($aa_aln,\%dnaseqs);


=head1 DESCRIPTION

This module contains utility methods for manipulating sequence
alignments ( L<Bio::Align::AlignI>) objects.

The B<aa_to_dna_aln> utility is essentially the same as the B<mrtrans>
program by Bill Pearson available at
ftp://ftp.virginia.edu/pub/fasta/other/mrtrans.shar.  Of course this
is a pure-perl implementation, but just to mention that if anything
seems odd you can check the alignments generated against Bill's
program.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

#' keep my emacs happy
# Let the code begin...


package Bio::Align::Utilities;
use vars qw(@ISA @EXPORT @EXPORT_OK $GAP $CODONGAP);
use strict;
use Carp;
use Bio::Root::Version;
require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw();
@EXPORT_OK = qw(aa_to_dna_aln);

BEGIN {
    use constant CODONSIZE => 3;
    $GAP = '-';
    $CODONGAP = $GAP x CODONSIZE;
}

=head2 aa_to_dna_aln

 Title   : aa_to_dna_aln
 Usage   : my $dnaaln = aa_to_dna_aln($aa_aln, \%seqs);
 Function: Will convert an AA alignment to DNA space given the 
           corresponding DNA sequences.  Note that this method expects 
           the DNA sequences to be in frame +1 (GFF frame 0) as it will
           start to project into coordinates starting at the first base of 
           the DNA sequence, if this alignment represents a different 
           frame for the cDNA you will need to edit the DNA sequences
           to remove the 1st or 2nd bases (and revcom if things should be).
 Returns : Bio::Align::AlignI object 
 Args    : 2 arguments, the alignment and a hashref.
           Alignment is a Bio::Align::AlignI of amino acid sequences. 
           The hash reference should have keys which are 
           the display_ids for the aa 
           sequences in the alignment and the values are a 
           Bio::PrimarySeqI object for the corresponding 
           spliced cDNA sequence. 

See also: L<Bio::Align::AlignI>, L<Bio::SimpleAlign>, L<Bio::PrimarySeq>

=cut

sub aa_to_dna_aln {
    my ($aln,$dnaseqs) = @_;
    unless( defined $aln && 
	    ref($aln) &&
	    $aln->isa('Bio::Align::AlignI') ) { 
	croak('Must provide a valid Bio::Align::AlignI object as the first argument to aa_to_dna_aln, see the documentation for proper usage and the method signature');
    }
    my $alnlen = $aln->length;
    #print "HSP length is $alnlen\n";
    my $dnaalign = new Bio::SimpleAlign;
    
    foreach my $seq ( $aln->each_seq ) {    
	my $aa_seqstr = $seq->seq();
	my $id = $seq->display_id;
	my $dnaseq = $dnaseqs->{$id} || $aln->throw("cannot find ".
						     $seq->display_id);
	my $start_offset = ($seq->start() - 1) * CODONSIZE;
	$dnaseq = $dnaseq->seq();
	my $dnalen = $dnaseqs->{$id}->length;
	my $nt_seqstr;
	my $j = 0;
	for( my $i = 0; $i < $alnlen; $i++ ) {
	    my $char = substr($aa_seqstr,$i + $start_offset,1);	    
	    if( $char eq '-' ) {
		$nt_seqstr .= $CODONGAP;
	    } else { 
		if( $j > $dnalen ) { 
		    $aln->warn("codons can't match up for $id, we've gone beyond the end of the DNA sequence\n");
		    next;
		}
		$nt_seqstr .= substr($dnaseq,$j,CODONSIZE);
		$j += CODONSIZE;
	    }
	}
	
        # funky looking math is to readjust to codon boundaries and deal
	# with fact that sequence start with 1
	my $newdna = new Bio::LocatableSeq(-display_id  => $id,
					   -alphabet    => 'dna',
					   -start       => $start_offset+1,
					   -end         => ($seq->end * 
							    CODONSIZE),
					   -strand      => 1,
					   -seq         => $nt_seqstr);    
	$dnaalign->add_seq($newdna);
    }
    return $dnaalign;
}


# This is the previous implementation of aa_to_dna_aln function
# which is ~98% slower.

sub OLD_aa_to_dna_aln {
    my ($aln,$dnaseqs) = @_;
    unless( defined $aln && 
	    ref($aln) &&
	    $aln->isa('Bio::Align::AlignI') ) { 
	croak('Must provide a valid Bio::Align::AlignI object as the first argument to aa_to_dna_aln, see the documentation for proper usage and the method signature');
    }
    my $alnlen = $aln->length;
    #print "HSP length is $alnlen\n";
    my $dnaalign = new Bio::SimpleAlign;
    foreach my $seq ( $aln->each_seq ) {    
	my $newseq;	    
	my $dnaseq = $dnaseqs->{$seq->display_id} || croak("cannot find ".
							 $seq->display_id);
	foreach my $pos ( 1..$alnlen ) {
	    my $loc = $seq->location_from_column($pos);
	    my $dna = ''; 
	    if( !defined $loc || $loc->location_type ne 'EXACT' ) {
		$dna = '---';
	    } else {
		# To readjust to codon boundaries
		# end needs to be +1 so we can just multiply by CODONSIZE 
		# to get this		    

		my ($start,$end) = ((($loc->start - 1)* CODONSIZE) +1,
				    ($loc->end)* CODONSIZE);
		
		if( $start <=0 || $end > $dnaseq->length() ) {
		    print STDERR "start is ", $loc->start, " end is ", $loc->end, " while dnaseq length is ", $dnaseq->length(), " and start/end projected are $start,$end \n";
		    warn("codons don't seem to be matching up for $start,$end");
		    $dna = '---';			    
		} else {
		    $dna = $dnaseq->subseq($start,$end);
		}
	    }
	    $newseq .= $dna;
	}
	# funky looking math is to readjust to codon boundaries and deal
	# with fact that sequence start with 1
	my $newdna = new Bio::LocatableSeq(-display_id  => $seq->id(),
					   -start => (($seq->start - 1) * 
						      CODONSIZE) + 1, 
					   -end   => ($seq->end * CODONSIZE),
					   -strand => $seq->strand,
					   -seq   => $newseq);    
	$dnaalign->add_seq($newdna);
    }
    return $dnaalign;
}

1;
