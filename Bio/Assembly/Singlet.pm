#
# BioPerl module for Bio::Assembly::Singlet
# 
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chad Matsalla <bioinformatics1 at dieselwurks.com>
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Assembly::Singlet - Perl module to hold and manipulate
                     singlets from sequence assembly contigs.

=head1 SYNOPSIS

    # Module loading
    use Bio::Assembly::IO;

    # Assembly loading methods
    $aio = Bio::Assembly::IO->new( -file   => 'test.ace.1',
                                   -format => 'phrap'      );

    $assembly = $aio->next_assembly;
    foreach $singlet ($assembly->all_singlets) {
      # do something
    }

    # OR, if you want to build the singlet yourself,

    use Bio::Assembly::Singlet;
    $singlet = Bio::Assembly::Singlet->new(
        -id     => 'Singlet1',
        -seqref => $seq
    );

=head1 DESCRIPTION

A singlet is a sequence that phrap was unable to align to any other sequences.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Chad S. Matsalla

bioinformatics1 at dieselwurks.com

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#'
package Bio::Assembly::Singlet;

use strict;

use Bio::SeqFeature::Collection;
use Bio::LocatableSeq;
use Bio::Seq::PrimaryQual;

use base qw(Bio::Assembly::Contig Bio::Root::Root Bio::Align::AlignI);

=head2 new

    Title   : new
    Usage   : $singlet = $io->new( -seqref => $seq )
    Function: Create a new singlet object
    Returns : A Bio::Assembly::Singlet object
    Args    : -seqref => Bio::Seq-compliant sequence object for the singlet

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($seqref) = $self->_rearrange([qw(SEQREF)], @args);
    $self->{'_seqref'} = undef;
    if (defined $seqref) {
        $self->seqref($seqref);   
    }
    return $self;
}

=head2 seqref

    Title   : seqref
    Usage   : $seqref = $singlet->seqref($seq);
    Function: Get/set the sequence to which this singlet refers
    Returns : A Bio::Seq-compliant object
    Args    : A Bio::Seq-compliant or Bio::Seq::Quality object

=cut

sub seqref {
    my ($self,$seq) = @_;
    if (defined $seq) { $self->_seq_to_singlet($seq) };
    return $self->{'_seqref'};
}

=head2 _seq_to_singlet

    Title   : _seq_to_singlet
    Usage   : $singlet->seqref($seq)
    Function: Transform a sequence into a singlet
    Returns : 1 for sucess
    Args    : A Bio::Seq-compliant object

=cut

sub _seq_to_singlet {
    my ($self, $seq) = @_;
    # Object type checking
    $self->throw("Unable to process non Bio::Seq-compliant object [".ref($seq)."]")
      unless ( defined $seq && ($seq->isa('Bio::PrimarySeqI') || $seq->isa('Bio::Seq::Quality')) );
    # Sanity check
    $self->throw("Unable to have more than one sequence reference in a singlet")
      if (defined $self->{'_seqref'});
    # From sequence to locatable sequence
    my $lseq = Bio::LocatableSeq->new(
        -id     => $seq->id,
        -seq    => $seq->seq,
        -strand => $seq->isa('Bio::LocatableSeq') ? $seq->strand : 1,
        -start  => 1,
        #-end   => we let Bio::LocatableSeq calculate it (Seq and LocatableSeq)
    );
    # Get End from $seq if $lseq can't figure it out (e.g. phrap output)
    if (not defined $lseq->end) {
        $lseq->end($seq->end);
    }
    # Add new sequence and its coordinates to the contig
    my $lcoord = Bio::SeqFeature::Generic->new( -start => $lseq->start,
                                                -end   => $lseq->end    );
    $self->set_seq_coord( $lcoord, $lseq );
    $self->{'_seqref'} = $lseq;
    # Creating consensus
    $self->set_consensus_sequence($lseq);
    if ($seq->isa("Bio::Seq::Quality")) {
        my $qual = Bio::Seq::PrimaryQual->new( -id   => $seq->id,
                                               -qual => $seq->qual );
        $self->set_consensus_quality($qual);
    }
    return 1;
}

1;
