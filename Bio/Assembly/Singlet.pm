# $Id$
#
# BioPerl module for Bio::Assembly::Singlet
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
    $aio = Bio::Assembly::IO->new(-file=>"test.ace.1",
                               -format=>'phrap');

    $assembly = $aio->next_assembly;
    foreach $singlet ($assembly->all_singlets) {
      # do something
    }

=head1 DESCRIPTION

A singlet is a sequence that phrap was unable to align to any other sequences.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

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
use Bio::Seq::PrimaryQual;
use Dumpvalue();
my $dumper = new Dumpvalue();
$dumper->veryCompact(1);
use base qw(Bio::Assembly::Contig Bio::Root::Root Bio::Align::AlignI);


sub new {
    my ($class,%ARG) = @_;
    my $self = $class->SUPER::new(%ARG);
    my $args = \%ARG;
    bless ($self,$class);
    if ($args->{'-seq'}) {
        $self->seq_to_singlet($args->{'-seq'});
    }
    return $self;
}

=head2 seq_to_singlet

    Title   : seq_to_singlet
    Usage   : my $singlet = $io->seq_to_singlet($seq)
    Function: Wrap the information for a singlet as a Bio::Assembly::Singlet
    Returns : A Bio::Assembly::Singlet object
    Args    : A Bio::Seq-compliant object

=cut

sub seq_to_singlet {
    my ($self, $seq) = @_;
    $self->seqref($seq);
    $self->strand(1);
    my $lseq = Bio::LocatableSeq->new(
        -seq   => $seq->seq(),
        -start => 1,
        -end   => $seq->length(),
        -id    => $seq->display_id() );
    $lseq->{chromatfilename} = $seq->{'chromatfilename'};
    $lseq->{phdfilename} = $seq->{'phdfilename'};
    $self->set_consensus_sequence($lseq);
    if (UNIVERSAL::isa($seq,"Bio::Seq::Quality")) {
        $self->set_consensus_quality($seq)
    } else {
        # print("seq_to_singlet: the sequence (".$seq->desc().") is not a Bio::Seq::quality. it is this ($seq)\n");
    }
    $self->add_seq($lseq);
}


=head2 id

    Title   : id
    Usage   : my $id = $singlet->id('chad matsalla')
    Function: 
    Returns : 
    Args    : 

=cut

sub id {
    my $self = shift;
    # print("Getting the id for this thing:\n");
    # $dumper->dumpValue($self->seqref());
    # print("This is the id: (".$self->seqref()->id().")\n");
    my $id = undef;
    if (defined($self->seqref())) {
      $id = $self->seqref()->id();
    } else {
      $self->warn("This singlet has no ID because no Bio::Seq-compliant ".
        "sequence is attached to it");
    }
    return $id;
}


=head2 seqref

    Title   : seqref
    Usage   : my $seqref = $singlet->seqref($seq);
    Function: Set the sequence to which this Singlet refers
    Returns : A Bio::Seq-compliant object
    Args    : 

=cut

sub seqref {
    my ($self,$seq) = @_;
    if ($seq) { $self->{'seqref'} = $seq; }
    return $self->{'seqref'};
}


=head2 chromatfilename

    Title   : chromatfilename
    Usage   : my $chromatfilename = $singlet->chromatfilename($newfilename);
    Function: Get the name of the chromatfile for this singlet
    Returns : A string.
    Args    : If a string is provided, the chromatfilename will be set to that value.

=cut

sub chromatfilename {
    my ($self,$name) = @_;
    if ($name) { $self->{'chromatfilename'} = $name; }
    return $self->{'chromatfilename'};
}

=head2 phdfilename

    Title   : phdfilename
    Usage   : my $phdfilename = $singlet->phdfilename($newfilename);
    Function: Get the name of the phdfile for this singlet
    Returns : A string.
    Args    : If a string is provided, the phdfilename will be set to that value.

=cut

sub phdfilename {
    my ($self,$name) = @_;
    if ($name) { $self->{phdfilename} = $name; }
    return $self->{'phdfilename'};
}


1;
