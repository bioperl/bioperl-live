#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Keith James
#
# Copyright Genome Research Ltd.
#
# You may distribute this module under the same terms as Perl itself

# POD documentation - main docs before the code


=encoding utf-8

=head1 NAME

Bio::Tools::Geneid - Results of one geneid run

=head1 SYNOPSIS

  use Bio::Tools::Geneid;
  my $gid = Bio::Tools::Geneid(-file => "geneid.out");

  while (my $gene = $gid->next_prediction)
  {
    my @transcripts = $gene->transcripts;
      foreach my $t (@transcripts)
      {
        my @exons = $t->exons;
        foreach my $e (@exons)
        {
          printf("Exon %d..%d\n", $e->start, $e->end);
        }
      }
  }

=head1 DESCRIPTION

This is the parser for the output of geneid by Enrique Blanco and
Roderic GuigE<243> (IMIM-UPF). See http://www1.imim.es/software/geneid. It
relies on native geneid output format internally and will work with
geneid versions 1.0 and 1.1. Currently this module supports only the
default mode of operation which is to predict exons and assemble an
optimal gene prediction.

It takes either a file handle or a file name and returns a
Bio::SeqFeature::Gene::GeneStructure object.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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

=head1 AUTHOR - Keith James

 Email: kdj@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Tools::Geneid;

use vars qw($SOURCE_TAG);
use strict;

use Bio::Tools::AnalysisResult;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Gene::Exon;
use Bio::SeqFeature::Gene::Transcript;
use Bio::SeqFeature::Gene::GeneStructure;

use base qw(Bio::Root::Root Bio::Root::IO);
$SOURCE_TAG = 'geneid';

=head2 new

 Title   : new
 Usage   : $obj->new(-file = "<geneid.out");
           $obj->new(-fh => \*GI);
 Function: Constructor for geneid wrapper. Takes either a file
         : or filehandle
 Returns : L<Bio::Tools::Geneid>

=cut

sub new
{
    my($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    $self->_initialize_io(@args);
    return $self;
}

=head2 next_prediction

 Title   : next_prediction
 Usage   : while($gene = $geneid->next_prediction)
           {
               # do something
           }
 Function: Returns the gene structure prediction of the geneid result
           file. Call this method repeatedly until FALSE is returned.
 Returns : A Bio::SeqFeature::Gene::GeneStructure object
 Args    : None

=cut

sub next_prediction
{
    my ($self) = @_;

    my ($gene, $transcript, $current_gene_id);
    my $transcript_score = 0;

    my ($gene_id, $exon_type, $exon_start, $exon_end, $exon_score,
        $exon_strand, $start_phase, $end_phase, $start_sig_score,
        $end_sig_score, $coding_pot_score, $homol_score);

    while (defined($_ = $self->_readline))
    {
        $self->debug($_);

        s/^\s+//;
        s/\s+$//;

        # We have a choice of geneid, gff or XML formats. The native
        # geneid format has more information than gff. However, we
        # then need to perform the hack of extracting the sequence ID
        # from the header of the embedded Fasta file which comes after
        # the exon data, as it is not stored elsewhere. Ack.
        # the new version of geneID includes the sequence ID in a slightly
        # different format and a new "or" statement was added below
        # also removed "unless defined $self->_target_id;" inorder to continue
        # generating new sequence IDs.
        
        if (/^>(\S+)\|GeneId/ or /^# Sequence (\S+)/)
        {
            my $target_id = $1;
            $self->_target_id($target_id);
            next;
        }

        next unless (/(Single|First|Internal|Terminal)/);

        my @fields = split(/\s+/, $_);

        # Grab gene_id from eol first as there are issues with
        # inconsistent whitespace in the AA coords field
        $gene_id = pop @fields;

        ($exon_type, $exon_start, $exon_end, $exon_score,
         $exon_strand, $start_phase, $end_phase, $start_sig_score,
         $end_sig_score, $coding_pot_score, $homol_score) = @fields[0..10];

        if (! defined $current_gene_id)
        {
            # Starting the requested prediction
            $current_gene_id = $gene_id;
            $transcript_score = $exon_score;

            $gene = Bio::SeqFeature::Gene::GeneStructure->new(-source =>
                                                              $SOURCE_TAG);
            $transcript = Bio::SeqFeature::Gene::Transcript->new(-source =>
                                                                 $SOURCE_TAG);

            $self->_add_exon($gene, $transcript, $exon_type, $exon_start, $exon_end, $exon_score,
                             $exon_strand, $start_phase, $end_phase, $start_sig_score,
                             $end_sig_score, $coding_pot_score, $homol_score);
        }
        elsif ($gene_id eq $current_gene_id)
        {
            # Still in requested prediction
            $transcript_score += $exon_score;

            $self->_add_exon($gene, $transcript, $exon_type, $exon_start, $exon_end, $exon_score,
                             $exon_strand, $start_phase, $end_phase, $start_sig_score,
                             $end_sig_score, $coding_pot_score, $homol_score);
        }
        else
        {
            # Found following prediction
            $self->_pushback($_);
            last;
        }
    }

    if (defined $gene)
    {
        $transcript->seq_id($self->_target_id);
        $transcript->score($transcript_score);
        $gene->add_transcript($transcript);
        $gene->seq_id($self->_target_id);

        foreach my $exon ($gene->exons)
        {
            $exon->seq_id($self->_target_id);
        }

        $self->_set_strand($gene);
    }

    return $gene;
}

=head2 _add_exon

 Title   : _add_exon
 Usage   : $obj->_add_exon($gene, $transcript, ... exon data ...)
 Function: Adds a new exon to both gene and transcript from the data
         : supplied as args
 Example :
 Returns : Nothing

=cut

sub _add_exon
{
    my ($self, $gene, $transcript, $exon_type, $exon_start, $exon_end,
        $exon_score, $exon_strand, $start_phase, $end_phase, $start_sig_score,
        $end_sig_score, $coding_pot_score, $homol_score) = @_;

    $exon_type =~ s/First/Initial/;

    my $strand = $exon_strand eq '+' ? 1 : -1;

    my $exon = Bio::SeqFeature::Gene::Exon->new(-source => $SOURCE_TAG,
                                                -start  => $exon_start,
                                                -end    => $exon_end,
                                                -strand => $strand,
                                                -score  => $exon_score);
    $exon->is_coding(1);
    $exon->add_tag_value("Type", $exon_type);
    $exon->add_tag_value('phase', $start_phase);
    $exon->add_tag_value('end_phase', $end_phase);
    $exon->add_tag_value('start_signal_score', $start_sig_score);
    $exon->add_tag_value('end_signal_score', $end_sig_score);
    $exon->add_tag_value('coding_potential_score', $coding_pot_score);
    $exon->add_tag_value('homology_score', $homol_score);

    $transcript->strand($strand) unless $transcript->strand != 0;
    $transcript->add_exon($exon, $exon_type);
}

=head2 _set_strand

 Title   : _set_strand
 Usage   : $obj->_set_strand($gene)
 Function: Sets the overall gene strand to the same strand as all
         : the exons if they are all on the same strand, or to strand 0
         : if the exons are on different strands.
 Example :
 Returns : Nothing

=cut

sub _set_strand
{
    my ($self, $gene) = @_;

    my $fwd = 0;
    my $rev = 0;

    my @exons = $gene->exons;
    foreach my $exon (@exons)
    {
        my $strand = $exon->strand;

        if ($strand == 1)
        {
            $fwd++;
        }
        elsif ($strand == -1)
        {
            $rev++;
        }
    }

    if ($#exons == $fwd)
    {
        $gene->strand(1);
    }
    elsif ($#exons == $rev)
    {
        $gene->strand(-1);
    }
    else
    {
        $gene->strand(0);
    }

    return $gene;
}

=head2 _target_id

 Title   : _target_id
 Usage   : $obj->_target_id
 Function: get/set for genomic sequence id
 Example :
 Returns : A target ID

=cut

sub _target_id
{
    my ($self,$val) = @_;
    if ($val)
    {
        $self->{'_target_id'} = $val;
    }

    return $self->{'_target_id'};
}

1;
