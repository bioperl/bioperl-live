
#
# BioPerl module for Bio::SeqFeature::Gene
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Gene - Object holding one particular gene, with Transcripts, Translations and Exons hanging off it

=head1 SYNOPSIS

    $gene  = Bio::SeqFeature::Gene->new();
    
    # build up a set of exons into an array
    foreach $exonpair ( split(/:/,"10,30:50,80:110,130") ) {
	($start,$end) = split(/,/,$exonpair);
	push(@exon,Bio::SeqFeature::Exon->new(-start => $start,-end => $end, -strand => 1));
    }

    # add them to a new Transcript
    $trans = Bio::SeqFeature::Transcript->new();
    $trans->add_Exon(@exon);	

    # assumme that the translation here is made from this transcript
    # starts at 20 in genomic coordinates and ends at 120
    $tls = $trans->create_Translation(20,120);
    $trans->Translation($tls);
 
    #alternative transcript with different translation
    pop(@exon);
    push(@exon,Bio::SeqFeature::Exon->new(-start => '90',-end => '120',-strand => 1));
    # add them to a new Transcript
    $alt1 = Bio::SeqFeature::Transcript->new();
    $alt1->add_Exon(@exon);	
    $alt_tls = $trans->create_Translation(20,100);
    $alt1->Translation($alt_tls);
    

    # alternative transcript with a different last exon
    # but the same translation as the first Transcript
    $alt2 = Bio::SeqFeature::Transcript->new();
    pop(@exon);
    push(@exon,Bio::SeqFeature::Exon->new(-start => '110',-end => '135',-strand => 1));
    $alt2->add_Exon(@exon);	

    # alt has the same translation as trans, even though a different
    # splicing structure
    $alt2->Translation($tls);

    # add things to the gene.

    $gene->add_Transcript($trans,$alt1,$alt2);
    # no need to add Translations. Translations must be added
    # associated with at least one Transcript

    # gene can now be added to an annseq. This is the only
    # valid way to associate a sequence with a gene.

    # making a new AnnSeq object from scratch. See Bio::AnnSeq
    $instream = Bio::SeqIO(-format =>'Fasta',-file =>'myseq.fa');
    $seq = $instream->next_seq();
    $annseq = Bio::AnnSeq->new();
    $annseq->seq($seq);
    
    $annseq->add_SeqFeature($gene);
    # accessing stuff.

    # SeqIO streams to write out sequences
    $pepstream = Bio::SeqIO->new(-format=> 'Fasta', -file => '>MyTranslations.fa');
    $cdnastream = Bio::SeqIO->new(-format=> 'Fasta', -file => '>MycDNA.fa');

    # gives you three transcripts.
    foreach $transcript ( $gene->each_Transcript() ) {
	if( $transcript->is_protein_coding() ) {
	    print "Transcript is protein coding\n";
	}
	foreach $exon ( $transcript->each_Exon() ) {
	    print "Exon ",$exon->start,":",$exon->end(), "\n";
	}
	$cdnastream->write_seq($transcript->seq());
    }

    # gives you only two translations.
    foreach $translation ( $gene->each_Translation() ) {
	$pepstream->write_seq($translation->seq());
    }


=head1 DESCRIPTION

This object represents genes on dna sequences. It basically provides
a coordinating object which can be used to hold transcripts and
translations, which are the two key objects that manage genes. Transcripts
manage the alternative processing products of a gene: translations manage
the alternative protein products of a gene. 

I need to write more docs about this!

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@sanger.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Gene;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

@ISA = qw(Bio::SeqFeature::Generic);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  $self->{'_transcript_hash'} = {};
  $self->{'_translation_hash'} = {};

  # set stuff in self from @args
  return $make; # success - we hope!
}

=head2 each_Transcript

 Title   : each_Transcript
 Usage   : foreach my $trans ( $gene->each_Transcript() )
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_Transcript{
   my ($self) = @_;
   return values %{$self->{'_transcript_hash'}};
}

=head2 each_Translation

 Title   : each_Translation
 Usage   : foreach my $tls ( $gene->each_Translation() )
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_Translation{
   my ($self) = @_;
   return values %{$self->{'_translation_hash'}};
}


=head2 add_Transcript

 Title   : add_Transcript
 Usage   : $gene->add_Transcript($trans)
 Function: adds a Transcript to this gene. If there is a translation 
           associated with a gene, adds this 
 Example :
 Returns : 
 Args    :


=cut

sub add_Transcript{
   my ($self) = shift @_;
   if( @_ == 0 ) {
       $self->warn("add_Transcript called with no Transcripts");
   }

   foreach my $trans ( @_ ) {
       if( ! ref $trans || !$trans->isa('Bio::SeqFeature::Transcript') ) {
	   $self->throw("add_Transcript must take a transcript, not a $trans");
       }
       my $key = "$trans";
       $self->{'_transcript_hash'}->{$key} = $trans;
       if( $trans->is_protein_coding ) {
	   $self->_add_Translation($trans->Translation);
       }
   }

}
=head2 attach_seq

 Title   : attach_seq
 Usage   : $gene->attach_seq
 Function: attaches sequence to genes/transcripts/translations.
 Example :
 Returns : 
 Args    :


=cut

sub attach_seq{
   my ($self,$seq) = @_;
   
   foreach my $trans ( $self->each_Transcript() ) {
       $trans->_attach_seq($seq);
   }
}


=head2 _add_Translation

 Title   : _add_Translation
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _add_Translation {
   my ($self) = shift @_;
   if( @_ == 0 ) {
       $self->warn("add_Translation called with no Translations");
   }


   foreach my $trans ( @_ ) {
       if( ! ref $trans || !$trans->isa('Bio::SeqFeature::Translation') ) {
	   $self->throw("add_Translation takes Translation, not a $trans");
       }
       my $key = "$trans";
       $self->{'_translation_hash'}->{$key} = $trans;
   }

}





