
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

    # building a gene

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

The object design is as follows:

A Gene is a collection of transcript and translation objects.

A Transcript is a particular RNA transcription and splicing pattern:
it is built from Exon objects. A Transcript might be associated with
a Translation or not.

A Translation is a particular protein encoding path in genomic DNA.
A particular Translation can be associated with more than one Transcript.
Translations can only be created from Transcripts.

An Exon represents a single start/end point in DNA which is spliced. When
an Translation is made from a Transcript, it immeadaitely promotes the
Exon into having protein coding properties, including frame.

If the gene object is attached to an AnnSeq object, then automatic
dumps of the processed sequences can occur. From Transcripts one can
get Bio::Seq objects representing the cDNA of the Transcript. From
Translations one can get Bio::Seq objects representing the protein
translations.

=head2 Interactions with SeqFeature attributes

Genes, Transcripts, Translations and Exons are all SeqFeatureI compliant objects.
The notional SeqFeature tree that they have is as follows:

                 Gene
                  |	     
          ________|__________  	 
          |    	       	    |
     Transcript	       Translation
      	  |    	       	       	  
      ____|_____
     | 	       	|    
   Exon	      Intron(created)   	
      		     
The Intron seq features are different from the rest as they are created 'on-the-fly'
when the Transcript object is asked for sub_SeqFeatures. This guarentees consistency
with the Exon objects.


A drawback of this system is that for transcripts which share many Exons and Introns, the
objects get duplicated in the sub_SeqFeature descent. Therefore "all_SeqFeatures" called
from aseq can get confused.

=head1 DEVELOPERS NOTES

These objects are pretty confusing, with some interacting arrays and hashes
Bascially we build on SeqFeature::Generic via inheritence which gives us
automatic attributes, an array of seq features for clients of the objects
to use and inheritence of the SeqFeatureI system. We set the primary and
source tags to the correct types in the intialisation systems.

Alot of effort has been put in to make the system have no circular references.
This means that clients need to access certain information through the Gene
object whereas it might have been more natural to access through the Transcript/Translation
object.
       
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
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

@ISA = qw(Bio::SeqFeature::Generic);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  $self->primary_tag("Gene");
  $self->source_tag("Bioperl");
  
  $self->{'_transcript_hash'} = {};
  $self->{'_translation_hash'} = {};

  # set stuff in self from @args
  return $make; # success - we hope!
}

=head2 sub_SeqFeature

 Title   : sub_SeqFeature
 Usage   : @sf = $gene->sub_SeqFeature();
 Function: Retrieves sequence features which are part of this
           gene - in particular the transcripts and the translations
           will be there, along with the exons, and any other features
           attached to the object
 Example : 
 Returns : Array of SeqFeatureI compliant objects
 Args    : None

 Don't use this to get out Transcripts or Translations specifically:
 use each_Transcript or each_Translation. Exons can be retrieved using
 each_Exon and each_coding_Exon


=cut

sub sub_SeqFeature{
   my ($self) = @_;
   my @out;

   push(@out,$self->SUPER::sub_SeqFeature());
   push(@out,$self->each_Transcript());
   push(@out,$self->each_Translation());
   
   return @out;
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
       if( !defined $self->start() ) {
	   $self->start($trans->start);
	   $self->end($trans->end);
       } else {
	   my ($start,$end) = $self->union($trans);
	   $self->start($start);
	   $self->end($end);
       }
       $self->strand($trans->strand);
   }

}


=head2 each_Exon

 Title   : each_Exon
 Usage   : @exons = $gene->each_Exon()
 Function: Provides all exons of all transcripts in a single array
 Example :
 Returns : 
 Args    :


=cut

sub each_Exon{
   my ($self,@args) = @_;
   my %th;

   foreach my $trans ( $self->each_Transcript) {
       foreach my $exon ( $trans->each_Exon ) {
	   my $key = "$exon";
	   $th{$key} = $exon;
       } 
   }
   return values %th;
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
   foreach my $tls   ( $self->each_Translation() ) {
       $tls->_attach_seq($seq);
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



1;






