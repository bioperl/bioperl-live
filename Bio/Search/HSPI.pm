# $Id$
#
# BioPerl module for Bio::Search::HSPI
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::HSPI - A High Scoring Pair found when doing a sequence alignment.

=head1 SYNOPSIS

    # get a Bio::Search::HSPI object somehow - perhaps from a SubjectI object
    my $hsp = $subject->next_hsp;

=head1 DESCRIPTION

Describe the interface here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::HSPI;
use vars qw(@ISA);

use Bio::SeqFeature::SimilarityPair;

use strict;

@ISA = qw(Bio::SeqFeature::SimilarityPair);

=head2 report_type

 Title   : report_type
 Usage   : my $r_type = $hsp->report_type
 Function: Obtain the report type for an HSP
 Returns : string
 Args    : none


=cut

sub report_type{
   my ($self,@args) = @_;
   $self->_abstractDeath('report_type');
}

=head2 P

 Title   : P
 Usage   : my $pvalue = $hsp->P();
 Function: Returns the P-value for this HSP (same as evalue)
 Returns : float or exponential (2e-10)
 Args    : none

=cut

sub P {
   my ($self) = @_;
   $self->evalue;
}

=head2 evalue

 Title   : evalue
 Usage   : my $evalue = $hsp->evalue();
 Function: Returns the e-value for this HSP
 Returns : float or exponential (2e-10)
 Args    : none

=cut

sub evalue {
   my ($self) = @_;
   $self->_abstractDeath('evalue');
}

=head2 percent_identical

 Title   : percent_identical
 Usage   : my $perc_id = $hsp->percent_identical();
 Function: Returns the Percent Identity (num positive / total HSP len)
           for this HSP
 Returns : Float in range 0 -> 100
 Args    : none

=cut

sub percent_identical {
   my ($self) = @_;
   $self->_abstractDeath('percent_identical');
}

=head2 positive

 Title    : positive
 Usage    : $hsp->positive();
 Function : returns the number of positive matches (symbols in the alignment
            with a positive score)
 Returns  : (int) number of positive matches in the alignment
 Args     : none

=cut

sub positive        {
    my ($self) = @_;
    $self->_abstractDeath('positive');
}

=head2 gaps

 Title    : gaps
 Usage    : $hsp->gaps();
 Function : returns the number of gaps or 0 if none
 Returns  : (int) number of gaps or 0 if none
 Args     : none

=cut

sub gaps        {
    my ($self) = @_;
    $self->_abstractDeath('gaps');
}

=head2 query_seq

 Title   : query_seq
 Usage   : my $qseq = $hsp->query_seq;
 Function: Retrieves the query sequence that is part of this HSP
 Returns : string
 Args    : none


=cut

sub query_seq{
   my ($self,@args) = @_;
   $self->_abstractDeath('query_seq');
}

=head2 subject_seq

 Title   : subject_seq
 Usage   : my $sseq = $hsp->subject_seq;
 Function: Retrieves the subject sequence that is part of this HSP
 Returns : string
 Args    : none


=cut

sub subject_seq{
   my ($self,@args) = @_;
   $self->_abstractDeath('subject_seq');
}

=head2 homology_seq

 Title   : homology_seq
 Usage   : my $homo_seq = $hsp->homology_seq;
 Function: Retrieves the homology sequence for this HSP
 Returns : string
 Args    : none

=cut

sub homology_seq{
   my ($self,@args) = @_;
   $self->_abstractDeath('homology_seq');
}

=head2 hsp_length

 Title   : hsp_length
 Usage   : my $len = $hsp->hsp_length
 Function: Returns the aggregate length of the HSP
           (which may be greater than either subject or query )
 Returns : integer
 Args    : none

=cut

sub hsp_length{
   my ($self,@args) = @_;
   $self->_abstractDeath('hsp_length');
}

=head2 Bio::SeqFeature::SimilarityPair methods

=cut

=head2 significance

 Title   : significance
 Usage   : $evalue = $obj->significance();
           $obj->significance($evalue);
 Function:
 Returns :
 Args    :

=head2 bits

 Title   : bits
 Usage   : $bits = $obj->bits();
           $obj->bits($value);
 Function:
 Returns :
 Args    :

=head2 score

 Title   : score
 Usage   : $score = $obj->score();
           $obj->score($value);
 Function:
 Returns :

=head2 Bio::SeqFeature::FeaturePair methods

=cut

=head2 query

 Title   : query
 Usage   : my $query = $hsp->query;
 Function: Access to the SeqFeature::Similarity
           for the query sequence that makes up this HSP
 Returns : Bio::SeqFeature::Similarity
 Args    : none

=head2 subject

 Title   : subject
 Usage   : my $subject = $hsp->subject
 Function: Access to the SeqFeature::Similarity
           for the subject (hit) sequence that makes up this HSP
 Returns : Bio::SeqFeature::Similarity
 Args    : none

=head2 Bio::SeqFeature::FeaturePair methods

=cut

=head2 score

 Title   : score
 Usage   : $score = $obj->score();
           $obj->score($value);
 Function:
 Returns :

=head2 feature1

 Title   : feature1
 Usage   : $f = $featpair->feature1
           $featpair->feature1($feature)
 Function: Get/set for the query feature
 Returns : Bio::SeqFeatureI
 Args    : Bio::SeqFeatureI

=head2 feature2

 Title   : feature2
 Usage   : $f = $featpair->feature2
           $featpair->feature2($feature)
 Function: Get/set for the hit feature
 Returns : Bio::SeqFeatureI
 Args    : Bio::SeqFeatureI

=head2 hseqname

 Title   : hseqname
 Usage   : $featpair->hseqname($newval)
 Function: Get/set method for the name of
           feature2.
 Returns : value of $feature2->seqname
 Args    : newvalue (optional)

=head2 hstart

 Title   : hstart
 Usage   : $start = $featpair->hstart
           $featpair->hstart(20)
 Function: Get/set on the start coordinate of feature2
 Returns : integer
 Args    : none

=head2 hend

 Title   : hend
 Usage   : $end = $featpair->hend
           $featpair->hend($end)
 Function: get/set on the end coordinate of feature2
 Returns : integer
 Args    : none

=head2 hstrand

 Title   : hstrand
 Usage   : $strand = $feat->strand()
           $feat->strand($strand)
 Function: get/set on strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : none

=head2 hscore

 Title   : hscore
 Usage   : $score = $feat->score()
           $feat->score($score)
 Function: get/set on score information
 Returns : float
 Args    : none if get, the new value if set

=head2 hframe

 Title   : hframe
 Usage   : $frame = $feat->frame()
           $feat->frame($frame)
 Function: get/set on frame information
 Returns : 0,1,2
 Args    : none if get, the new value if set


=head2 hprimary_tag

 Title   : hprimary_tag
 Usage   : $ptag = $featpair->hprimary_tag
 Function: Get/set on the primary_tag of feature2
 Returns : 0,1,2
 Args    : none if get, the new value if set

=head2 hsource_tag

 Title   : hsource_tag
 Usage   : $tag = $feat->hsource_tag()
           $feat->source_tag('genscan');
 Function: Returns the source tag for a feature,
           eg, 'genscan'
 Returns : a string
 Args    : none

=head2 invert

 Title   : invert
 Usage   : $tag = $feat->invert
 Function: Swaps feature1 and feature2 around
 Returns : Nothing
 Args    : none

=head2 Bio::SeqFeatureI methods

=cut

=head2 start

 Title   : start
 Usage   : $start = $featpair->start
           $featpair->start(20)
 Function: Get/set on the start coordinate of feature1
 Returns : integer
 Args    : [optional] beginning of feature

=head2 end

 Title   : end
 Usage   : $end = $featpair->end
           $featpair->end($end)
 Function: get/set on the end coordinate of feature1
 Returns : integer
 Args    : [optional] ending point of feature

=head2 strand

 Title   : strand
 Usage   : $strand = $feat->strand()
           $feat->strand($strand)
 Function: get/set on strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : [optional] strand information to set

=head2 location

 Title   : location
 Usage   : $location = $featpair->location
           $featpair->location($location)
 Function: Get/set location object (using feature1)
 Returns : Bio::LocationI object
 Args    : [optional] LocationI to store

=head2 primary_tag

 Title   : primary_tag
 Usage   : $ptag = $featpair->primary_tag
 Function: get/set on the primary_tag of feature1
 Returns : 0,1,2
 Args    : none if get, the new value if set

=head2 source_tag

 Title   : source_tag
 Usage   : $tag = $feat->source_tag()
           $feat->source_tag('genscan');
 Function: Returns the source tag for a feature,
           eg, 'genscan'
 Returns : a string
 Args    : none

=head2 seqname

 Title   : seqname
 Usage   : $obj->seqname($newval)
 Function: There are many cases when you make a feature that you
           do know the sequence name, but do not know its actual
           sequence. This is an attribute such that you can store
           the seqname.

           This attribute should *not* be used in GFF dumping, as
           that should come from the collection in which the seq
           feature was found.
 Returns : value of seqname
 Args    : newvalue (optional)

1;
