# BioPerl module for Bio::PrimedSeq
#
# Cared for by Chad Matsalla <bioinformatics1@dieselwurks.com>
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 Bio::Seq::PrimedSeq

Bio::Seq::PrimedSeq - A representation of a sequence and two primers flanking a
target region for amplification

=head1 SYNOPSIS

  # create a sequence
  my $sequence = "ctagctagctagctagctagctagctagctgatcgtagctagctagct";
  # create left and right primer seqfeatures
  # unfortunately, I haven't created constructors for these yet.
  my $left = Bio::SeqFeature::Primer();
  my $right = Bio::SeqFeature::Primer();
  # now create the PrimedSeq
  $primedseq = new Bio::Seq::PrimedSeq(
                                       -seq => $sequence,
                                       -display_id => "chads_fantastic_sequence",
                                       -LEFT_PRIMER => $left,
                                       -RIGHT_PRIMER => $right,
                                       -TARGET => '513,26'
                                       -PRIMER_PRODUCT_SIZE_RANGE => '100-500'
                                       -PRIMER_FILE_FLAG => '0'
                                       -PRIMER_LIBERAL_BASE => '1'
                                       -PRIMER_NUM_RETURN => '1'
                                       -PRIMER_FIRST_BASE_INDEX => '1'
                                       -PRIMER_EXPLAIN_FLAG => '1'
                                       -PRIMER_PRODUCT_SIZE => '185'
                                       );
  # get the amplified region
  my $amplified_sequence = $primed_seq->get_amplified_sequence();

=head1 DESCRIPTION

This module is a slightly glorified capsule containg a primed seqBuence. It was
created to address the fact that a primer is more the a seqfeature and there
need to be ways to represent the primer-sequence complex and the behaviors and
attributes that are associated with the complex.

=head1 FEEDBACK

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::PrimedSeq;
use vars qw(@ISA);
use strict;

use Bio::RangeI;

@ISA = qw(Bio::Seq);


=head2 new

 Title   : new()
 Usage   : $primed_sequence = new Bio::SeqFeature::Primer( -seq => $sequence,
                                                            -left_primer => $left_primer,
                                                            -right_primer => $right_primer);
 Function: A constructor for an object representing a primed sequence 
 Returns : A Bio::Seq::PrimedSeq object
 Args    : 
     -seq => a Bio::Seq object
     -left_primer => a Bio::SeqFeature::Primer object
     -right_primer => a Bio::SeqFeature::Primer object
     Many other parameters can be included including all of the output
     parameters from the primer3 program.
Developer Notes: This is incomplete and doesn't work. As of ISMB2002 I am working on it.


=cut

sub new {
     my($class,@args) = @_;
     my %arguments = @args;
     my $self = $class->SUPER::new(@args);
          # these are the absolute minimum components required to make
          # a primedseq
     my $newkey;
     foreach my $key (sort keys %arguments) {
          ($newkey = $key) =~ s/-//;
          $self->{$newkey} = $arguments{$key};
          push @{$self->{arguments}},$newkey;
     }
          # and now the insurance- make sure that things are ok
     if (!$self->{target_sequence} || !$self->{left_primer} || !$self->{right_primer} ) {
          $self->throw("You must provide a target_sequence, left_primer, and right_primer to create this object.");
     }
     if (ref($self->{target_sequence}) ne "Bio::Seq") {
          $self->throw("The target_sequence must be a Bio::Seq to create this object.");
     }
     if (ref($self->{left_primer}) ne "Bio::SeqFeature::Primer" || ref($self->{right_primer}) ne "Bio::SeqFeature::Primer") {
          $self->throw("You must provide a left_primer and right_primer, both as Bio::SeqFeature::Primer to create this object.");
     }
     return $self;
}


=head2 get_left_primer

 Title   : get_left_primer();
 Usage   : $left_primer = $primedseq->get_left_primer();
 Function: A getter for the left primer in thie PrimedSeq object.
 Returns : A Bio::SeqFeature::Primer object
 Args    : None.

=cut

sub get_left_primer() {
     my $self = shift;




}












=head2 Bio::RangeI methods

List of interfaces inherited from Bio::RangeI (see L<Bio::RangeI>
for details).

=head2 start

 Title   : start
 Usage   : $start = $feat->start
 Function: Returns the start coordinate of the feature
 Returns : integer
 Args    : none
Developer Notes:
          This is entirely dependent on the sequence to which this primer is attached!
          I think that there could be trouble if one takes this primer from sequence 1
          and naively place it on sequence 2 without updating this
          ** This is incomplete at this time.
=cut

sub start() {
     my $self = shift;


}




=head2 end

 Title   : end
 Usage   : $end = $feat->end
 Function: Returns the end coordinate of the feature
 Returns : integer
 Args    : none
Developer Notes:
          ** This is incomplete at this time.
=cut

sub end() {
     my $self = shift;


}

=head2 strand

 Title   : strand
 Usage   : $strand = $feat->strand()
 Function: Returns strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : none
Developer Notes:
          ** This is incomplete at this time.


=cut

sub strand() {
     my $self = shift;
}


=head2 SeqFeatureI specific methods

New method interfaces.

=head2 sub_SeqFeature

 Title   : sub_SeqFeature
 Usage   : @feats = $feat->sub_SeqFeature();
 Function: Returns an array of sub Sequence Features
 Returns : An array
 Args    : none

=cut

sub sub_SeqFeature{
   my ($self,@args) = @_;

   $self->throw_not_implemented();
}

=head2 display_id

 Title   : display_id
 Usage   : $name = $feat->display_id()
 Function: Returns the human-readable ID of the
           feature for displays.
 Returns : a string
 Args    : none

=cut

sub display_id { 
  my ($self,@args) = @_;
  $self->throw_not_implemented();
}

=head2 primary_tag

 Title   : primary_tag
 Usage   : $tag = $feat->primary_tag()
 Function: Returns the primary tag for a feature,
           eg 'exon'
 Returns : a string 
 Args    : none


=cut

sub primary_tag{
   my ($self,@args) = @_;

   $self->throw_not_implemented();

}

=head2 source_tag

 Title   : source_tag
 Usage   : $tag = $feat->source_tag()
 Function: Returns the source tag for a feature,
           eg, 'genscan' 
 Returns : a string 
 Args    : none


=cut

sub source_tag{
   my ($self,@args) = @_;

   $self->throw_not_implemented();
}

=head2 has_tag

 Title   : has_tag
 Usage   : $tag_exists = $self->has_tag('some_tag')
 Function: 
 Returns : TRUE if the specified tag exists, and FALSE otherwise
 Args    :


=cut

sub has_tag{
   my ($self,@args) = @_;

   $self->throw_not_implemented();

}

=head2 each_tag_value

 Title   : each_tag_value
 Usage   : @values = $self->each_tag_value('some_tag')
 Function: 
 Returns : An array comprising the values of the specified tag.
 Args    :


=cut

sub each_tag_value {
   my ($self,@args) = @_;

   $self->throw_not_implemented();
}

=head2 all_tags

 Title   : all_tags
 Usage   : @tags = $feat->all_tags()
 Function: gives all tags for this feature
 Returns : an array of strings
 Args    : none


=cut

sub all_tags{
   my ($self,@args) = @_;

   $self->throw_not_implemented();
}

=head2 gff_string

 Title   : gff_string
 Usage   : $str = $feat->gff_string;
           $str = $feat->gff_string($gff_formatter);
 Function: Provides the feature information in GFF format.

           The implementation provided here returns GFF2 by default. If you
           want a different version, supply an object implementing a method
           gff_string() accepting a SeqFeatureI object as argument. E.g., to
           obtain GFF1 format, do the following:

                my $gffio = Bio::Tools::GFF->new(-gff_version => 1);
                $gff1str = $feat->gff_string($gff1io);

 Returns : A string
 Args    : Optionally, an object implementing gff_string().


=cut

sub gff_string{
   my ($self,$formatter) = @_;

   $formatter = $self->_static_gff_formatter unless $formatter;
   return $formatter->gff_string($self);
}

my $static_gff_formatter = undef;

=head2 _static_gff_formatter

 Title   : _static_gff_formatter
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _static_gff_formatter{
   my ($self,@args) = @_;

   if( !defined $static_gff_formatter ) {
       $static_gff_formatter = Bio::Tools::GFF->new('-gff_version' => 2);
   }
   return $static_gff_formatter;
}



=head1 RangeI methods

These methods are inherited from RangeI and can be used
directly from a SeqFeatureI interface. Remember that a 
SeqFeature is-a RangeI, and so wherever you see RangeI you
can use a feature ($r in the below documentation).

=head2 overlaps

  Title   : overlaps
  Usage   : if($feat->overlaps($r)) { do stuff }
            if($feat->overlaps(200)) { do stuff }
  Function: tests if $feat overlaps $r
  Args    : a RangeI to test for overlap with, or a point
  Returns : true if the Range overlaps with the feature, false otherwise


=head2 contains

  Title   : contains
  Usage   : if($feat->contains($r) { do stuff }
  Function: tests whether $feat totally contains $r
  Args    : a RangeI to test for being contained
  Returns : true if the argument is totaly contained within this range


=head2 equals

  Title   : equals
  Usage   : if($feat->equals($r))
  Function: test whether $feat has the same start, end, strand as $r
  Args    : a RangeI to test for equality
  Returns : true if they are describing the same range


=head1 Geometrical methods

These methods do things to the geometry of ranges, and return
triplets (start, stop, strand) from which new ranges could be built.

=head2 intersection

  Title   : intersection
  Usage   : ($start, $stop, $strand) = $feat->intersection($r)
  Function: gives the range that is contained by both ranges
  Args    : a RangeI to compare this one to
  Returns : nothing if they do not overlap, or the range that they do overlap

=head2 union

  Title   : union
  Usage   : ($start, $stop, $strand) = $feat->union($r);
          : ($start, $stop, $strand) = Bio::RangeI->union(@ranges);
  Function: finds the minimal range that contains all of the ranges
  Args    : a range or list of ranges to find the union of
  Returns : the range containing all of the ranges

=cut

=head2 location

 Title   : location
 Usage   : my $location = $seqfeature->location()
 Function: returns a location object suitable for identifying location 
	   of feature on sequence or parent feature  
 Returns : Bio::LocationI object
 Args    : none


=cut

sub location {
   my ($self) = @_;

   $self->throw_not_implemented();
}


1;
