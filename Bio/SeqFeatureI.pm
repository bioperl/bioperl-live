
#
# BioPerl module for Bio::SeqFeatureI
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeatureI - Abstract interface of a Sequence Feature

=head1 SYNOPSIS

    # get a seqfeature somehow, eg,

    foreach $feat ( $annseq->all_seqfeatures() ) {
            print "Feature from ", $feat->start, "to ", 
	          $feat->end, " Primary tag  ", $feat->primary_tag, 
	          ", produced by ", $feat->source_tag() "\n";

            if( $feat->strand == 0 ) {
		print "Feature applicable to either strand\n";
            } else {
                print "Feature on strand ", $feat->strand,"\n"; # -1,1
            }

            foreach $tag ( $feat->all_tags() ) {
		print "Feature has tag ", $tag, "with values, ",
		      join(' ',@feat->each_tag_value($tag)), "\n";
            }
	    print "new feature\n" if $feat->has_tag('new');
	    # features can have sub features
	    my @subfeat = $feat->sub_SeqFeature();
	}

=head1 DESCRIPTION

This interface is the functions one can expect for any Sequence Feature, whatever
its implemtation or whether it is a more complex type (eg, a Gene). This object
doesn't actually provide any implemention, it just provides the definitions
of what methods one can call. See Bio::SeqFeature::Generic for a good standard
implementation of this object

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
  http://bio.perl.org/bioperl-bugs/

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut
#'

# Let the code begin...


package Bio::SeqFeatureI;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::RangeI;
use Carp;

@ISA = qw(Bio::RangeI);

# utility method
#
# Prints out a method like:
# Abstract method stop defined in interface Bio::SeqFeatureI not implemented by package You::BadFeature
sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller)[1];
  
  confess "Abstract method '$caller' defined in interface Bio::SeqFeatureI not implemented by pacakge $package";
}

=head2 start

 Title   : start
 Usage   : $start = $feat->start
 Function: Returns the start coordinate of the feature
 Returns : integer
 Args    : none


=cut

sub start{
   my ($self,@args) = @_;

   $self->_abstractDeath();
}

=head2 end

 Title   : end
 Usage   : $end = $feat->end
 Function: Returns the end coordinate of the feature
 Returns : integer
 Args    : none


=cut

sub end{
   my ($self,@args) = @_;

   $self->_abstractDeath();

}


=head2 strand

 Title   : strand
 Usage   : $strand = $feat->strand()
 Function: Returns strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : none


=cut

sub strand{
   my ($self,@args) = @_;

   $self->_abstractDeath();

}

=head2 sub_SeqFeature

 Title   : sub_SeqFeature
 Usage   : @feats = $feat->sub_SeqFeature();
 Function: Returns an array of sub Sequence Features
 Returns : An array
 Args    : none


=cut

sub sub_SeqFeature{
   my ($self,@args) = @_;

   $self->_abstractDeath();
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

   $self->_abstractDeath();

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

   $self->_abstractDeath();

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

   $self->_abstractDeath();

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

   $self->_abstractDeath();
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

   $self->_abstractDeath();

}

=head2 gff_string

 Title   : gff_string
 Usage   : $str = $feat->gff_string
 Function: provides the feature information in GFF
           version 2 format.
 Returns : A string
 Args    : None


=cut

sub gff_string{
   my ($feat) = @_;
   my ($str,$score,$frame,$name);

   if( $feat->can('score') ) {
       $score = $feat->score();
   }
   $score = '.' unless defined $score;

   if( $feat->can('frame') ) {
       $frame = $feat->frame();
   }
   $frame = '.' unless defined $frame;

   my $strand;
   if( $feat->strand == 1 ) {
       $strand = '+';
   } elsif ( $feat->strand == -1 ) {
       $strand = '-';
   } else {
       $strand = '.';
   }
   
   if( $feat->can('seqname') ) {
       $name = $feat->seqname();
       $name ||= 'SEQ';
   } else {
       $name = 'SEQ';
   }


   $str = join("\t",
                 $name,
		 $feat->source_tag(),
		 $feat->primary_tag(),
		 $feat->start(),
		 $feat->end(),
		 $score,
		 $strand,
		 $frame);

   foreach my $tag ( $feat->all_tags ) {
       foreach my $value ( $feat->each_tag_value($tag) ) {
	   $str .= " $tag=$value";
       }
   }


   return $str;
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

=head2

  Title   : intersection
  Usage   : ($start, $stop, $strand) = $feat->intersection($r)
  Function: gives the range that is contained by both ranges
  Args    : a RangeI to compare this one to
  Returns : nothing if they don't overlap, or the range that they do overlap


=head2 union

  Title   : union
  Usage   : ($start, $stop, $strand) = $feat->union($r);
          : ($start, $stop, $strand) = Bio::RangeI->union(@ranges);
  Function: finds the minimal range that contains all of the ranges
  Args    : a range or list of ranges to find the union of
  Returns : the range containing all of the ranges

=cut

1;

