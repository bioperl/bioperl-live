
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
            print "Feature from ", $feat->start, "to ", $feat->end, " Primary tag  " $feat->primary_tag, 
	    " From", $feat->source_tag() "\n";

            if( $feat->strand == 0 ) {
		print "Feature applicable to either strand\n";
            } else {
                print "Feature on strand ", $feat->strand,"\n"; # -1,1
            }

            foreach $tag ( $feat->all_tags() ) {
		print "Feature has tag ",$tag,"with value," $feat->has_tag($tag), "\n";
            }
	}

=head1 DESCRIPTION

This interface is the functions one can expect for any Sequence Feature, whatever
its implemtation or whether it is a more complex type (eg, a Gene). This object
doesn't actually provide any implemention, it just provides the definitions
of what methods one can call. See Bio::SeqFeature::Generic for a good standard
implementation of this object

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeatureI;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::RangeI;
use Carp;

@ISA = qw(Bio::RangeI Exporter);

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
 Usage   : $value = $self->has_tag('some_tag')
 Function: Returns the value of the tag (undef if 
           none)
 Returns : 
 Args    :


=cut

sub has_tag{
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




