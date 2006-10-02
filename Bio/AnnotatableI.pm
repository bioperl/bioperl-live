# $Id$
#
# BioPerl module for Bio::AnnotatableI
#
# Cared for by Hilmar Lapp <hlapp at gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::AnnotatableI - the base interface an annotatable object must implement

=head1 SYNOPSIS

    use Bio::SeqIO;
    # get an annotatable object somehow: for example, Bio::SeqI objects
    # are annotatable
    my $seqio = Bio::SeqIO->new(-fh => \*STDIN, -format => 'genbank');
    while (my $seq = $seqio->next_seq()) {
        # $seq is-a Bio::AnnotatableI, hence:
        my $ann_coll = $seq->annotation();
        # $ann_coll is-a Bio::AnnotationCollectionI, hence:
        my @all_anns = $ann_coll->get_Annotations();
        # do something with the annotation objects
    }

=head1 DESCRIPTION

This is the base interface that all annotatable objects must implement. A 
good example is Bio::Seq which is an AnnotableI object.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR

 Hilmar Lapp E<lt>hlapp@gmx.netE<gt>
 Allen Day E<lt>allenday@ucla.eduE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::AnnotatableI;
use strict;
use Carp;

use Bio::Annotation::Comment;
use Bio::Annotation::DBLink;
#use Bio::Annotation::OntologyTerm;
use Bio::Annotation::Reference;
use Bio::Annotation::SimpleValue;

our %tagclass = (
  comment        => 'Bio::Annotation::Comment',
  dblink         => 'Bio::Annotation::DBLink',
  description    => 'Bio::Annotation::SimpleValue',
  gene_name      => 'Bio::Annotation::SimpleValue',
  ontology_term  => 'Bio::Annotation::OntologyTerm',
  reference      => 'Bio::Annotation::Reference',
  __DEFAULT__    => 'Bio::Annotation::SimpleValue',
);

our %tag2text = (
  'Bio::Annotation::Comment'        => 'text',
  'Bio::Annotation::DBLink'         => 'primary_id',
  'Bio::Annotation::SimpleValue'    => 'value',
  'Bio::Annotation::SimpleValue'    => 'value',
  'Bio::Annotation::OntologyTerm'   => 'name',
  'Bio::Annotation::Reference'      => 'title',
  __DEFAULT__                       => 'value',

);

use base qw(Bio::Root::RootI);

=head2 annotation

 Title   : annotation
 Usage   : $obj->annotation($newval)
 Function: Get the annotation collection for this annotatable object.
 Example : 
 Returns : a Bio::AnnotationCollectionI implementing object, or undef
 Args    : on set, new value (a Bio::AnnotationCollectionI
           implementing object, optional) (an implementation may not
           support changing the annotation collection)

See L<Bio::AnnotationCollectionI>

=cut

sub annotation{
  shift->throw_not_implemented();
}


=head1 "*_tag_*" METHODS

The methods below allow mapping of the old "get_tag_values()"-style
annotation access to Bio::AnnotationCollectionI.  These need not be
implemented in a Bio::AnnotationCollectionI compliant class, as they
are built on top of the methods.

B<DEPRECATED>: DO NOT USE THESE FOR FUTURE DEVELOPMENT.

=cut

=head2 has_tag

 Usage   : $count = $obj->has_tag($tag)
 Function: returns the number of annotations corresponding to $tag
 Returns : an integer
 Args    : tag name
 Note    : DEPRECATED

Use L</get_Annotations> instead.

=cut

sub has_tag {
  my ($self,$tag) = @_;
  #uncomment in 1.6
  #$self->deprecated('has_tag() is deprecated, use get_Annotations()');

  return scalar($self->annotation->get_Annotations($tag));
}

=head2 add_tag_value

 Usage   : See add_Annotation
 Function:
 Returns : 
 Args    : DEPRECATED

See L<Bio::AnnotationCollectionI::add_Annotation>

=cut

sub add_tag_value {
  my ($self,$tag,@vals) = @_;

  #uncomment in 1.6
  #$self->deprecated('add_tag_value() is deprecated, use add_Annotation()');

  foreach my $val (@vals){
    my $class = $tagclass{$tag}   || $tagclass{__DEFAULT__};
    my $slot  = $tag2text{$class};

    my $a = $class->new();
    $a->$slot($val);

    $self->annotation->add_Annotation($tag,$a);
  }

  return 1;
  #return $self->annotation->add_Annotation(@args);
}


=head2 get_Annotations

 Usage   : my $parent   = $obj->get_Annotations('Parent');
           my @parents = $obj->get_Annotations('Parent');
 Function: a wrapper around Bio::Annotation::Collection::get_Annotations().
 Returns : returns annotations as
           Bio::Annotation::Collection::get_Annotations() does, but
           additionally returns a single scalar in scalar context
           instead of list context so that if an annotation tag
           contains only a single value, you can do:

           $parent = $feature->get_Annotations('Parent');

           instead of:

           ($parent) = ($feature->get_Annotations('Parent'))[0];

           if the 'Parent' tag has multiple values and is called in a
           scalar context, the number of annotations is returned.

 Args    : an annotation tag name.

=cut

sub get_Annotations {
    my $self = shift;

    my @annotations = $self->annotation->get_Annotations(@_);

    if(wantarray){
        return @annotations;
    } elsif(scalar(@annotations) == 1){
        return $annotations[0];
    } else {
        return scalar(@annotations);
    }
}

=head2 get_tag_values

 Usage   : @annotations = $obj->get_tag_values($tag)
 Function: returns annotations corresponding to $tag
 Returns : a list of scalars
 Args    : tag name
 Note    : DEPRECATED

This method is essentially L</get_Annotations>, use it instead.

=cut

sub get_tag_values {
    my ($self,$tag) = @_;
    
    #uncomment in 1.6
    #$self->deprecated('get_tag_values() is deprecatedk, use get_Annotations()');

    if(!$tagclass{$tag} && $self->annotation->get_Annotations($tag)){
        #new tag, haven't seen it yet but it exists.  add to registry
        my($proto) = $self->annotation->get_Annotations($tag);
        # we can only register if there's a method known for obtaining the value
        if (exists($tag2text{ref($proto)})) {
            $tagclass{$tag} = ref($proto);
        }
    }

    my $slot  = $tag2text{ $tagclass{$tag} || $tagclass{__DEFAULT__} };
    
    return map { $_->$slot } $self->annotation->get_Annotations($tag);
}

=head2 get_tagset_values

 Usage   : @annotations = $obj->get_tagset_values($tag1,$tag2)
 Function: returns annotations corresponding to a list of tags.
           this is a convenience method equivalent to multiple calls
           to get_tag_values with each tag in the list.
 Returns : a list of Bio::AnnotationI objects.
 Args    : a list of tag names
 Note    : DEPRECATED

See L<Bio::AnnotationCollectionI::get_Annotations>

=cut

sub get_tagset_values {
  my ($self,@tags) = @_;

  #uncomment in 1.6
  #$self->deprecated('get_tagset_values() is deprecated, use get_Annotations()');

  my @r = ();
  foreach my $tag (@tags){
    my $slot  = $tag2text{ $tagclass{$tag} || $tagclass{__DEFAULT__} };
    push @r, map { $_->$slot } $self->annotation->get_Annotations($tag);
  }
  return @r;
}

=head2 get_all_tags

 Usage   : @tags = $obj->get_all_tags()
 Function: returns a list of annotation tag names.
 Returns : a list of tag names
 Args    : none
 Note    : DEPRECATED

See L<Bio::AnnotationCollectionI::get_all_annotation_keys>

=cut

sub get_all_tags {
  my ($self,@args) = @_;

  #uncomment in 1.6
  #$self->deprecated('get_all_tags() is deprecated, use get_all_annotation_keys()');

  return $self->annotation->get_all_annotation_keys(@args);
}

=head2 remove_tag

 Usage   : See remove_Annotations().
 Function:
 Returns : 
 Args    : DEPRECATED
 Note    : Contrary to what the name suggests, this method removes
           all annotations corresponding to $tag, not just a
           single anntoation.

See L<Bio::AnnotationCollectionI::remove_Annotations>

=cut

sub remove_tag {
  my ($self,@args) = @_;

  #uncomment in 1.6
  #$self->deprecated('remove_tag() is deprecated, use remove_Annotations()');

  return $self->annotation->remove_Annotations(@args);
}

1;
