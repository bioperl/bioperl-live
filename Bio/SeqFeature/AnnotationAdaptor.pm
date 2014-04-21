#
# BioPerl module for Bio::SeqFeature::AnnotationAdaptor
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Hilmar Lapp <hlapp at gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

#
# (c) Hilmar Lapp, hlapp at gmx.net, 2002.
# (c) GNF, Genomics Institute of the Novartis Research Foundation, 2002.
#
# You may distribute this module under the same terms as perl itself.
# Refer to the Perl Artistic License (see the license accompanying this
# software package, or see http://www.perl.com/language/misc/Artistic.html)
# for the terms under which you may use, modify, and redistribute this module.
# 
# THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
# MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::AnnotationAdaptor - integrates SeqFeatureIs annotation

=head1 SYNOPSIS

   use Bio::SeqFeature::Generic;
   use Bio::SeqFeature::AnnotationAdaptor;

   # obtain a SeqFeatureI implementing object somehow
   my $feat = Bio::SeqFeature::Generic->new(-start => 10, -end => 20);

   # add tag/value annotation
   $feat->add_tag_value("mytag", "value of tag mytag");
   $feat->add_tag_value("mytag", "another value of tag mytag");

   # Bio::SeqFeature::Generic also provides annotation(), which returns a
   # Bio::AnnotationCollectionI compliant object
   $feat->annotation->add_Annotation("dbxref", $dblink);

   # to integrate tag/value annotation with AnnotationCollectionI
   # annotation, use this adaptor, which also implements 
   # Bio::AnnotationCollectionI
   my $anncoll = Bio::SeqFeature::AnnotationAdaptor->new(-feature => $feat);

   # this will now return tag/value pairs as 
   # Bio::Annotation::SimpleValue objects
   my @anns = $anncoll->get_Annotations("mytag");
   # other added before annotation is available too
   my @dblinks = $anncoll->get_Annotations("dbxref");

   # also supports transparent adding of tag/value pairs in 
   # Bio::AnnotationI flavor
   my $tagval = Bio::Annotation::SimpleValue->new(-value => "some value",
                                                  -tagname => "some tag");
   $anncoll->add_Annotation($tagval);
   # this is now also available from the feature's tag/value system
   my @vals = $feat->get_tag_values("some tag");

=head1 DESCRIPTION

L<Bio::SeqFeatureI> defines light-weight annotation of features
through tag/value pairs. Conversely, L<Bio::AnnotationCollectionI>
together with L<Bio::AnnotationI> defines an annotation bag, which is
better typed, but more heavy-weight because it contains every single
piece of annotation as objects. The frequently used base
implementation of Bio::SeqFeatureI, Bio::SeqFeature::Generic, defines
an additional slot for AnnotationCollectionI-compliant annotation.

This adaptor provides a L<Bio::AnnotationCollectionI> compliant,
unified, and integrated view on the annotation of L<Bio::SeqFeatureI>
objects, including tag/value pairs, and annotation through the
annotation() method, if the object supports it. Code using this
adaptor does not need to worry about the different ways of possibly
annotating a SeqFeatureI object, but can instead assume that it
strictly follows the AnnotationCollectionI scheme. The price to pay is
that retrieving and adding annotation will always use objects instead
of light-weight tag/value pairs.

In other words, this adaptor allows us to keep the best of both
worlds. If you create tens of thousands of feature objects, and your
only annotation is tag/value pairs, you are best off using the
features' native tag/value system. If you create a smaller number of
features, but with rich and typed annotation mixed with tag/value
pairs, this adaptor may be for you. Since its implementation is by
double-composition, you only need to create one instance of the
adaptor. In order to transparently annotate a feature object, set the
feature using the feature() method. Every annotation you add will be
added to the feature object, and hence will not be lost when you set
feature() to the next object.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Hilmar Lapp

Email hlapp at gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


#' Let the code begin...


package Bio::SeqFeature::AnnotationAdaptor;
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Annotation::SimpleValue;

use base qw(Bio::Root::Root Bio::AnnotationCollectionI Bio::AnnotatableI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::SeqFeature::AnnotationAdaptor->new();
 Function: Builds a new Bio::SeqFeature::AnnotationAdaptor object 
 Returns : an instance of Bio::SeqFeature::AnnotationAdaptor
 Args    : Named parameters
            -feature    the Bio::SeqFeatureI implementing object to adapt
                        (mandatory to be passed here, or set via feature()
                        before calling other methods)
            -annotation the Bio::AnnotationCollectionI implementing object
                        for storing richer annotation (this will default to
                        the $feature->annotation() if it supports it)
            -tagvalue_factory the object factory to use for creating tag/value
                        pair representing objects


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  my ($feat,$anncoll,$fact) =
	$self->_rearrange([qw(FEATURE
                          ANNOTATION
                          TAGVALUE_FACTORY)], @args);

  $self->feature($feat) if $feat;
  $self->annotation($anncoll) if $feat;
  $self->tagvalue_object_factory($fact) if $fact;

  return $self;
}

=head2 feature

 Title   : feature
 Usage   : $obj->feature($newval)
 Function: Get/set the feature that this object adapts to an
           AnnotationCollectionI.
 Example : 
 Returns : value of feature (a Bio::SeqFeatureI compliant object)
 Args    : new value (a Bio::SeqFeatureI compliant object, optional)


=cut

sub feature{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'feature'} = $value;
    }
    return $self->{'feature'};
}

=head2 annotation

 Title   : annotation
 Usage   : $obj->annotation($newval)
 Function: Get/set the AnnotationCollectionI implementing object used by
           this adaptor to store additional annotation that cannot be stored
           by the SeqFeatureI itself.

           If requested before having been set, the value will default to the
           annotation object of the feature if it has one.
 Example : 
 Returns : value of annotation (a Bio::AnnotationCollectionI compliant object)
 Args    : new value (a Bio::AnnotationCollectionI compliant object, optional)


=cut

sub annotation{
    my ($self,$value) = @_;

    if( defined $value) {
        $self->{'annotation'} = $value;
    }
    if((! exists($self->{'annotation'})) &&
       $self->feature()->can('annotation')) {
        return $self->feature()->annotation();
    }
    return $self->{'annotation'};
}

=head1 AnnotationCollectionI implementing methods

=cut

=head2 get_all_annotation_keys

 Title   : get_all_annotation_keys
 Usage   : $ac->get_all_annotation_keys()
 Function: gives back a list of annotation keys, which are simple text strings
 Returns : list of strings
 Args    : none

=cut

sub get_all_annotation_keys{
    my ($self) = @_;
    my @keys = ();
    
    # get the tags from the feature object
    if ($self->feature()->can('get_all_tags')) {
        push(@keys, $self->feature()->get_all_tags());
    } else {
        push(@keys, $self->feature()->all_tags());
    }
    # ask the annotation implementation in addition, while avoiding duplicates
    if($self->annotation()) {
	push(@keys,
	     grep { ! $self->feature->has_tag($_); }
	          $self->annotation()->get_all_annotation_keys());
    }
    # done
    return @keys;
}


=head2 get_Annotations

 Title   : get_Annotations
 Usage   : my @annotations = $collection->get_Annotations('key')
 Function: Retrieves all the Bio::AnnotationI objects for a specific key
 Returns : list of Bio::AnnotationI - empty if no objects stored for a key
 Args    : string which is key for annotations

=cut

sub get_Annotations{
    my ($self, @keys) = @_;
    my @anns = ();

    # we need a annotation object factory
    my $fact = $self->tagvalue_object_factory();

    # get all tags if no keys have been provided
    @keys = $self->feature->all_tags() unless @keys;

    # build object for each value for each tag
    foreach my $key (@keys) {
      # protect against keys that aren't tags
      next unless $self->feature->has_tag($key);
      # add each tag/value pair as a SimpleValue object
      foreach my $val ($self->feature()->get_tag_values($key)) {
       my $ann;
       if($fact) {
          $ann = $fact->create_object(-value => $val, -tagname => $key);
       } else {
          $ann = Bio::Annotation::SimpleValue->new(-value => $val,
                                                   -tagname => $key);
       }
       push(@anns, $ann);
      }
    }

    # add what is in the annotation implementation if any
    if($self->annotation()) {
      push(@anns, $self->annotation->get_Annotations(@keys));
    }

    # done
    return @anns;
}

=head2 get_num_of_annotations

 Title   : get_num_of_annotations
 Usage   : my $count = $collection->get_num_of_annotations()
 Function: Returns the count of all annotations stored in this collection 
 Returns : integer
 Args    : none


=cut

sub get_num_of_annotations{
  my ($self) = @_;

  # first, count the number of tags on the feature
  my $num_anns = 0;

  foreach ($self->feature()->all_tags()) {
	$num_anns += scalar( $self->feature()->get_tag_values($_));
  }

  # add from the annotation implementation if any
  if($self->annotation()) {
 	$num_anns += $self->annotation()->get_num_of_annotations();
  }

  # done
  return $num_anns;
}

=head1 Implementation specific functions - to allow adding

=cut

=head2 add_Annotation

 Title   : add_Annotation
 Usage   : $self->add_Annotation('reference',$object);
           $self->add_Annotation($object,'Bio::MyInterface::DiseaseI');
           $self->add_Annotation($object);
           $self->add_Annotation('disease',$object,'Bio::MyInterface::DiseaseI');
 Function: Adds an annotation for a specific key.

           If the key is omitted, the object to be added must provide a value
           via its tagname().

           If the archetype is provided, this and future objects added under
           that tag have to comply with the archetype and will be rejected
           otherwise.

           This implementation will add all Bio::Annotation::SimpleValue
           objects to the adapted features as tag/value pairs. Caveat: this
           may potentially result in information loss if a derived object
           is supplied.

 Returns : none
 Args    : annotation key ('disease', 'dblink', ...)
           object to store (must be Bio::AnnotationI compliant)
           [optional] object archetype to map future storage of object 
                      of these types to

=cut

sub add_Annotation{
    my ($self,$key,$object,$archetype) = @_;
   
    # if there's no key we use the tagname() as key
    if(ref($key) && $key->isa("Bio::AnnotationI") &&
       (! ($object && ref($object)))) {
	$archetype = $object if $object;
	$object = $key;
	$key = $object->tagname();
	$key = $key->name() if $key && ref($key); # OntologyTermI
	$self->throw("Annotation object must have a tagname if key omitted")
	    unless $key;
    }
    
    if( !defined $object ) {
	$self->throw("Must have at least key and object in add_Annotation");
    }
    
    if( ! (ref($object) && $object->isa("Bio::AnnotationI")) ) {
	$self->throw("object must be a Bio::AnnotationI compliant object, otherwise we wont add it!");
    }
    
    # ready to add -- if it's a SimpleValue, we add to the feature's tags,
    # otherwise we'll add to the annotation collection implementation

    if($object->isa("Bio::Annotation::SimpleValue") &&
       $self->feature()->can('add_tag_value')) {
	return $self->feature()->add_tag_value($key, $object->value());
    } else {
	my $anncoll = $self->annotation();
	if(! $anncoll) {
	    $anncoll = Bio::Annotation::Collection->new();
	    $self->annotation($anncoll);
	}
	if($anncoll->can('add_Annotation')) {
	    return $anncoll->add_Annotation($key,$object,$archetype);
	}
	$self->throw("Annotation implementation does not allow adding!");
    }
}

=head2 remove_Annotations

 Title   : remove_Annotations
 Usage   :
 Function: Remove the annotations for the specified key from this
           collection.

           If the key happens to be a tag, then the tag is removed
           from the feature.

 Example :
 Returns : an array Bio::AnnotationI compliant objects which were stored
           under the given key(s)
 Args    : the key(s) (tag name(s), one or more strings) for which to
           remove annotations (optional; if none given, flushes all
           annotations)


=cut

sub remove_Annotations{
    my ($self, @keys) = @_;

    # set to all keys if none are supplied
    @keys = $self->get_all_annotation_keys() unless @keys;
    # collect existing annotation
    my @anns = $self->get_Annotations(@keys);
    # flush
    foreach my $key (@keys) {
	# delete the tag if it is one
	$self->feature->remove_tag($key) if $self->feature->has_tag($key);
	# and delegate to the annotation implementation 
	my $anncoll = $self->annotation();
	if($anncoll && $anncoll->can('remove_Annotations')) {
	    $anncoll->remove_Annotations($key);
	} elsif($anncoll) {
	    $self->warn("Annotation bundle implementation ".ref($anncoll).
			" does not allow remove!");
	}
    }
    return @anns;
}

=head1 Additional methods

=cut

=head2 tagvalue_object_factory

 Title   : tagvalue_object_factory
 Usage   : $obj->tagval_object_factory($newval)
 Function: Get/set the object factory to use for creating objects that
           represent tag/value pairs (e.g.,
           Bio::Annotation::SimpleValue).

           The object to be created is expected to follow
           Bio::Annotation::SimpleValue in terms of supported
           arguments at creation time, and the methods.

 Example : 
 Returns : A Bio::Factory::ObjectFactoryI compliant object
 Args    : new value (a Bio::Factory::ObjectFactoryI compliant object, 
           optional)


=cut

sub tagvalue_object_factory{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'tagval_object_factory'} = $value;
    }
    return $self->{'tagval_object_factory'};
}

1;
