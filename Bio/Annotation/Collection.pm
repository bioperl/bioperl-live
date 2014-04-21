#
# BioPerl module for Bio::Annotation::Collection.pm
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Annotation::Collection - Default Perl implementation of 
AnnotationCollectionI

=head1 SYNOPSIS

   # get an AnnotationCollectionI somehow, eg

   $ac = $seq->annotation();

   foreach $key ( $ac->get_all_annotation_keys() ) {
       @values = $ac->get_Annotations($key);
       foreach $value ( @values ) {
          # value is an Bio::AnnotationI, and defines a "as_text" method
          print "Annotation ",$key," stringified value ",$value->as_text,"\n";

          # also defined hash_tree method, which allows data orientated
          # access into this object
          $hash = $value->hash_tree();
       }
   }

=head1 DESCRIPTION

Bioperl implementation for Bio::AnnotationCollectionI 

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Annotation::Collection;

use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Annotation::TypeManager;
use Bio::Annotation::SimpleValue;


use base qw(Bio::Root::Root Bio::AnnotationCollectionI Bio::AnnotationI);


=head2 new

 Title   : new
 Usage   : $coll = Bio::Annotation::Collection->new()
 Function: Makes a new Annotation::Collection object. 
 Returns : Bio::Annotation::Collection
 Args    : none

=cut

sub new{
   my ($class,@args) = @_;

   my $self = $class->SUPER::new(@args);

   $self->{'_annotation'} = {};
   $self->_typemap(Bio::Annotation::TypeManager->new());

   return $self;
}


=head1 L<Bio::AnnotationCollectionI> implementing methods

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
   return keys %{$self->{'_annotation'}};
}

=head2 get_Annotations

 Title   : get_Annotations
 Usage   : my @annotations = $collection->get_Annotations('key')
 Function: Retrieves all the Bio::AnnotationI objects for one or more
           specific key(s).

           If no key is given, returns all annotation objects.

           The returned objects will have their tagname() attribute set to
           the key under which they were attached, unless the tagname was
           already set.

 Returns : list of Bio::AnnotationI - empty if no objects stored for a key
 Args    : keys (list of strings) for annotations (optional)

=cut

sub get_Annotations{
    my ($self,@keys) = @_;

    my @anns = ();
    @keys = $self->get_all_annotation_keys() unless @keys;
    foreach my $key (@keys) {
      if(exists($self->{'_annotation'}->{$key})) {
        push(@anns,
            map {
            $_->tagname($key) if ! $_->tagname(); $_;
            } @{$self->{'_annotation'}->{$key}});
      }
    }
    return @anns;
}


=head2 get_nested_Annotations

 Title   : get_nested_Annotations
 Usage   : my @annotations = $collection->get_nested_Annotations(
                                '-key' => \@keys,
                                '-recursive => 1);
 Function: Retrieves all the Bio::AnnotationI objects for one or more
           specific key(s). If -recursive is set to true, traverses the nested 
           annotation collections recursively and returns all annotations 
           matching the key(s).

           If no key is given, returns all annotation objects.

           The returned objects will have their tagname() attribute set to
           the key under which they were attached, unless the tagname was
           already set.

 Returns : list of Bio::AnnotationI - empty if no objects stored for a key
 Args    : -keys      => arrayref of keys to search for (optional)
           -recursive => boolean, whether or not to recursively traverse the 
            nested annotations and return annotations with matching keys.

=cut

sub get_nested_Annotations {
  my ($self, @args) = @_;
  my ($keys, $recursive) = $self->_rearrange([qw(KEYS RECURSIVE)], @args);
  $self->verbose(1);
  
  my @anns = ();
  # if not recursive behave exactly like get_Annotations()
  if (!$recursive) {
	  my @keys = $keys? @$keys : $self->get_all_annotation_keys();
    foreach my $key (@keys) {
      if(exists($self->{'_annotation'}->{$key})) {
        push(@anns,
            map {
            $_->tagname($key) if ! $_->tagname(); $_;
            } @{$self->{'_annotation'}->{$key}});
      }
    }
  }
  # if recursive search for keys recursively
  else {
    my @allkeys = $self->get_all_annotation_keys();
    foreach my $key (@allkeys) {
      my $keymatch = 0;
      foreach my $searchkey (@$keys) {
        if ($key eq $searchkey) { $keymatch = 1;}
      }
      if ($keymatch) {
        if(exists($self->{'_annotation'}->{$key})) {
          push(@anns,
              map {
              $_->tagname($key) if ! $_->tagname(); $_;
              } @{$self->{'_annotation'}->{$key}});
        }
      }
      else {
        my @annotations = @{$self->{'_annotation'}->{$key}};
        foreach (@annotations) {
          if ($_->isa("Bio::AnnotationCollectionI")) {
            push (@anns, 
                  $_->get_nested_Annotations('-keys' => $keys, '-recursive' => 1)
                 );
          }
        }
      }
    }
  }
  return @anns;
}

=head2 get_all_Annotations

 Title   : get_all_Annotations
 Usage   :
 Function: Similar to get_Annotations, but traverses and flattens nested
           annotation collections. This means that collections in the
           tree will be replaced by their components.

           Keys will not be passed on to nested collections. I.e., if the
           tag name of a nested collection matches the key, it will be
           flattened in its entirety.

           Hence, for un-nested annotation collections this will be identical
           to get_Annotations.
 Example :
 Returns : an array of L<Bio::AnnotationI> compliant objects
 Args    : keys (list of strings) for annotations (optional)


=cut

sub get_all_Annotations{
    my ($self,@keys) = @_;

    return map {
	$_->isa("Bio::AnnotationCollectionI") ?
	    $_->get_all_Annotations() : $_;
    } $self->get_Annotations(@keys);
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
   my $count = 0;
   map { $count += scalar @$_ } values %{$self->{'_annotation'}};
   return $count;
}

=head1 Implementation specific functions - mainly for adding

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

 Returns : none
 Args    : annotation key ('disease', 'dblink', ...)
           object to store (must be Bio::AnnotationI compliant)
           [optional] object archetype to map future storage of object 
                      of these types to

=cut

sub add_Annotation{
   my ($self,$key,$object,$archetype) = @_;
   
   # if there's no key we use the tagname() as key
   if(ref($key) && $key->isa("Bio::AnnotationI") && (!ref($object))) {
       $archetype = $object if defined($object);
       $object = $key;
       $key = $object->tagname();
       $key = $key->name() if ref($key); # OntologyTermI
       $self->throw("Annotation object must have a tagname if key omitted")
	   unless $key;
   }

   if( !defined $object ) {
       $self->throw("Must have at least key and object in add_Annotation");
   }

   if( !ref $object ) {
       $self->throw("Must add an object. Use Bio::Annotation::{Comment,SimpleValue,OntologyTerm} for simple text additions");
   }

   if( !$object->isa("Bio::AnnotationI") ) {
       $self->throw("object must be AnnotationI compliant, otherwise we won't add it!");
   }

   # ok, now we are ready! If we don't have an archetype, set it
   # from the type of the object

   if( !defined $archetype ) {
       $archetype = ref $object;
   }

   # check typemap, storing if needed.
   my $stored_map = $self->_typemap->type_for_key($key);

   if( defined $stored_map ) {
       # check validity, irregardless of archetype. A little cheeky
       # this means isa stuff is executed correctly

       if( !$self->_typemap()->is_valid($key,$object) ) {
	   $self->throw("Object $object was not valid with key $key. ".
         "If you were adding new keys in, perhaps you want to make use\n".
         "of the archetype method to allow registration to a more basic type");
       }
   } else {
       $self->_typemap->_add_type_map($key,$archetype);
   }

   # we are ok to store

   if( !defined $self->{'_annotation'}->{$key} ) {
       $self->{'_annotation'}->{$key} = [];
   }

   push(@{$self->{'_annotation'}->{$key}},$object);

   return 1;
}

=head2 remove_Annotations

 Title   : remove_Annotations
 Usage   :
 Function: Remove the annotations for the specified key from this collection.
 Example :
 Returns : an array Bio::AnnotationI compliant objects which were stored
           under the given key(s)
 Args    : the key(s) (tag name(s), one or more strings) for which to
           remove annotations (optional; if none given, flushes all
           annotations)


=cut

sub remove_Annotations{
    my ($self, @keys) = @_;

    @keys = $self->get_all_annotation_keys() unless @keys;
    my @anns = $self->get_Annotations(@keys);
    # flush
    foreach my $key (@keys) {
      delete $self->{'_annotation'}->{$key};
      delete $self->{'_typemap'}->{'_type'}->{$key};
    }
    return @anns;
}

=head2 flatten_Annotations

 Title   : flatten_Annotations
 Usage   :
 Function: Flattens part or all of the annotations in this collection.

           This is a convenience method for getting the flattened
           annotation for the given keys, removing the annotation for
           those keys, and adding back the flattened array.

           This should not change anything for un-nested collections.
 Example :
 Returns : an array Bio::AnnotationI compliant objects which were stored
           under the given key(s)
 Args    : list of keys (strings) the annotation for which to flatten,
           defaults to all keys if not given


=cut

sub flatten_Annotations{
    my ($self,@keys) = @_;

    my @anns = $self->get_all_Annotations(@keys);
    my @origanns = $self->remove_Annotations(@keys);
    foreach (@anns) {
	$self->add_Annotation($_);
    }
    return @origanns;
}

=head1 Bio::AnnotationI methods implementations

   This is to allow nested annotation: you can use a collection as an
   annotation object for an annotation collection.

=cut

=head2 as_text

 Title   : as_text
 Usage   :
 Function: See L<Bio::AnnotationI>
 Example :
 Returns : a string
 Args    : none


=cut

sub as_text{
    my $self = shift;

    my $txt = "Collection consisting of ";
    my @texts = ();
    foreach my $ann ($self->get_Annotations()) {
	push(@texts, $ann->as_text());
    }
    if(@texts) {
	$txt .= join(", ", map { '['.$_.']'; } @texts);
    } else {
	$txt .= "no elements";
    }
    return $txt;
}

=head2 display_text

 Title   : display_text
 Usage   : my $str = $ann->display_text();
 Function: returns a string. Unlike as_text(), this method returns a string
           formatted as would be expected for te specific implementation.

           One can pass a callback as an argument which allows custom text
           generation; the callback is passed the current instance and any text
           returned
 Example :
 Returns : a string
 Args    : [optional] callback

=cut

{
   # this just calls the default display_text output for
   # any AnnotationI
  my $DEFAULT_CB = sub {
    my $obj = shift;
    my $txt;
    foreach my $ann ($obj->get_Annotations()) {
      $txt .= $ann->display_text()."\n";
    }
    return $txt;
    };

  sub display_text {
    my ($self, $cb) = @_;
    $cb ||= $DEFAULT_CB;
    $self->throw("") if ref $cb ne 'CODE';
    return $cb->($self);
  }
}


=head2 hash_tree

 Title   : hash_tree
 Usage   :
 Function: See L<Bio::AnnotationI>
 Example :
 Returns : a hash reference
 Args    : none


=cut

sub hash_tree{
    my $self = shift;
    my $tree = {};

    foreach my $key ($self->get_all_annotation_keys()) {
	# all contained objects will support hash_tree() 
	# (they are AnnotationIs)
	$tree->{$key} = [$self->get_Annotations($key)];
    }
    return $tree;
}

=head2 tagname

 Title   : tagname
 Usage   : $obj->tagname($newval)
 Function: Get/set the tagname for this annotation value.

           Setting this is optional. If set, it obviates the need to
           provide a tag to Bio::AnnotationCollectionI when adding
           this object. When obtaining an AnnotationI object from the
           collection, the collection will set the value to the tag
           under which it was stored unless the object has a tag
           stored already.

 Example : 
 Returns : value of tagname (a scalar)
 Args    : new value (a scalar, optional)


=cut

sub tagname{
    my $self = shift;

    return $self->{'tagname'} = shift if @_;
    return $self->{'tagname'};
}


=head1 Backward compatible functions

Functions put in for backward compatibility with old
Bio::Annotation.pm stuff

=cut

=head2 description

 Title   : description
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub description{
   my ($self,$value) = @_;

   $self->deprecated("Using old style annotation call on new Annotation::Collection object");

   if( defined $value ) {
       my $val = Bio::Annotation::SimpleValue->new();
       $val->value($value);
       $self->add_Annotation('description',$val);
   }

   my ($desc) = $self->get_Annotations('description');
   
   # If no description tag exists, do not attempt to call value on undef:
   return $desc ? $desc->value : undef;
}


=head2 add_gene_name

 Title   : add_gene_name
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_gene_name{
   my ($self,$value) = @_;

   $self->deprecated("Old style add_gene_name called on new style Annotation::Collection");

   my $val = Bio::Annotation::SimpleValue->new();
   $val->value($value);
   $self->add_Annotation('gene_name',$val);
}

=head2 each_gene_name

 Title   : each_gene_name
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_gene_name{
   my ($self) = @_;

   $self->deprecated("Old style each_gene_name called on new style Annotation::Collection");

   my @out;
   my @gene = $self->get_Annotations('gene_name');

   foreach my $g ( @gene ) {
       push(@out,$g->value);
   }

   return @out;
}

=head2 add_Reference

 Title   : add_Reference
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Reference{
   my ($self, @values) = @_;

   $self->deprecated("add_Reference (old style Annotation) on new style Annotation::Collection");
   
   # Allow multiple (or no) references to be passed, as per old method
   foreach my $value (@values) {
       $self->add_Annotation('reference',$value);
   }
}

=head2 each_Reference

 Title   : each_Reference
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_Reference{
   my ($self) = @_;

   $self->deprecated("each_Reference (old style Annotation) on new style Annotation::Collection");
   
   return $self->get_Annotations('reference');
}


=head2 add_Comment

 Title   : add_Comment
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Comment{
   my ($self,$value) = @_;

   $self->deprecated("add_Comment (old style Annotation) on new style Annotation::Collection");

   $self->add_Annotation('comment',$value);

}

=head2 each_Comment

 Title   : each_Comment
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_Comment{
   my ($self) = @_;

   $self->deprecated("each_Comment (old style Annotation) on new style Annotation::Collection");
   
   return $self->get_Annotations('comment');
}



=head2 add_DBLink

 Title   : add_DBLink
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_DBLink{
   my ($self,$value) = @_;

   $self->deprecated("add_DBLink (old style Annotation) on new style Annotation::Collection");

   $self->add_Annotation('dblink',$value);

}

=head2 each_DBLink

 Title   : each_DBLink
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_DBLink{
   my ($self) = @_;

   $self->deprecated("each_DBLink (old style Annotation) on new style Annotation::Collection - use get_Annotations('dblink')");
   
   return $self->get_Annotations('dblink');
}



=head1 Implementation management functions

=cut

=head2 _typemap

 Title   : _typemap
 Usage   : $obj->_typemap($newval)
 Function: 
 Example : 
 Returns : value of _typemap
 Args    : newvalue (optional)


=cut

sub _typemap{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_typemap'} = $value;
    }
    return $self->{'_typemap'};

}

1;
