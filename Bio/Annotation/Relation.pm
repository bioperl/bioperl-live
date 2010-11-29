# $Id: Relation.pm 14708 2008-06-10 00:08:17Z heikki $
#
# BioPerl module for Bio::Annotation::Relation
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by bioperl <bioperl-l@bioperl.org>
#
# Copyright bioperl
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Annotation::Relation - Describe a pairwise relationship (pairwise) with
other BioPerl objects

=head1 SYNOPSIS

   use Bio::Annotation::Relation;
   use Bio::Annotation::Collection;

   my $col = Bio::Annotation::Collection->new();
   my $sv = Bio::Annotation::Relation->new(-type => "paralogy"
                                           -to => "someSeqI");
   $col->add_Annotation('tagname', $sv);

=head1 DESCRIPTION

Object which presents a simple but defined relation between two BioPerl objects
of the same type (designated as 'from' and 'to'). This relationship is given a
defined type (required) as well as an optional confidence type and confidence
score, and may indicate directionality. 

The best analogy is when describing the relations in a graph, where this
represents the typed (and possibly scored) edge between two nodes in that graph.
This implementation currently recognizes directionality between relationships
with the default assumption the relationship is bidirectional (undirected).   

B<NOTE:> This should not be used as a direct replacement for already-defined
classes that are used to describe certain types of graphs within BioPerl such as
L<Bio::Tree::Tree>, but instead should be used to describe additional relations
of note within such a graph, such as orthologs in a phylogenetic tree. Do not be
surprised if this class is replaced or reimplemented with something that is a
bit more ontology friendly or utilizes something like Graph to define such
relationships.

With this current implementation, relations can defined as one of two types:

=over2

=item * From 'self' to another instance

The object that this annotation is stored in is related to another instance. An
example would be to describe the relationship one node in a tree has with
another, from the perspective of the B<node>.

In this case, one would not need to set 'from'; it should be assumed that an
undefined 'from' indicates the relationship is between 'self' (the current node)
and another node in the graph.

=item * Between two instances

The object that this annotation is stored in contains two other instances that
are related to one another in some way. An example would be a annotated tree
object, which contains two nodes that are related to one another. The key
difference is this is from the perspective of the B<tree>, which contains the
nodes of interest,

In this case, one would set 'from' as the instance the relationship is from
(Node A) and 'to' as the instance 'from' has a relationship to (Node B).

=back

=head1 NOTE

Care must be taken to avoid circular references.  In particular, do not
set 'from' in cases such as the first one above with an actual 'self' instance.

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

=head1 AUTHOR  - Mira Han

Email mirhan@indiana.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::Annotation::Relation;
use strict;

# Object preamble - inherits from Bio::Root::Root

use base qw(Bio::Root::Root Bio::AnnotationI);

# we shouldn't have to do this, but there is no canonically defined id method
# for all bioperl objs (yet)
 
our %ID_METHOD = (
    'Bio::SeqFeatureI'  => 'display_name',
    'Bio::Tree::NodeI'  => 'id',
    'Bio::SeqI'         => 'object_id'
);

=head2 new

 Title   : new
 Usage   : my $sv = Bio::Annotation::Relation->new();
 Function: Instantiate a new Relation object
 Returns : Bio::Annotation::Relation object
 Args    : -type     => $type of relation [required]
           -class    => Class/interface of instances this relation ties together
                        [required]
           -from     => $obj relation is from [optional], see notes
           -to       => $obj which 'from' is related to [required]
           -is_directed => is this relationship directed [optional, default = 0]
           -id_method => method called from the class [optional]
           -tagname  => $tag to initialize the tagname [optional]
           -tag_term => ontology term representation of the tag [optional]

=cut

sub new{
    my ($class,@args) = @_;
 
    my $self = $class->SUPER::new(@args);
 
    my ($type, $from, $to, $is_directed, $tag, $term, $relation_class, $method) =
        $self->_rearrange([qw(TYPE FROM TO IS_DIRECTED
                           TAGNAME TAG_TERM CLASS
                           ID_METHOD)], @args);
    
    $is_directed ||= 0;
    $self->is_directed($is_directed);

    $self->throw("Must define a relation class") unless defined $relation_class;
    $self->relation_class($relation_class);

    # set the term first
    defined $term   && $self->tag_term($term);
    defined $type   && $self->type($type);
    defined $from   && $self->from($from);
    defined $to     && $self->to($to);
    defined $tag    && $self->tagname($tag);
    defined $method && $self->id_method($method);
    
    return $self;
}


=head1 AnnotationI implementing functions

=cut

=head2 as_text

 Title   : as_text
 Usage   : my $text = $obj->as_text
 Function: return the string "Value: $v" where $v is the value
 Returns : string
 Args    : none

=cut

sub as_text{
    my ($self) = @_;
    my $method = $self->id_method;
    my $dir = $self->is_directed ? '->' : '<->';
    my $from = $self->from ? $self->from->$method : 'self';
    return "$from $dir ".$self->to->$method."(".$self->type.")";
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
    my $DEFAULT_CB = sub { return $_[0]->as_text };
  
    sub display_text {
        my ($self, $cb) = @_;
        $cb ||= $DEFAULT_CB;
        $self->throw("Callback must be a code reference") if ref $cb ne 'CODE';
        return $cb->($self);
    }
}

=head2 hash_tree

 Title   : hash_tree
 Usage   : my $hashtree = $value->hash_tree
 Function: For supporting the AnnotationI interface just returns the value
           as a hashref with the key 'value' pointing to the value
 Returns : hashrf
 Args    : none


=cut

sub hash_tree{
    my $self = shift;

    my $h = {};
    $h->{'type'} = $self->type;
    $h->{'to'} = $self->to;
    return $h;
}

=head2 tagname

 Title   : tagname
 Usage   : $obj->tagname($newval)
 Function: Get/set the tagname for this annotation value.

           Setting this is optional. If set, it obviates the need to
           provide a tag to AnnotationCollection when adding this
           object.

 Example :
 Returns : value of tagname (a scalar)
 Args    : new value (a scalar, optional)


=cut

sub tagname{
    my $self = shift;

    # check for presence of an ontology term
    if($self->{'_tag_term'}) {
	# keep a copy in case the term is removed later
	$self->{'tagname'} = $_[0] if @_;
	# delegate to the ontology term object
	return $self->tag_term->name(@_);
    }
    return $self->{'tagname'} = shift if @_;
    return $self->{'tagname'};
}


=head1 Specific accessors for Relation

=cut

=head2 type 

 Title   : type 
 Usage   : $obj->type($newval)
 Function: Get/Set the type
 Returns : type of relation
 Args    : newtype (optional)


=cut

sub type{
   my ($self,$type) = @_;

   if( defined $type) {
      $self->{'type'} = $type;
    }
    return $self->{'type'};
}

=head2 from

 Title   : from
 Usage   : $obj->from($newval)
 Function: Get/Set the object which $self is in relation from
 Returns : the object which the relation is from
 Args    : new target object (optional), class type (optional)
 Note    : If this is unset, this should be assumed to point to 'self'.
           When an instance is passed, the class name of the instance is
           used for setting the class type

=cut

sub from {
    my ($self,$from) = @_;
 
    if( defined $from) {
        my $c = $self->relation_class;
        $self->throw("Object type not $c") unless $from->isa($c);
        $self->{'from'} = $from;
    }
    return $self->{'from'};
}

=head2 to

 Title   : to
 Usage   : $obj->to($newval)
 Function: Get/Set the object which $self is in relation to
 Returns : the object which the relation applies to
 Args    : new target object (optional), class type (optional)

=cut

sub to {
    my ($self,$to) = @_;
 
    if( defined $to) {
        my $c = $self->relation_class;
        $self->throw("Object type not $c") unless $to->isa($c);
        $self->{'to'} = $to;
    }
    return $self->{'to'};
}

=head2 is_directed

 Title   : is_directed
 Usage   : $obj->is_directed($newval)
 Function: Boolean, indicates whether this relationship is directed or not
 Returns : Boolean value (0 or 1)
 Args    : Boolean value (Default is 0)


=cut

sub is_directed {
    my ($self,$is_directed) = @_;

    if( defined $is_directed) {
       $self->{'is_directed'} = $is_directed ? 1 : 0;
    }
    return $self->{'is_directed'};
}

=head2 confidence

 Title   : confidence
 Usage   : $self->confidence($newval)
 Function: Gives the confidence value.
 Example :
 Returns : value of confidence
 Args    : newvalue (optional)


=cut

sub confidence{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'confidence'} = $value;
    }
    return $self->{'confidence'};
}

=head2 confidence_type

 Title   : confidence_type
 Usage   : $self->confidence_type($newtype)
 Function: Gives the confidence type.
 Example :
 Returns : type of confidence
 Args    : newtype (optional)

=cut

sub confidence_type{
   my ($self,$type) = @_;
   if( defined $type) {
      $self->{'confidence_type'} = $type;
    }
    return $self->{'confidence_type'};
}

=head2 tag_term

 Title   : tag_term
 Usage   : $obj->tag_term($newval)
 Function: Get/set the L<Bio::Ontology::TermI> object representing
           the tag name.

           This is so you can specifically relate the tag of this
           annotation to an entry in an ontology. You may want to do
           this to associate an identifier with the tag, or a
           particular category, such that you can better match the tag
           against a controlled vocabulary.

           This accessor will return undef if it has never been set
           before in order to allow this annotation to stay
           light-weight if an ontology term representation of the tag
           is not needed. Once it is set to a valid value, tagname()
           will actually delegate to the name() of this term.

 Example :
 Returns : a L<Bio::Ontology::TermI> compliant object, or undef
 Args    : on set, new value (a L<Bio::Ontology::TermI> compliant
           object or undef, optional)

=cut

sub tag_term{
    my $self = shift;
    return $self->{'_tag_term'} = shift if @_;
    return $self->{'_tag_term'};
}

=head2 confidence_type

 Title   : confidence_type
 Usage   : $self->confidence_type($newtype)
 Function: Gives the confidence type.
 Example :
 Returns : type of confidence
 Args    : newtype (optional)

=cut

sub relation_class {
    my $self = shift;
    $self->{relation_class} = shift if @_;
    $self->{relation_class};
}

=head2 id_method

 Title   : id_method
 Usage   : $self->id_method($newtype)
 Function: Gives the id method called on the contained instances when
           generating text output.
 Example :
 Returns : type of confidence
 Args    : newtype (optional)

=cut

sub id_method {
    my $self = shift;
    $self->{id_method} = shift if @_;
    $self->{id_method} ? $self->{id_method} :
        exists $ID_METHOD{$self->relation_class} ?
            $ID_METHOD{$self->relation_class} : 
            'object_id';  
}

1;
