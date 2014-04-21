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

Bio::Annotation::Relation - Relationship (pairwise) with other objects SeqI and NodeI;

=head1 SYNOPSIS

   use Bio::Annotation::Relation;
   use Bio::Annotation::Collection;

   my $col = Bio::Annotation::Collection->new();
   my $sv = Bio::Annotation::Relation->new(-type => "paralogy" -to => "someSeqI");
   $col->add_Annotation('tagname', $sv);

=head1 DESCRIPTION

Scalar value annotation object

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

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Annotation::Relation;
use strict;

# Object preamble - inherits from Bio::Root::Root

use base qw(Bio::Root::Root Bio::AnnotationI);

=head2 new

 Title   : new
 Usage   : my $sv = Bio::Annotation::Relation->new();
 Function: Instantiate a new Relation object
 Returns : Bio::Annotation::Relation object
 Args    : -type    => $type of relation [optional]
           -to     => $obj which $self is in relation to [optional]
           -tagname  => $tag to initialize the tagname [optional]
           -tag_term => ontology term representation of the tag [optional]

=cut

sub new{
   my ($class,@args) = @_;

   my $self = $class->SUPER::new(@args);

   my ($type, $to, $tag, $term) =
       $self->_rearrange([qw(TYPE TO TAGNAME TAG_TERM)], @args);

   # set the term first
   defined $term   && $self->tag_term($term);
   defined $type && $self->type($type);
   defined $to  && $self->to($to);
   defined $tag    && $self->tagname($tag);

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

   return $self->type." to  ".$self->to->id;
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
  my $DEFAULT_CB = sub { return $_[0]->type." to  ".$_[0]->to->id };
  #my $DEFAULT_CB = sub { $_[0]->value};

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

=head2 to

 Title   : to
 Usage   : $obj->to($newval)
 Function: Get/Set the object which $self is in relation to
 Returns : the object which the relation applies to
 Args    : new target object (optional)


=cut

sub to{
   my ($self,$to) = @_;

   if( defined $to) {
      $self->{'to'} = $to;
    }
    return $self->{'to'};
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

1;
