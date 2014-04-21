# $Id: TagTree.pm 11693 2007-09-17 20:54:04Z cjfields $
#
# BioPerl module for Bio::Annotation::TagTree
#
# Cared for Chris Fields
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

Bio::Annotation::TagTree - AnnotationI with tree-like hierarchal key-value
relationships ('structured tags') that can be represented as simple text.

=head1 SYNOPSIS

   use Bio::Annotation::TagTree;
   use Bio::Annotation::Collection;

   my $col = Bio::Annotation::Collection->new();

   # data structure can be an array reference with a data structure
   # corresponding to that defined by Data::Stag:

   my $sv = Bio::Annotation::TagTree->new(-tagname => 'mytag1',
                                          -value => $data_structure);
   $col->add_Annotation($sv);

   # regular text passed is parsed based on the tagformat().
   my $sv2 = Bio::Annotation::TagTree->new(-tagname => 'mytag2',
                                          -tagformat => 'xml',
                                          -value => $xmltext);
   $col->add_Annotation($sv2);

=head1 DESCRIPTION

This takes tagged data values and stores them in a hierarchal structured
element-value hierarchy (complements of Chris Mungall's Data::Stag module). Data
can then be represented as text using a variety of output formats (indention,
itext, xml, spxr). Furthermore, the data structure can be queried using various
means. See L<Data::Stag> for details.

Data passed in using value() or the '-value' parameter upon instantiation
can either be:

1) an array reference corresponding to the data structure for Data::Stag;

2) a text string in 'xml', 'itext', 'spxr', or 'indent' format. The default
format is 'xml'; this can be changed using tagformat() prior to using value() or
by passing in the proper format using '-tagformat' upon instantiation;

3) another Bio::Annotation::TagTree or Data::Stag node instance.  In both cases
a deep copy (duplicate) of the instance is generated.

Beyond checking for an array reference no format guessing occurs (so, for
roundtrip tests ensure that the IO formats correspond). For now, we recommend
when using text input to set tagformat() to one of these formats prior to data
loading to ensure the proper Data::Stag parser is selected. After data loading,
the tagformat() can be changed to change the text string format returned by
value(). (this may be rectified in the future)

This Annotation type is fully BioSQL compatible and could be considered a
temporary replacement for nested Bio::Annotation::Collections, at least until
BioSQL and bioperl-db can support nested annotation collections.

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
or the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR 

Chris Fields

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Annotation::TagTree;
use strict;

# Object preamble - inherits from Bio::Root::Root

use base qw(Bio::Annotation::SimpleValue);
use Data::Stag;

=head2 new

 Title   : new
 Usage   : my $sv = Bio::Annotation::TagTree->new();
 Function: Instantiate a new TagTree object
 Returns : Bio::Annotation::TagTree object
 Args    : -value => $value to initialize the object data field [optional]
           -tagname => $tag to initialize the tagname [optional]
           -tagformat => format for output [optional]
                      (types 'xml', 'itext', 'sxpr', 'indent', default = 'itext')
           -node => Data::Stag node or Bio::Annotation::TagTree instance

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new();
    my ( $node, $value, $tag, $format, $verbose ) = $self->_rearrange(
        [
            qw(
              NODE
              VALUE
              TAGNAME
              TAGFORMAT
              VERBOSE)
        ],
        @args
    );
    $self->throw("Cant use both node and value; mutually exclusive")
      if defined $node && defined $value;
    defined $tag && $self->tagname($tag);
    $format ||= 'itext';
    $self->tagformat($format);
    defined $value   && $self->value($value);
    defined $node    && $self->node($node);
    defined $verbose && $self->verbose($verbose);
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

sub as_text {
    my ($self) = @_;
    return "TagTree: " . $self->value;
}

=head2 display_text

 Title   : display_text
 Usage   : my $str = $ann->display_text();
 Function: returns a string. Unlike as_text(), this method returns a string
           formatted as would be expected for the specific implementation.

           One can pass a callback as an argument which allows custom text
           generation; the callback is passed the current instance and any text
           returned
 Example :
 Returns : a string
 Args    : [optional] callback

=cut

{
    my $DEFAULT_CB = sub { $_[0]->value || '' };

    sub display_text {
        my ( $self, $cb ) = @_;
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
           Maybe reimplement using Data::Stag::hash()?
 Returns : hashrf
 Args    : none

=cut

sub hash_tree {
    my ($self) = @_;
    my $h = {};
    $h->{'value'} = $self->value;
}

=head2 tagname

 Title   : tagname
 Usage   : $obj->tagname($newval)
 Function: Get/set the tagname for this annotation value.

           Setting this is optional. If set, it obviates the need to provide
           a tag to AnnotationCollection when adding this object.
 Example :
 Returns : value of tagname (a scalar)
 Args    : new value (a scalar, optional)

=cut

sub tagname {
    my ( $self, $value ) = @_;
    if ( defined $value ) {
        $self->{'tagname'} = $value;
    }
    return $self->{'tagname'};
}

=head1 Specific accessors for TagTree

=cut

=head2 value

 Title   : value
 Usage   : $obj->value($newval)
 Function: Get/set the value for this annotation.
 Returns : value of value
 Args    : newvalue (optional)

=cut

sub value {
    my ( $self, $value ) = @_;

    # set mode? This resets the entire tagged database
    my $format = $self->tagformat;
    if ($value) {
        if ( ref $value ) {
            if ( ref $value eq 'ARRAY' ) {

                # note the tagname() is not used here; it is only used for
                # storing this AnnotationI in the annotation collection
                eval { $self->{db} = Data::Stag->nodify($value) };
            }
            else {

                # assuming this is blessed; passing on to node() and copy
                $self->node( $value, 'copy' );
            }
        }
        else {

            # not trying to guess here for now; we go by the tagformat() setting
            my $h = Data::Stag->getformathandler($format);
            eval { $self->{db} = Data::Stag->from( $format . 'str', $value ) };
        }
        $self->throw("Data::Stag error:\n$@") if $@;
    }

    # get mode?
    # How do we return a data structure?
    # for now, we use the output (if there is a Data::Stag node present)
    # may need to run an eval {} to catch Data::Stag output errors
    $self->node->$format;
}

=head2 tagformat

 Title   : tagformat
 Usage   : $obj->tagformat($newval)
 Function: Get/set the output tag format for this annotation.
 Returns : value of tagformat
 Args    : newvalue (optional) - format for the data passed into value
           must be of values 'xml', 'indent', 'sxpr', 'itext', 'perl'

=cut

my %IS_VALID_FORMAT = map { $_ => 1 } qw(xml indent sxpr itext);

sub tagformat {
    my ( $self, $value ) = @_;
    if ( defined $value ) {
        $self->throw( "$value is not a valid format; valid format types:\n"
              . join( ',', map { "'$_'" } keys %IS_VALID_FORMAT ) )
          if !exists $IS_VALID_FORMAT{$value};
        $self->{'tagformat'} = $value;
    }
    return $self->{'tagformat'};
}

=head2 node

 Title   : node
 Usage   : $obj->node()
 Function: Get/set the topmost Data::Stag node used for this annotation.  
 Returns : Data::Stag node implementation
           (default is Data::Stag::StagImpl)
 Args    : (optional) Data::Stag node implementation
           (optional)'copy' => flag to create a copy of the node

=cut

sub node {
    my ( $self, $value, $copy ) = @_;
    if ( defined $value && ref $value ) {
        $self->{'db'} =
          $value->isa('Data::Stag::StagI')
          ? ( $copy && $copy eq 'copy' ? $value->duplicate : $value )
          : $value->isa('Bio::Annotation::TagTree') ? ( $copy
              && $copy eq 'copy' ? $value->node->duplicate : $value->node )
          : $self->throw(
            'Object must be Data::Stag::StagI or Bio::Annotation::TagTree');
    }
    
    # lazily create Data::Stag instance if not present
    if (!$self->{'db'}) {
        $self->{'db'} = Data::Stag->new();
    }
    return $self->{'db'};
}

=head2 Data::Stag convenience methods

Because Data::Stag uses blessed arrays and the core Bioperl class uses blessed
hashes, TagTree uses an internal instance of a Data::Stag node for data storage.
Therefore the following methods actually delegate to the Data:::Stag internal
instance.

For consistency (since one could recursively check child nodes), methods retain
the same names as Data::Stag. Also, no 'magic' (AUTOLOAD'ed) methods are
employed, simply b/c full-fledged Data::Stag functionality can be attained by
grabbing the Data::Stag instance using node().

=head2 element

 Title   : element
 Usage   :
 Function: Returns the element name (key name) for this node
 Example :
 Returns : scalar
 Args    : none

=cut

sub element {
    my $self = shift;
    return $self->node->element;
}

=head2 data

 Title   : data
 Usage   :
 Function: Returns the data structure (array ref) for this node
 Example :
 Returns : array ref
 Args    : none

=cut

sub data {
    my $self = shift;
    return $self->node->data;
}

=head2 children

 Title   : children
 Usage   :
 Function: Get the top-level array of Data::Stag nodes or (if the top level is
           a terminal node) a scalar value.

           This is similar to StructuredValue's get_values() method, with the
           key difference being instead of array refs and scalars you get either
           Data::Stag nodes or the value for this particular node.

           For consistency (since one could recursively check nodes),
           we use the same method name as Data::Stag children().
 Example :
 Returns : an array
 Args    : none

=cut

sub children {
    my $self = shift;
    return $self->node->children;
}

=head2 subnodes

 Title   : subnodes
 Usage   :
 Function: Get the top-level array of Data::Stag nodes.  Unlike children(),
           this only returns an array of nodes (if this is a terminal node,
           no value is returned)
 Example :
 Returns : an array of nodes
 Args    : none

=cut

sub subnodes {
    my $self = shift;
    return $self->node->subnodes;
}

=head2 get

 Title   : get
 Usage   : 
 Function: Returns the nodes or value for the named element or path
 Example : 
 Returns : returns array of nodes or a scalar (if node is terminal)
           dependent on wantarray
 Args    : none

=cut

sub get {
    my ( $self, @vals ) = @_;
    return $self->node->get(@vals);
}

=head2 find

 Title   : find
 Usage   : 
 Function: Recursively searches for and returns the nodes or values for the
           named element or path
 Example : 
 Returns : returns array of nodes or scalars (for terminal nodes)
 Args    : none

=cut

sub find {
    my ( $self, @vals ) = @_;
    return $self->node->find(@vals);
}

=head2 findnode

 Title   : findnode
 Usage   : 
 Function: Recursively searches for and returns a list of nodes
           of the given element path
 Example : 
 Returns : returns array of nodes
 Args    : none

=cut

sub findnode {
    my ( $self, @vals ) = @_;
    return $self->node->findnode(@vals);
}

=head2 findval

 Title   : findval
 Usage   : 
 Function: 
 Example : 
 Returns : returns array of nodes or values
 Args    : none

=cut

sub findval {
    my ( $self, @vals ) = @_;
    return $self->node->findval(@vals);
}

=head2 addchild

 Title   : addchild
 Usage   : $struct->addchild(['name' => [['foo'=> 'bar1']]]);
 Function: add new child node to the current node.  One can pass in a node, TagTree,
           or data structure; for instance, in the above, this would translate
           to (in XML):

           <name>
             <foo>bar1</foo>
           </name>

 Returns : node
 Args    : first arg = element name
           all other args are added as tag-value pairs

=cut

sub addchild {
    my ( $self, @vals ) = @_;

    # check for element tag first (if no element, must be empty Data::Stag node)
    if ( !$self->element ) {

        # try to do the right thing; if more than one element, wrap in array ref
        @vals > 1 ? $self->value( \@vals ) : $self->value( $vals[0] );
        return $self->{db};
    }
    elsif ( !$self->node->ntnodes ) {

        # if this is a terminal node, can't add to it (use set?)
        $self->throw("Can't add child to node; only terminal node is present!");
    }
    else {
        return $self->node->addchild(@vals);
    }
}

=head2 add

 Title   : add
 Usage   : $struct->add('foo', 'bar1', 'bar2', 'bar3');
 Function: add tag-value nodes to the current node.  In the above, this would
           translate to (in XML):
           <foo>bar1</foo>
           <foo>bar2</foo>
           <foo>bar3</foo>
 Returns : 
 Args    : first arg = element name
           all other args are added as tag-value pairs

=cut

sub add {
    my ( $self, @vals ) = @_;

    # check for empty object and die for now
    if ( !$self->node->element ) {
        $self->throw("Can't add to terminal element!");
    }
    return $self->node->add(@vals);
}

=head2 set

 Title   : set
 Usage   : $struct->set('foo','bar');
 Function: sets a single tag-value pair in the current node.  Note this
           differs from add() in that this replaces any data already present
 Returns : node
 Args    : first arg = element name
           all other args are added as tag-value pairs

=cut

sub set {
    my ( $self, @vals ) = @_;

    # check for empty object
    if ( !$self->node->element ) {
        $self->throw("Can't add to tree; empty tree!");
    }
    return $self->node->set(@vals);
}

=head2 unset

 Title   : unset
 Usage   : $struct->unset('foo');
 Function: unsets all key-value pairs of the passed element from the
           current node
 Returns : node
 Args    : element name

=cut

sub unset {
    my ( $self, @vals ) = @_;
    return $self->node->unset(@vals);
}

=head2 free

 Title   : free
 Usage   : $struct->free
 Function: removes all data from the current node
 Returns : 
 Args    : 

=cut

sub free {
    my ($self) = @_;
    return $self->node->free;
}

=head2 hash

 Title   : hash
 Usage   : $struct->hash;
 Function: turns the tag-value tree into a hash, all data values are array refs
 Returns : hash
 Args    : first arg = element name
           all other args are added as tag-value pairs

=cut

sub hash {
    my ($self) = @_;
    return $self->node->hash;
}

=head2 pairs

 Title   : pairs
 Usage   : $struct->pairs;
 Function: turns the tag-value tree into a hash, all data values are scalar
 Returns : hash
 Args    : first arg = element name
           all other args are added as tag-value pairs, note that duplicates
           will be lost

=cut

sub pairs {
    my ($self) = @_;
    return $self->node->pairs;
}

=head2 qmatch

 Title    : qmatch
 Usage    : @persons = $s->qmatch('person', ('name'=>'fred'));
 Function : returns all elements in the node tree which match the
            element name and the key-value pair
 Returns  : Array of nodes
 Args     : return-element str, match-element str, match-value str

=cut

sub qmatch {
    my ( $self, @vals ) = @_;
    return $self->node->qmatch(@vals);
}

=head2 tnodes

 Title    : tnodes
 Usage    : @termini = $s->tnodes;
 Function : returns all terminal nodes below this node
 Returns  : Array of nodes
 Args     : return-element str, match-element str, match-value str

=cut

sub tnodes {
    my ($self) = @_;
    return $self->node->tnodes;
}

=head2 ntnodes

 Title    : ntnodes
 Usage    : @termini = $s->ntnodes;
 Function : returns all nonterminal nodes below this node
 Returns  : Array of nodes
 Args     : return-element str, match-element str, match-value str

=cut

sub ntnodes {
    my ($self) = @_;
    return $self->node->ntnodes;
}

=head2 StructureValue-like methods

=cut

=head2 get_all_values

 Title    : get_all_values
 Usage    : @termini = $s->get_all_values;
 Function : returns all terminal node values
 Returns  : Array of values
 Args     : return-element str, match-element str, match-value str

This is meant to emulate the values one would get from StructureValue's
get_all_values() method. Note, however, using this method dissociates the
tag-value relationship (i.e. you only get the value list, no elements)

=cut

sub get_all_values {
    my $self = shift;
    my @kids = $self->children;
    my @vals;
    while ( my $val = shift @kids ) {
        ( ref $val ) ? push @kids, $val->children : push @vals, $val;
    }
    return @vals;
}

1;
