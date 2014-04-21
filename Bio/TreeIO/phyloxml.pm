# $Id: phyloxml.pm 11507 2007-06-23 01:37:45Z jason $
#
# BioPerl module for Bio::TreeIO::phyloxml
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Mira Han <mirhan@indiana.edu>
#
# Copyright Mira Han
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::TreeIO::phyloxml - TreeIO implementation for parsing PhyloXML format.

=head1 SYNOPSIS

  # do not use this module directly
  use Bio::TreeIO;
  my $treeio = Bio::TreeIO->new(-format => 'phyloxml',
                                -file => 'tree.dnd');
  my $tree = $treeio->next_tree;

=head1 DESCRIPTION

This module handles parsing and writing of phyloXML format.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted viax the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Mira Han

Email mirhan@indiana.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::TreeIO::phyloxml;
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Tree::Tree;
use Bio::Tree::AnnotatableNode;
use Bio::Annotation::SimpleValue;
use Bio::Annotation::Relation;
use XML::LibXML;
use XML::LibXML::Reader;
use base qw(Bio::TreeIO);


sub _initialize 
{
  my($self, %args) = @_;
  $args{-treetype} ||= 'Bio::Tree::Tree';
  $args{-nodetype} ||= 'Bio::Tree::AnnotatableNode';
  $self->SUPER::_initialize(%args);

  # phyloxml TreeIO does not use SAX, 
  # therefore no need to attach EventHandler
  # instead we will define a reader that is a pull-parser of libXML
  if ($self->mode eq 'r') {
    if ($self->_fh) {
      $self->{'_reader'} = XML::LibXML::Reader->new( 
                        IO => $self->_fh,
                        no_blanks => 1
                        );
    }
    if (!$self->{'_reader'}) {
      $self->throw("XML::LibXML::Reader not initialized");
    }
  }
  elsif ($self->mode eq 'w') {
    # print default lines
    $self->_print('<?xml version="1.0" encoding="UTF-8"?>',"\n");
    $self->_print('<phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns="http://www.phyloxml.org" xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd">');
  }

  $self->treetype($args{-treetype});
  $self->nodetype($args{-nodetype});
  $self->{'_lastitem'} = {}; # holds open items and the attribute hash
  $self->_init_func();
}

sub _init_func
{
  my ($self) = @_;
  my %start_elements = (
    'phylogeny' => \&element_phylogeny,
    'clade' => \&element_clade, 
    'sequence_relation' => \&element_relation, 
    'clade_relation' => \&element_relation, 
  );
  $self->{'_start_elements'} = \%start_elements;
  my %end_elements = (
    'phylogeny' => \&end_element_phylogeny,
    'clade' => \&end_element_clade, 
    'sequence_relation' => \&end_element_relation, 
    'clade_relation' => \&end_element_relation, 
  );
  $self->{'_end_elements'} = \%end_elements;
}

sub DESTROY {
  my $self = shift;
  if ($self->mode eq 'w') {
    $self->_print('</phyloxml>');
    $self->flush if $self->_flush_on_write && defined $self->_fh;
  }
  $self->SUPER::DESTROY;
}

=head2 next_tree

 Title   : next_tree
 Usage   : my $tree = $treeio->next_tree
 Function: Gets the next tree in the stream
 Returns : Bio::Tree::TreeI
 Args    : none

=cut

sub next_tree 
{
  my ($self) = @_;
  my $reader = $self->{'_reader'};
  my $tree;
  while ($reader->read) 
  {
    if ($reader->nodeType == XML_READER_TYPE_END_ELEMENT) 
    {
      if ($reader->name eq 'phylogeny') 
      {
        $tree = $self->end_element_phylogeny();
        last;
      }
    }
    $self->processXMLNode;
  }
  return $tree;
}

=head2 add_attribute

 Title   : add_phyloXML_annotation
 Usage   : my $node = $treeio->add_phyloXML_annotation(-obj=>$node, -attr=>"id_source = \"A\"")
 Function: add attributes to an object 
 Returns : the node that we added annotations to
 Args    : -obj   => object that will have the Annotation. (Bio::Tree::AnnotatableNode)
           -attr  => string in the form "A = B", where A is the attribute name and B is the attribute value

=cut

sub add_attribute
{
  my ($self, @args) = @_;
  my ($obj, $attr) = $self->_rearrange([qw(OBJ ATTR)], @args);

  if ($attr) { 
    $attr = '<dummy '.$attr.'/>';
  }
  
  my $oldreader = $self->{'_reader'};   # save reader
  $self->{'_reader'} = XML::LibXML::Reader->new( 
                string => $attr,
                no_blanks => 1
                );      
  my $reader = $self->{'_reader'};
  $self->{'_currentannotation'} = []; # holds annotationcollection 
  $self->{'_currenttext'} = '';
  #$self->{'_id_link'} = {};

  # pretend we saw a <clade> element 
  $self->{'_lastitem'}->{'dummy'}++;
  push @{$self->{'_lastitem'}->{'current'}}, { 'dummy'=>{}};  # current holds current element and empty hash for its attributes

  # push object to annotate
  push @{$self->{'_currentitems'}}, $obj;

  # read attributes of element
  while ($reader->read) 
  {
    #$self->processXMLNode;
    $self->processAttribute($self->current_attr);
  }

  # if there is id_source add sequence to _id_link
  if (exists $self->current_attr->{'id_source'}) { 
    my $idsrc = $self->current_attr->{'id_source'}; 
    $self->{'_id_link'}->{$idsrc} = $obj;
  }

  # check idref
  my $idref = '';
  if (exists $self->current_attr->{'id_ref'}) { 
    $idref = $self->current_attr->{'id_ref'}; 
  }

  my $srcbyidref = '';
  $srcbyidref = $self->{'_id_link'}->{$idref};

  # exception when id_ref is defined but id_src is not, or vice versa.
  if ($idref xor $srcbyidref) {
    $self->throw("id_ref and id_src incompatible: $idref, $srcbyidref");
  }

  # if attribute exists then add Annotation::Collection with tag '_attr'
  my $newac = $obj->annotation;
  if ( scalar keys %{$self->current_attr} ) {
    my $newattr = Bio::Annotation::Collection->new();
    foreach my $tag (keys %{$self->current_attr}) {
      my $sv = Bio::Annotation::SimpleValue->new(
          -value => $self->current_attr->{$tag}
          );
      $newattr->add_Annotation($tag, $sv);
    }
    $newac->add_Annotation('_attr', $newattr);
  }

  # pop from temporary list
  pop @{$self->{'_currentitems'}};
  $self->{'_lastitem'}->{ $reader->name }-- if $reader->name;
  pop @{$self->{'_lastitem'}->{'current'}};

  $self->{'_reader'} = $oldreader;  # restore reader
  return $obj;

}

=head2 add_phyloXML_annotation

 Title   : add_phyloXML_annotation
 Usage   : my $node = $treeio->add_phyloXML_annotation(-obj=>$node, -xml=>$xmlstring)
           my $tree = $treeio->add_phyloXML_annotation('-obj'=>$tree, '-xml'=>'<sequence_relation id_ref_0="A" id_ref_1="B" type="orthology"/>')

 Function: add annotations to a node in the phyloXML format string
 Returns : the node that we added annotations to
 Args    : -obj   => object that will have the Annotation. (Bio::Tree::AnnotatableNode)
           -xml  => string in phyloXML format that describes the annotation for the node

=cut

sub add_phyloXML_annotation
{
  my ($self, @args) = @_;
  my ($obj, $xml_string) = $self->_rearrange([qw(OBJ XML)], @args);
  
  $xml_string = '<phyloxml>'.$xml_string.'</phyloxml>';
  $self->debug( $xml_string );

  my $oldreader = $self->{'_reader'};   # save reader
  $self->{'_reader'} = XML::LibXML::Reader->new( 
                string => $xml_string,
                no_blanks => 1
                );
  my $reader = $self->{'_reader'};
  #$self->{'_currentannotation'} = []; # holds annotationcollection 
  #$self->{'_currenttext'} = '';
  #$self->{'_id_link'} = {};

  # pretend we saw a <clade> element 
  $self->{'_lastitem'}->{'clade'}++;
  push @{$self->{'_lastitem'}->{'current'}}, { 'clade'=>{}};  # current holds current element and empty hash for its attributes
  # our object to annotate (nodeI) 
  # push into temporary list
  push @{$self->{'_currentitems'}}, $obj;

  $reader->read;    #read away the first element 'phyloxml'
  while ($reader->read) 
  {
    $self->processXMLNode;
  }

  # pop from temporary list
  pop @{$self->{'_currentitems'}};
  $self->{'_lastitem'}->{ $reader->name }-- if $reader->name;
  pop @{$self->{'_lastitem'}->{'current'}};
  
  $self->{'_reader'} = $oldreader;  # restore reader
  return $obj;
}


=head2 write_tree

 Title   : write_tree
 Usage   : $treeio->write_tree($tree);
 Function: Write a tree out to data stream in phyloxml format
 Returns : none
 Args    : Bio::Tree::TreeI object

=cut

sub write_tree
{
  my ($self, @trees) = @_;
  foreach my $tree (@trees) {
    my $root = $tree->get_root_node;
    $self->_print("<phylogeny");
    my @tags = $tree->get_all_tags();
    my $attr_str = '';
    foreach my $tag (@tags) {
      my @values = $tree->get_tag_values($tag);
      foreach (@values) {
        $attr_str .= " ".$tag."=\"".$_."\"";
      }
    }
    # check if rooted
    my ($b_rooted) = $tree->get_tag_values('rooted');
    if ($b_rooted) {
      $attr_str .= " rooted=\"true\"";
    }
    else {
      if($tree->is_binary($tree->get_root_node)) {
        $attr_str .= " rooted=\"true\"";
      }
      else {
        $attr_str .= " rooted=\"false\"";
      }
    }
    $self->_print($attr_str); 
    $self->_print(">");
    if ($root->isa('Bio::Tree::AnnotatableNode')) {
      $self->_print($self->_write_tree_Helper_annotatableNode($root));
    }
    else {
      $self->_print($self->_write_tree_Helper_generic($root));
    }

    # print clade relations
    while (my $str = pop (@{$self->{'_tree_attr'}->{'clade_relation'}})) {
      $self->_print($str);
    }
    # print sequence relations
    while (my $str = pop (@{$self->{'_tree_attr'}->{'sequence_relation'}})) {
      $self->_print($str);
    }
    $self->_print("</phylogeny>");
  }
  $self->flush if $self->_flush_on_write && defined $self->_fh;
  return;
}

=head2 _write_tree_Helper_annotatableNode

 Title   : _write_tree_Helper_annotatableNode
 Usage   : internal method used by write_tree, not to be used directly
 Function: recursive helper function of write_tree for the annotatableNodes. 
           translates annotations into xml elements.
 Returns : string describing the node
 Args    : Bio::Node::AnnotatableNode object, string

=cut

sub _write_tree_Helper_annotatableNode
{
  my ($self, $node, $str) = @_;     # this self is a Bio::Tree::phyloxml
  
  my $ac = $node->annotation;

  # if clade_relation exists
  my @relations = $ac->get_Annotations('clade_relation');
  foreach (@relations) {
    my $clade_rel = $self->_relation_to_string($node, $_, '');
    # set as tree attr
    push (@{$self->{'_tree_attr'}->{'clade_relation'}}, $clade_rel);
  }

  # start <clade>
  $str .= '<clade';
  my ($attr) = $ac->get_Annotations('_attr'); # check id_source
    if ($attr) { 
      my ($id_source) = $attr->get_Annotations('id_source');
      if ($id_source) {
        $str .= " id_source=\"".$id_source->value."\"";
      }
    }
  $str .= ">";

  # print all descendent nodes
  foreach my $child ( $node->each_Descendent() ) {
    $str = $self->_write_tree_Helper_annotatableNode($child, $str);
  }

  # print all annotations
  $str = print_annotation( $node, $str, $ac );

  # print all sequences
  if ($node->has_sequence) {
    foreach my $seq (@{$node->sequence}) {
      # if sequence_relation exists
      my @relations = $seq->annotation->get_Annotations('sequence_relation');
      foreach (@relations) {
        my $sequence_rel = $self->_relation_to_string($seq, $_, '');
        # set as tree attr
        push (@{$self->{'_tree_attr'}->{'sequence_relation'}}, $sequence_rel);
      }
      $str = print_seq_annotation( $node, $str, $seq );
    }
  }

  $str .= "</clade>";

  return $str;
}

=head2 _write_tree_Helper_generic

 Title   : _write_tree_Helper_generic
 Usage   : internal method used by write_tree, not to be used directly
 Function: recursive helper function of write_tree for generic NodesI. 
           all tags are translated into property elements.
 Returns : string describing the node
 Args    : Bio::Node::NodeI object, string

=cut

sub _write_tree_Helper_generic
{
  my ($self, $node, $str) = @_;     # this self is a Bio::Tree::phyloxml
  
  # start <clade>
  $str .= '<clade>';

  # print all descendent nodes
  foreach my $child ( $node->each_Descendent() ) {
    $str = $self->_write_tree_Helper_generic($child, $str);
  }

  # print all tags
  my @tags = $node->get_all_tags();
  foreach my $tag (@tags) {
    my @values = $node->get_tag_values($tag);
    foreach my $val (@values) {
      $str .= "<property datatype=\"xsd:string\" ref=\"tag:$tag\" applies_to=\"clade\">";
      $str .=$val;
      $str .= "</property>";
    }
  }

  # print NodeI features
  if ($node->id) {
    $str .= "<name>";
    $str .= $node->id;
    $str .= "</name>";
  }
  if ($node->branch_length) {
    $str .= "<branch_length>";
    $str .= $node->branch_length;
    $str .= "</branch_length>";
  }
  if ($node->bootstrap) {
    $str .= "<confidence type = \"bootstrap\">";
    $str .= $node->bootstrap;
    $str .= "</confidence>";
  }

  $str .= "</clade>";
  return $str;
}

=head2 _relation_to_string

 Title   : _relation_to_string
 Usage   : internal method used by write_tree, not to be used directly
 Function: internal function used by write_tree to translate Annotation::Relation objects into xml elements. 
 Returns : string describing the node
 Args    : Bio::Node::AnnotatableNode (or Bio::SeqI) object that contains the Annotation::Relation, 
           the Annotation::Relation object, 
           the string

=cut

# It may be more appropriate to make Annotation::Relation have 
# a to_string callback function, 
# and have this subroutine set as the callback when we are in 
# phyloXML context.  
# I've put it here for now, since write_tree is the only place it is used.

sub _relation_to_string {
  my ($self, $obj, $rel, $str) = @_;

  my @attr = $obj->annotation->get_Annotations('_attr'); # check id_source
  if (@attr) { 
    my @id_source = $attr[0]->get_Annotations('id_source');
  }
  my ($id_ref_0) = $obj->annotation->get_nested_Annotations(
                                      '-keys' => ['id_source'],
                                      '-recursive' => 1); 
  my ($id_ref_1) = $rel->to->annotation->get_nested_Annotations( 
                                      '-keys' => ['id_source'],
                                      '-recursive' => 1); 

  my $confidence = $rel->confidence();
  my $confidence_type = $rel->confidence_type(); 
  $str .= "<";
  $str .= $rel->tagname;
  $str .= " id_ref_0=\"".$id_ref_0->value."\"";
  $str .= " id_ref_1=\"".$id_ref_1->value."\"";
  $str .= " type=\"".$rel->type."\"";
  if ($confidence) {
    $str .= " ><confidence";
    if ($confidence_type) {
      $str .= " type=\"".$confidence_type."\"";
    }
    $str .= ">";
    $str .= $confidence;
    $str .= "</confidence>";
    $str .= "</";
    $str .= $rel->tagname;
    $str .= ">";
  }
  else {
    $str .= "/>";
  }
  return $str;
}

=head2 read_annotation

 Title   : read_annotation
 Usage   : $treeio->read_annotation(-obj=>$node, -path=>$path, -attr=>1);
 Function: read text value (or attribute value) of the annotations corresponding to the element path 
 Returns : list of text values of the annotations matching the path
 Args    : -obj   => object that contains the Annotation. (Bio::Tree::AnnotatableNode or Bio::SeqI)
           -path  => path of the nested elements
           -attr  => Boolean value to indicate whether to get the attribute of the element or the text value. 
                    (default is 0, meaning text value is returned)

=cut

# It may be more appropriate to make a separate Annotation::phyloXML object
# and have this subroutine within that object so it can handle the 
# reading and writing of the values and attributes.
# but since tagTree is a temporary stub and I didn't want to make 
# a redundant object I've put it here for now.

sub read_annotation
{
  my ($self, @args) = @_;
  my ($obj, $path, $attr) = $self->_rearrange([qw(OBJ PATH ATTR)], @args);
  my $ac = $obj->annotation;
  if ($attr) {
    my @elements = split ('/', $path);
    my $final = pop @elements;
    push (@elements, '_attr');
    push (@elements, $final);
    $path = join ('/', @elements);
    return $self->_read_annotation_attr_Helper( [$ac], $path);
  } 
  else {
    return $self->_read_annotation_text_Helper( [$ac], $path);
  }
}

sub _read_annotation_text_Helper 
{
  my ($self, $acs, $path) = @_;
  my @elements = split ('/', $path);
  my $key = shift @elements;
  my @nextacs = ();
  foreach my $ac (@$acs) {
    foreach my $ann ($ac->get_Annotations($key)) {
      if ($ann->isa('Bio::AnnotationCollectionI')) {push (@nextacs, $ann)}
    }
  }
  if (@elements == 0) {
    my @values = ();
    my @texts = map {$_->get_Annotations('_text')} @nextacs;
    foreach (@texts) {
      $_ && push (@values, $_->value);
    }
    return @values;
  }
  else {
    $path = join ('/', @elements);
    return $self->_read_annotation_text_Helper( \@nextacs, $path);
  }
}

sub _read_annotation_attr_Helper 
{
  my ($self, $acs, $path) = @_;
  my @elements = split ('/', $path);
  my $key = shift @elements;
  my @nextacs = ();
  foreach my $ac (@$acs) {
    foreach my $ann ($ac->get_Annotations($key)) {
      if ($ann->isa('Bio::AnnotationCollectionI')) {push (@nextacs, $ann)}
    }
  }
  if (@elements == 1) {
    my $attrname = $elements[0];
    my @sv = map {$_->get_Annotations($attrname)} @nextacs;
    return map {$_->value} @sv;
  }
  else {
    $path = join ('/', @elements);
    return $self->_read_annotation_attr_Helper( \@nextacs, $path);
  }
}

=head1 Methods for parsing the XML document

=cut

=head2 processXMLNode

 Title   : processXMLNode
 Usage   : $treeio->processXMLNode
 Function: read the XML node and process according to the node type
 Returns : none
 Args    : none

=cut

sub processXMLNode 
{
  my ($self) = @_;
  my $reader = $self->{'_reader'};
  my $nodetype = $reader->nodeType;
  if ( $nodetype == XML_READER_TYPE_ELEMENT) 
  {
    $self->{'_lastitem'}->{$reader->name}++;
    push @{$self->{'_lastitem'}->{'current'}}, { $reader->name=>{}};  # current holds current element and empty hash for its attributes

    if (exists $self->{'_start_elements'}->{$reader->name}) {
      my $method = $self->{'_start_elements'}->{$reader->name};
      $self->$method();
    }
    else {
      $self->element_default();
    }
    if ($reader->isEmptyElement) {
      # element is complete
      # set nodetype so it can jump and
      # do procedures for XML_READER_TYPE_END_ELEMENT 
      $nodetype = XML_READER_TYPE_END_ELEMENT; 
    }
  }
  if ($nodetype == XML_READER_TYPE_TEXT)
  {
    $self->{'_currenttext'} = $reader->value;
  } 
  if ($nodetype == XML_READER_TYPE_END_ELEMENT)
  {
    if (exists $self->{'_end_elements'}->{$reader->name}) {
      my $method = $self->{'_end_elements'}->{$reader->name};
      $self->$method();
    }
    else {
      $self->end_element_default();
    }
    $self->{'_lastitem'}->{ $reader->name }--;
    pop @{$self->{'_lastitem'}->{'current'}};
    $self->{'_currenttext'} = '';
  }
}


=head2 processAttribute

 Title   : processAttribute
 Usage   : $treeio->processAttribute(\%hash_for_attribute);
 Function: reads the attributes of the current element into a hash
 Returns : none
 Args    : hash reference where the attributes will be stored.

=cut

sub processAttribute
{
  my ($self, $data) = @_;
  my $reader = $self->{'_reader'};

  # several ways of reading attributes:
  # read all attributes:
  if ($reader-> moveToFirstAttribute) {
    do {
       $data->{$reader->name()} = $reader->value;
    } while ($reader-> moveToNextAttribute);
    $reader-> moveToElement;
  }
}


=head2 element_phylogeny

 Title   : element_phylogeny
 Usage   : $treeio->element_phylogeny
 Function: initialize the parsing of a tree
 Returns : none 
 Args    : none

=cut

sub element_phylogeny 
{
  my ($self) = @_;   
  $self->{'_currentitems'} = []; # holds nodes while parsing current level
  $self->{'_currentnodes'} = []; # holds nodes while constructing tree
  $self->{'_currentannotation'} = []; # holds annotationcollection 
  $self->{'_currenttext'} = '';
  $self->{'_levelcnt'} = [];
  $self->{'_id_link'} = {};
  $self->{'_tree_attr'} = $self->current_attr;
  $self->processAttribute($self->current_attr);
  return; 
}

=head2 end_element_phylogeny

 Title   : end_element_phylogeny
 Usage   : $treeio->end_element_phylogeny
 Function: ends the parsing of a tree building a Tree::TreeI object.
 Returns : Tree::TreeI
 Args    : none

=cut

sub end_element_phylogeny
{
  my ($self) = @_;

  my $root;
  # if there is more than one node in _currentnodes
  # aggregate the nodes into trees basically ad-hoc.
  if ( @{$self->{'_currentnodes'}} > 1) 
  {
    $root = $self->nodetype->new(  
                                  -id => '',
                                  tostring => \&node_to_string,
                                );
    while ( @{$self->{'_currentnodes'}} ) {
      my ($node) = ( shift @{$self->{'_currentnodes'}});
      $root->add_Descendent($node);
    }
  }
  # if there is only one node in _currentnodes 
  # that node is root.
  elsif ( @{$self->{'_currentnodes'}} == 1) 
  {
    $root = shift @{$self->{'_currentnodes'}};
  }

  my $tree = $self->treetype->new(
    -root => $root,
    -id => $self->current_attr->{'name'},
    %{$self->current_attr}
  );
  foreach my $tag ( keys %{$self->current_attr} ) {
    $tree->add_tag_value( $tag, $self->current_attr->{$tag} );
  }
  return $tree;
}

=head2 element_clade

 Title   : element_clade
 Usage   : $treeio->element_clade
 Function: initialize the parsing of a node
           creates a new AnnotatableNode with annotations
 Returns : none 
 Args    : none

=cut

sub element_clade
{
  my ($self) = @_;
  my $reader = $self->{'_reader'};
  my %clade_attr = ();    # doesn't use current attribute in order to save memory 
  $self->processAttribute(\%clade_attr);
  # create a node (Annotatable Node)
  my $tnode = $self->nodetype->new(  
                                    -id => '', 
                                    tostring => \&node_to_string,
                                    %clade_attr,
                                  );
  # add all attributes as annotation collection with tag '_attr'
  my $ac = $tnode->annotation;
  my $newattr = Bio::Annotation::Collection->new();
  foreach my $tag (keys %clade_attr) {
    my $sv = Bio::Annotation::SimpleValue->new(
              -value => $clade_attr{$tag}
             );
    $newattr->add_Annotation($tag, $sv);
  }
  $ac->add_Annotation('_attr', $newattr);
  
  # if there is id_source add clade to _id_link
  if (exists $clade_attr{'id_source'}) {
    $self->{'_id_link'}->{$clade_attr{'id_source'}} = $tnode;
  }
  # push into temporary list
  push @{$self->{'_currentitems'}}, $tnode;
}

=head2 end_element_clade

 Title   : end_element_clade
 Usage   : $treeio->end_element_clade
 Function: ends the parsing of a node
 Returns : none 
 Args    : none

=cut

sub end_element_clade
{
  my ($self) = @_;
  my $reader = $self->{'_reader'};

  my $curcount = scalar @{$self->{'_currentnodes'}};
  my $level   = $reader->depth() - 2;
  my $childcnt = $self->{'_levelcnt'}->[$level+1] || 0; 

  # pop from temporary list
  my $tnode = pop @{$self->{'_currentitems'}};
  if ( $childcnt > 0) {
    if( $childcnt > $curcount) 
    {
      $self->throw("something wrong with event construction treelevel ".
        "$level is recorded as having $childcnt nodes  ".
        "but current nodes at this level is $curcount\n");
    }
    my @childnodes = splice( @{$self->{'_currentnodes'}}, - $childcnt);
    for ( @childnodes ) {
      $tnode->add_Descendent($_);
    }
    $self->{'_levelcnt'}->[$level+1] = 0;
  }
  push @{$self->{'_currentnodes'}}, $tnode;
  $self->{'_levelcnt'}->[$level]++;

}

=head2 element_relation

 Title   : element_relation
 Usage   : $treeio->element_relation
 Function: starts the parsing of clade relation & sequence relation
 Returns : none 
 Args    : none

=cut

sub element_relation
{
  my ($self) = @_;
  $self->processAttribute($self->current_attr);
  my $relationtype = $self->current_attr->{'type'};
  my $id_ref_0 = $self->current_attr->{'id_ref_0'};
  my $id_ref_1 = $self->current_attr->{'id_ref_1'};
  
  my @srcbyidref = ();
  $srcbyidref[0] = $self->{'_id_link'}->{$id_ref_0};
  $srcbyidref[1] = $self->{'_id_link'}->{$id_ref_1};
  
  # exception when id_ref is defined but id_src is not, or vice versa.
  if ( ($id_ref_0 xor $srcbyidref[0])||($id_ref_1 xor $srcbyidref[1]) ) {
    $self->throw("id_ref and id_src incompatible: $id_ref_0, $id_ref_1, ", $srcbyidref[0], $srcbyidref[1]);
  }

  # set id_ref_0 
  my $ac0 = $srcbyidref[0]->annotation;
  my $newann = Bio::Annotation::Relation->new(
                    '-type' => $relationtype,
                    '-to' => $srcbyidref[1],
                    '-tagname' => $self->current_element
                    );
  $ac0->add_Annotation($self->current_element, $newann);
  # set id_ref_1 
  my $ac1 = $srcbyidref[1]->annotation;
  $newann = Bio::Annotation::Relation->new(
                    '-type' => $relationtype,
                    '-to' => $srcbyidref[0],
                    '-tagname' => $self->current_element
                    );
  $ac1->add_Annotation($self->current_element, $newann);
  push (@{$self->{'_currentannotation'}}, $newann);
}

=head2 end_element_relation

 Title   : end_element_relation
 Usage   : $treeio->end_element_relation
 Function: ends the parsing of clade relation & sequence relation
 Returns : none 
 Args    : none

=cut

sub end_element_relation
{
  my ($self) = @_;
  my $ac = pop (@{$self->{'_currentannotation'}});
}


=head2 element_default

 Title   : element_default
 Usage   : $treeio->element_default
 Function: starts the parsing of all other elements
 Returns : none 
 Args    : none

=cut

sub element_default
{
  my ($self) = @_;
  my $reader = $self->{'_reader'};
  my $current = $self->current_element();
  my $prev = $self->prev_element();
  
  # read attributes of element
  $self->processAttribute($self->current_attr);

  # check idref
  my $idref = '';
  if (exists $self->current_attr->{'id_ref'}) { 
    $idref = $self->current_attr->{'id_ref'}; 
  }
 
  my $srcbyidref = '';
  $srcbyidref = $self->{'_id_link'}->{$idref};

  # exception when id_ref is defined but id_src is not, or vice versa.
  if ($idref xor $srcbyidref) {
    $self->throw("id_ref and id_src incompatible: $idref, $srcbyidref");
  }
  
  # we are annotating a Node
  # set _currentannotation
  if ( ($srcbyidref && $srcbyidref->isa($self->nodetype)) || ((!$srcbyidref) && $prev eq 'clade') ) {
      # find node to annotate
      my $tnode;
      if ($srcbyidref) {
        $tnode = $srcbyidref;
      }
      else {
        $tnode = $self->{'_currentitems'}->[-1];
      }
      my $ac = $tnode->annotation();
      # add the new anncollection with the current element as key
      my $newann = Bio::Annotation::Collection->new();
      $ac->add_Annotation($current, $newann);
      # push to current annotation
      push (@{$self->{'_currentannotation'}}, $newann);
  }
  # we are within sequence_relation or clade_relation
  elsif ($prev eq 'clade_relation' || $prev eq 'sequence_relation') {
    # do nothing?
  }
  # we are already within an annotation
  else {
    my $ac = $self->{'_currentannotation'}->[-1];
    if ($ac) {
      # add the new anncollection with the current element as key
      my $newann = Bio::Annotation::Collection->new();
      $ac->add_Annotation($current, $newann);
      push (@{$self->{'_currentannotation'}}, $newann);
    }
  }
}


=head2 end_element_default

 Title   : end_element_default
 Usage   : $treeio->end_element_default
 Function: ends the parsing of all other elements
 Returns : none 
 Args    : none

=cut

sub end_element_default
{
  my ($self) = @_;
  my $reader = $self->{'_reader'};
  my $current = $self->current_element();
  my $prev = $self->prev_element();
  
  # check idsrc
  my $idsrc = $self->current_attr->{'id_source'};

  # check idref
  my $idref = '';
  if (exists $self->current_attr->{'id_ref'}) { 
    $idref = $self->current_attr->{'id_ref'}; 
    delete $self->current_attr->{'id_ref'};
  }
 
  my $srcbyidref = '';
  $srcbyidref = $self->{'_id_link'}->{$idref};

  # exception when id_ref is defined but id_src is not, or vice versa.
  if ($idref xor $srcbyidref) {
    $self->throw("id_ref and id_src incompatible: $idref, $srcbyidref");
  }
 
  # we are annotating a Tree
  if ((!$srcbyidref) && $prev eq 'phylogeny') {
    # annotate Tree via tree attribute
    $self->prev_attr->{$current} = $self->{'_currenttext'};
  }
  # we are within sequence_relation or clade_relation
  elsif ($prev eq 'clade_relation' || $prev eq 'sequence_relation') {
    my $ann_relation = $self->{'_currentannotation'}->[-1];
    # we are here only with <confidence>
    if ($current eq 'confidence') {
      if (exists $self->current_attr->{'type'}) {
        $ann_relation->confidence_type($self->current_attr->{'type'});
      }
      $ann_relation->confidence($self->{'_currenttext'});
    }
    else {
      $self->throw($current, " is not allowed within <*_relation>");
    }
  }
  # we are annotating a Node
  elsif (( $srcbyidref && $srcbyidref->isa($self->nodetype)) || ((!$srcbyidref) && $prev eq 'clade'))  
  {
    # pop from current annotation
    my $ac = pop (@{$self->{'_currentannotation'}});
    $self->annotateNode( $current, $ac);

    # additional setups for compatibility with NodeI
    my $tnode;
    if ($srcbyidref) {
      $tnode = $srcbyidref;
    }
    else {
      $tnode = $self->{'_currentitems'}->[-1];
    }
    if ($current eq 'name') {
      $tnode->id($self->{'_currenttext'});
    }
    elsif ($current eq 'branch_length') {
      $tnode->branch_length($self->{'_currenttext'});
    }
    elsif ($current eq 'confidence') {
      if ((exists $self->current_attr->{'type'}) && ($self->current_attr->{'type'} eq 'bootstrap')) {
        $tnode->bootstrap($self->{'_currenttext'}); # this needs to change (adds 'B' annotation)
      }
    }
    elsif ($current eq 'sequence') {
      # if annotation is <sequence> 
      # transform the Bio::Annotation object into a Bio::Seq object
      my $str = '';
      # retrieve the sequence 
      if (my ($molseq) = $ac->get_Annotations('mol_seq')) {
        my ($strac) = $molseq->get_Annotations('_text');
        $str = $strac->value();
      }
      # create Seq object with sequence
      my $newseq = Bio::Seq->new( -seq => $str, 
          -annotation=>$ac, 
          -nowarnonempty=>1);
      $tnode->sequence($newseq);
      $ac->remove_Annotations('mol_seq');
      $tnode->annotation->remove_Annotations($current);
      # if there is id_source add sequence to _id_link
      if ($idsrc) {
        $self->{'_id_link'}->{$idsrc} = $newseq;
      }
    }
    elsif ($idsrc && $current eq 'taxonomy') {
      # if there is id_source add sequence to _id_link
      $self->{'_id_link'}->{$idsrc} = $ac;
    }
  }
  # we are within a default Annotation
  else {
    my $ac = pop (@{$self->{'_currentannotation'}});
    if ($ac) {
      $self->annotateNode( $current, $ac);
    }
  }
}


=head2 annotateNode

 Title   : annotateNode
 Usage   : $treeio->annotateNode($element, $ac)
 Function: adds text value and attributes to the AnnotationCollection 
           that has element name as key. If there are nested elements, 
           optional AnnotationCollections are added recursively, 
           with the nested element name as key.
           The structure of each AnnotationCollection is 
           'element' => AnnotationCollection {
               '_text' => SimpleValue (text value)
               '_attr' => AnnotationCollection { 
                   attribute1 => SimpleValue (attribute value 1)
                   attribute2 => SimpleValue (attribute value 2)
                   ...
               } 
               ['nested element' => AnnotationCollection ]
           }
 Returns : none 
 Args    : none

=cut

sub annotateNode 
{
  my ($self, $element,  $newac) = @_;
  # if attribute exists then add Annotation::Collection with tag '_attr'
  if ( scalar keys %{$self->current_attr} ) {
    my $newattr = Bio::Annotation::Collection->new();
    foreach my $tag (keys %{$self->current_attr}) {
      my $sv = Bio::Annotation::SimpleValue->new(
                -value => $self->current_attr->{$tag}
               );
      $newattr->add_Annotation($tag, $sv);
    }
    $newac->add_Annotation('_attr', $newattr);
  }
  # if text exists add text as SimpleValue with tag '_text'
  if ( $self->{'_currenttext'} ) {
    my $newvalue = Bio::Annotation::SimpleValue->new( -value => $self->{'_currenttext'} );
    $newac->add_Annotation('_text', $newvalue);
  }
}


=head1 Methods for exploring the document

=cut

=head2 current_attr

 Title   : current_attr
 Usage   : $attr_hash = $treeio->current_attr;
 Function: returns the attribute hash for current item
 Returns : reference of the attribute hash
 Args    : none

=cut

sub current_attr {
  my ($self) = @_;

  return 0 if ! defined $self->{'_lastitem'} ||
    ! defined $self->{'_lastitem'}->{'current'}->[-1];
  my @keys = keys %{$self->{'_lastitem'}->{'current'}->[-1]};
  (@keys == 1) || die "there should be only one key for each hash";
  return $self->{'_lastitem'}->{'current'}->[-1]->{$keys[0]};
}

=head2 prev_attr

 Title   : prev_attr
 Usage   : $hash_ref = $treeio->prev_attr
 Function: returns the attribute hash for previous item
 Returns : reference of the attribute hash
 Args    : none

=cut

sub prev_attr {
  my ($self) = @_;

  return 0 if ! defined $self->{'_lastitem'} ||
    ! defined $self->{'_lastitem'}->{'current'}->[-2];
  my @keys = keys %{$self->{'_lastitem'}->{'current'}->[-2]};
  (@keys == 1) || die "there should be only one key for each hash";
  return $self->{'_lastitem'}->{'current'}->[-2]->{$keys[0]};
}

=head2 current_element

 Title   : current_element
 Usage   : $element = $treeio->current_element
 Function: returns the name of the current element
 Returns : string (element name)
 Args    : none

=cut

sub current_element {
  my ($self) = @_;

  return 0 if ! defined $self->{'_lastitem'} ||
    ! defined $self->{'_lastitem'}->{'current'}->[-1];
  my @keys = keys %{$self->{'_lastitem'}->{'current'}->[-1]};
  (@keys == 1) || die "there should be only one key for each hash";
  return $keys[0];
}

=head2 prev_element

 Title   : prev_element
 Usage   : $element = $treeio->current_element
 Function: returns the name of the previous element
 Returns : string (element name)
 Args    : none

=cut

sub prev_element {
  my ($self) = @_;

  return 0 if ! defined $self->{'_lastitem'} ||
    ! defined $self->{'_lastitem'}->{'current'}->[-2];
  my @keys = keys %{$self->{'_lastitem'}->{'current'}->[-2]};
  (@keys == 1) || die "there should be only one key for each hash";
  return $keys[0];
}


=head2 treetype

 Title   : treetype
 Usage   : $obj->treetype($newval)
 Function: returns the tree type (default is Bio::Tree::Tree)
 Returns : value of treetype
 Args    : newvalue (optional)


=cut

sub treetype{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'treetype'} = $value; 
    }
    return $self->{'treetype'};
}

=head2 nodetype

 Title   : nodetype
 Usage   : $obj->nodetype($newval)
 Function: returns the node type (default is Bio::Node::AnnotatableNode)
 Returns : value of nodetype
 Args    : newvalue (optional)

=cut

sub nodetype{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'nodetype'} = $value;
    }
    return $self->{'nodetype'};
}


=head1 Methods for implementing to_string callback for AnnotatableNode

=cut

=head2 node_to_string

 Title   : node_to_string
 Usage   : $annotatablenode->to_string_callback(\&node_to_string)
 Function: set as callback in AnnotatableNode, prints the node information in string 
 Returns : string of node information
 Args    : none

=cut

# this function is similar to _write_tree_Helper_annotatableNode, 
# but it is not recursive
sub node_to_string 
{
  my ($self) = @_;     # this self is a Bio::Tree::AnnotatableNode
                       # not a Bio::TreeIO::phyloxml
  my $str = '';
  my $ac = $self->annotation;

  # start <clade>
  $str .= '<clade';
  my @attr = $ac->get_Annotations('_attr'); # check id_source
  if (@attr) { 
    my @id_source = $attr[0]->get_Annotations('id_source');
    if (@id_source) {
      $str .= " id_source=\"".$id_source[0]->value."\"";
    }
  }
  $str .= '>';

  # print all annotations
  $str = print_annotation( $self, $str, $ac );
  # print all sequences
  if ($self->has_sequence) {
    foreach my $seq (@{$self->sequence}) {
      $str = print_seq_annotation( $self, $str, $seq );
    }
  }
  
  $str .= '</clade>';
  return $str;
}

=head2 print_annotation

 Title   : print_annotation
 Usage   : $str = $annotatablenode->print_annotation($str, $annotationcollection)
 Function: prints the annotationCollection in a phyloXML format.
 Returns : string of annotation information
 Args    : string to which the Annotation should be concatenated to,
           annotationCollection that holds the Annotations

=cut

# Again, it may be more appropriate to make a separate Annotation::phyloXML object
# and have this subroutine within that object so it can handle the 
# reading and writing of the values and attributes.
# especially since this function is used both by 
# Bio::TreeIO::phyloxml (through write_tree) and 
# Bio::Node::AnnotatableNode (through node_to_string).
# but since tagTree is a temporary stub and I didn't want to make 
# a redundant object I've put it here for now.

sub print_annotation 
{
  my ($self, $str, $ac) = @_; 
 
  my @all_anns = $ac->get_Annotations();
  foreach my $ann (@all_anns) {
    my $key = $ann->tagname;
    if ($key eq '_attr') { next; } # attributes are already printed in the previous level 
    if  ($ann->isa('Bio::Annotation::SimpleValue')) 
    {
      if ($key eq '_text') {
        $str .= $ann->value;
      }
      else {
        $str .= "<$key>";
        $str .= $ann->value;
        $str .= "</$key>";
      }
    }
    elsif ($ann->isa('Bio::Annotation::Collection')) 
    {
      my @attrs = $ann->get_Annotations('_attr');
      if (@attrs) {   # if there is a attribute collection
        $str .= "<$key";
        $str = print_attr($self, $str, $attrs[0]);
        $str .= ">";
      }
      else {
        $str .= "<$key>";
      }
      $str = print_annotation($self, $str, $ann);
      $str .= "</$key>";
    }
  } 
  return $str;
}

=head2 print_attr

 Title   : print_attr
 Usage   : $str = $annotatablenode->print_attr($str, $annotationcollection)
 Function: prints the annotationCollection in a phyloXML format.
 Returns : string of attributes
 Args    : string to which the Annotation should be concatenated to,
           AnnotationCollection that holds the attributes

=cut

# Again, it may be more appropriate to make a separate Annotation::phyloXML object
# and have this subroutine within that object so it can handle the 
# reading and writing of the values and attributes.
# especially since this function is used both by 
# Bio::TreeIO::phyloxml and Bio::Node::AnnotatableNode 
# (through print_annotation).
# but since tagTree is a temporary stub and I didn't want to make 
# a redundant object I've put it here for now.

sub print_attr
{
  my ($self, $str, $ac) = @_; 
  my @all_attrs = $ac->get_Annotations();
  foreach my $attr (@all_attrs) {
    if  (!$attr->isa('Bio::Annotation::SimpleValue')) {
      $self->throw("attribute should be a SimpleValue");
    }
    $str .= ' ';
    $str .= $attr->tagname;
    $str .= '=';
    $str .= '"'.$attr->value.'"';
  }
  return $str;
} 

=head2 print_sequence_annotation

 Title   : print_sequence_annotation
 Usage   : $str = $node->print_seq_annotation( $str, $seq );
 Function: prints the Bio::Seq object associated with the node 
           in a phyloXML format.
 Returns : string that describes the sequence
 Args    : string to which the Annotation should be concatenated to,
           Seq object to print in phyloXML

=cut

# Again, it may be more appropriate to make a separate Annotation::phyloXML object
# and have this subroutine within that object so it can handle the 
# reading and writing of the values and attributes.
# especially since this function is used both by 
# Bio::TreeIO::phyloxml (through write_tree) and 
# Bio::Node::AnnotatableNode (through node_to_string).
# but since tagTree is a temporary stub and I didn't want to make 
# a redundant object I've put it here for now.


sub print_seq_annotation 
{
  my ($self, $str, $seq) = @_; 
  
  $str .= "<sequence";
  my ($attr) = $seq->annotation->get_Annotations('_attr'); # check id_source
  if ($attr) { 
    my ($id_source) = $attr->get_Annotations('id_source');
    if ($id_source) {
      $str .= " id_source=\"".$id_source->value."\"";
    }
  }
  $str .= ">";

  my @all_anns = $seq->annotation->get_Annotations();
  foreach my $ann (@all_anns) {
    my $key = $ann->tagname;
    if ($key eq '_attr') { next; } # attributes are already printed in the previous level 
    if  ($ann->isa('Bio::Annotation::SimpleValue')) 
    {
      if ($key eq '_text') {
        $str .= $ann->value;
      }
      else {
        $str .= "<$key>";
        $str .= $ann->value;
        $str .= "</$key>";
      }
    }
    elsif ($ann->isa('Bio::Annotation::Collection')) 
    {
      my @attrs = $ann->get_Annotations('_attr');
      if (@attrs) {   # if there is a attribute collection
        $str .= "<$key";
        $str = print_attr($self, $str, $attrs[0]);
        $str .= ">";
      }
      else {
        $str .= "<$key>";
      }
      $str = print_annotation($self, $str, $ann);
      $str .= "</$key>";
    }
  }
  # print mol_seq 
  if ($seq->seq()) {
    $str .= "<mol_seq>";
    $str .= $seq->seq();
    $str .= "</mol_seq>";
  }

  $str .= "</sequence>";
  return $str;
}

1;
