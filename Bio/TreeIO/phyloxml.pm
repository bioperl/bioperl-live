# $Id: phyloxml.pm 11507 2007-06-23 01:37:45Z jason $
#
# BioPerl module for Bio::TreeIO::phyloxml
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

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted viax the
web:

  http://bugzilla.open-bio.org/

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
  if ($self->{'_file'}) {
    $self->{'_reader'} = XML::LibXML::Reader->new( 
                        location => $self->{'_file'},
                        no_blanks => 1
                        );
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
  );
  $self->{'_start_elements'} = \%start_elements;
  my %end_elements = (
    'phylogeny' => \&end_element_phylogeny,
    'clade' => \&end_element_clade, 
  );
  $self->{'_end_elements'} = \%end_elements;
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
    processNode($self);
  }
  return $tree;
}


=head2 write_tree

 Title   : write_tree
 Usage   : $treeio->write_tree($tree);
 Function: Write a tree out to data stream in phyloxml format
 Returns : none
 Args    : Bio::Tree::TreeI object

=cut

sub write_tree{
  my ($self, @trees) = @_;
  foreach my $tree (@trees) {
    my $root = $tree->get_root_node;
    $self->_print("<phylogeny>");
    $self->_print($self->_write_tree_Helper($root));
    $self->_print("</phylogeny>");
    $self->_print("\n");
  }
  $self->flush if $self->_flush_on_write && defined $self->_fh;
  return;
}


sub _write_tree_Helper 
{
  my ($self, $node, $str) = @_;     # this self is a Bio::Tree::phyloxml
  if (ref($node) ne 'Bio::Tree::AnnotatableNode') {
    $self->throw( "node must be a Bio::Tree::AnnotatableNode" );
  }
  my $ac = $node->annotation;
  my $seq = $node->sequence;

  # start <clade>
  $str .= '<clade';
  my @attr = $ac->get_Annotations('_attr'); # check id_source
  if (@attr) { 
    my @id_source = $attr[0]->get_Annotations('id_source');
    if (@id_source) {
      $str .= " id_source=\"".$id_source[0]->value."\"";
    }
  }
  $str .= ">";

  # print all descendent nodes
  foreach my $child ( $node->each_Descendent() ) {
    $str = $self->_write_tree_Helper($child, $str);
  }

  # print all annotations
  $str = print_annotation( $node, $str, $ac );
  # print all sequences
  if ($seq) {
    $str = print_seq_annotation( $self, $str, $seq );
  }
  
  $str .= "</clade>";
  return $str;
}


=head2 processNode

 Title   : processNode
 Usage   : 
 Function: 
 Returns : none
 Args    : 

=cut

sub processNode 
{
  my ($self) = @_;
  my $reader = $self->{'_reader'};
  if ($reader->nodeType == XML_READER_TYPE_ELEMENT) 
  {
    $self->{'_lastitem'}->{$reader->name}++;
    push @{$self->{'_lastitem'}->{'current'}}, { $reader->name=>{}};  # current holds current element and empty hash for its attributes

    if (exists $self->{'_start_elements'}->{$reader->name}) {
      my $method = $self->{'_start_elements'}->{$reader->name};
      $self->$method();
    }
    else {
      $self->element_annotation();
    }
    if ($reader->isEmptyElement) {
      # do procedures for XML_READER_TYPE_END_ELEMENT since element is complete
      
      if (exists $self->{'_end_elements'}->{$reader->name}) {
        my $method = $self->{'_end_elements'}->{$reader->name};
        $self->$method();
      }
      else {
        $self->end_element_annotation();
      }
      $self->{'_lastitem'}->{ $reader->name }--;
      pop @{$self->{'_lastitem'}->{'current'}};
    }
  }
  if ($reader->nodeType == XML_READER_TYPE_TEXT)
  {
    $self->{'_currenttext'} = $reader->value;
  } 
  if ($reader->nodeType == XML_READER_TYPE_END_ELEMENT)
  {
    if (exists $self->{'_end_elements'}->{$reader->name}) {
      my $method = $self->{'_end_elements'}->{$reader->name};
      $self->$method();
    }
    else {
      $self->end_element_annotation();
    }
    $self->{'_lastitem'}->{ $reader->name }--;
    pop @{$self->{'_lastitem'}->{'current'}};
    $self->{'_currenttext'} = '';
  }
}


=head2 processAttribute

 Title   : processAttribute
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

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
  # back at the element
  # ...

  # read a specific attribute:
  #print "Attribute b: ",$reader-> getAttribute('b'),"\n";
}


=head2 element_phylogeny

 Title   : element_phylogeny
 Usage   : $handler->element_phylogeny
 Function: Begins a Tree event cycle
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

  $self->processAttribute($self->current_attr);
  return; 
}

sub end_element_phylogeny
{
  my ($self) = @_;

  my $root;
  # if there is more than one node in _currentnodes
  # aggregate the nodes into trees basically ad-hoc.
  if ( @{$self->{'_currentnodes'}} > 1) 
  {
    $root = $self->nodetype->new( -verbose => $self->verbose, 
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
    -verbose => $self->verbose, 
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
 Usage   : $->element_clade
 Function: Begins a clade cycle
 Returns : none 
 Args    : none

=cut

sub element_clade
{
  my ($self) = @_;
  my $reader = $self->{'_reader'};
  my %data = ();    # doesn't use current attribute in order to save memory 
  $self->processAttribute(\%data);
  # create a node (Annotatable Node)
  my $tnode = $self->nodetype->new( -verbose => $self->verbose, 
                                    -id => '', 
                                    tostring => \&node_to_string,
                                    %data,
                                  );
  # add all attributes as tags (Annotation::SimpleValue)
  foreach my $tag ( keys %data ) {
    $tnode->add_tag_value( $tag, $data{$tag} );
  }
  # if there is id_source add clade to _id_link
  if (exists $data{'id_source'}) {
    $self->{'_id_link'}->{$data{'id_source'}} = $tnode;
  }
  # push into temporary list
  push @{$self->{'_currentitems'}}, $tnode;
}

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


=head2 element_annotation

 Title   : element_annotation
 Usage   : $->element_annotation
 Function: 
 Returns : none 
 Args    : none

=cut

sub element_annotation
{
  my ($self) = @_;
  my $reader = $self->{'_reader'};
  my $current = $self->current_element();
  my $prev = $self->prev_element();
  
  # read attributes of element
  $self->processAttribute($self->current_attr);

  # check idref
  my @idrefs = ();
  map { if ($_ =~ /^id_ref/) {push @idrefs, $self->current_attr->{$_};} } keys %{$self->current_attr};
 
  my @srcbyidrefs = ();
  foreach my $idref (@idrefs) { push @srcbyidrefs, $self->{'_id_link'}->{$idref}; }

  # exception when id_ref is defined but id_src is not, or vice versa.
  if (@idrefs xor @srcbyidrefs) {
    $self->throw("id_ref and id_src incompatible: @idrefs, @srcbyidrefs");
  }
  
  # we are annotating a Node
  # set _currentannotation
  if ( (@srcbyidrefs && $srcbyidrefs[0]->isa($self->nodetype)) || ((@srcbyidrefs == 0) && $prev eq 'clade') ) {
      # find node to annotate
      my $tnode;
      if (@srcbyidrefs) {
        $tnode = $srcbyidrefs[0];
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


=head2 end_element_annotation

 Title   : end_element_annotation
 Usage   : $->end_element_annotation
 Function: 
 Returns : none 
 Args    : none

=cut

sub end_element_annotation
{
  my ($self) = @_;
  my $reader = $self->{'_reader'};
  my $current = $self->current_element();
  my $prev = $self->prev_element();
  
  # check idsrc
  my $idsrc = $self->current_attr->{'id_source'};

  # check idref
  my @idrefs = ();
  map { if ($_ =~ /^id_ref/) {
    push @idrefs, $self->current_attr->{$_};
    delete $self->current_attr->{$_};
  } } keys %{$self->current_attr};
 
  my @srcbyidrefs = ();
  foreach my $idref (@idrefs) { push @srcbyidrefs, $self->{'_id_link'}->{$idref}; }

  # exception when id_src is defined but id_ref is not, or vice versa.
  if (@idrefs xor @srcbyidrefs) {
    $self->throw("id_ref and id_src incompatible: @idrefs, @srcbyidrefs");
  }

  # we are annotating a Tree
  if ( $prev eq 'phylogeny') {
    if (@srcbyidrefs) {
      # if annotation regards nodes
      if ($srcbyidrefs[0]->isa($self->nodetype))  {
        # goto case node
      }
      # if annotation regards sequences
      elsif ($srcbyidrefs[0]->isa("Bio::SeqI"))  {
        # add code to implement sequence_relation among Bio::SeqI's
      }
    }
    else {
      # annotate Tree via tree attribute
      $self->prev_attr->{$current} = $self->{'_currenttext'};
    }
  }
  # we are annotating a Node
  if (( @srcbyidrefs && $srcbyidrefs[0]->isa($self->nodetype)) || ((@srcbyidrefs == 0) && $prev eq 'clade'))  
  {
    # pop from current annotation
    my $ac = pop (@{$self->{'_currentannotation'}});
    $self->annotateNode( $current, $ac);
    # additional setups for compatibility with NodeI
      my $tnode;
      if (@srcbyidrefs) {
        $tnode = $srcbyidrefs[0];
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
        my $str = '';
        if (my @molseq = $ac->get_Annotations('mol_seq')) {
          my @strac = $molseq[0]->get_Annotations('_text');
          $str = $strac[0]->value();
        }
        my $newseq = Bio::Seq->new(-seq => $str, -annotation=>$ac);
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
  elsif ($prev eq 'clade_relation') {
  }
  # we are within an Annotation
  else {
    my $ac = pop (@{$self->{'_currentannotation'}});
    if ($ac) {
      $self->annotateNode( $current, $ac);
    }
  }
}


=head2 annotateNode

 Title   : annotateNode
 Usage   : $->annotateNode( $element, $ac)
 Function: 
 Returns : none 
 Args    : none

=cut

sub annotateNode 
{
  my ($self, $element,  $newac) = @_;
  # if attribute exists then add Annotation::Collection
  if ( scalar keys %{$self->current_attr} ) {
    my $newattr = Bio::Annotation::Collection->new();
    foreach my $tag (keys %{$self->current_attr}) {
      my $sv = new Bio::Annotation::SimpleValue(
                -value => $self->current_attr->{$tag}
               );
      $newattr->add_Annotation($tag, $sv);
    }
    $newac->add_Annotation('_attr', $newattr);
  }
  # if text exists add text as SimpleValue
  if ( $self->{'_currenttext'} ) {
    my $newvalue = new Bio::Annotation::SimpleValue( -value => $self->{'_currenttext'} );
    $newac->add_Annotation('_text', $newvalue);
  }
}


=head2 element_id

 Title   : element_id
 Usage   : $->element_id
 Function: identifier element used by phylogeny, clade, taxonomy
 Returns : none 
 Args    : none

=cut

sub element_id
{
  my ($self) = @_;
  my $reader = $self->{'_reader'};
}



=head2 current_attr

 Title   : current_attr
 Usage   :
 Function: returns the attribute hash for current item
 Example :
 Returns : 
 Args    :

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
 Usage   :
 Function: returns the attribute hash for previous item
 Example :
 Returns : 
 Args    :

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
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

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
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub prev_element {
  my ($self) = @_;

  return 0 if ! defined $self->{'_lastitem'} ||
    ! defined $self->{'_lastitem'}->{'current'}->[-2];
  my @keys = keys %{$self->{'_lastitem'}->{'current'}->[-2]};
  (@keys == 1) || die "there should be only one key for each hash";
  return $keys[0];
}

=head2 in_element

 Title   : in_element
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub in_element{
  my ($self,$e) = @_;

  return 0 if ! defined $self->{'_lastitem'} ||
    ! defined $self->{'_lastitem'}->{'current'}->[-1];
  my @keys = keys %{$self->{'_lastitem'}->{'current'}->[-1]};
  (@keys == 1) || die "there should be only one key for each hash";
  return ($e eq $keys[0]);
}

=head2 within_element

 Title   : within_element
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub within_element{
  my ($self,$e) = @_;
  return $self->{'_lastitem'}->{$e};
}

=head2 treetype

 Title   : treetype
 Usage   : $obj->treetype($newval)
 Function: 
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
 Function: 
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
 Args    : 

=cut


sub node_to_string 
{
  my ($self) = @_;     # this self is a Bio::Tree::AnnotatableNode
                       # not a Bio::TreeIO::phyloxml
  my $str = '';
  my $ac = $self->annotation;
  my $seq = $self->sequence;

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
  if ($seq) {
    $str = print_seq_annotation( $self, $str, $seq );
  }
  
  $str .= '</clade>';
  return $str;
}

sub print_annotation 
{
  my ($self, $str, $ac) = @_; 
 
  my @all_anns = $ac->get_Annotations();
  foreach my $ann (@all_anns) {
    my $key = $ann->tagname;
    if ($key eq '_attr') { next; } # attributes are already printed in the previous level 
    if  (ref($ann) eq 'Bio::Annotation::SimpleValue') 
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
    elsif (ref($ann) eq 'Bio::Annotation::Collection') 
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

sub print_attr
{
  my ($self, $str, $ac) = @_; 
  my @all_attrs = $ac->get_Annotations();
  foreach my $attr (@all_attrs) {
    if  (ref($attr) ne 'Bio::Annotation::SimpleValue') {
      $self->throw("attribute should be a SimpleValue");
    }
    $str .= ' ';
    $str .= $attr->tagname;
    $str .= '=';
    $str .= $attr->value;
  }
  return $str;
} 

sub print_seq_annotation 
{
  my ($self, $str, $seq) = @_; 
  
  $str .= "<sequence>";
  my @all_anns = $seq->annotation->get_Annotations();
  foreach my $ann (@all_anns) {
    my $key = $ann->tagname;
    if ($key eq '_attr') { next; } # attributes are already printed in the previous level 
    if  (ref($ann) eq 'Bio::Annotation::SimpleValue') 
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
    elsif (ref($ann) eq 'Bio::Annotation::Collection') 
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

