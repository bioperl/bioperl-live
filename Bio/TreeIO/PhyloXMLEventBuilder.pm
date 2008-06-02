# $Id: PhyloXMLEventBuilder.pm 11480 2007-06-14 14:16:21Z sendu $
#
# BioPerl module for Bio::TreeIO::PhyloXMLEventBuilder
#
# Cared for by Mira Han <mirhan@indiana.edu>
#
# Copyright Mira Han
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::TreeIO::PhyloXMLEventBuilder - Build Bio::Tree::Tree's and 
  Bio::Tree::Node's from Events 

=head1 SYNOPSIS

# internal use only

=head1 DESCRIPTION

This object will take events and build a Bio::Tree::TreeI compliant
object makde up of Bio::Tree::NodeI objects.

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

=head1 AUTHOR - Mira Han

Email mirhan@indiana.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::TreeIO::PhyloXMLEventBuilder;
use strict;

use Bio::Tree::Tree;
use Bio::Tree::Node;

use base qw(Bio::Root::Root Bio::Event::EventHandlerI Bio::TreeIO::TreeEventBuilder);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::TreeIO::PhyloXMLEventBuilder->new();
 Function: Builds a new Bio::TreeIO::PhyloXMLEventBuilder object 
 Returns : Bio::TreeIO::PhyloXMLEventBuilder
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($treetype, $nodetype) = $self->_rearrange([qw(TREETYPE 
						    NODETYPE)], @args);
  $treetype ||= 'Bio::Tree::Tree';
  $nodetype ||= 'Bio::Tree::NodePhyloXML';

  eval { 
      $self->_load_module($treetype);
      $self->_load_module($nodetype);
  };

  if( $@ ) {
      $self->throw("Could not load module $treetype or $nodetype. \n$@\n")
  }
  $self->treetype($treetype);
  $self->nodetype($nodetype);
  $self->{'_treelevel'} = 0;
  $self->debug("Creating obj PhyloXMLEventBuilder\n");
  return $self;
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


=head2 SAX methods

=cut

=head2 start_document

 Title   : start_document
 Usage   : $handler->start_document
 Function: Begins a Tree event cycle
 Returns : none 
 Args    : none

=cut

sub start_document {
   my ($self) = @_;   
   $self->{'_lastitem'} = {};
   $self->{'_currentitems'} = [];
   $self->{'_currentnodes'} = [];
   $self->debug("Starting Document\n");
   return;
}

=head2 end_document

 Title   : end_document
 Usage   : my @trees = $parser->end_document
 Function: Finishes a Phylogeny cycle
 Returns : An array  Bio::Tree::TreeI
 Args    : none

=cut

sub end_document {
  my ($self,$label) = @_; 
  $self->debug("Ending Document\n");
  my $root = $self->nodetype->new(
      -id => $label,
      -verbose => $self->verbose);
# aggregate the nodes into trees basically ad-hoc.
  while ( @{$self->{'_currentnodes'}} ) {	
    my ($node) = ( shift @{$self->{'_currentnodes'}});
    $root->add_Descendent($node);
  }

  $self->debug("Root node is " . $root->to_string()."\n");
  if( $self->verbose > 0 ) { 
    foreach my $node ( $root->get_Descendents  ) {
      $self->debug("node is ". $node->to_string(). "\n");
    }
  }
  my $tree = $self->treetype->new(-verbose => $self->verbose,
      -root => $root);
  return $tree;       
}

=head2 start_element

 Title   : start_element
 Usage   :
 Function:
 Example :
 Returns : 
 Args    : $element => hashref with key 'Name'

=cut

sub start_element
{
  my ($self,$element) =@_;
  $self->{'_lastitem'}->{$element->{'Name'}}++;   

  $self->debug("starting element: $element->{Name}\n");   
  push @{$self->{'_lastitem'}->{'current'}},$element->{'Name'};

  my %data;

  if( $element->{'Name'} eq 'clade' ) 
  {
    elementAttribute($self, $element, \%data);
    push @{$self->{'_currentitems'}}, \%data; 
  }
  elsif ( $element->{'Name'} eq 'phylogeny' ) 
  {
    $self->{'_treelevel'}++;
  }
}

=head2 end_element

 Title   : end_element
 Usage   : 
 Function:
 Returns : none
 Args    : $data => hashref with key 'Name'

=cut

sub end_element{
  my ($self,$data) = @_;   

  $self->debug("end of element: $data->{Name}\n");
# this is the stack where we push/pop items from it
  my $curcount = scalar @{$self->{'_currentnodes'}};
  my $level   = $self->{'_treelevel'};
  my $levelct = $self->{'_nodect'}->[$self->{'_treelevel'}+1] || 0;

  if( $data->{'Name'} eq 'clade' ) 
  {
    my $tnode;
    my $node = pop @{$self->{'_currentitems'}};	   
    $tnode = $self->nodetype->new( -verbose => $self->verbose,
                                  %{$node});       
    $self->debug( "new node will be ".$tnode->to_string."\n");
    if ( !$node->{'-leaf'} && $levelct > 0) {
      $self->debug(join(',', map { $_->to_string } 
            @{$self->{'_currentnodes'}}). "\n");
      $self->throw("something wrong with event construction treelevel ".
          "$level is recorded as having $levelct nodes  ".
          "but current nodes at this level is $curcount\n")
        if( $levelct > $curcount);	
      for ( splice( @{$self->{'_currentnodes'}}, - $levelct)) {
        $self->debug("adding desc: " . $_->to_string . "\n");
        $tnode->add_Descendent($_);
      }
      $self->{'_nodect'}->[$self->{'_treelevel'}+1] = 0;
    }
    push @{$self->{'_currentnodes'}}, $tnode;
    $self->{'_nodect'}->[$self->{'_treelevel'}]++;

    $self->debug ("added node: nodes in stack is $curcount, treelevel: $level, nodect: $levelct\n");

  } 
  elsif(  $data->{'Name'} eq 'phylogeny' ) 
  { 
    $self->debug("end of tree: nodes in stack is $curcount\n");
    $self->{'_treelevel'}--;
  }

  $self->{'_lastitem'}->{ $data->{'Name'} }--; 

  pop @{$self->{'_lastitem'}->{'current'}};
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
   return ($e eq $self->{'_lastitem'}->{'current'}->[-1]);

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

=head2 characters

 Title   : characters
 Usage   : $handler->characters($text);
 Function: Processes characters 
 Returns : none
 Args    : text string


=cut

sub characters{
  my ($self,$ch) = @_;
  if ( $self->within_element('clade') ) 
  {
    my $hash = pop @{$self->{'_currentitems'}};
    if( $self->in_element('name')  ) 
    {
      $hash->{'-id'} = $ch->{'Data'}; # change name to id since Node.pm uses id
    } 
    push @{$self->{'_currentitems'}}, $hash;
  }
  elsif ( $self->within_element('phylogeny') )
  {
    my $hash = pop @{$self->{'_currentitems'}};
    if( $self->in_element('name')  ) 
    {
      # name of the phylogeny
    }
    if( $self->in_element('description')  ) 
    {
      # description of the phylogeny
    }
  }
}


sub elementAttribute
{
  my ( $self, $element, $data) = @_;
  my $attr = $element->{'Attributes'};
  if( $element->{'Name'} eq 'clade' ) 
  {
    if (exists $attr->{'{}distance'}) 
    {
      $data->{'-branch_length'} = $attr->{'{}distance'}->{'Value'};
    }
  }
  return $data;
}

1;
