# $Id$
#
# BioPerl module for Bio::TreeIO::TreeEventBuilder
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::TreeIO::TreeEventBuilder - Build Bio::Tree::Tree\'s and Bio::Tree::Node\'s from Events 

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

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::TreeIO::TreeEventBuilder;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::Event::EventHandlerI;
use Bio::Tree::Tree;
use Bio::Tree::Node;

@ISA = qw(Bio::Root::Root Bio::Event::EventHandlerI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::TreeIO::TreeEventBuilder();
 Function: Builds a new Bio::TreeIO::TreeEventBuilder object 
 Returns : Bio::TreeIO::TreeEventBuilder
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  
  return $self;
}


=head2 SAX methods

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
    my ($self) = @_;    
    
     my $root = new Bio::Tree::Node();
    
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
    my $tree = new Bio::Tree::Tree(-root => $root);
    return $tree;       
}

=head2 start_element

 Title   : start_element
 Usage   :
 Function:
 Example :
 Returns : 
 Args    : $data => hashref with key 'Name'

=cut

sub start_element{
   my ($self,$data) =@_;
   $self->{'_lastitem'}->{$data->{'Name'}}++;   
   push @{$self->{'_lastitem'}->{'current'}},$data->{'Name'};

   my %data;
   
   if( $data->{'Name'} eq 'node' ) {
       push @{$self->{'_currentitems'}}, \%data; 
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
   if( $data->{'Name'} eq 'node' ) {
       my $tnode;
       my $node = pop @{$self->{'_currentitems'}};	   
       if( $node->{'-id'} ) { 
	   $tnode = new Bio::Tree::Node( %{$node});
       } else {
	   my ($right,$left) = ( pop @{$self->{'_currentnodes'}},
				 pop @{$self->{'_currentnodes'}},
				 );
	   $tnode = new Bio::Tree::Node(%{$node});
	   $tnode->add_Descendent($left);
	   $tnode->add_Descendent($right);	   
       }
       push @{$self->{'_currentnodes'}}, $tnode;       
       $self->debug ("added node: nodes in stack is ". scalar @{$self->{'_currentnodes'}}. "\n");       
   } elsif(  $data->{'Name'} eq 'tree' ) { 
       $self->debug("end of tree: nodes in stack is ". scalar @{$self->{'_currentnodes'}}. "\n");

   }
   $self->{'_lastitem'}->{$data->{'Name'}}--; 
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
   return ( $self->{'_lastitem'}->{$e} );
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
   if( $self->within_element('node') ) {
       my $hash = pop @{$self->{'_currentitems'}};
       if( $self->in_element('bootstrap') ) {
	   $hash->{'-bootstrap'} = $ch;
       } elsif( $self->in_element('branch_length') ) {
	   $hash->{'-branch_length'} = $ch;
       } elsif( $self->in_element('id')  ) {
	   $hash->{'-id'} = $ch;
       } elsif( $self->in_element('description') ) {
	   $hash->{'-desc'} = $ch;
       }
       push @{$self->{'_currentitems'}}, $hash;
   } 
}


1;
