# $Id$
#
# BioPerl module for Bio::Tree::NodeI
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::NodeI - Interface describing a Tree Node

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the interface here

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


package Bio::Tree::NodeI;
use vars qw(@ISA);
use strict;
use Carp;


sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller)[1];
  confess "Abstract method '$caller' defined in interfaceBio::Tree::NodeI not implemented by pacakge $package. Not your fault - author of $package should be blamed!";
}


=head2 get_Descendents

 Title   : get_Descendents
 Usage   : my @nodes = $node->get_Descendents;
 Function: Recursively fetches all the child nodes
 Returns : Array of Bio::Tree::NodeI objects
 Args    : none

=cut

sub get_Descendents{
   my ($self) = @_;
   my @children;
   if( my $left = $self->get_Left_Descendent ) { 
       push @children, ($left, $left->get_Descendents); 
   } 
   if( my $right = $self->get_Right_Descendent ) { 
       push @children, ($right, $right->get_Descendents); 
   }
   return @children;
}

=head2 is_Leaf

 Title   : is_Leaf
 Usage   : if( $node->is_Leaf ) 
 Function: Get Leaf status
 Returns : boolean
 Args    : none

=cut

sub is_Leaf{
    my ($self) = @_;
    return if( ! $self->get_Right_Descendent &&
	       ! $self->get_Left_Descendent );
}


=head2 descendent_count

 Title   : descendent_count
 Usage   : my $count = $node->descendent_count;
 Function: Counts the number of descendents a node has 
           (and all of their subnodes)
 Returns : integer
 Args    : none

=cut

sub descendent_count{
   my ($self) = @_;
   my $count = 0;
   
   # how to avoid the possiblility that a node has 2 places in a tree...
   if( my $left = $self->get_Left_Descendent ) { 
       $count += $left->descendent_count;
   }
   
   if( my $right = $self->get_Right_Descendent ) { 
       $count += $right->descendent_count;
   }   
   return $count;
}

=head2 branch_length

 Title   : branch_length
 Usage   : $obj->branch_length($newval)
 Function: 
 Example : 
 Returns : value of branch_length
 Args    : newvalue (optional)


=cut

sub branch_length{
    my ($self)= @_;
    $self->_abstractDeath;

}


=head2 to_string

 Title   : to_string
 Usage   : my $str = $node->to_string()
 Function: For debugging, provide a node as a string
 Returns : string
 Args    : none


=cut

sub to_string{
   my ($self) = @_;
   return sprintf("BL:%s",defined $self->branch_length ? 
		  $self->branch_length : ' ');
}


1;
