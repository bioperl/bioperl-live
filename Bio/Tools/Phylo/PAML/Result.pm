# $Id$
#
# BioPerl module for Bio::Tools::Phylo::PAML::Result
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich, Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Phylo::PAML::Result - A PAML result set object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

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

=head1 AUTHOR - Jason Stajich, Aaron Mackey

Email jason@bioperl.org
Email amackey@virginia.edu

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Phylo::PAML::Result;
use vars qw(@ISA);
use strict;


use Bio::Root::Root;
use Bio::Factory::TreeFactoryI;

@ISA = qw(Bio::Root::Root Bio::Factory::TreeFactoryI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tools::Phylo::PAML::Result();
 Function: Builds a new Bio::Tools::Phylo::PAML::Result object 
 Returns : Bio::Tools::Phylo::PAML::Result
 Args    : -trees => array reference of L<Bio::Tree::TreeI> objects
           -MLmatrix => ML matrix
           .... MORE ARGUMENTS LISTED HERE BY AARON AND JASON 

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  my ($trees) = $self->_rearrange([qw(TREES)],@args);
  
  if( defined $trees ) {
      if(ref($trees) !~ /ARRAY/i ) { 
	  $self->warn("Must have provided a valid array reference to initialize trees");
      } else { 
	  foreach my $t ( @$trees ) {
	      $self->add_tree($t);
	  }
      }
  }
  $self->{'_treeiterator'} = 0;

  return $self;
}

=head2 next_tree

 Title   : next_tree
 Usage   : my $tree = $factory->next_tree;
 Function: Get the next tree from the factory
 Returns : L<Bio::Tree::TreeI>
 Args    : none

=cut

sub next_tree{
   my ($self,@args) = @_;
   return $self->{'_trees'}->[$self->{'_treeiterator'}++] || undef;
}

=head2 rewind_tree

 Title   : rewind_tree_iterator
 Usage   : $result->rewind_tree()
 Function: Rewinds the tree iterator so that next_tree can be 
           called again from the beginning
 Returns : none
 Args    : none

=cut

sub rewind_tree_iterator {
    shift->{'_treeiterator'} = 0;
}

=head2 add_tree

 Title   : add_tree
 Usage   : $result->add_tree($tree);
 Function: Adds a tree 
 Returns : integer which is the number of trees stored
 Args    : L<Bio::Tree::TreeI>

=cut

sub add_tree{
   my ($self,$tree) = @_;
   if( $tree && ref($tree) && $tree->isa('Bio::Tree::TreeI') ) {
       push @{$self->{'_trees'}},$tree;
   }
   return scalar @{$self->{'_trees'}};
}

1;
