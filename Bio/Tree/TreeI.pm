# $Id$
#
# BioPerl module for Bio::Tree::TreeI
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::TreeI - A Tree object suitable for lots of things, designed originally for Phylogenetic Trees.

=head1 SYNOPSIS

# get a Bio::Tree::TreeI somehow
# like from a TreeIO
    my $treeio = new Bio::TreeIO(-format => 'newick', -file => 'treefile.dnd');
    my $tree = $treeio->next_tree;
    my @nodes = $tree->get_nodes;
    my $root = $tree->get_root_node;

=head1 DESCRIPTION

This object holds a pointer to the Root of a Tree which is a
Bio::Tree::NodeI.

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

Elia Stupka, elia@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tree::TreeI;
use vars qw(@ISA);
use strict;

use Carp;

sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller)[1];
  
  confess "Abstract method '$caller' defined in interface Bio::Tree::TreeI not implemented by pacakge $package. Not your fault - author of $package should be blamed!";
}

=head2 get_nodes

 Title   : get_nodes
 Usage   : my @nodes = $tree->get_nodes()
 Function: Return list of Tree::NodeI objects
 Returns : array of Tree::NodeI objects
 Args    : order => 'b|breadth' first order or 'd|depth' first order

=cut

sub get_nodes{
   my ($self) = @_;
   $self->_abstractDeath;
}

=head2 get_root_node

 Title   : get_root_node
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_root_nodes{
   my ($self,@args) = @_;
   $self->_abstractDeath;
}

1;
