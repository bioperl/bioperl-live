#
# BioPerl module for Bio::Tree::TreeI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::TreeI - A Tree object suitable for lots of things, designed
  originally for Phylogenetic Trees.

=head1 SYNOPSIS

  # get a Bio::Tree::TreeI somehow
  # like from a TreeIO
  my $treeio = Bio::TreeIO->new(-format => 'newick', -file => 'treefile.dnd');
  my $tree   = $treeio->next_tree;
  my @nodes  = $tree->get_nodes;
  my @leaves = $tree->get_leaf_nodes;
  my $root   = $tree->get_root_node;

=head1 DESCRIPTION

This object holds a pointer to the Root of a Tree which is a
Bio::Tree::NodeI.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 CONTRIBUTORS

Aaron Mackey, amackey@virginia.edu
Elia Stupka,  elia@fugu-sg.org
Sendu Bala,   bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tree::TreeI;
use strict;

use base qw(Bio::Root::RootI);

sub is_rooted {
   shift->throw_not_implemented();
}

sub root {
   shift->throw_not_implemented();
}

sub id{
   shift->throw_not_implemented();
}

sub score{
   shift->throw_not_implemented();
}

1;
