#
# BioPerl module for Bio::Tree::AnnotatableTree
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chris Fields <cjfields - at - bioperl dot org>
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::AnnotatableTree - An Implementation of TreeI interface that is
Bio::AnnotatableI

=head1 SYNOPSIS

  # Do things that Bio::Tree::Tree can do...
  my $tree = Bio::TreeIO->new(-format => 'newick', -file => 'foo.txt')
    ->next_tree;
  my $root = $tree->get_root_node;
  
  # get Tree-level annotation collection
  my $ac = $tree->annotation;
  
  # do annotation-collection-y things
  ...

=head1 DESCRIPTION

This is a simple sub-class of Bio::Tree::Tree that holds tree-level annotation
data.  See Bio::TreeIO::phyloxml for an example, see L<Bio::Tree::TreeI> and
L<Bio::AnnotatableI> for specific details about those interface methods.

BEWARE: May be converged into Bio::Tree::Tree at some point. I'm mainly creating
a subclass so as not to collide with on-going refactors to Bio::Tree code by
other developers. 

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

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Chris Fields

Email cjfields at bioperl dot org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::Tree::AnnotatableTree;
use strict;

# Object preamble - inherits from Bio::Root::Root

use base qw(Bio::Tree::Tree Bio::AnnotatableI);

1;
