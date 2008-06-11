# BioPerl module for Bio::Tree::TreePhyloXML
#
# Cared for by Mira Han <mirhan@indiana.edu>
#
# Copyright Mira Han
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::TreePhyloXML - A Tree with support for PhyloXML tags.

=head1 SYNOPSIS

    # like from a TreeIO
    my $treeio = Bio::TreeIO->new(-format => 'phyloxml', -file => 'treefile.xml');
    my $tree = $treeio->next_tree;
    my @nodes = $tree->get_nodes;
    my $root = $tree->get_root_node;


=head1 DESCRIPTION

This object subclasses Tree to handle phyloxml specific variables.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Mira Han

Email mirhan@indiana.edu

=head1 CONTRIBUTORS

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tree::TreePhyloXML;
use strict;

use base qw(Bio::Tree::Tree);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tree::TreePhyloXML->new();
 Function: Builds a new Bio::Tree::TreePhyloXML object 
 Returns : Bio::Tree::TreePhyloXML
 Args    : -root     => L<Bio::Tree::NodeI> object which is the root
             OR
           -node     => L<Bio::Tree::NodeI> object from which the root will be
                        determined

           -nodelete => boolean, whether or not to try and cleanup all
                                 the nodes when this this tree goes out
                                 of scope.
           -id       => optional tree ID
           -score    => optional tree score value
           -rooted    => optional boolean, rooted
           -rerootable    => optional boolean, rerootable
           -type    => optional i.e. 'gene tree' 
           -branch_len_unit => optional 

=cut

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  $self->debug("new TreePhyloXML\n");
  my ($root, $rooted,$rerootable,$type,$branch_len_unit)= $self->_rearrange([qw(ROOT ROOTED REROOTABLE TYPE BRANCH_LEN_UNIT)], @args);
  # need to check if root is a NodePhyloXML object?
  $self->{'_rooted'} = $rooted;
  $self->{'_rerootable'} = $rerootable;
  $self->{'_type'} = $type;
  $self->{'_branch_len_unit'} = $branch_len_unit;
  return $self;
}

1;
