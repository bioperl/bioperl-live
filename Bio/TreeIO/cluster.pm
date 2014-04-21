#
# BioPerl module for Bio::TreeIO::cluster
#
# Contributed by Guillaume Rousse <Guillaume-dot-Rousse-at-inria-dot-fr>
#
# Copyright INRIA
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::TreeIO::cluster - A TreeIO driver module for parsing Algorithm::Cluster::treecluster() output

=head1 SYNOPSIS

  # do not use this module directly
  use Bio::TreeIO;
  use Algorithm::Cluster;
  my ($result, $linkdist) = Algorithm::Cluster::treecluster(
    distances => $matrix
  );
  my $treeio = Bio::TreeIO->new(
    -format   => 'cluster',
    -result   =>  $result,
    -linkdist =>  $linkdist,
    -labels   =>  $labels
  );
  my $tree = $treeio->next_tree;

=head1 DESCRIPTION

This is a driver module for parsing Algorithm::Cluster::treecluster() output.

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

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Guillaume Rousse

Email Guillaume-dot-Rousse-at-inria-dot-fr

=head1 CONTRIBUTORS

Jason Stajich - jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::TreeIO::cluster;
use strict;

use Bio::Event::EventGeneratorI;
use IO::String;

use base qw(Bio::TreeIO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::TreeIO::cluster->new();
 Function: Builds a new Bio::TreeIO::cluster object for reading Algorithm::Cluster::treecluster output
 Returns : Bio::TreeIO::cluster
 Args    :-result   => Algorithm::Cluster result
          -linkdist => distance between links
          -labels   => node labels

=cut

sub _initialize {
  my $self = shift;
  ($self->{_result},$self->{_linkdist},
   $self->{_labels}) = $self->_rearrange([qw
					  (RESULT LINKDIST LABELS)],
					 @_);
  $self->SUPER::_initialize(@_);
}

=head2 next_tree

 Title   : next_tree
 Usage   : my $tree = $treeio->next_tree
 Function: Gets the next tree in the stream
 Returns : Bio::Tree::TreeI
 Args    : none


=cut

sub next_tree {
    my ($self) = @_;
    if( ! $self->{_result} ){
	$self->warn("Must provide value 'result' and 'linkdist' and 'labels' when initializing a TreeIO::cluster object");
	return;
    }
    $self->_eventHandler->start_document();

    # build tree from the root
    $self->_eventHandler->start_element({Name => 'tree'});
    $self->_recurse(-1, 0);
    $self->_recurse(-1, 1);
    $self->_eventHandler->end_element({Name => 'tree'});

    return $self->_eventHandler->end_document;
}

sub _recurse {
    my ($self, $line, $column) = @_;

    my $id  = $self->{_result}->[$line]->[$column];
    if ($id >= 0) {
	# leaf
	$self->debug("leaf $id\n");
	$self->debug("distance $self->{_linkdist}->[$line]\n");
	$self->debug("label $self->{_labels}->[$id]\n");
	$self->_eventHandler->start_element({Name => 'node'});
	$self->_eventHandler->start_element({Name => 'branch_length'});
	$self->_eventHandler->characters($self->{_linkdist}->[$line]);
	$self->_eventHandler->end_element({Name => 'branch_length'});
	$self->_eventHandler->start_element({Name => 'id'});
	$self->_eventHandler->characters($self->{_labels}->[$id]);
	$self->_eventHandler->end_element({Name => 'id'});
	$self->_eventHandler->start_element({Name => 'leaf'});
	$self->_eventHandler->characters(1);
	$self->_eventHandler->end_element({Name => 'leaf'});
	$self->_eventHandler->end_element({Name => 'node'});
    } else {
	# internal node
	$self->debug("internal node $id\n");
	$self->debug("distance $self->{_linkdist}->[$line]\n");
	$self->_eventHandler->start_element({Name => 'node'});
	$self->_eventHandler->start_element({Name => 'branch_length'});
	$self->_eventHandler->characters($self->{_linkdist}->[$line]);
	$self->_eventHandler->end_element({Name => 'branch_length'});
	$self->_eventHandler->start_element({Name => 'leaf'});
	$self->_eventHandler->characters(0);
	$self->_eventHandler->end_element({Name => 'leaf'});
	$self->_eventHandler->start_element({Name => 'tree'});
	my $child_id = - ($id + 1);
	$self->_recurse($child_id, 0);
	$self->_recurse($child_id, 1);
	$self->_eventHandler->end_element({Name => 'tree'});
	$self->_eventHandler->end_element({Name => 'node'});

    }
}

=head2 write_tree

 Title   : write_tree
 Usage   :
 Function: Sorry not possible with this format
 Returns : none
 Args    : none


=cut

sub write_tree{
    $_[0]->throw("Sorry the format 'cluster' can only be used as an input format");
}

1;
