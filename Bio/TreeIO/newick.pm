#
# BioPerl module for Bio::TreeIO::newick
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

Bio::TreeIO::newick

=head1 SYNOPSIS

  # do not use this module directly
  use Bio::TreeIO;

  my $treeio = Bio::TreeIO->new(-format => 'newick', 
                               -file => 't/data/LOAD_Ccd1.dnd');
  my $tree = $treeio->next_tree;

=head1 DESCRIPTION

This module handles parsing and writing of Newick/PHYLIP/New Hampshire format.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::TreeIO::newick;
use strict;

use Bio::Event::EventGeneratorI;

use base qw(Bio::TreeIO Bio::TreeIO::NewickParser);

=head2 new

Title   : new
Args    : -print_count     => boolean  default is false
           -bootstrap_style => set the bootstrap style (one of nobranchlength,
							molphy, traditional)
           -order_by => set the order by sort method 

See L<Bio::Node::Node::each_Descendent()>

=cut

sub _initialize {
    my $self = shift;
    $self->SUPER::_initialize(@_);
    my ( $print_count ) = $self->_rearrange(
        [
            qw(PRINT_COUNT)
        ],
        @_
    );
    $self->print_tree_count( $print_count || 0 );
    return;
}

=head2 next_tree

Title   : next_tree
Usage   : my $tree = $treeio->next_tree
Function: Gets the next tree in the stream
Returns : L<Bio::Tree::TreeI>
Args    : none

=cut

sub next_tree {
    my ($self) = @_;
    local $/ = ";\n";
    return unless $_ = $self->_readline;

    s/[\r\n]//gs;
    my $score;
    my $despace = sub { my $dirty = shift; $dirty =~ s/\s+//gs; return $dirty };
    my $dequote = sub {
        my $dirty = shift;
        $dirty =~ s/^"?\s*(.+?)\s*"?$/$1/;
        return $dirty;
    };
s/([^"]*)(".+?")([^"]*)/$despace->($1) . $dequote->($2) . $despace->($3)/egsx;

    if (s/^\s*\[([^\]]+)\]//) {
        my $match = $1;
        $match =~ s/\s//g;
        $match =~ s/lh\=//;
        if ( $match =~ /([-\d\.+]+)/ ) {
            $score = $1;
        }
    }

    $self->_eventHandler->start_document;

    # Call the parse_newick method as defined in NewickParser.pm
    $self->parse_newick($_);

    my $tree = $self->_eventHandler->end_document;

    # Add the tree score afterwards if it exists.
    if (defined $tree) {
      $tree->score($score);
      return $tree;
    }
}

# Returns the default set of parsing & writing parameters for the Newick format.
sub get_default_params {
  my $self = shift;

  return {
    newline_each_node => 0,
    order_by => '', # ???
    bootstrap_style => 'traditional', # Can be 'traditional', 'molphy', 'nobranchlength'
    internal_node_id => 'id', # Can be 'id' or 'bootstrap'
    
    no_branch_lengths => 0,
    no_bootstrap_values => 0,
    no_internal_node_labels => 0
  };
}


=head2 write_tree

Title   : write_tree
Usage   : $treeio->write_tree($tree);
Function: Write a tree out to data stream in newick/phylip format
Returns : none
Args    : L<Bio::Tree::TreeI> object

=cut

sub write_tree {
    my ( $self, @trees ) = @_;
    if ( $self->print_tree_count ) {
        $self->_print( sprintf( " %d\n", scalar @trees ) );
    }

    my $params = $self->get_params;

    foreach my $tree (@trees) {
      if (  !defined $tree
            || ref($tree) =~ /ARRAY/i
            || !$tree->isa('Bio::Tree::TreeI') )
      {
        $self->throw(
                     "Calling write_tree with non Bio::Tree::TreeI object\n");
      }
      my @data = $self->_write_tree_Helper( $tree->get_root_node, $params);
      $self->_print( join( ',', @data ).";" );
    }
    
    $self->flush if $self->_flush_on_write && defined $self->_fh;
    return;
}

sub _write_tree_Helper {
  my $self = shift;
    my ( $node, $params ) = @_;
    my @data;

    foreach my $n ( $node->each_Descendent($params->{order_by}) ) {
        push @data, $self->_write_tree_Helper( $n, $params );
    }

    my $label = $self->_node_as_string($node,$params);

    if ( scalar(@data) >= 1) {
      $data[0] = "(" . $data[0];
      $data[-1] .= ")";
      $data[-1] .= $label;
    } else {
      push @data, $label;
    }

    return @data;    
}

sub _node_as_string {
  my $self = shift;
  my $node = shift;
  my $params = shift;

  my $label_stringbuffer = '';

  if ($params->{no_bootstrap_values} != 1 &&
      !$node->is_Leaf && 
      defined $node->bootstrap &&
      $params->{bootstrap_style} eq 'traditional' &&
      $params->{internal_node_id} eq 'bootstrap') {
    # If we're an internal node and we're using 'traditional' bootstrap style,
    # we output the bootstrap instead of any label.
    my $bootstrap = $node->bootstrap;
    $label_stringbuffer .= $bootstrap if (defined $bootstrap);
  } elsif ($params->{no_internal_node_labels} != 1) {
    my $id = $node->id;
    $label_stringbuffer .= $id  if( defined $id );
  }

  if ($params->{no_branch_lengths} != 1) {
    my $blen  = $node->branch_length;
    $label_stringbuffer .= ":". $blen if (defined $blen);
  }  

  if ($params->{bootstrap_style} eq 'molphy') {
    my $bootstrap = $node->bootstrap;
    $label_stringbuffer .= "[$bootstrap]" if (defined $bootstrap);
  }

  if ($params->{newline_each_node} == 1) {
    $label_stringbuffer .= "\n";
  }

  return $label_stringbuffer;
}


=head2 print_tree_count

Title   : print_tree_count
Usage   : $obj->print_tree_count($newval)
Function: Get/Set flag for printing out the tree count (paml,protml way)
Returns : value of print_tree_count (a scalar)
Args    : on set, new value (a scalar or undef, optional)

=cut

sub print_tree_count {
    my $self = shift;
    return $self->{'_print_tree_count'} = shift if @_;
    return $self->{'_print_tree_count'} || 0;
}

=head2 bootstrap_style

Title   : bootstrap_style
Usage   : $obj->bootstrap_style($newval)
Function: A description of how bootstraps and branch lengths are
           written, as the ID part of the internal node or else in []
           in the branch length (Molphy-like; I am sure there is a
           better name for this but am not sure where to go for some
           sort of format documentation)

           If no branch lengths are requested then no bootstraps are usually
           written (unless someone REALLY wants this functionality...)

           Can take on strings which contain the possible values of
           'nobranchlength'   --> don't draw any branch lengths - this
                                  is helpful if you don't want to have to 
                                  go through and delete branch len on all nodes
           'molphy' --> draw bootstraps (100) like
                                  (A:0.11,B:0.22):0.33[100];
           'traditional' --> draw bootstraps (100) like
                                  (A:0.11,B:0.22)100:0.33;
Returns : value of bootstrap_style (a scalar)
Args    : on set, new value (a scalar or undef, optional)

=cut

sub bootstrap_style {
    my $self = shift;
    my $val  = shift;
    if ( defined $val ) {

        if ( $val !~ /^nobranchlength|molphy|traditional/i ) {
            $self->warn(
"requested an unknown bootstrap style $val, expect one of nobranchlength,molphy,traditional, not updating value.\n"
            );
        }
        else {
            $self->{'_bootstrap_style'} = $val;
        }
    }
    return $self->{'_bootstrap_style'} || 'traditional';
}

=head2 order_by

Title   : order_by
Usage   : $obj->order_by($newval)
Function: Allow node order to be specified (typically "alpha")
           See L<Bio::Node::Node::each_Descendent()>
Returns : value of order_by (a scalar)
Args    : on set, new value (a scalar or undef, optional)

=cut

sub order_by {
    my $self = shift;

    return $self->{'order_by'} = shift if @_;
    return $self->{'order_by'};
}

1;
