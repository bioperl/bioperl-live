# $Id$
#
# BioPerl module for Bio::TreeIO::tabtree
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::TreeIO::tabtree - A simple output format which displays a tree as an ASCII drawing

=head1 SYNOPSIS

  use Bio::TreeIO;
  my $in = new Bio::TreeIO(-file => 'input', -format => 'newick');
  my $out = new Bio::TreeIO(-file => '>output', -format => 'tabtree');

  while( my $tree = $in->next_tree ) {
      $out->write_tree($tree);
  }

=head1 DESCRIPTION

This is a made up format just for outputting trees as an ASCII drawing.

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


package Bio::TreeIO::tabtree;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::TreeIO;


@ISA = qw(Bio::TreeIO );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::TreeIO::tabtree();
 Function: Builds a new Bio::TreeIO::tabtree object 
 Returns : Bio::TreeIO::tabtree
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

}

=head2 write_tree

 Title   : write_tree
 Usage   : $treeio->write_tree($tree);
 Function: Write a tree out to data stream in newick/phylip format
 Returns : none
 Args    : Bio::Tree::TreeI object

=cut

sub write_tree{
   my ($self,$tree) = @_;      
   my @data = _write_tree_Helper($tree->get_root_node,0);
   $self->_print(join("\n", @data), "\n");   
   return;
}

sub _write_tree_Helper {
    my ($node,$depth) = @_;
    return () if (!defined $node);

    my @data;
    my @d = $node->each_Descendent();
    my $half = ( int ( (scalar @d) / 2)) || 1;
    if( @d ) {
	my $c = 0;
	foreach my $n ( @d ) {
	    if( $c++ == $half ) {
		push @data,sprintf("%s%s","\t"x$depth,
				   defined $node->id && $node->id ? 
				   $node->id : 'O',
				   defined $node->branch_length ? ":" .
				   $node->branch_length : '');
	    }   
	    push @data, _write_tree_Helper($n,$depth+1);
	}
    } else {
	push @data,sprintf("%s%s","\t"x$depth,
			   defined $node->id && $node->id ? 
			   $node->id : 'O',
			   defined $node->branch_length ? ":" .
			   $node->branch_length : '');
    }

    return @data;
}

=head2 next_tree

 Title   : next_tree
 Usage   : 
 Function: Sorry not possible with this format
 Returns : none
 Args    : none


=cut

sub next_tree{
    $_[0]->throw("Sorry the format 'tabtree' can only be used as an output format at this time");
}

1;
