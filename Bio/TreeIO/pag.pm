#
# BioPerl module for Bio::TreeIO::pag
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::TreeIO::pag - Bio::TreeIO driver for Pagel format

=head1 SYNOPSIS

  use Bio::TreeIO;
  my $in = Bio::TreeIO->new(-format => 'nexus',
                            -file   => 't/data/adh.mb_tree.nexus');

  my $out = Bio::TreeIO->new(-format => 'pag');
  while( my $tree = $in->next_tree ) {
    $out->write_tree($tree);
  }

=head1 DESCRIPTION

Convert a Bio::TreeIO to Pagel format.
More information here http://www.evolution.reading.ac.uk/index.html

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

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::TreeIO::pag;
use strict;

our $TaxonNameLen = 10;

use base qw(Bio::TreeIO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::TreeIO::pag->new();
 Function: Builds a new Bio::TreeIO::pag object 
 Returns : an instance of Bio::TreeIO::pag
 Args    : -file/-fh for filename or filehandles
           -name_length for minimum name length (default = 10)

=cut

sub _initialize {
    my $self = shift;
    $self->SUPER::_initialize(@_);
    my ( $name_length ) = $self->_rearrange(
        [
            qw(NAME_LENGTH)
        ],
        @_
    );
    $self->name_length( defined $name_length ? $name_length : $TaxonNameLen );
}

=head2 write_tree

 Title   : write_tree
 Usage   :
 Function: Write a tree out in Pagel format
           Some options are only appropriate for bayesianmultistate and
           the simpler output is only proper for discrete
 Returns : none
 Args    : -no_outgroups => (number)
           -print_header => 0/1 (leave 0 for discrete, 1 for bayesianms)
           -special_node => special node - not sure what they wanted to do here
           -keep_outgroup => 0/1 (keep the outgroup node in the output)
           -outgroup_ancestor => Bio::Tree::Node (if we want to exclude or include the outgroup this is what we operate on)
           -tree_no       => a tree number label - only useful for BayesianMultistate


=cut

sub write_tree {
    my ($self,$tree,@args) = @_;
    my ($keep_outgroup,
	$print_header,
	$no_outgroups,
	$special_node, 
	$outgroup_ancestor,
	$tree_no) = (0,0,1);
    my $name_len = $self->name_length;
    if( @args ) {
	($no_outgroups,
	 $print_header,
	 $special_node, 
	 $outgroup_ancestor,
	 $tree_no,
	 $keep_outgroup) = $self->_rearrange([qw(
                         NO_OUTGROUPS
						 PRINT_HEADER
						 SPECIAL_NODE
						 OUTGROUP_ANCESTOR
						 TREE_NO
						 KEEP_OUTGROUP
                         NAME_LENGTH)],@args);
    }
    my $newname_base = 1;

    my $root = $tree->get_root_node;
    my $eps = 0.0001;
    my (%chars,%names);
    my @nodes = $tree->get_nodes;
    my $species_ct;
    my $traitct;
    for my $node ( @nodes ) {
	if ((defined $special_node) && ($node eq $special_node)) {
	    my $no_of_tree_nodes = scalar(@nodes);
	    my $node_name = sprintf("N%d",$no_of_tree_nodes+1);
	    $names{$node->internal_id} = $node_name;

	} elsif ($node->is_Leaf) {
	    $species_ct++;

	    my $node_name = $node->id;
	    if( length($node_name)> $name_len ) {
		$self->warn( "Found a taxon name longer than $name_len letters, \n",
			     "name will be abbreviated.\n");
		$node_name = substr($node_name, 0,$name_len);
	    } else { 
		# $node_name = sprintf("%-".$TaxonNameLen."s",$node_name);
	    }
	    $names{$node->internal_id} = $node_name;
	    my @tags = sort $node->get_all_tags;
	    my @charstates = map { ($node->get_tag_values($_))[0] } @tags;
	    $traitct = scalar @charstates unless defined $traitct;
	    $chars{$node->internal_id} = [@charstates];
	} else {
	    $names{$node->internal_id} = sprintf("N%d", $newname_base++);
	}
    }

    # generate PAG representation
    if( $print_header ) { 
	if ($keep_outgroup) {
	    $self->_print(sprintf("%d %d\n",$species_ct,$traitct));
	} else {
	    $self->_print( sprintf("%d %d\n",$species_ct-$no_outgroups,$traitct));
	}
    }

    my @ancestors = ();
    if ($keep_outgroup) {
        push @ancestors, $root;
    } else {
	push @ancestors, ( $root, $outgroup_ancestor);
    }
    my @rest;
    foreach my $node (@nodes) {
        my $i = 0;
        foreach my $anc (@ancestors) {
            if ($anc && $node eq $anc) { $i = 1; last }
        }
        unless ($i > 0) {       # root not given in PAG
            my $current_name = $names{$node->internal_id};
	    my $branch_length_to_output;
            if ($node->branch_length < $eps) {
                my $msg_nodename = $current_name;
                $msg_nodename =~ s/\s+$//;
                warn( "TREE $tree_no, node \"$msg_nodename\": branch too ",
		      "short (", $node->branch_length, "): increasing length to ",
		      "$eps\n");
                $branch_length_to_output = $eps;
            } else {
                $branch_length_to_output = $node->branch_length;
            }
	    my @line = ( $current_name,
			 $names{$node->ancestor->internal_id},
			 $branch_length_to_output);
	    
	    if ($node->is_Leaf) {		
		push @line, @{$chars{$node->internal_id}};
		$self->_print(join(',', @line),"\n");
	    } else { 
		push @rest, \@line;
	    }
        }
    }
    for ( @rest ) { 
	$self->_print(join(',', @$_),"\n");
    }
}

=head2 next_tree

 Title   : next_tree
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub next_tree{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 name_length

 Title   : name_length
 Usage   : $self->name_length(20);
 Function: set mininum taxon name length
 Returns : integer (length of name)
 Args    : integer

=cut

sub name_length {
    my ($self, $val) = @_;
    return $self->{'name_len'} = $val if $val;
    return $self->{'name_len'};
}

1;
