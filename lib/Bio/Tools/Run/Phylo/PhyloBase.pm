#
# BioPerl module for Bio::Tools::Run::Phylo::PhyloBase
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Sendu Bala
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::Phylo::PhyloBase- base module for phylo wrappers

=head1 SYNOPSIS

  use base qw(Bio::Tools::Run::Phylo::PhyloBase);

=head1 DESCRIPTION

For use by Bio::Tools::Run::Phylo modules as a base in place of
Bio::Tools::Run::WrapperBase.

This is based on WrapperBase but provides additional phylo-related private
helper subs.

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

  http://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Tools::Run::Phylo::PhyloBase;

use strict;

use Bio::AlignIO;
use Bio::TreeIO;

use base qw(Bio::Tools::Run::WrapperBase);


=head2 _alignment

 Title   : _alignment
 Usage   : $aln = $obj->_alignment()
 Function: Get/set an alignment object, generating one from a file if desired.
 Returns : Bio::Align::AlignI (probably a Bio::SimpleAlign)
 Args    : none to get
           OR filename & input format of the alignment file (latter defaults to
           guess) to set from file
           OR Bio::Align::AlignI to set

=cut

sub _alignment {
    my ($self, $thing, $format) = @_;
    
    if (ref($thing) && $thing->isa('Bio::Align::AlignI')) {
        $self->{_align_obj} = $thing;
    }
    elsif ($thing && -e $thing) {
        my $align_in = Bio::AlignIO->new(-verbose => $self->verbose, -file => $thing, $format ? (-format => $format) : ());
        my $aln = $align_in->next_aln || $self->throw("Alignment file '$thing' had no alignment!");
        $align_in->close();
        $self->{_align_obj} = $aln;
    }
    
    return $self->{_align_obj};
}

=head2 _write_alignment

 Title   : _write_alignment
 Usage   : $obj->_write_alignment()
 Function: Writes the alignment object returned by _alignment() out in the
           desired format to a temp file.
 Returns : filename
 Args    : string to describe format (default 'fasta'), any other options to pass
           to AlignIO

=cut

sub _write_alignment {
    my ($self, $format, @options) = @_;
    my $align = $self->_alignment || $self->throw("_write_alignment called when _alignment had not been set");
    $format ||= 'fasta';
    
    my ($tfh, $tempfile) = $self->io->tempfile(-dir => $self->tempdir);
    
    my $out = Bio::AlignIO->new(-verbose => $self->verbose, '-fh' => $tfh, '-format' => $format, @options);
    $align->set_displayname_flat;
    $out->write_aln($align);
    
    $out->close();
    $out = undef;
    close($tfh);
    undef $tfh;
    
    return $tempfile;
}

=head2 _tree

 Title   : _tree
 Usage   : $tree = $obj->_tree()
 Function: Get/set a tree object, generating one from a file/database if desired
 Returns : Bio::Tree::TreeI
 Args    : none to get, OR to set:
           OR filename & input format of the tree file (latter defaults to
           guess) to set from file
           OR Bio::Tree::TreeI
           OR Bio::DB::Taxonomy when _alignment() has been set and where
           sequences in the alignment have ids matching species in the taxonomy
           database

=cut

sub _tree {
    my ($self, $thing, $format) = @_;
    
    if ($thing) {
        my $tree;
        if (ref($thing) && $thing->isa('Bio::Tree::TreeI')) {
            $tree = $thing;
        }
        elsif (ref($thing) && $thing->isa('Bio::DB::Taxonomy')) {
            # get all the alignment sequence names
            my @species_names = $self->_get_seq_names;
            
            $tree = $thing->get_tree(@species_names);
            
            # convert node ids to their seq_ids for correct output with TreeIO
            foreach my $node ($tree->get_nodes) {
                my $seq_id = $node->name('supplied');
                $seq_id = $seq_id ? shift @{$seq_id} : ($node->node_name ? $node->node_name : $node->id);
                
                $node->id($seq_id);
            }
        }
        elsif (-e $thing) {
            my $tree_in = Bio::TreeIO->new(-verbose => $self->verbose, -file => $thing, $format ? (-format => $format) : ());
            $tree = $tree_in->next_tree || $self->throw("Tree file '$thing' had no tree!");
            $tree_in->close;
        }
        
        $self->{_tree_obj} = $tree || $self->throw("'$thing' supplied but unable to generate a tree from it");
    }
    
    return $self->{_tree_obj};
}

=head2 _write_tree

 Title   : _write_tree
 Usage   : $obj->_write_tree()
 Function: Writes the tree object returned by _tree() out in the desired format
           to a temp file.
 Returns : filename
 Args    : string to describe format (default 'newick')

=cut

sub _write_tree {
    my ($self, $format) = @_;
    my $tree = $self->_tree || $self->throw("_write_tree called when _tree had not been set");
    $format ||= 'newick';
    
    my ($tfh, $tempfile) = $self->io->tempfile(-dir => $self->tempdir);
    
    my $out = Bio::TreeIO->new(-verbose => $self->verbose, -fh => $tfh, -format => $format);
    $out->write_tree($tree);
    
    $out->close();
    $out = undef;
    close($tfh);
    undef $tfh;
    
    return $tempfile;
}

=head2 _get_seq_names

 Title   : _get_seq_names
 Usage   : @names = $obj->_get_seq_names()
 Function: Get all the sequence names (from id()) of the sequenes in the
           alignment.  _alignment() must be set prior to calling this.
 Returns : list of strings (seq ids)
 Args    : none

=cut

sub _get_seq_names {
    my $self = shift;
    my $aln = $self->_alignment || $self->throw("_get_seq_names called when _alignment had not been set");
    
    my @names;
    foreach my $seq ($aln->each_seq) {
        push(@names, $seq->id);
    }
    
    return @names;
}

=head2 _get_node_names

 Title   : _get_node_names
 Usage   : @names = $obj->_get_node_names()
 Function: Get all the node names (from id()) of the nodes in the tree.
           _tree must be set prior to calling this.
 Returns : list of strings (node ids)
 Args    : none

=cut

sub _get_node_names {
    my $self = shift;
    my $tree = $self->_tree || $self->throw("_get_node_names called when _tree had not been set");
    
    my @names;
    foreach my $node ($tree->get_leaf_nodes) {
        push(@names, $node->id);
    }
    
    return @names;
}

=head2 _check_names

 Title   : _check_names
 Usage   : if ($obj->_check_names) { ... }
 Function: Determine if all sequences in the alignment file have a corresponding
           node in the tree file. _alignment() and _tree() must be set
           prior to calling this.
 Returns : boolean (will also warn about the specific problems when returning
           false)
 Args    : none

=cut

sub _check_names {
    my $self = shift;
    
    my @seq_names = $self->_get_seq_names;
    my %node_names = map { $_ => 1 } $self->_get_node_names;
    
    # (not interested in tree nodes that don't map to sequence, since we
    #  expect the tree to have internal nodes not represented by sequence)
    foreach my $name (@seq_names) {
        $self->{_unmapped}{$name} = 1 unless defined $node_names{$name};
    }
    
    if (defined($self->{_unmapped})) {
        my $count = scalar(keys %{$self->{_unmapped}});
        my $unmapped = join(", ", keys %{$self->{_unmapped}});
        $self->warn("$count unmapped ids between the supplied alignment and tree: $unmapped");
        return 0;
    }
    
    return 1;
}

1;
