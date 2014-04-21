#
# BioPerl module for Bio::TreeIO::nexus
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-open-bio-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::TreeIO::nexus - A TreeIO driver module for parsing Nexus tree output from PAUP

=head1 SYNOPSIS

  use Bio::TreeIO;
  my $in = Bio::TreeIO->new(-file => 't/data/cat_tre.tre');
  while( my $tree = $in->next_tree ) {
  }

=head1 DESCRIPTION

This is a driver module for parsing PAUP Nexus tree format which
basically is just a remapping of trees.

=head2 Comments

The nexus format allows node comments that are placed inside square
brackets. Usually the comments (implemented as tags for nodes) are
used to give a name for an internal node or record the bootstap value,
but other uses are possible.

The FigTree program by Andrew Rambaut adds various rendering
parameters inside comments and flags these comments by starting them
with '&!'. The parameters implemented here are 'label' and 'color'.

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

Email jason-at-open-bio-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::TreeIO::nexus;
use strict;

use Bio::Event::EventGeneratorI;
use IO::String;

use base qw(Bio::TreeIO);

=head2 new

 Title   : new
 Args    : -header    => boolean  default is true 
                         print/do not print #NEXUS header
           -translate => boolean default is true
                         print/do not print Node Id translation to a number

=cut

sub _initialize {
    my $self = shift;
    $self->SUPER::_initialize(@_);
    my ( $hdr, $trans ) = $self->_rearrange(
        [
            qw(HEADER
              TRANSLATE)
        ],
        @_
    );
    $self->header( defined $hdr           ? $hdr   : 1 );
    $self->translate_node( defined $trans ? $trans : 1 );
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
    unless ( $self->{'_parsed'} ) {
        $self->_parse;
    }
    return $self->{'_trees'}->[ $self->{'_treeiter'}++ ];
}

sub rewind {
    shift->{'_treeiter'} = 0;
}

sub _parse {
    my ($self) = @_;

    $self->{'_parsed'}   = 1;
    $self->{'_treeiter'} = 0;

    while ( defined( $_ = $self->_readline ) ) {
        next if /^\s+$/;
        last;
    }
    return unless ( defined $_ );

    unless (/^\#NEXUS/i) {
        $self->warn("File does not start with #NEXUS");    #'
        return;
    }

    my $line;
    while ( defined( $_ = $self->_readline ) ) {
        $line .= $_;
    }
    my @sections = split( /#NEXUS/i, $line );
    for my $s (@sections) {
        my %translate;        
        if ( $self->verbose > 0 ) {
            while ( $s =~ s/(\[[^\]]+\])// ) {
                $self->debug("removing comment $1\n");
            }
        }
        else {
            $s =~ s/(\[[^\]]+\])//g;
        }
        
        if ( $s =~ /begin trees;(.+)(end;)?/si ) {
            my $trees = $1;
            if ( $trees =~ s/\s+translate\s+([^;]+);//i ) {
                my @trans;
                my $tr = $1;

                while ($tr =~ m{\s*([^,\s]+?\s+(?:'[^']+'|[^'\s]+)),?}gc) {
	            push @trans, $1;
                }
                for my $n ( @trans ) {
                    if ($n  =~ /^\s*(\S+)\s+(.+)$/) {
                        my ($id,$tag) = ($1,$2);
                        $tag =~ s/[\s,]+$//;  # remove the extra spaces of the last taxon
                        $translate{$id} = $tag;
                    }                    
                }
            }
            else {
                $self->debug("no translate in: $trees\n");
            }
            while ($trees =~ /\s+tree\s+\*?\s*(\S+)\s*\=
             \s*(?:\[\S+\])?\s*([^\;]+;)/igx)
            {
                my ( $tree_name, $tree_str ) = ( $1, $2 );

                # MrBayes does not print colons for node label
                # $tree_str =~ s/\)(\d*\.\d+)\)/:$1/g;
                my $buf    = IO::String->new($tree_str);
                my $treeio = Bio::TreeIO->new(
                    -format => 'newick',
                    -fh     => $buf
                );
                my $tree = $treeio->next_tree;
                foreach my $node ( grep { $_->is_Leaf } $tree->get_nodes ) {
                    my $id     = $node->id;
                    my $lookup = $translate{$id};
                    $node->id( $lookup || $id );
                }
                $tree->id($tree_name) if defined $tree_name;
                push @{ $self->{'_trees'} }, $tree;
            }
        }
        else {
            $self->debug("begin_trees failed: $s\n");
        }
    }
    if ( !@sections ) {
        $self->debug("warn no sections: $line\n");
    }
}

=head2 write_tree

 Title   : write_tree
 Usage   : $treeio->write_tree($tree);
 Function: Writes a tree onto the stream
 Returns : none
 Args    : Bio::Tree::TreeI


=cut

sub write_tree {
    my ( $self, @trees ) = @_;
    if ( $self->header ) {
        $self->_print("#NEXUS\n\n");
    }
    my $translate = $self->translate_node;
    my $time      = localtime();
    $self->_print( sprintf( "Begin trees; [Treefile created %s]\n", $time ) );

    my ( $first, $nodecter, %node2num ) = ( 0, 1 );
    foreach my $tree (@trees) {

        if (   $first == 0
            && $translate )
        {
            $self->_print("\tTranslate\n");
            $self->_print(
                join(
                    ",\n",
                    map {
                        $node2num{ $_->id } = $nodecter;
                        sprintf( "\t\t%d %s", $nodecter++, $_->id )
                      }
                      grep { $_->is_Leaf } $tree->get_nodes
                ),
                "\n;\n"
            );
        }
        my @data = _write_tree_Helper( $tree->get_root_node, \%node2num );
        if ( $data[-1] !~ /\)$/ ) {
            $data[0] = "(" . $data[0];
            $data[-1] .= ")";
        }

        # by default all trees in bioperl are currently rooted
        # something we'll try and fix one day....
        $self->_print(
            sprintf(
                "\t tree %s = [&%s] %s;\n",
                ( $tree->id || sprintf( "Bioperl_%d", $first + 1 ) ),
                ( $tree->get_root_node ) ? 'R' : 'U',
                join( ',', @data )
            )
        );
        $first++;
    }
    $self->_print("End;\n");
    $self->flush if $self->_flush_on_write && defined $self->_fh;
    return;
}

sub _write_tree_Helper {
    my ( $node, $node2num ) = @_;
    return () if ( !defined $node );
    my @data;

    foreach my $n ( $node->each_Descendent() ) {
        push @data, _write_tree_Helper( $n, $node2num );
    }
    if ( @data > 1 ) {		# internal node
        $data[0] = "(" . $data[0];
        $data[-1] .= ")";

        # FigTree comments start
        my $comment_flag;
        $comment_flag = 0
            if ( $node->has_tag('color') or  $node->has_tag('label') );
    
        $data[-1] .= '[&!' if defined $comment_flag;
    
            if ( $node->has_tag('color')) {
            my $color = $node->get_tag_values('color');
                $data[-1] .= "color=$color";
            $comment_flag++;
            }
            if ( $node->has_tag('label')) {
            my $label = $node->get_tag_values('label');
            $data[-1] .= ',' if $comment_flag;
                $data[-1] .= 'label="'. $label. '"';
            }
        $data[-1] .= ']' if defined $comment_flag;
        # FigTree comments end
        
        # let's explicitly write out the bootstrap if we've got it
        my $b;

        my $bl = $node->branch_length;
        if ( !defined $bl ) {
        }
        elsif ( $bl =~ /\#/ ) {
            $data[-1] .= $bl;
        }
        else {
            $data[-1] .= ":$bl";
        }
        if ( defined( $b = $node->bootstrap ) ) {
            $data[-1] .= sprintf( "[%s]", $b );
        }
        elsif ( defined( $b = $node->id ) ) {
            $b = $node2num->{$b} if ( $node2num->{$b} );    # translate node2num
            $data[-1] .= sprintf( "[%s]", $b ) if defined $b;
        }

    }
    else {			# leaf node
        if ( defined $node->id || defined $node->branch_length ) {
            my $id = defined $node->id ? $node->id : '';
            if ( length($id) && $node2num->{$id} ) {
                $id = $node2num->{$id};
            }
	    if ( $node->has_tag('color')) {
		my ($color) = $node->get_tag_values('color');
		$id .= "[&!color=$color\]";
	    }
            push @data,
              sprintf( "%s%s",
                $id,
                defined $node->branch_length
                ? ":" . $node->branch_length
                : '' );
        }
    }
    return @data;
}

=head2 header

 Title   : header
 Usage   : $obj->header($newval)
 Function: 
 Example : 
 Returns : value of header (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub header {
    my $self = shift;

    return $self->{'header'} = shift if @_;
    return $self->{'header'};
}

=head2 translate_node

 Title   : translate_node
 Usage   : $obj->translate_node($newval)
 Function: 
 Example : 
 Returns : value of translate_node (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub translate_node {
    my $self = shift;

    return $self->{'translate_node'} = shift if @_;
    return $self->{'translate_node'};
}

1;
