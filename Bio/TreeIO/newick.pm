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

Bio::TreeIO::newick - TreeIO implementation for parsing 
  Newick/New Hampshire/PHYLIP format.

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

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::TreeIO::newick;
use vars qw($DefaultBootstrapStyle);
use strict;

use Bio::Event::EventGeneratorI;

#initialize some package variables, could use 'our' but fails in perl < 5.6

$DefaultBootstrapStyle = 'traditional';
use base qw(Bio::TreeIO);


=head2 new

 Title   : new
 Args    : -print_count     => boolean  default is false
           -bootstrap_style => set the bootstrap style (one of nobranchlength,
							molphy, traditional)
           -order_by        => set the order by sort method 
                               (see L<Bio::Node::Node::each_Descendent()> )

=cut

sub _initialize { 
    my $self = shift;
    $self->SUPER::_initialize(@_);
    my ($print_count,$style,$order_by) = $self->_rearrange([qw(PRINT_COUNT 
							       BOOTSTRAP_STYLE
							       ORDER_BY)],
					  @_);
    $self->print_tree_count($print_count || 0);
    $self->bootstrap_style($style || $DefaultBootstrapStyle);
    $self->order_by($order_by) if defined $order_by;
    return;
}


=head2 next_tree

 Title   : next_tree
 Usage   : my $tree = $treeio->next_tree
 Function: Gets the next tree in the stream
 Returns : L<Bio::Tree::TreeI>
 Args    : none


=cut

sub next_tree{
   my ($self) = @_;
   local $/ = ";\n";
   return unless $_ = $self->_readline;
   s/[\r\n]//gs;
   my $score;
   my $despace = sub {my $dirty = shift; $dirty =~ s/\s+//gs; return $dirty};
   my $dequote = sub {my $dirty = shift; $dirty =~ s/^"?\s*(.+?)\s*"?$/$1/; return $dirty};
   s/([^"]*)(".+?")([^"]*)/$despace->($1) . $dequote->($2) . $despace->($3)/egsx;
   if( s/^\s*\[([^\]]+)\]// ) {
       my $match = $1;
       $match =~ s/\s//g;
       $match =~ s/lh\=//;
       if( $match =~ /([-\d\.+]+)/ ) {
	   $score = $1;
       }
   }

   $self->debug("entry is $_\n");
#   my $empty = chr(20);
 
   # replace empty labels with a tag
#   s/\(,/\($empty,/ig;
#   s/,,/,$empty,/ig;
#   s/,,/,/ig;
#   s/,\)/,$empty\)/ig;
#   s/\"/\'/ig;

   my $chars = '';
   $self->_eventHandler->start_document;
   my ($prev_event,$lastevent,$id) = ('','','');
   foreach my $ch ( split(//,$_) ) {
       if( $ch eq ';' ) {
	   my $tree = $self->_eventHandler->end_document($chars);
	   $tree->score($score) if defined $score;
	   if( $self->internal_node_id eq 'bootstrap' ) {
	       $tree->move_id_to_bootstrap;
	   }
	   return $tree;
       } elsif( $ch eq '(' ) {
	   $chars = '';
	   $self->_eventHandler->start_element( {'Name' => 'tree'} );
       } elsif($ch eq ')' ) {
	   if( length($chars) ) {
	       if( $lastevent eq ':' ) {
		   $self->_eventHandler->start_element( { 'Name' => 'branch_length'});
		   $self->_eventHandler->characters($chars);
		   $self->_eventHandler->end_element( {'Name' => 'branch_length'});
		   $lastevent = $prev_event;
	       } else { 
		   $self->debug("internal node, id with no branchlength is $chars\n");
		   $self->_eventHandler->start_element( { 'Name' => 'node' } );
		   $self->_eventHandler->start_element( { 'Name' => 'id' } );
		   $self->_eventHandler->characters($chars);
		   $self->_eventHandler->end_element( { 'Name' => 'id' } );
		   $id = $chars;
	       }
	       my $leafstatus = 0;
	       if( $lastevent ne ')' ) {
		   $leafstatus = 1;
	       }

	       $self->_eventHandler->start_element({'Name' => 'leaf'});
	       $self->_eventHandler->characters($leafstatus);
	       $self->_eventHandler->end_element({'Name' => 'leaf'});
	       $id = '';
	   } else {
	       $self->_eventHandler->start_element( {'Name' => 'node'} );
	   }

 	   $self->_eventHandler->end_element( {'Name' => 'node'} );
	   $self->_eventHandler->end_element( {'Name' => 'tree'} );
	   $chars = '';
       } elsif ( $ch eq ',' ) {
	   if( length($chars) ) {
	       if( $lastevent eq ':' ) {
		   $self->_eventHandler->start_element( { 'Name' => 'branch_length'});
		   $self->_eventHandler->characters($chars);
		   $self->_eventHandler->end_element( {'Name' => 'branch_length'});
		   $lastevent = $prev_event;
		   $chars = '';		   
	       } else { 
		   $self->debug("leaf id with no branchlength is $chars\n");
		   $self->_eventHandler->start_element( { 'Name' => 'node' } );
		   $self->_eventHandler->start_element( { 'Name' => 'id' } );
		   $self->_eventHandler->characters($chars);
		   $self->_eventHandler->end_element( { 'Name' => 'id' } );
		   $id = $chars;
	       }
	   } else {
	       $self->_eventHandler->start_element( { 'Name' => 'node' } );
	   }
	   my $leafstatus = 0;
	   if( $lastevent ne ')' ) {
	       $leafstatus = 1;
	   }
	   $self->_eventHandler->start_element({'Name' => 'leaf'});
	   $self->_eventHandler->characters($leafstatus);
	   $self->_eventHandler->end_element({'Name' => 'leaf'});
	   $self->_eventHandler->end_element( {'Name' => 'node'} );
	   $chars = '';
	   $id    = '';
       } elsif( $ch eq ':' ) {
	   $self->debug("id with a branchlength coming is $chars\n");
	   $self->_eventHandler->start_element( { 'Name' => 'node' } );
	   $self->_eventHandler->start_element( { 'Name' => 'id' } );	   
	   $self->_eventHandler->characters($chars);
	   $self->_eventHandler->end_element( { 'Name' => 'id' } );	   
	   $id = $chars;
	   $chars = '';
       } else { 	   
	   $chars .= $ch;
	   next;
       }
       $prev_event = $lastevent;
       $lastevent = $ch;
   }
   my $tree = $self->_eventHandler->end_document($chars);
   return $tree if $tree;
   return;
}

=head2 write_tree

 Title   : write_tree
 Usage   : $treeio->write_tree($tree);
 Function: Write a tree out to data stream in newick/phylip format
 Returns : none
 Args    : L<Bio::Tree::TreeI> object

=cut

sub write_tree{
   my ($self,@trees) = @_;  
   my $orderby = $self->order_by;
   my $bootstrap_style = $self->bootstrap_style;
   if( $self->print_tree_count ){ 
       $self->_print(sprintf(" %d\n",scalar @trees));
   }
   my $nl = $self->newline_each_node;
   foreach my $tree( @trees ) {
       
       if( ! defined $tree || ref($tree) =~ /ARRAY/i ||
	   ! $tree->isa('Bio::Tree::TreeI') ) {
	   $self->throw("Calling write_tree with non Bio::Tree::TreeI object\n");
       }
       my @data = _write_tree_Helper($tree->get_root_node,
				     $bootstrap_style,
				     $orderby,
				     $nl);
       if( $nl ) {
	   chomp($data[-1]);# remove last newline
	   $self->_print(join(",\n", @data), ";\n");
       } else {
	   $self->_print(join(',', @data), ";\n");
       }
   }
   $self->flush if $self->_flush_on_write && defined $self->_fh;
   return;
}

sub _write_tree_Helper {
    my ($node,$style,$orderby,$nl) = @_;
    $style = '' unless defined $style;
    return () if (!defined $node);

    my @data;
    foreach my $n ( $node->each_Descendent($orderby) ) {
	push @data, _write_tree_Helper($n,$style,$orderby,$nl);
    }
    
    # let's explicitly write out the bootstrap if we've got it
    my $id = $node->id_output;
    my $bs = $node->bootstrap; # bs better not have any spaces?
    $bs =~ s/\s+//g if defined $bs;
    my $bl = $node->branch_length;
    if( @data ) {
	if( $nl ) {
	    $data[0] = "(\n" . $data[0];
	    $data[-1] .= ")\n";
	} else {
	    $data[0] = "(" . $data[0];
	    $data[-1] .= ")";
	}

	if( $node->is_Leaf ) { 
	    $node->debug("node is a leaf!  This is unexpected...");

	    $id ||= '';
	    if( ! defined $bl || ! length($bl) ||
		($style && $style =~ /nobranchlength/i) ) {
		$data[-1] .= $id;
	    } elsif( defined $bl && length($bl) ) { 
		$data[-1] .= "$id:$bl";
	    } else { 
		$data[-1] .= $id;
	    }
	} else { 
	    if( ! defined $bl || ! length($bl) ||
		($style && $style =~ /nobranchlength/i) ) {
		
		if( defined $id || defined $bs ) {
		    $data[-1] .= defined $bs ? $bs : $id;
		}
	    } elsif( $style =~ /molphy/i ) {
		if( defined $id ) {
		    $data[-1] .= $id;
		}
		if( $bl =~ /\#/) {
		    $data[-1] .= $bl;
		} else { 
		    $data[-1] .= ":$bl";
		}
		if( defined $bs ) { 
		    $data[-1] .= "[$bs]";
		}
	    } else {
		# traditional style of 
		# ((A:1,B:2)81:3);   where 3 is internal node branch length
		#                    and 81 is bootstrap/node label
		if( defined $bs || defined $id ) {
		    $data[-1] .= defined $bs ? "$bs:$bl" : "$id:$bl";
		} elsif( $bl =~ /\#/ ) {
		    $data[-1] .= $bl;
		} else { 
		    $data[-1] .= ":$bl"; 
		}
	    }
	}
    } elsif( defined $id || defined $bl ) {
	my $str;
	$id ||= '';
	if( ! defined $bl || ! length($bl) ||
	    ($style && $style =~ /nobranchlength/i) ) {
	    $str = $id;
	} elsif( defined $bl && length($bl) ) { 
	    $str = "$id:$bl";
	} else { 
	    $str = $id;
	}
	push @data, $str;
    }
    return @data;
}

=head2 print_tree_count

 Title   : print_tree_count
 Usage   : $obj->print_tree_count($newval)
 Function: Get/Set flag for printing out the tree count (paml,protml way)
 Returns : value of print_tree_count (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub print_tree_count{
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

sub bootstrap_style{
    my $self = shift;
    my $val = shift;
    if( defined $val ) {

	if( $val !~ /^nobranchlength|molphy|traditional/i ) {
	    $self->warn("requested an unknown bootstrap style $val, expect one of nobranchlength,molphy,traditional, not updating value.  Default is $DefaultBootstrapStyle\n");
	} else { 
	    $self->{'_bootstrap_style'} = $val;
	}
    }
    return $self->{'_bootstrap_style'} || $DefaultBootstrapStyle;
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
