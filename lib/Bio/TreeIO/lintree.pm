#
# BioPerl module for Bio::TreeIO::lintree
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

Bio::TreeIO::lintree - Parser for lintree output trees

=head1 SYNOPSIS

  # do not use directly, use through Bio::TreeIO
  use Bio::TreeIO;
  my $treeio = Bio::TreeIO->new(-format => 'lintree',
                               -file   => 't/data/crab.nj');
  my $tree = $treeio->next_tree;

=head1 DESCRIPTION

Parser for the lintree output which looks like this

  13 sequences     1000 bootstraping
1 A-salina
2 C-vittat
3 C-sp.
4 L-aequit
5 P-camtsc
6 E-tenuim
7 L-splend
8 P-bernha
9 P-acadia
10 P-p(NE)
11 P-p(GU)
12 P-l(NE)
13 P-l(GU)
 14 and   2        0.098857      1000
 14 and   3        0.127932      1000
 15 and   1        0.197471      1000
 15 and  14        0.029273       874
 16 and  10        0.011732      1000
 16 and  11        0.004529      1000
 17 and  12        0.002258      1000
 17 and  13        0.000428      1000
 18 and  16        0.017512      1000
 18 and  17        0.010824       998
 19 and   4        0.006534      1000
 19 and   5        0.006992      1000
 20 and  15        0.070461      1000
 20 and  18        0.030579       998
 21 and   8        0.003339      1000
 21 and   9        0.002042      1000
 22 and   6        0.011142      1000
 22 and  21        0.010693       983
 23 and  20        0.020714       996
 23 and  19        0.020350      1000
 24 and  23        0.008665       826
 24 and  22        0.013457       972
 24 and   7        0.025598      1000

See http://www.bio.psu.edu/People/Faculty/Nei/Lab/software.htm for access
to the program and N Takezaki, A Rzhetsky, and M Nei, "Phylogenetic test
of the molecular clock and linearized trees." Mol Biol Evol 12(5):823-33.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 CONTRIBUTORS

Ideas and discussion from:
 Alan Christoffels
 Avril Coghlan

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::TreeIO::lintree;
use vars qw(%Defaults);
use strict;


use base qw(Bio::TreeIO);
$Defaults{'NodeType'} = "Bio::Tree::Node";

=head2 new

 Title   : new
 Usage   : my $obj = Bio::TreeIO::lintree->new();
 Function: Builds a new Bio::TreeIO::lintree object 
 Returns : an instance of Bio::TreeIO::lintree
 Args    : -nodetype => Node type to create [default Bio::Tree::Node]


=cut

sub _initialize { 
    my ($self,@args) = @_;
    $self->SUPER::_initialize(@args);
    my ($nodetype) = $self->_rearrange([qw(NODETYPE)],@args);
    $nodetype ||= $Defaults{'NodeType'};
    $self->nodetype($nodetype);
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
    my $seentop = 0;
    my ($tipcount,%data,@nodes) = (0);
    my $nodetype = $self->nodetype;   

    while( defined( $_ = $self->_readline) ) {
	if( /^\s*(\d+)\s+sequences/ox ) {
	    if( $seentop ) { 
		$self->_pushback($_);
		last;
	    }
	    $tipcount = $1;
	    $seentop = 1;
	} elsif( /^(\d+)\s+(\S+)\s*$/ox ) {
	    # deal with setting an outgroup
	    unless( defined $data{'outgroup'} ) {
		$data{'outgroup'} = [$1,$2];
	    }
	    $nodes[$1 - 1] = { '-id' => $2 }; 
	} elsif( m/^\s*(\d+)\s+and\s+(\d+)\s+(\-?\d+\.\d+)(?:\s+(\d+))?/ox ) {
	    my ($node,$descend,$blength,$bootstrap) = ( $1, $2, $3, $4 );
	    # need to -- descend and node because
	    # array is 0 based
	    $node--;$descend--;
	    $nodes[$descend]->{'-branch_length'} = $blength;
	    $nodes[$descend]->{'-bootstrap'}     = $bootstrap; #? here
	    $nodes[$node]->{'-id'} = $node+1;
	    push @{$nodes[$node]->{'-d'}}, $descend;
	    
	} elsif( /\s+(\S+)\-distance was used\./ox ) {
	    $data{'method'} = $1;
	} elsif( /\s*seed=(\d+)/ox ) {
	    $data{'seed'} = $1;
	} elsif( m/^outgroup:\s+(\d+)\s+(\S+)/ox ) {
	    $data{'outgroup'} = [$1,$2];
	}
    }
    if( @nodes ) {
	my @treenodes;
	foreach my $n ( @nodes ) { 	
	    push @treenodes, $nodetype->new(%{$n});
	}
	
	foreach my $tn ( @treenodes ) {
	    my $n = shift @nodes;
	    for my $ptr ( @{ $n->{'-d'} || [] } ) {
		$tn->add_Descendent($treenodes[$ptr]);
	    }
	}
	my $T = Bio::Tree::Tree->new(-root => (pop @treenodes) );
	if( $data{'outgroup'} ) {
	    my ($outgroup) = $treenodes[$data{'outgroup'}->[0]];
	    if( ! defined $outgroup) {
		$self->warn("cannot find '". $data{'outgroup'}->[1]. "'\n");
	    } else { 
		$T->reroot($outgroup->ancestor);
	    }
	}
	return $T;
    }
    return; # if there are no more trees, return undef
	
}

=head2 nodetype

 Title   : nodetype
 Usage   : $obj->nodetype($newval)
 Function: 
 Example : 
 Returns : value of nodetype (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub nodetype{
   my ($self,$value) = @_;
   if( defined $value) {
       eval "require $value";
       if( $@ ) { $self->throw("$@: Unrecognized Node type for ".ref($self). 
			       "'$value'");}
       
       my $a = bless {},$value;
       unless( $a->isa('Bio::Tree::NodeI')  ) {
	   $self->throw("Must provide a valid Bio::Tree::NodeI or child class to SeqFactory Not $value");
       }
      $self->{'nodetype'} = $value;
    }
   return $self->{'nodetype'};
}

1;
