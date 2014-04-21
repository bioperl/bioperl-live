#
# BioPerl module for Bio::Tree::DistanceFactory
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::DistanceFactory - Construct a tree using distance based methods

=head1 SYNOPSIS

  use Bio::Tree::DistanceFactory;
  use Bio::AlignIO;
  use Bio::Align::DNAStatistics;
  my $tfactory = Bio::Tree::DistanceFactory->new(-method => "NJ");
  my $stats    = Bio::Align::DNAStatistics->new();

  my $alnin    = Bio::AlignIO->new(-format => 'clustalw',
                                   -file   => 'file.aln');
  my $aln = $alnin->next_aln;
  # Of course matrix can come from a different place
  # like PHYLIP if you prefer, Bio::Matrix::IO should be able
  # to parse many things
  my $jcmatrix = $stats->distance(-align => $aln, 
                                  -method => 'Jukes-Cantor');
  my $tree = $tfactory->make_tree($jcmatrix);


=head1 DESCRIPTION

This is a factory which will construct a phylogenetic tree based on
the pairwise sequence distances for a set of sequences.  Currently
UPGMA (Sokal and Michener 1958) and NJ (Saitou and Nei 1987) tree
construction methods are implemented.

=head1 REFERENCES

Eddy SR, Durbin R, Krogh A, Mitchison G, (1998) "Biological Sequence Analysis",
Cambridge Univ Press, Cambridge, UK.

Howe K, Bateman A, Durbin R, (2002) "QuickTree: building huge
Neighbour-Joining trees of protein sequences." Bioinformatics
18(11):1546-1547.

Saitou N and Nei M, (1987) "The neighbor-joining method: a new method
for reconstructing phylogenetic trees." Mol Biol Evol 4(4):406-25.

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
of the bugs and their resolution. Bug reports can be submitted the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Tree::DistanceFactory;
use vars qw($DefaultMethod $Precision);
use strict;

# some defaults
$DefaultMethod = 'UPGMA';
$Precision = 5;

use Bio::Tree::Node;
use Bio::Tree::Tree;

use base qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tree::DistanceFactory->new();
 Function: Builds a new Bio::Tree::DistanceFactory object 
 Returns : an instance of Bio::Tree::DistanceFactory
 Args    : -method => 'NJ' or 'UPGMA'


=cut

sub new {
  my($class,@args) = @_;  
  my $self = $class->SUPER::new(@args);

  my ($method) = $self->_rearrange([qw(METHOD)],
				   @args);
  $self->method($method || $DefaultMethod);
  return $self;
}

=head2 make_tree

 Title   : make_tree
 Usage   : my $tree = $disttreefact->make_tree($matrix);
 Function: Build a Tree based on a distance matrix
 Returns : L<Bio::Tree::TreeI>
 Args    : L<Bio::Matrix::MatrixI> object


=cut

sub make_tree{
   my ($self,$matrix) = @_;
   if( ! defined $matrix || !ref($matrix) || 
       ! $matrix->isa('Bio::Matrix::MatrixI') ) {
       $self->warn("Need to provide a valid Bio::Matrix::MatrixI object to make_tree");
       return;
   }

   my $method = uc ($self->method);
   if( $method =~ /NJ/i ) {
       return $self->_nj($matrix);
   } elsif( $method =~ /UPGMA/i ) {
       return $self->_upgma($matrix);
   } else { 
       $self->warn("Unknown tree construction method '$method'.  Cannot run.");
       return;
   }
   
}


=head2 _nj

 Title   : _nj
 Usage   : my $tree = $disttreefact->_nj($matrix);
 Function: Construct a tree based on distance matrix using the 
           Neighbor Joining algorithm (Saitou and Nei, 1987)
           Implementation based on Kevin Howe's Quicktree implementation
           and uses his tricks (some based on Bill Bruno's work) to eliminate
           negative branch lengths
 Returns : L<Bio::Tree::TreeI>
 Args    : L<Bio::Matrix::MatrixI> object

=cut

sub _nj {
   my ($self,$distmat) = @_;

   # we assume type checking of $aln has already been done
   # client shouldn't be calling this directly anyways, using the
   # make_tree method is preferred
   
   # so that we can trim the number of digits shown as the branch length
   my $precisionstr = "%.$Precision"."f";

   my @names =  $distmat->column_names;
   my $N = scalar @names;
   my ($i,$j,$m,@nodes,$mat,@r);
   my $L = $N;

   if( $N < 2 ) {
       $self->warn("Can only perform NJ treebuilding on sets of 2 or more species\n");
       return;
   } elsif( $N == 2 ) {
       $i = 0;
       my $d = sprintf($precisionstr,
		       $distmat->get_entry($names[0],$names[1]) / 2);
       my $root = Bio::Tree::Node->new();
       for my $nm ( @names ) {
	   $root->add_Descendents( Bio::Tree::Node->new(-id => $nm,
							-branch_length => $d));
       }
       return Bio::Tree::Tree(-root => $root);
   }
   my $c = 0;
   
   for ( $i = 0; $i < $N; $i++ ) {
       push @nodes, Bio::Tree::Node->new(-id => $names[$i]);
       my $ri = 0;
       for( $j = 0; $j < $N; $j++ ) {
	   $mat->[$i][$j] = $distmat->get_entry($names[$i],$names[$j]);
	   $ri += $mat->[$i][$j];
       }
       $r[$i] = $ri / ($L -2);
   }
   
   for( my $nodecount = 0; $nodecount < $N-3; $nodecount++) {
       my ($mini,$minj,$min);
       for($i = 0; $i < $N; $i++ ) {
	   next unless defined $nodes[$i];
	   for( $j = 0; $j < $i; $j++ ) {
	       next unless defined $nodes[$j];
	       my $dist = $mat->[$i][$j] - ($r[$i] + $r[$j]);
	       if( ! defined $min ||
		   $dist <= $min) {
		   ($mini,$minj,$min) = ($i,$j,$dist);
	       }
	   }
       }
       my $dij    = $mat->[$mini][$minj];
       my $dist_i = ($dij + $r[$mini] - $r[$minj]) / 2;
       my $dist_j = $dij - $dist_i;
       
       # deal with negative branch lengths
       # per code in K.Howe's quicktree
       if( $dist_i < 0 ) {
	   $dist_i = 0;
	   $dist_j = $dij;
	   $dist_j = 0 if( $dist_j < 0 );
       } elsif( $dist_j < 0 ) { 
	   $dist_j = 0;
	   $dist_i = $dij;
	   $dist_i = 0 if( $dist_i < 0 );
       }
       
       $nodes[$mini]->branch_length(sprintf($precisionstr,$dist_i));
       $nodes[$minj]->branch_length(sprintf($precisionstr,$dist_j));
       
       my $newnode = Bio::Tree::Node->new(-descendents => [ $nodes[$mini],
							    $nodes[$minj] ]);

       $nodes[$mini] = $newnode;
       delete $nodes[$minj];
       
       # update the distance matrix
       $r[$mini] = 0;
       my ($dmi,$dmj);
       for( $m = 0; $m < $N; $m++ ) {	   
	   next unless defined $nodes[$m];
	   if( $m != $mini ) {
	       $dmj = $mat->[$m][$minj];
	       
	       my ($row,$col);
	       ($row,$col) = ($m,$mini);
	       $dmi = $mat->[$row][$col];
	       
	       # from K.Howe's notes in quicktree
	       # we can actually adjust r[m] here, by using the form:
	       # rm = ((rm * numseqs) - dmi - dmj + dmk) / (numseqs-1)

	       # Note: in Bill Bruno's method for negative branch
	       # elimination, then if either dist_i is positive and
	       # dist_j is 0, or dist_i is zero and dist_j is positive
	       # (after adjustment) then the matrix entry is formed
	       # from the distance to the node in question (m) to the
	       # node with the zero branch length (whichever it was).
	       # I think my code already has the same effect; this is
	       # certainly true if dij is equal to dist_i + dist_j,
	       # which it should have been fixed to

	       my $dmk = $mat->[$row][$col] = $mat->[$col][$row] = 
		   ($dmi + $dmj - $dij) / 2;
	       
	       # If we don't want to try and correct negative brlens
	       # this is essentially what is in Edddy et al, BSA book.
	       # $r[$m] = (($r[$m] * $L) - $dmi - $dmj + $dmk) / ($L-1);
	       # 
	       $r[$m] = (($r[$m] * ($L - 2)) - $dmi - $dmj + 
			 $mat->[$row][$col]) / ( $L - 3);
	       $r[$mini] += $dmk;
	   }
       }
       $L--;
       $r[$mini] /= $L - 2;
   }
   
   # should be 3 nodes left
   my (@leftovernodes,@leftovers);
   for( my $k = 0; $k < $N; $k++ ) {
       if( defined $nodes[$k] ) {
	   push @leftovers, $k;
	   push @leftovernodes, $nodes[$k];
       }
   }
   my ($l_0,$l_1,$l_2) = @leftovers;
   
   my $dist_i = ( $mat->[$l_1][$l_0] + $mat->[$l_2][$l_0] -
		  $mat->[$l_2][$l_1] ) / 2;
   
   my $dist_j = ( $mat->[$l_1][$l_0] - $dist_i);
   my $dist_k = ( $mat->[$l_2][$l_0] - $dist_i);

   # This is Kev's code to get rid of negative branch lengths
   if( $dist_i < 0 ) { 
       $dist_i = 0;
       $dist_j = $mat->[$l_1][$l_0];
       $dist_k = $mat->[$l_2][$l_0];
       if( $dist_j < 0 ) { 
	   $dist_j = 0;
	   $dist_k = ( $mat->[$l_2][$l_0] + $mat->[$l_2][$l_1] ) / 2;
	   $dist_k = 0 if( $dist_k < 0 );
       } elsif( $dist_k < 0 ) {
	   $dist_k = 0;
	   $dist_j = ($mat->[$l_1][$l_0] + $mat->[$l_2][$l_1]) / 2;
	   $dist_j = 0 if( $dist_j < 0 );
       }
   } elsif( $dist_j < 0 ) {
       $dist_j = 0;
       $dist_i = $mat->[$l_1][$l_0];
       $dist_k = $mat->[$l_2][$l_1];
       if( $dist_i < 0 ) { 
	   $dist_i = 0;
	   $dist_k = ( $mat->[$l_2][$l_0] + $mat->[$l_2][$l_1]) / 2;
	   $dist_k = 0 if( $dist_k  < 0 );
       } elsif( $dist_k < 0 ) { 
	   $dist_k = 0;
	   $dist_i = ( $mat->[$l_1][$l_0] + $mat->[$l_2][$l_0]) / 2;
	   $dist_i = 0 if( $dist_i < 0 );
       }
   } elsif( $dist_k < 0 ) {
       $dist_k = 0;
       $dist_i = $mat->[$l_2][$l_0];
       $dist_j = $mat->[$l_2][$l_1];
       if( $dist_i < 0 ) { 
	   $dist_i = 0;
	   $dist_j = ( $mat->[$l_1][$l_0] + $mat->[$l_2][$l_1] ) / 2;
	   $dist_j = 0 if $dist_j < 0;
       } elsif( $dist_j < 0  ) {
	   $dist_j = 0;
	   $dist_i = ($mat->[$l_1][$l_0] + $mat->[$l_2][$l_0]) / 2;
	   $dist_i = 0 if $dist_i < 0;
       }
   }
   $leftovernodes[0]->branch_length(sprintf($precisionstr,$dist_i));
   $leftovernodes[1]->branch_length(sprintf($precisionstr,$dist_j));
   $leftovernodes[2]->branch_length(sprintf($precisionstr,$dist_k));

   Bio::Tree::Tree->new(-root => Bio::Tree::Node->new
			(-descendents => \@leftovernodes));
}

=head2 _upgma

 Title   : _upgma
 Usage   : my $tree = $disttreefact->_upgma($matrix);
 Function: Construct a tree based on alignment using UPGMA
 Returns : L<Bio::Tree::TreeI>
 Args    : L<Bio::Matrix::MatrixI> object


=cut

sub _upgma{
   my ($self,$distmat) = @_;
   # we assume type checking of $matrix has already been done
   # client shouldn't be calling this directly anyways, using the
   # make_tree method is preferred
   
   # algorithm, from Eddy, Durbin, Krogh, Mitchison, 1998
   # originally by Sokal and Michener 1956

   my $precisionstr = "%.$Precision"."f";
   
   my ($i,$j,$x,$y,@dmat,@orig,@nodes);

   my @names = $distmat->column_names;
   my $c = 0;
   my @clusters = map { 
       my $r = { 'id'        => $c,
		 'height'    => 0,
		 'contains'  => [$c],
	     };
       $c++;
       $r;
   } @names;

   my $K = scalar @clusters;
   my (@mins,$min);
   for ( $i = 0; $i < $K; $i++ ) {
       for( $j = $i+1; $j < $K; $j++ ) {
	   my $d =  $distmat->get_entry($names[$i],$names[$j]);
	   # get Min here on first time around, save 1 cycle
	   $dmat[$j][$i] = $dmat[$i][$j] = $d;
	   $orig[$i][$j] = $orig[$j][$i] = $d;
	   if ( ! defined $min || $d <= $min ) {
	       if( defined $min && $min == $d ) { 
		   push @mins, [$i,$j];
	       } else { 
		   @mins = [$i,$j];
		   $min  = $d;
	       }
	   }
       }
   }
   # distance between each cluster is avg distance
   # between pairs of sequences from each cluster
   while( $K > 1 ) {       
       # fencepost - we already have found the $min
       # so very first time loop is executed we can skip checking
       unless( defined $min ) {
	   for($i = 0; $i < $K; $i++ ) {
	       for( $j = $i+1; $j < $K; $j++ ) {
		   my $dij = $dmat[$i][$j];
		   if( ! defined $min ||
		       $dij <= $min) {
		       if( defined $min &&
			   $min == $dij ) { 
			   push @mins, [$i,$j];
		       } else { 
			   @mins = [ $i,$j ];
			   $min = $dij;
		       }
		   }
	       }
	   }
       }
       # randomly break ties
       ($x,$y) = @{ $mins[int(rand(scalar @mins))] };   

       # now we are going to join clusters x and y, make a new cluster

       my $node = Bio::Tree::Node->new();   
       my @subids;
       for my $cid ( $x,$y ) {
	   my $nid = $clusters[$cid]->{'id'};
	   if( ! defined $nodes[$nid] ) {
	       $nodes[$nid] = Bio::Tree::Node->new(-id => $names[$nid]);
	   }
	   $nodes[$nid]->branch_length
	       (sprintf($precisionstr,$min/2 - $clusters[$cid]->{'height'}));
	   $node->add_Descendent($nodes[$nid]);
	   push @subids, @{ $clusters[$cid]->{'contains'} };
       }
       my $cluster = { 'id'       => $c++,
		       'height'   => $min / 2,
		       'contains' => [@subids],
		   };

       $K--; # we are going to drop the last node so go ahead and decrement K
       $nodes[$cluster->{'id'}] = $node;
       if ( $y != $K ) {
	   $clusters[$y] = $clusters[$K];
	   $dmat[$y] = $dmat[$K];
	   for ( $i = 0; $i < $K; $i++ ) {
	       $dmat[$i][$y] = $dmat[$y][$i];
	   }
       }
       delete $clusters[$K];
       $clusters[$x] = $cluster;
       # now recalculate @dmat
       for( $i = 0; $i < $K; $i++ ) {	   
	   if( $i != $x) {
	       $dmat[$i][$x] = $dmat[$x][$i] = 
		   &_upgma_distance($clusters[$i],$clusters[$x],\@orig);
	   } else { 
	       $dmat[$i][$i] = 0;
	   }
       }
       # reset so next loop iteration
       # we will find minimum distance
       @mins = ();
       $min = undef;
   }
   Bio::Tree::Tree->new(-root => $nodes[-1]);
}

# calculate avg distance between clusters - be they
# single sequences or the combination of multiple seqences
# $cluster_i and $cluster_j are the clusters to operate on
# and $distances is a matrix (arrayref of arrayrefs) of pairwise 
# differences indexed on the sequence ids - 
# so $distances->[0][1] is the distance between sequences 0 and 1

sub _upgma_distance { 
    my ($cluster_i, $cluster_j, $distances) = @_;
    my $ilen = scalar @{ $cluster_i->{'contains'} };
    my $jlen = scalar @{ $cluster_j->{'contains'} };
    my ($d,$count);
    for( my $i = 0; $i < $ilen; $i++ ) {
	my $i_id = $cluster_i->{'contains'}->[$i];
	for( my $j = 0; $j < $jlen; $j++) {	    
	    my $j_id = $cluster_j->{'contains'}->[$j];
	    if( ! defined $distances->[$i_id][$j_id] ) {
		warn("no value for $i_id $j_id\n");
	    } else { 
		$d += $distances->[$i_id][$j_id];
	    }
	    $count++;
	}
    }
    return $d / $count;
}

=head2 method

 Title   : method
 Usage   : $obj->method($newval)
 Function: 
 Example : 
 Returns : value of method (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub method{
    my $self = shift;
    return $self->{'_method'} = shift if @_;
    return $self->{'_method'};
}


=head2 check_additivity

 Title     : check_additivity
 Usage     : if( $distance->check_additivity($matrix) ) {
             }
 Function  : See if matrix obeys additivity principal
 Returns   : boolean
 Args      : Bio::Matrix::MatrixI 
 References: Based on a Java implementation by
             Peter Sestoft, sestoft@dina.kvl.dk 1999-12-07 version 0.3
             http://www.dina.kvl.dk/~sestoft/bsa.html
             which in turn is based on algorithms described in 
             R. Durbin, S. Eddy, A. Krogh, G. Mitchison. 
             Biological Sequence Analysis CUP 1998, Chapter 7.

=cut

sub check_additivity{
   my ($self,$matrix) = @_;
   my @names = $matrix->column_names;
   my $len = scalar @names;
   return unless $len >= 4;
   # look at all sets of 4
   for( my $i = 0; $i < $len; $i++ ) { 
       for( my $j = $i+1; $j< $len; $j++) {
	   for( my $k = $j+1; $k < $len; $k ++ ) {
	       for( my $m = $k +1; $m < $len; $m++ ) {
		   my $DijDkm = $matrix->get_entry($names[$i],$names[$j]) + 
		       $matrix->get_entry($names[$k],$names[$m]);
		   my $DikDjm = $matrix->get_entry($names[$i],$names[$k]) + 
		       $matrix->get_entry($names[$j],$names[$m]);
		   my $DimDjk = $matrix->get_entry($names[$i],$names[$m]) + 
		       $matrix->get_entry($names[$j],$names[$k]);
		   if( !( ( $DijDkm == $DikDjm && $DijDkm >= $DimDjk)
			  || ( $DijDkm == $DimDjk && $DijDkm >= $DikDjm)
			  || ( $DikDjm == $DimDjk && $DikDjm >= $DijDkm) )) {
		       return 0;
		   }
	       }
	   }
       } 
   }
   return 1;
}

1;
