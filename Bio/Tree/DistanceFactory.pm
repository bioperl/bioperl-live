# $Id$
#
# BioPerl module for Bio::Tree::DistanceFactory
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
  my $tfactory = Bio::Tree::DistanceFactory->new(-method => "UPGMA");
  my $stats    = Bio::Align::DNAStatistics->new();

  my $alnin    = Bio::AlignIO->new(-format => 'clustalw',
                                   -file   => 'file.aln');
  my $aln = $alnin->next_aln;
  my $jcmatrix = $stats->distance(-align => $aln, 
                                  -method => 'Jukes-Cantor');
  my $tree = $tfactory->make_tree($jcmatrix);


=head1 DESCRIPTION

This is a factory which will construct a phylogenetic tree based on
the pairwise sequence distances for a set of sequences.  Currently
UPGMA (Sokal and Rolf) and NJ (Saitou and Nei) tree construction
methods are implemented.


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

Email jason-at-bioperl.org

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Tree::DistanceFactory;
use vars qw(@ISA $DefaultMethod);
use strict;

$DefaultMethod = 'UPGMA';

use Bio::Root::Root;
use Bio::Tree::Node;
use Bio::Tree::Tree;

@ISA = qw(Bio::Root::Root );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tree::DistanceFactory();
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
       return undef;
   }

   my $method = uc ($self->method);
   if( $method eq 'NJ' ) {
       return $self->_nj($matrix);
   } elsif( $method eq 'UPGMA' ) {
       return $self->_upgma($matrix);
   } else { 
       $self->warn("Unknown tree construction method '$method'.  Cannot run.");
       return undef;
   }
   
}


=head2 _nj

 Title   : _nj
 Usage   : my $tree = $disttreefact->_nj($matrix);
 Function: Construct a tree based on distance matrix using the 
           Neighbor Joining algorithm (Saitou and Nei, 1987)
 Returns : L<Bio::Tree::TreeI>
 Args    : L<Bio::Matrix::MatrixI> object


=cut

sub _nj {
   my ($self,$matrix) = @_;

   $self->throw("Not currently implemented - this is dev code, be patient!");

   # we assume type checking of $aln has already been done
   # client shouldn't be calling this directly anyways, using the
   # make_tree method is preferred
   
   # This implementation is highly influenced by the PHYLIP
   # neighbor.c code
   my @names =  @{$matrix->names || []};
   my $sp_count = scalar @names;
   if( $sp_count < 3 ) {
       $self->warn("Can only perform NJ treebuilding on sets of 3 or more species\n");
       return undef;
   }
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
   my ($i,$j,$x,$y,@dmat,@orig,@nodes);

   my @names = @{$distmat->names || []};
   my $c = 0;
   my @clusters = map { 
       my $r = { 'id'        => $c,
		 'height'    => 0,
		 'contains'  => [$c],
	     };
       $c++;
       $r;
   } @names;
   my $n = scalar @clusters;
   my $max = 2*$n - 1;
   my $K = $n;
   my (@mins,$min);
   for ( $i = 0; $i < $n; $i++ ) {
       for( $j = $i+1; $j < $n; $j++ ) {
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
	   $nodes[$nid]->branch_length($min/2 - $clusters[$cid]->{'height'});
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
   my @names = @{$matrix->names || {}};
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
