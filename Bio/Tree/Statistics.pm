# $Id$
#
# BioPerl module for Bio::Tree::Statistics
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::Statistics - Calculate certain statistics for a coalescent population

=head1 SYNOPSIS

  use Bio::Tree::Statistics;
  use Bio::Tree::Tree;
  use Bio::Tree::RandomFactory;
  my $factory = new Bio::Tree::RandomFactory( -sample_size => 6,
					      -maxcount => 50);


  my $tree = $factory->next_tree;
  
  $factory->add_Mutations($tree,20);

  my $stats = new Bio::Tree::Statistics();
  my $pi = $stats->pi($tree);
  my $D = $stats->tajima_d($tree);

=head1 DESCRIPTION

This object is a place to accumulate routines for calculating various
statistics from the coalescent simulation.  This module will likely
move to Bio::PopGen as it is more appropriate there.

Currently implemented:
 Fu and Li's D  (fu_and_li_D)
 Fu and Li's D* (fu_and_li_D_star)
 Fu and Li's F  (fu_and_li_F)
 Tajima's D     (tajima_D)
 theta          (theta)
 pi --- not currently working. (pi)

References forthcoming.

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
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 CONTRIBUTORS

Matt Hahn E<lt>matthew.hahn@duke.dukeE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tree::Statistics;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;

@ISA = qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tree::Statistics();
 Function: Builds a new Bio::Tree::Statistics object 
 Returns : Bio::Tree::Statistics
 Args    :


=cut

=head2 fu_and_li_D

 Title   : fu_and_li_D
 Usage   : my $D = $statistics->fu_an_li_D($tree,$nummut);
 Function:
           For this we assume that the tree is made up of
           Bio::Tree::AlleleNode's which contain markers and alleles
           each marker is a 'mutation' 
 Returns : Fu and Li's D statistic for this Tree
 Args    : $tree - Bio::Tree::TreeI which contains Bio::Tree::AlleleNodes

=cut

sub fu_and_li_D{
   my ($self,$tree) = @_;
   
   # for this we assume that the tree is made up of
   # allele nodes which contain markers and alleles
   # each marker is a 'mutation' 
   my @nodes = $tree->get_nodes();
   my $muttotal =0;
   my $tipmutcount = 0;
   my $sample_size = 0;
   foreach my $n ( @nodes ) {
   
       if( ! $n->isa('Bio::Tree::AlleleNode') ) {
	   $self->warn("Cannot run fu_and_li_D without AlleleNodes");
	   return 0;
       }
       $muttotal += scalar $n->get_marker_names;
       
       if ($n->is_Leaf() ) {
	   $sample_size++;	   
	   $tipmutcount += $n->get_marker_names();
       }
   }
   
   if( $muttotal <= 0 ) { 
       $self->warn("mutation total was not > 0, cannot calculate a Fu and Li D");
       return 0;
   }

   my $a = 0;
   for(my $k= 1; $k < $sample_size; $k++ ) {
        $a += ( 1 / $k );
    }
   
   my $b = 0;
    for(my $k= 1; $k < $sample_size; $k++ ) {
        $b += ( 1 / $k**2 );
    }
 
    my $c = 2 * ( ( ( $sample_size * $a ) - (2 * ( $sample_size -1 ))) /
                  ( ( $sample_size - 1) * ( $sample_size - 2 ) ) );
 
    my $v = 1 + ( ( $a**2 / ( $b + $a**2 ) ) * ( $c - ( ( $sample_size + 1) /
                                                        ( $sample_size - 1) ) ));
 
    my $u = $a - 1 - $v;
    my $D = ( $muttotal - (  $a * $tipmutcount) ) /
            ( sqrt ( ($u * $muttotal) + ( $v * $muttotal**2) ) );
 
    return $D;
}

# fu_and_li_D*

sub fu_and_li_D_star {
    my ($self,$tree) = @_;
    my ($sample_size, $seg_sites, $singletons);
    my %markers;
    foreach my $n ( $tree->get_nodes() ) {
	if( ! $n->isa('Bio::Tree::AlleleNode') ) {
	    $self->warn("Cannot run fu_and_li_D_star without AlleleNodes");
	    return 0;
	}
	foreach my $m ( $n->get_marker_names ) { $markers{$m}++ }
	$sample_size++ if $n->is_Leaf();
    }
    $seg_sites = scalar keys %markers;

    my $a = 0;
    for(my $k= 1; $k < $sample_size; $k++ ) {
	$a += ( 1 / $k );
    }

    my $a1 = 0;
    for(my $k= 1; $k <= $sample_size; $k++ ) {
	$a1 += ( 1 / $k );
    }

    my $b = 0;
    for(my $k= 1; $k < $sample_size; $k++ ) {
	$b += ( 1 / $k**2 );
    }

    my $c = 2 * ( ( ( $sample_size * $a ) - (2 * ( $sample_size -1 ))) / 
		  ( ( $sample_size - 1) * ( $sample_size - 2 ) ) );

    my $d = $c + ( ($sample_size -2) / ($sample_size - 1)**2 ) +
	( 2 / ($sample_size -1) * 
	  ( (3/2) - ( (2*$a1 - 3) / ($sample_size -2) ) - 
	    ( 1/ $sample_size) ) 
	  );
    my $v_star = ( ( ($sample_size/($sample_size-1) )**2)*$b + (($a**2)*$d) -
		 (2*( ($sample_size*$a*($a+1)) )/(($sample_size-1)**2)) )  /
		   (($a**2) + $b);

    my $u_star = ( ($sample_size/($sample_size-1))*
		   ($a - ($sample_size/
			  ($sample_size-1)))) - $v_star;

    my $D_star = ( (($sample_size/($sample_size-1))*$seg_sites) -
		   ($a*$seg_sites) ) / 
		   ( sqrt( ($u_star*$seg_sites) + ($v_star*($seg_sites**2)) ));
    return $D_star;
}

sub fu_and_li_F {
    my ($self,$tree) = @_;
    
    my ($sample_size, $seg_sites, $pi, $ext_mutations);
    my %markers;
    foreach my $n ( $tree->get_nodes() ) {
	if( ! $n->isa('Bio::Tree::AlleleNode') ) {
	    $self->warn("Cannot run fu_and_li_F without AlleleNodes");
	    return 0;
	}
	my @names = $n->get_marker_names;
	foreach my $m ( @names ) { $markers{$m}++ }
	if ($n->is_Leaf() ) {
	    $sample_size++;
	    $ext_mutations += scalar @names;
	}
    }
    $seg_sites = keys %markers;
    $pi = $self->pi($tree);

    my $a = 0;
    for(my $k= 1; $k < $sample_size; $k++ ) {
	$a += ( 1 / $k );
    }

    my $a1 = 0;
    for(my $k= 1; $k <= $sample_size; $k++ ) {
	$a1 += ( 1 / $k );
    }

    my $b = 0;
    for(my $k= 1; $k < $sample_size; $k++ ) {
	$b += ( 1 / $k**2 );
    }

    my $c = 2 * ( ( ( $sample_size * $a ) - (2 * ( $sample_size -1 ))) / 
		  ( ( $sample_size - 1) * ( $sample_size - 2 ) ) );

    my $v_F = ( $c + ( (2*(($sample_size**2)+$sample_size+3)) / 
		       ( (9*$sample_size)*($sample_size-1) ) ) -
		(2/($sample_size-1)) ) / ( ($a**2)+$b );

    my $u_F = ( 1 + ( ($sample_size+1)/(3*($sample_size-1)) )-
		( 4*( ($sample_size+1)/(($sample_size-1)**2) ))*
		($a1 - ((2*$sample_size)/($sample_size+1))) ) /
		($sample_size - $v_F);

    my $F = ($pi - $ext_mutations) / ( sqrt( ($u_F*$seg_sites) +
					     ($v_F*($seg_sites**2)) ) );

    return $F;
}


sub tajima_D {
    my ($self,$tree) = @_;
    my @compare;
    # we are calculating pi - all pairwise differences between 
    # tips
    my $pi = $self->pi($tree);
    my $a1 = 0; 
    my ($muttotal,$sample_size,$ext_mutations);
    my @nodes = $tree->get_nodes();
    foreach my $n ( @nodes ) {
       if( ! $n->isa('Bio::Tree::AlleleNode') ) {
	   $self->warn("Cannot run tajima_D without AlleleNodes");
	   return 0;
       }       
       $muttotal += $n->get_marker_names();
       if ($n->is_Leaf() ) {
	   $sample_size++;	   
	   $ext_mutations += $n->get_marker_names();
       }
   }
    for(my $k= 1; $k < $sample_size; $k++ ) {
	$a1 += ( 1 / $k );
    }

     my $a2 = 0;
     for(my $k= 1; $k < $sample_size; $k++ ) {
	 $a2 += ( 1 / $k**2 );
     }

    my $b1 = ( $sample_size + 1 ) / ( 3* ( $sample_size - 1) );
    my $b2 = ( 2 * ( $sample_size ** 2 + $sample_size + 3) ) / 
	     ( ( 9 * $sample_size) * ( $sample_size - 1) );
    my $c1 = $b1 - ( 1 / $a1 );
    my $c2 = $b2 - ( ( $sample_size + 2 ) /
		     ( $a1 * $sample_size))+( $a2 / $a1 ** 2);
    my $e1 = $c1 / $a1;
    my $e2 = $c2 / ( $a1**2 + $a2 );
    
    my $D = ( $pi - ( $muttotal / $a1 ) ) / 
	sqrt ( ($e1 * $muttotal) + (( $e2 * $muttotal) * ( $muttotal - 1)));

    return $D;
}

sub pi {
    my ($self, $tree) = @_;
    my ($count,$total,$sample_size) = (0,0,0);
    
    my %markers;
    my @nodes = $tree->get_nodes;
    foreach my $n ( @nodes ) {
	if( ! $n->isa('Bio::Tree::AlleleNode') ) {
	    $self->warn("Cannot run pi without AlleleNodes");
	    return 0;
	}
	# this could be a bit expensive...
	for my $m ( $n->get_marker_names ) { $markers{$m}++ }
	$sample_size++ if $n->is_Leaf();
    }
   
    my @marker_names = keys %markers;

    for( my $i = 0; $i < $sample_size; $i++ ) {
	for( my $j=$i+1; $j < $sample_size; $j++) {	    
	    $total++;
	    for my $name ( @marker_names ) {
		if( $nodes[$i]->has_marker($name) &&
		    $nodes[$j]->has_marker($name) ) {
		    $self->debug(join('',"alleles for $i $name ", 
				      $nodes[$i]->get_alleles($name), "\n"));
		    $self->debug(join('', "\t"x$j,"alleles for $j $name ", 
				      $nodes[$j]->get_alleles($name), "\n"));
		}
	    }
	    #( $nodes[$i]->get_alleles($name)->[0] != 
	    #  $nodes[$j]->get_alleles($name)->[0] ));
	}
    }
    my $pi = $count / $total;
    $self->debug( "total=$total count=$count pi=$pi\n");

    return $pi;
}

sub theta {
    my ($self,$tree) = @_;
    my $sample_size = 0;
    my @nodes = $tree->get_nodes();
    my %markers;
    foreach my $n ( @nodes ) {
	if( ! $n->isa('Bio::Tree::AlleleNode') ) {
	    $self->warn("Cannot run theta without AlleleNodes");
	    return 0;
	}
	for my $m ( $n->get_marker_names ) { $markers{$m}++ }
	if ($n->is_Leaf() ) {
	    # this could be a bit expensive...
	    $sample_size++;
	}
    }
    my $seg_sites = scalar keys %markers;
    my $a1 = 0; 
    for(my $k= 1; $k < $sample_size; $k++ ) {
	$a1 += ( 1 / $k );
    }    
    return $seg_sites / $a1;
}

1;
