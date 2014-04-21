#
# BioPerl module for Bio::PopGen::Statistics
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

Bio::PopGen::Statistics - Population Genetics statistical tests  

=head1 SYNOPSIS

  use Bio::PopGen::Statistics;
  use Bio::AlignIO;
  use Bio::PopGen::IO;
  use Bio::PopGen::Simulation::Coalescent;

  my $sim = Bio::PopGen::Simulation::Coalescent->new( -sample_size => 12);

  my $tree = $sim->next_tree;

  $sim->add_Mutations($tree,20);

  my $stats = Bio::PopGen::Statistics->new();
  my $individuals = [ $tree->get_leaf_nodes];
  my $pi = $stats->pi($individuals);
  my $D  = $stats->tajima_D($individuals);

  # Alternatively to do this on input data from
  # See the tests in t/PopGen.t for more examples
  my $parser = Bio::PopGen::IO->new(-format => 'prettybase',
                                   -file   => 't/data/popstats.prettybase');
  my $pop = $parser->next_population;
  # Note that you can also call the stats as a class method if you like
  # the only reason to instantiate it (as above) is if you want
  # to set the verbosity for debugging
  $pi     = Bio::PopGen::Statistics->pi($pop);
  $theta  = Bio::PopGen::Statistics->theta($pop);

  # Pi and Theta also take additional arguments,
  # see the documentation for more information

  use Bio::PopGen::Utilities;
  use Bio::AlignIO;

  my $in = Bio::AlignIO->new(-file   => 't/data/t7.aln',
                            -format => 'clustalw');
  my $aln = $in->next_aln;
  # get a population, each sequence is an individual and 
  # for the default case, every site which is not monomorphic
  # is a 'marker'.  Each individual will have a 'genotype' for the
  # site which will be the specific base in the alignment at that
  # site

  my $pop = Bio::PopGen::Utilities->aln_to_population(-alignment => $aln);


=head1 DESCRIPTION

This object is intended to provide implementations some standard
population genetics statistics about alleles in populations.

This module was previously named Bio::Tree::Statistics.

This object is a place to accumulate routines for calculating various
statistics from the coalescent simulation, marker/allele, or from
aligned sequence data given that you can calculate alleles, number of
segregating sites.

Currently implemented:
 Fu and Li's D    (fu_and_li_D)
 Fu and Li's D*   (fu_and_li_D_star)
 Fu and Li's F    (fu_and_li_F)
 Fu and Li's F*   (fu_and_li_F_star)
 Tajima's D       (tajima_D)
 Watterson's theta (theta)
 pi               (pi) - number of pairwise differences
 composite_LD     (composite_LD)
 McDonald-Kreitman (mcdonald_kreitman or MK)

Count based methods also exist in case you have already calculated the
key statistics (seg sites, num individuals, etc) and just want to
compute the statistic.

In all cases where a the method expects an arrayref of
L<Bio::PopGen::IndividualI> objects and L<Bio::PopGen::PopulationI>
object will also work.

=head2 REFERENCES

Fu Y.X and Li W.H. (1993) "Statistical Tests of Neutrality of
Mutations." Genetics 133:693-709.

Fu Y.X. (1996) "New Statistical Tests of Neutrality for DNA samples
from a Population." Genetics 143:557-570.

McDonald J, Kreitman M.

Tajima F. (1989) "Statistical method for testing the neutral mutation
hypothesis by DNA polymorphism." Genetics 123:585-595.


=head2 CITING THIS WORK

Please see this reference for use of this implementation.

Stajich JE and Hahn MW "Disentangling the Effects of Demography and Selection in Human History." (2005) Mol Biol Evol 22(1):63-73. 

If you use these Bio::PopGen modules please cite the Bioperl
publication (see FAQ) and the above reference.


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

=head1 AUTHOR - Jason Stajich, Matthew Hahn

Email jason-at-bioperl-dot-org
Email matthew-dot-hahn-at-duke-dot-edu

McDonald-Kreitman implementation based on work by Alisha Holloway at
UC Davis.


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::PopGen::Statistics;
use strict;
use constant { 
    in_label => 'ingroup',
    out_label => 'outgroup',
    non_syn   => 'non_synonymous',
    syn       => 'synonymous',
    default_codon_table => 1, # Standard Codon table
};

use Bio::MolEvol::CodonModel;
use List::Util qw(sum);

use base qw(Bio::Root::Root);
our $codon_table => default_codon_table;
our $has_twotailed => 0;
BEGIN {
    eval { require Text::NSP::Measures::2D::Fisher2::twotailed };
    if( $@ ) { $has_twotailed = 0; }
    else { $has_twotailed = 1; }
}






=head2 new

 Title   : new
 Usage   : my $obj = Bio::PopGen::Statistics->new();
 Function: Builds a new Bio::PopGen::Statistics object 
 Returns : an instance of Bio::PopGen::Statistics
 Args    : none


=cut


=head2 fu_and_li_D

 Title   : fu_and_li_D
 Usage   : my $D = $statistics->fu_and_li_D(\@ingroup,\@outgroup);
	    OR
	   my $D = $statistics->fu_and_li_D(\@ingroup,$extmutations);
 Function: Fu and Li D statistic for a list of individuals
           given an outgroup and the number of external mutations
           (either provided or calculated from list of outgroup individuals)
 Returns : decimal
 Args    : $individuals - array reference which contains ingroup individuals 
           (L<Bio::PopGen::Individual> or derived classes)
           $extmutations - number of external mutations OR
           arrayref of outgroup individuals

=cut

sub fu_and_li_D { 
    my ($self,$ingroup,$outgroup) = @_;

    my ($seg_sites,$n,$ancestral,$derived) = (0,0,0,0);
    if( ref($ingroup) =~ /ARRAY/i ) {
	$n = scalar @$ingroup;
	# pi - all pairwise differences 
	$seg_sites   = $self->segregating_sites_count($ingroup);
    } elsif( ref($ingroup) && 
	     $ingroup->isa('Bio::PopGen::PopulationI')) {
	$n = $ingroup->get_number_individuals;
	$seg_sites   = $self->segregating_sites_count($ingroup);
    } else { 
	$self->throw("expected an array reference of a list of Bio::PopGen::IndividualI OR a Bio::PopGen::PopulationI object to fu_and_li_D");
	return 0;
    }
    
    if( $seg_sites <= 0 ) { 
	$self->warn("mutation total was not > 0, cannot calculate a Fu and Li D");
	return 0;
    }

    if( ! defined $outgroup ) {
	$self->warn("Need to provide either an array ref to the outgroup individuals or the number of external mutations");
	return 0;
    } elsif( ref($outgroup) ) {
	($ancestral,$derived) = $self->derived_mutations($ingroup,$outgroup);
	$ancestral = 0 unless defined $ancestral;
    } else { 
	$ancestral = $outgroup;
    }
   
    return $self->fu_and_li_D_counts($n,$seg_sites,
				     $ancestral,$derived);
}

=head2 fu_and_li_D_counts

 Title   : fu_li_D_counts
 Usage   : my $D = $statistics->fu_and_li_D_counts($samps,$sites,
                                                   $external);
 Function: Fu and Li D statistic for the raw counts of the number
           of samples, sites, external and internal mutations
 Returns : decimal number
 Args    : number of samples (N)
           number of segregating sites (n)
           number of external mutations (n_e)

=cut


sub fu_and_li_D_counts {
    my ($self,$n,$seg_sites, $external_mut) = @_;
    my $a_n = 0;
    if( $n <= 3 ) {
	$self->warn("n is $n, too small, must be > 3\n");
	return;
    }
    for(my $k= 1; $k < $n; $k++ ) {
	$a_n += ( 1 / $k );
    }
    my $b = 0;
    for(my $k= 1; $k < $n; $k++ ) {
        $b += ( 1 / $k**2 );
    }

    my $c = 2 * ( ( ( $n * $a_n ) - (2 * ( $n -1 ))) /
                  ( ( $n - 1) * ( $n - 2 ) ) );

    my $v = 1 + ( ( $a_n**2 / ( $b + $a_n**2 ) ) * 
		  ( $c - ( ( $n + 1) /
			   ( $n - 1) ) ));
    
    my $u = $a_n - 1 - $v;

    ($seg_sites - $a_n * $external_mut) / 
	sqrt( ($u * $seg_sites) + ($v * $seg_sites*$seg_sites));
    
}


=head2 fu_and_li_D_star

 Title   : fu_and_li_D_star
 Usage   : my $D = $statistics->fu_an_li_D_star(\@individuals);
 Function: Fu and Li's D* statistic for a set of samples
            Without an outgroup
 Returns : decimal number
 Args    : array ref of L<Bio::PopGen::IndividualI> objects
           OR
           L<Bio::PopGen::PopulationI> object

=cut

#'
# fu_and_li_D*

sub fu_and_li_D_star {
    my ($self,$individuals) = @_;

    my ($seg_sites,$n,$singletons);
    if( ref($individuals) =~ /ARRAY/i ) {
	$n = scalar @$individuals;
	$seg_sites   = $self->segregating_sites_count($individuals);
	$singletons  = $self->singleton_count($individuals);
    } elsif( ref($individuals) && 
	     $individuals->isa('Bio::PopGen::PopulationI')) {
	my $pop = $individuals;
	$n = $pop->get_number_individuals;
	$seg_sites   = $self->segregating_sites_count($pop);
	$singletons  = $self->singleton_count($pop);
    } else { 
	$self->throw("expected an array reference of a list of Bio::PopGen::IndividualI OR a Bio::PopGen::PopulationI object to fu_and_li_D_star");
	return 0;
    }

    return $self->fu_and_li_D_star_counts($n,$seg_sites, $singletons);
}

=head2 fu_and_li_D_star_counts

 Title   : fu_li_D_star_counts
 Usage   : my $D = $statistics->fu_and_li_D_star_counts($samps,$sites,
                                                        $singletons);

 Function: Fu and Li D statistic for the raw counts of the number
           of samples, sites, external and internal mutations
 Returns : decimal number
 Args    : number of samples (N)
           number of segregating sites (n)
           singletons (n_s)

=cut


sub fu_and_li_D_star_counts {
    my ($self,$n,$seg_sites, $singletons) = @_;
    my $a_n;
    for(my $k = 1; $k < $n; $k++ ) {
	$a_n += ( 1 / $k );
    }

    my $a1 = $a_n + 1 / $n;

    my $b = 0;
    for(my $k= 1; $k < $n; $k++ ) {
        $b += ( 1 / $k**2 );
    }

    my $c = 2 * ( ( ( $n * $a_n ) - (2 * ( $n -1 ))) /
                  ( ( $n - 1) * ( $n - 2 ) ) );

    my $d = $c + ($n -2) / ($n - 1)**2 +
	2 / ($n -1) * 
	( 1.5 - ( (2*$a1 - 3) / ($n -2) ) - 
	  1 / $n ); 
    
    my $v_star = ( ( ($n/($n-1) )**2)*$b + (($a_n**2)*$d) -
		 (2*( ($n*$a_n*($a_n+1)) )/(($n-1)**2)) )  /
		   (($a_n**2) + $b);

    my $u_star = ( ($n/($n-1))*
		   ($a_n - ($n/
			  ($n-1)))) - $v_star;


    return (($n / ($n - 1)) * $seg_sites - 
	    $a_n * $singletons) / 
	    sqrt( ($u_star * $seg_sites) + ($v_star * $seg_sites*$seg_sites));
}


=head2 fu_and_li_F

 Title   : fu_and_li_F
 Usage   : my $F = Bio::PopGen::Statistics->fu_and_li_F(\@ingroup,$ext_muts);
 Function: Calculate Fu and Li's F on an ingroup with either the set of 
           outgroup individuals, or the number of external mutations
 Returns : decimal number
 Args    : array ref of L<Bio::PopGen::IndividualI> objects for the ingroup
           OR a L<Bio::PopGen::PopulationI> object
           number of external mutations OR list of individuals for the outgroup

=cut

#'

sub fu_and_li_F {
    my ($self,$ingroup,$outgroup) = @_;
    my ($seg_sites,$pi,$n,$external,$internal);
    if( ref($ingroup) =~ /ARRAY/i ) {
	$n = scalar @$ingroup;
	# pi - all pairwise differences 
	$pi          = $self->pi($ingroup);  
	$seg_sites   = $self->segregating_sites_count($ingroup);
    } elsif( ref($ingroup) && 
	     $ingroup->isa('Bio::PopGen::PopulationI')) {
	$n = $ingroup->get_number_individuals;
	$pi          = $self->pi($ingroup);
	$seg_sites   = $self->segregating_sites_count($ingroup);
    } else { 
	$self->throw("expected an array reference of a list of Bio::PopGen::IndividualI OR a Bio::PopGen::PopulationI object to Fu and Li's F");
	return 0;
    }
    
    if( ! defined $outgroup ) {
	$self->warn("Need to provide either an array ref to the outgroup individuals or the number of external mutations");
	return 0;
    } elsif( ref($outgroup) ) {
	($external,$internal) = $self->derived_mutations($ingroup,$outgroup);
    } else { 
	$external = $outgroup;
    }
    $self->fu_and_li_F_counts($n,$pi,$seg_sites,$external);
}

=head2 fu_and_li_F_counts

 Title   : fu_li_F_counts
 Usage   : my $F = $statistics->fu_and_li_F_counts($samps,$pi,
                                                   $sites,
                                                   $external);
 Function: Fu and Li F statistic for the raw counts of the number
           of samples, sites, external and internal mutations
 Returns : decimal number
 Args    : number of samples (N)
           average pairwise differences (pi)
           number of segregating sites (n)
           external mutations (n_e)

=cut


sub fu_and_li_F_counts {
    my ($self,$n,$pi,$seg_sites, $external) = @_;
    my $a_n = 0;
    for(my $k= 1; $k < $n; $k++ ) {
	$a_n += ( 1 / $k );
    }

    my $a1 = $a_n + (1 / $n );

    my $b = 0;
    for(my $k= 1; $k < $n; $k++ ) {
	$b += ( 1 / $k**2 );
    }

    my $c = 2 * ( ( ( $n * $a_n ) - (2 * ( $n -1 ))) / 
		  ( ( $n - 1) * ( $n - 2 ) ) );

    my $v_F = ( $c + ( (2*(($n**2)+$n+3)) / 
		       ( (9*$n)*($n-1) ) ) -
		(2/($n-1)) ) / ( ($a_n**2)+$b );

    my $u_F = ( 1 + ( ($n+1)/(3*($n-1)) )-
		( 4*( ($n+1)/(($n-1)**2) ))*
		($a1 - ((2*$n)/($n+1))) ) /
		$a_n - $v_F;

    # warn("$v_F vf $u_F uf n = $n\n");
    my $F = ($pi - $external) / ( sqrt( ($u_F*$seg_sites) +
					($v_F*($seg_sites**2)) ) );

    return $F;
}

=head2 fu_and_li_F_star

 Title   : fu_and_li_F_star
 Usage   : my $F = Bio::PopGen::Statistics->fu_and_li_F_star(\@ingroup);
 Function: Calculate Fu and Li's F* on an ingroup without an outgroup
           It uses count of singleton alleles instead 
 Returns : decimal number
 Args    : array ref of L<Bio::PopGen::IndividualI> objects for the ingroup
           OR
           L<Bio::PopGen::PopulationI> object

=cut

#' keep my emacs happy

sub fu_and_li_F_star {
    my ($self,$individuals) = @_;

    my ($seg_sites,$pi,$n,$singletons);
    if( ref($individuals) =~ /ARRAY/i ) {
	$n = scalar @$individuals;
	# pi - all pairwise differences 
	$pi          = $self->pi($individuals);  
	$seg_sites   = $self->segregating_sites_count($individuals);
	$singletons  = $self->singleton_count($individuals);
    } elsif( ref($individuals) && 
	     $individuals->isa('Bio::PopGen::PopulationI')) {
	my $pop = $individuals;
	$n = $pop->get_number_individuals;
	$pi          = $self->pi($pop);
	$seg_sites   = $self->segregating_sites_count($pop);
	$singletons  = $self->singleton_count($pop);
    } else { 
	$self->throw("expected an array reference of a list of Bio::PopGen::IndividualI OR a Bio::PopGen::PopulationI object to fu_and_li_F_star");
	return 0;
    }
    return $self->fu_and_li_F_star_counts($n,
					  $pi,
					  $seg_sites,
					  $singletons);
} 

=head2 fu_and_li_F_star_counts

 Title   : fu_li_F_star_counts
 Usage   : my $F = $statistics->fu_and_li_F_star_counts($samps,
                                                   $pi,$sites,
                                                   $singletons);
 Function: Fu and Li F statistic for the raw counts of the number
           of samples, sites, external and internal mutations
 Returns : decimal number
 Args    : number of samples (N)
           average pairwise differences (pi)
           number of segregating sites (n)
           singleton  mutations (n_s)

=cut


sub fu_and_li_F_star_counts {
    my ($self,$n,$pi,$seg_sites, $singletons) = @_;
    if( $n <= 1 ) {
	$self->warn("N must be > 1\n");
	return;
    }
    if( $n == 2) { 
	return 0;
    } 

    my $a_n = 0;
    

    my $b = 0;
    for(my $k= 1; $k < $n; $k++ ) {
	$b += (1 / ($k**2));
	$a_n += ( 1 / $k );     # Eq (2)
    }
    my $a1 = $a_n + (1 / $n );

    # warn("a_n is $a_n a1 is $a1 n is $n b is $b\n");

    # From Simonsen et al (1995) instead of Fu and Li 1993
    my $v_F_star = ( (( 2 * $n ** 3 + 110 * $n**2 - (255 * $n) + 153)/
		      (9 * ($n ** 2) * ( $n - 1))) +
		     ((2 * ($n - 1) * $a_n ) / $n ** 2) -
		     (8 * $b / $n) ) / 
		     ( ($a_n ** 2) + $b );
    
    my $u_F_star = ((( (4* ($n**2)) + (19 * $n) + 3 - (12 * ($n + 1)* $a1)) /
		    (3 * $n * ( $n - 1))) / $a_n) - $v_F_star;

    # warn("vf* = $v_F_star uf* = $u_F_star n = $n\n");
    my $F_star = ( $pi - ($singletons*( ( $n-1) / $n)) ) /
	sqrt ( $u_F_star*$seg_sites + $v_F_star*$seg_sites**2);
    return $F_star;
}

=head2 tajima_D

 Title   : tajima_D
 Usage   : my $D = Bio::PopGen::Statistics->tajima_D(\@samples);
 Function: Calculate Tajima's D on a set of samples 
 Returns : decimal number
 Args    : array ref of L<Bio::PopGen::IndividualI> objects
           OR 
           L<Bio::PopGen::PopulationI> object


=cut

#'

sub tajima_D {
    my ($self,$individuals) = @_;
    my ($seg_sites,$pi,$n);

    if( ref($individuals) =~ /ARRAY/i ) {
	$n = scalar @$individuals;
	# pi - all pairwise differences 
	$pi          = $self->pi($individuals);  
	$seg_sites = $self->segregating_sites_count($individuals);

    } elsif( ref($individuals) && 
	     $individuals->isa('Bio::PopGen::PopulationI')) {
	my $pop = $individuals;
	$n = $pop->get_number_individuals;
	$pi          = $self->pi($pop);
	$seg_sites = $self->segregating_sites_count($pop);
    } else { 
	$self->throw("expected an array reference of a list of Bio::PopGen::IndividualI OR a Bio::PopGen::PopulationI object to tajima_D");
	return 0;
    }
    $self->tajima_D_counts($n,$seg_sites,$pi);
}

=head2 tajima_D_counts

 Title   : tajima_D_counts
 Usage   : my $D = $statistics->tajima_D_counts($samps,$sites,$pi);
 Function: Tajima's D statistic for the raw counts of the number
           of samples, sites, and avg pairwise distances (pi)
 Returns : decimal number
 Args    : number of samples (N)
           number of segregating sites (n)
           average pairwise differences (pi)

=cut

#'

sub tajima_D_counts {
    my ($self,$n,$seg_sites,$pi) = @_;
    my $a1 = 0; 
    for(my $k= 1; $k < $n; $k++ ) {
	$a1 += ( 1 / $k );
    }

     my $a2 = 0;
     for(my $k= 1; $k < $n; $k++ ) {
	 $a2 += ( 1 / $k**2 );
     }
    
    my $b1 = ( $n + 1 ) / ( 3* ( $n - 1) );
    my $b2 = ( 2 * ( $n ** 2 + $n + 3) ) / 
	     ( ( 9 * $n) * ( $n - 1) );
    my $c1 = $b1 - ( 1 / $a1 );
    my $c2 = $b2 - ( ( $n + 2 ) /
		     ( $a1 * $n))+( $a2 / $a1 ** 2);
    my $e1 = $c1 / $a1;
    my $e2 = $c2 / ( $a1**2 + $a2 );
    
    my $denom = sqrt ( ($e1 * $seg_sites) + (( $e2 * $seg_sites) * ( $seg_sites - 1)));
    return if $denom == 0;
    my $D = ( $pi - ( $seg_sites / $a1 ) ) / $denom;
    return $D;
}


=head2 pi

 Title   : pi
 Usage   : my $pi = Bio::PopGen::Statistics->pi(\@inds)
 Function: Calculate pi (average number of pairwise differences) given
           a list of individuals which have the same number of markers
           (also called sites) as available from the get_Genotypes()
           call in L<Bio::PopGen::IndividualI>
 Returns : decimal number
 Args    : Arg1= array ref of L<Bio::PopGen::IndividualI> objects
             which have markers/mutations.  We expect all individuals to
             have a marker - we will deal with missing data as a special case.
           OR
           Arg1= L<Bio::PopGen::PopulationI> object.  In the event that
                 only allele frequency data is available, storing it in
                 Population object will make this available.
           num sites [optional], an optional second argument (integer)
             which is the number of sites, then pi returned is pi/site.

=cut

sub pi {
    my ($self,$individuals,$numsites) = @_;
    my (%data,%marker_total,@marker_names,$n);

    if( ref($individuals) =~ /ARRAY/i ) {
	# one possible argument is an arrayref of Bio::PopGen::IndividualI objs
	@marker_names = $individuals->[0]->get_marker_names;
	$n = scalar @$individuals;

	# Here we are calculating the allele frequencies
	foreach my $ind ( @$individuals ) {
	    if( ! $ind->isa('Bio::PopGen::IndividualI') ) {
		$self->warn("Expected an arrayref of Bio::PopGen::IndividualI objects, this is a ".ref($ind)."\n");
		return 0;
	    }
	    foreach my $m ( @marker_names ) {
		foreach my $allele (map { $_->get_Alleles} 
				    $ind->get_Genotypes($m) ) {
		    $data{$m}->{$allele}++;
		    $marker_total{$m}++;
		}
	    }
	}
#	while( my ($marker,$count) =  each %marker_total ) {
#	    foreach my $c ( values %{$data{$marker}} ) {
#		$c /= $count;
#	    }
#	}
	# %data will contain allele frequencies for each marker, allele
    } elsif( ref($individuals) &&
	     $individuals->isa('Bio::PopGen::PopulationI') ) {
	my $pop = $individuals;
	$n = $pop->get_number_individuals;
	foreach my $marker( $pop->get_Markers ) {
	    push @marker_names, $marker->name;
	    #$data{$marker->name} = {$marker->get_Allele_Frequencies};
	    my @genotypes = $pop->get_Genotypes(-marker => $marker->name);
	    for my $al ( map { $_->get_Alleles} @genotypes ) {
	      $data{$marker->name}->{$al}++;
	      $marker_total{$marker->name}++;
	   }
	}
    } else {
	$self->throw("expected an array reference of a list of Bio::PopGen::IndividualI to pi");
    }
    # based on Kevin Thornton's code:
    # http://molpopgen.org/software/libsequence/doc/html/PolySNP_8cc-source.html#l00152
    # For now we assume that all individuals have the same markers
    my ($diffcount,$totalcompare) = (0,0);
    my $pi = 0;
    while ( my ($marker,$markerdat) = each %data ) {
      my $sampsize = $marker_total{$marker};
      my $ssh = 0;
      my @alleles = keys %$markerdat;
      if ( $sampsize > 1 ) {
	my $denom = $sampsize * ($sampsize - 1.0);
	foreach my $al ( @alleles ) {
	  $ssh += ($markerdat->{$al} * ($markerdat->{$al} - 1)) / $denom;
	}
	$pi += 1.0 - $ssh;
      }
    }
    $self->debug( "pi=$pi\n");
    if( $numsites ) {
	return $pi / $numsites;
    } else {
	return $pi;
    }
}


=head2 theta

 Title   : theta
 Usage   : my $theta = Bio::PopGen::Statistics->theta($sampsize,$segsites);
 Function: Calculates Watterson's theta from the sample size 
           and the number of segregating sites.
           Providing the third parameter, total number of sites will
           return theta per site.
           This is also known as K-hat = K / a_n   
 Returns : decimal number 
 Args    : sample size (integer),
           num segregating sites (integer)
           total sites (integer) [optional] (to calculate theta per site)
           OR
           provide an arrayref of the L<Bio::PopGen::IndividualI> objects
           total sites (integer) [optional] (to calculate theta per site)
           OR
           provide an L<Bio::PopGen::PopulationI> object
           total sites (integer)[optional]

=cut

#'

sub theta {
    my $self = shift;
    my ( $n, $seg_sites,$totalsites) = @_;
    if( ref($n) =~ /ARRAY/i ) {
	my $samps = $n;
	$totalsites = $seg_sites; # only 2 arguments if one is an array
	my %data;
	my @marker_names = $samps->[0]->get_marker_names;
	# we need to calculate number of polymorphic sites
	$seg_sites = $self->segregating_sites_count($samps);
	$n = scalar @$samps;

    } elsif(ref($n) &&
	    $n->isa('Bio::PopGen::PopulationI') ) {
	# This will handle the case when we pass in a PopulationI object
	my $pop = $n;
	$totalsites = $seg_sites; # shift the arguments over by one
	$n = $pop->haploid_population->get_number_individuals;
	$seg_sites = $self->segregating_sites_count($pop);
    }
    my $a1 = 0; 
    for(my $k= 1; $k < $n; $k++ ) {
	$a1 += ( 1 / $k );
    }    
    if( $totalsites ) { # 0 and undef are the same can't divide by them
	$seg_sites /= $totalsites;
    }
    if( $a1 == 0 ) { 
	return 0;
    } 
    return $seg_sites / $a1;
}

=head2 singleton_count

 Title   : singleton_count
 Usage   : my ($singletons) = Bio::PopGen::Statistics->singleton_count(\@inds)
 Function: Calculate the number of mutations/alleles which only occur once in
           a list of individuals for all sites/markers
 Returns : (integer) number of alleles which only occur once (integer)
 Args    : arrayref of L<Bio::PopGen::IndividualI> objects
           OR
           L<Bio::PopGen::PopulationI> object

=cut

sub singleton_count {
    my ($self,$individuals) = @_;

    my @inds;
    if( ref($individuals) =~ /ARRAY/ ) {
	@inds = @$individuals;
    } elsif( ref($individuals) && 
	     $individuals->isa('Bio::PopGen::PopulationI') ) {
	my $pop = $individuals;
	@inds = $pop->get_Individuals();
	unless( @inds ) { 
	    $self->warn("Need to provide a population which has individuals loaded, not just a population with allele frequencies");
	    return 0;
	}
    } else {
	$self->warn("Expected either a PopulationI object or an arrayref of IndividualI objects");
	return 0;
    }
    # find number of sites where a particular allele is only seen once

    my ($singleton_allele_ct,%sites) = (0);
    # first collect all the alleles into a hash structure
    
    foreach my $n ( @inds ) {
	if( ! $n->isa('Bio::PopGen::IndividualI') ) {
	    $self->warn("Expected an arrayref of Bio::PopGen::IndividualI objects, this is a ".ref($n)."\n");
	    return 0;
	}
	foreach my $g ( $n->get_Genotypes ) {
	    my ($nm,@alleles) = ($g->marker_name, $g->get_Alleles);
	    foreach my $allele (@alleles ) {
		$sites{$nm}->{$allele}++;
	    }
	}
    }
    foreach my $site ( values %sites ) { # don't really care what the name is
	foreach my $allelect ( values %$site ) { # 
            # find the sites which have an allele with only 1 copy
 	    $singleton_allele_ct++ if( $allelect == 1 );
	}
    }
    return $singleton_allele_ct;
}

# Yes I know that singleton_count and segregating_sites_count are
# basically processing the same data so calling them both is
# redundant, something I want to fix later but want to make things
# correct and simple first

=head2 segregating_sites_count

 Title   : segregating_sites_count
 Usage   : my $segsites = Bio::PopGen::Statistics->segregating_sites_count
 Function: Gets the number of segregating sites (number of polymorphic sites)
 Returns : (integer) number of segregating sites
 Args    : arrayref of L<Bio::PopGen::IndividualI> objects 
           OR
           L<Bio::PopGen::PopulationI> object

=cut

# perhaps we'll change this in the future 
# to return the actual segregating sites
# so one can use this to pull in the names of those sites.
# Would be trivial if it is useful.

sub segregating_sites_count {
   my ($self,$individuals) = @_;
   my $type = ref($individuals);
   my $seg_sites = 0;
   if( $type =~ /ARRAY/i ) {
       my %sites;
       foreach my $n ( @$individuals ) {
	   if( ! $n->isa('Bio::PopGen::IndividualI') ) {
	       $self->warn("Expected an arrayref of Bio::PopGen::IndividualI objects, this is a ".ref($n)."\n");
	       return 0;
	   }
	   foreach my $g ( $n->get_Genotypes ) {
	       my ($nm,@alleles) = ($g->marker_name, $g->get_Alleles);
	       foreach my $allele (@alleles ) {
		   $sites{$nm}->{$allele}++;
	       }
	   }
       }
       foreach my $site ( values %sites ) { # use values b/c we don't 
	                                    # really care what the name is
	   # find the sites which >1 allele
	   $seg_sites++ if( keys %$site > 1 );
       }
   } elsif( $type && $individuals->isa('Bio::PopGen::PopulationI') ) {
       foreach my $marker ( $individuals->haploid_population->get_Markers ) {  
	   my @alleles = $marker->get_Alleles;	    
	   $seg_sites++ if ( scalar @alleles > 1 );
       }
   } else { 
       $self->warn("segregating_sites_count expects either a PopulationI object or a list of IndividualI objects");
       return 0;
   } 
   return $seg_sites;
}


=head2 heterozygosity

 Title   : heterozygosity
 Usage   : my $het = Bio::PopGen::Statistics->heterozygosity($sampsize,$freq1);
 Function: Calculate the heterozgosity for a sample set for a set of alleles
 Returns : decimal number
 Args    : sample size (integer)
           frequency of one allele (fraction - must be less than 1)
           [optional] frequency of another allele - this is only needed
                      in a non-binary allele system

Note     : p^2 + 2pq + q^2

=cut


sub heterozygosity {
    my ($self,$samp_size, $freq1,$freq2) = @_;
    if( ! $freq2 ) { $freq2 = 1 - $freq1 }
    if( $freq1 > 1 || $freq2 > 1 ) { 
	$self->warn("heterozygosity expects frequencies to be less than 1");
    }
    my $sum = ($freq1**2) + (($freq2)**2);
    my $h = ( $samp_size*(1- $sum) ) / ($samp_size - 1) ;
    return $h;
}


=head2 derived_mutations

 Title   : derived_mutations
 Usage   : my $ext = Bio::PopGen::Statistics->derived_mutations($ingroup,$outgroup);
 Function: Calculate the number of alleles or (mutations) which are ancestral
           and the number which are derived (occurred only on the tips)
 Returns : array of 2 items - number of external and internal derived 
           mutation
 Args    : ingroup - L<Bio::PopGen::IndividualI>s arrayref OR 
                     L<Bio::PopGen::PopulationI>
           outgroup- L<Bio::PopGen::IndividualI>s arrayref OR 
                     L<Bio::PopGen::PopulationI> OR
                     a single L<Bio::PopGen::IndividualI>

=cut

sub derived_mutations {
   my ($self,$ingroup,$outgroup) = @_;
   my (%indata,%outdata,@marker_names);

   # basically we have to do some type checking
   # if that perl were typed...
   my ($itype,$otype) = (ref($ingroup),ref($outgroup));

   return $outgroup unless( $otype ); # we expect arrayrefs or objects, nums
                                      # are already the value we 
                                      # are searching for
   # pick apart the ingroup
   # get the data
   if( ref($ingroup) =~ /ARRAY/i ) {
       if( ! ref($ingroup->[0]) ||
	   ! $ingroup->[0]->isa('Bio::PopGen::IndividualI') ) {
	   $self->warn("Expected an arrayref of Bio::PopGen::IndividualI objects or a Population for ingroup in external_mutations");
	   return 0;
       }
       # we assume that all individuals have the same markers 
       # i.e. that they are aligned
       @marker_names = $ingroup->[0]->get_marker_names;
       for my $ind ( @$ingroup ) {
	   for my $m ( @marker_names ) {
	       for my $allele ( map { $_->get_Alleles }
				    $ind->get_Genotypes($m) ) {
		   $indata{$m}->{$allele}++;
	       }
	   }
       }	   
   } elsif( ref($ingroup) && $ingroup->isa('Bio::PopGen::PopulationI') ) {
       @marker_names = $ingroup->get_marker_names;
       for my $ind ( $ingroup->haploid_population->get_Individuals() ) {
	   for my $m ( @marker_names ) {
	       for my $allele ( map { $_->get_Alleles} 
				    $ind->get_Genotypes($m) ) {
		   $indata{$m}->{$allele}++;
	       }
	   }
       }
   } else { 
       $self->warn("Need an arrayref of Bio::PopGen::IndividualI objs or a Bio::PopGen::Population for ingroup in external_mutations");
       return 0;
   }
    
   if( $otype =~ /ARRAY/i ) {
       if( ! ref($outgroup->[0]) ||
	   ! $outgroup->[0]->isa('Bio::PopGen::IndividualI') ) {
	   $self->warn("Expected an arrayref of Bio::PopGen::IndividualI objects or a Population for outgroup in external_mutations");
	   return 0;
       }
       for my $ind ( @$outgroup ) {
	   for my $m ( @marker_names ) {
	       for my $allele ( map { $_->get_Alleles }
				$ind->get_Genotypes($m) ) {
		   $outdata{$m}->{$allele}++;
	       }
	   }
       }
   
   } elsif( $otype->isa('Bio::PopGen::PopulationI') ) {
       for my $ind ( $outgroup->haploid_population->get_Individuals() ) {
	   for my $m ( @marker_names ) {
	       for my $allele ( map { $_->get_Alleles} 
				    $ind->get_Genotypes($m) ) {
		   $outdata{$m}->{$allele}++;
	       }
	   }
       }
   } else {
       $self->warn("Need an arrayref of Bio::PopGen::IndividualI objs or a Bio::PopGen::Population for outgroup in external_mutations");
       return 0;
   }
   
   # derived mutations are defined as 
   # 
   # ingroup  (G A T)
   # outgroup (A)
   # derived mutations are G and T, A is the external mutation
   
   # ingroup  (A T)
   # outgroup (C)
   # derived mutations A,T no external/ancestral mutations
   
   # ingroup  (G A T)
   # outgroup (A T)
   # cannot determine
  
   my ($internal,$external);
   foreach my $marker ( @marker_names ) {
       my @outalleles = keys %{$outdata{$marker}};
       my @in_alleles = keys %{$indata{$marker}};
       next if( @outalleles > 1 || @in_alleles == 1);
       for my $allele ( @in_alleles ) {
	   if( ! exists $outdata{$marker}->{$allele} ) { 
	       if( $indata{$marker}->{$allele} == 1 ) { 
		   $external++;
	       } else { 
		   $internal++;
	       }
	   }
       }
   }
   return ($external, $internal);
}


=head2 composite_LD

 Title   : composite_LD
 Usage   : %matrix = Bio::PopGen::Statistics->composite_LD($population);
 Function: Calculate the Linkage Disequilibrium 
           This is for calculating LD for unphased data. 
           Other methods will be appropriate for phased haplotype data.

 Returns : Hash of Hashes - first key is site 1,second key is site 2
           and value is LD for those two sites.
           my $LDarrayref = $matrix{$site1}->{$site2};
           my ($ldval, $chisquared) = @$LDarrayref;
 Args    : L<Bio::PopGen::PopulationI> or arrayref of 
           L<Bio::PopGen::IndividualI>s 
 Reference: Weir B.S. (1996) "Genetic Data Analysis II", 
                      Sinauer, Sunderlanm MA.

=cut

sub composite_LD {
    my ($self,$pop) = @_;
    if( ref($pop) =~ /ARRAY/i ) {
	if( ref($pop->[0]) && $pop->[0]->isa('Bio::PopGen::IndividualI') ) {
	    $pop = Bio::PopGen::Population->new(-individuals => @$pop);
	} else { 
	    $self->warn("composite_LD expects a Bio::PopGen::PopulationI or an arrayref of Bio::PopGen::IndividualI objects");
	    return ();
	}
    } elsif( ! ref($pop) || ! $pop->isa('Bio::PopGen::PopulationI') ) {
	$self->warn("composite_LD expects a Bio::PopGen::PopulationI or an arrayref of Bio::PopGen::IndividualI objects");
	return ();
    }

    my @marker_names = $pop->get_marker_names;
    my @inds = $pop->get_Individuals;
    my $num_inds = scalar @inds;
    my (%lookup);
    # calculate allele frequencies for each marker from the population
    # use the built-in get_Marker to get the allele freqs
    # we still need to calculate the genotype frequencies
    foreach my $marker_name ( @marker_names ) {	
	my(%allelef);

	foreach my $ind ( @inds ) {
	    my ($genotype) = $ind->get_Genotypes(-marker => $marker_name);
	    if( ! defined $genotype ) { 
		$self->warn("no genotype for marker $marker_name for individual ". $ind->unique_id. "\n");
		next;
	    }
	    my @alleles  = sort $genotype->get_Alleles;
	    next if( scalar @alleles != 2);
	    my $genostr  = join(',', @alleles);
            $allelef{$alleles[0]}++;
            $allelef{$alleles[1]}++;
	}

	# we should check for cases where there > 2 alleles or
	# only 1 allele and throw out those markers.
	my @alleles      = sort keys %allelef;
	my $allele_count = scalar @alleles;
	# test if site is polymorphic
	if( $allele_count != 2) { 
	    # only really warn if we're seeing multi-allele
	    $self->warn("Skipping $marker_name because it has $allele_count alleles (".join(',',@alleles)."), \ncomposite_LD will currently only work for biallelic markers") if $allele_count > 2;
	    next;		# skip this marker
	}

	# Need to do something here to detect alleles which aren't 
	# a single character
	if( length($alleles[0]) != 1 ||
	    length($alleles[1]) != 1 ) {
	    $self->warn("An individual has an allele which is not a single base, this is currently not supported in composite_LD - consider recoding the allele as a single character");
	    next;
	}

	# fix the call for allele 1 (A or B) and 
	# allele 2 (a or b) in terms of how we'll do the 
	# N square from Weir p.126
	$self->debug( "$alleles[0] is 1, $alleles[1] is 2 for $marker_name\n");
	$lookup{$marker_name}->{'1'} = $alleles[0];
	$lookup{$marker_name}->{'2'} = $alleles[1];
    }

    @marker_names = sort keys %lookup;
    my $site_count   = scalar @marker_names;
    # where the final data will be stored
    my %stats_for_sites;

    # standard way of generating pairwise combos
    # LD is done by comparing all the pairwise site (marker)
    # combinations and keeping track of the genotype and 
    # pairwise genotype (ie genotypes of the 2 sites) frequencies
    for( my $i = 0; $i < $site_count - 1; $i++ ) {
	my $site1 = $marker_names[$i];

	for( my $j = $i+1; $j < $site_count ; $j++) { 	 
	    my (%genotypes, %total_genotype_count,$total_pairwisegeno_count,
		%pairwise_genotypes);
	 
	    my $site2 = $marker_names[$j];
	    my (%allele_count,%allele_freqs) = (0,0);
	    foreach my $ind ( @inds ) {
		# build string of genotype at site 1
		my ($genotype1) = $ind->get_Genotypes(-marker => $site1);
		my @alleles1  = sort $genotype1->get_Alleles;

                # if an individual has only one available allele
		# (has a blank or N for one of the chromosomes)
		# we don't want to use it in our calculation

		next unless( scalar @alleles1 == 2);
		my $genostr1  = join(',', @alleles1);

		# build string of genotype at site 2
		my ($genotype2) = $ind->get_Genotypes(-marker => $site2);
		my @alleles2  = sort $genotype2->get_Alleles;
		my $genostr2  = join(',', @alleles2);
		
		next unless( scalar @alleles2 == 2);
		for (@alleles1) {
		    $allele_count{$site1}++;
		    $allele_freqs{$site1}->{$_}++;
		}
		$genotypes{$site1}->{$genostr1}++;
		$total_genotype_count{$site1}++;

		for (@alleles2) {
		    $allele_count{$site2}++;
		    $allele_freqs{$site2}->{$_}++;
		}
		$genotypes{$site2}->{$genostr2}++;
		$total_genotype_count{$site2}++;

		# We are using the $site1,$site2 to signify
		# a unique key
		$pairwise_genotypes{"$genostr1,$genostr2"}++;
		# some individuals 
		$total_pairwisegeno_count++;
	    }
	    for my $site ( %allele_freqs ) {
		for my $al ( keys %{ $allele_freqs{$site} } ) {
		    $allele_freqs{$site}->{$al} /= $allele_count{$site};
		}
	    }
	    my $n = $total_pairwisegeno_count;	# number of pairs of comparisons
	    # 'A' and 'B' are two loci or in our case site1 and site2  
	    my $allele1_site1 = $lookup{$site1}->{'1'};	# this is the BigA allele
	    my $allele1_site2 = $lookup{$site2}->{'1'};	# this is the BigB allele
	    my $allele2_site1 = $lookup{$site1}->{'2'};	# this is the LittleA allele
	    my $allele2_site2 = $lookup{$site2}->{'2'};	# this is the LittleB allele
	    # AABB
	    my $N1genostr = join(",",( $allele1_site1, $allele1_site1,
				       $allele1_site2, $allele1_site2));
	    $self->debug(" [$site1,$site2](AABB) N1genostr=$N1genostr\n");
	    # AABb
	    my $N2genostr = join(",",( $allele1_site1, $allele1_site1,
				       $allele1_site2, $allele2_site2));
	    $self->debug(" [$site1,$site2](AABb) N2genostr=$N2genostr\n");
	    # AaBB
	    my $N4genostr = join(",",( $allele1_site1, $allele2_site1,
				       $allele1_site2, $allele1_site2));
	    $self->debug(" [$site1,$site2](AaBB) N4genostr=$N4genostr\n");
	    # AaBb
	    my $N5genostr = join(",",( $allele1_site1, $allele2_site1,
				       $allele1_site2, $allele2_site2));
	    $self->debug(" [$site1,$site2](AaBb) N5genostr=$N5genostr\n");
	    # count of AABB in 
	    my $n1 = $pairwise_genotypes{$N1genostr} || 0;
	    # count of AABb in 
	    my $n2 = $pairwise_genotypes{$N2genostr} || 0;
	    # count of AaBB in 
	    my $n4 = $pairwise_genotypes{$N4genostr} || 0;
	    # count of AaBb in 
	    my $n5 = $pairwise_genotypes{$N5genostr} || 0;

	    my $homozA_site1 = join(",", ($allele1_site1,$allele1_site1));
	    my $homozB_site2 = join(",", ($allele1_site2,$allele1_site2));
	    my $p_AA = ($genotypes{$site1}->{$homozA_site1} || 0) / $n;
	    my $p_BB = ($genotypes{$site2}->{$homozB_site2} || 0) / $n;
	    my $p_A  = $allele_freqs{$site1}->{$allele1_site1} || 0;	# an individual allele freq
	    my $p_a  =  1 - $p_A;

	    my $p_B  = $allele_freqs{$site2}->{$allele1_site2} || 0;	# an individual allele freq
	    my $p_b  =  1 - $p_B;

	    # variance of allele frequencies
	    my $pi_A = $p_A * $p_a;
	    my $pi_B = $p_B * $p_b;

	    # hardy weinberg
	    my $D_A  = $p_AA - $p_A**2;
	    my $D_B  = $p_BB - $p_B**2;
	    my $n_AB = 2*$n1 + $n2 + $n4 + 0.5 * $n5;
	    $self->debug("n_AB=$n_AB -- n1=$n1, n2=$n2 n4=$n4 n5=$n5\n");

	    my $delta_AB = (1 / $n ) * ( $n_AB ) - ( 2 * $p_A * $p_B );
	    $self->debug("delta_AB=$delta_AB -- n=$n, n_AB=$n_AB p_A=$p_A, p_B=$p_B\n");
	    $self->debug(sprintf(" (%d * %.4f) / ( %.2f + %.2f) * ( %.2f + %.2f) \n",
				 $n,$delta_AB**2, $pi_A, $D_A, $pi_B, $D_B));
	    
	    my $chisquared;
	    eval { $chisquared = ( $n * ($delta_AB**2) ) / 
		       ( ( $pi_A + $D_A) * ( $pi_B + $D_B) );
	       };
	    if( $@ ) {
		$self->debug("Skipping the site because the denom is 0.\nsite1=$site1, site2=$site2 : pi_A=$pi_A, pi_B=$pi_B D_A=$D_A, D_B=$D_B\n");
		next;
	    }
	    # this will be an upper triangular matrix
	    $stats_for_sites{$site1}->{$site2} = [$delta_AB,$chisquared];
	}
    }
    return %stats_for_sites;
}

=head2 mcdonald_kreitman

 Title   : mcdonald_kreitman
 Usage   : $Fstat = mcdonald_kreitman($ingroup, $outgroup);
 Function: Calculates McDonald-Kreitman statistic based on a set of ingroup
           individuals and an outgroup by computing the number of 
           differences at synonymous and non-synonymous sites
           for intraspecific comparisons and with the outgroup 
 Returns : 2x2 table, followed by a hash reference indicating any 
           warning messages about the status of the alleles or codons 
 Args    : -ingroup    => L<Bio::PopGen::Population> object or 
                          arrayref of L<Bio::PopGen::Individual>s 
           -outgroup   => L<Bio::PopGen::Population> object or 
                          arrayef of L<Bio::PopGen::Individual>s
           -polarized  => Boolean, to indicate if this should be 
                          a polarized test. Must provide two individuals 
                          as outgroups.

=cut

sub mcdonald_kreitman {
    my ($self,@args) = @_;
    my ($ingroup, $outgroup,$polarized) = 
	$self->_rearrange([qw(INGROUP OUTGROUP POLARIZED)],@args);
    my $verbose = $self->verbose;
    my $outgroup_count;
    my $gapchar = '\-';
    if( ref($outgroup) =~ /ARRAY/i ) {
	$outgroup_count = scalar @$outgroup;
    } elsif( UNIVERSAL::isa($outgroup,'Bio::PopGen::PopulationI') ) {
	$outgroup_count = $outgroup->get_number_individuals;
    } else {
	$self->throw("Expected an ArrayRef of Individuals OR a Bio::PopGen::PopulationI");
    }
	
    if( $polarized ) {
	if( $outgroup_count < 2 ) {
	    $self->throw("Need 2 outgroups with polarized option\n");
	}
    } elsif( $outgroup_count > 1 ) {
	$self->warn(sprintf("%s outgroup sequences provided, but only first will be used",$outgroup_count ));
    } elsif( $outgroup_count == 0 ) {
	$self->throw("No outgroup sequence provided");
    }
    
    my $codon_path = Bio::MolEvol::CodonModel->codon_path;
    
    my (%marker_names,%unique,@inds);
    for my $p ( $ingroup, $outgroup)  {
	if( ref($p) =~ /ARRAY/i ) {
	    push @inds, @$p;
	} else {
	    push @inds, $p->get_Individuals;
	}
    }
    for my $i ( @inds ) {
	if( $unique{$i->unique_id}++ ) {
	    $self->warn("Individual ". $i->unique_id. " is seen more than once in the ingroup or outgroup set\n");
	}
	for my $n ( $i->get_marker_names ) {
	    $marker_names{$n}++;
	}
    }

    my @marker_names = keys %marker_names;
    if( $marker_names[0] =~ /^(Site|Codon)/ ) {
	# sort by site or codon number and do it in 
	# a schwartzian transformation baby!
	@marker_names = map { $_->[1] } 
	sort { $a->[0] <=> $b->[0] }
	map { [$_ =~ /^(?:Codon|Site)-(\d+)/, $_] } @marker_names;
    }


    my $num_inds = scalar @inds;
    my %vals = ( 'ingroup'  => $ingroup,
		 'outgroup' => $outgroup,		 
		 );

    # Make the Codon Table type a parameter!
    my $table = Bio::Tools::CodonTable->new(-id => $codon_table);
    my @vt = qw(outgroup ingroup);
    my %changes;
    my %status;
    my %two_by_two = ( 'fixed_N' => 0,
		       'fixed_S' => 0,
		       'poly_N'  => 0,
		       'poly_S'  => 0);

    for my $codon ( @marker_names ) {
	my (%codonvals);
	my %all_alleles;
	for my $t ( @vt ) {
	    my $outcount = 1;
	    for my $ind ( @{$vals{$t}} ) {
		my @alleles = $ind->get_Genotypes($codon)->get_Alleles;
		if( @alleles > 2 ) {
		    warn("Codon $codon saw ", scalar @alleles, " alleles for ind ", $ind->unique_id, "\n");
		    die;
		} else {
		    my ($allele) = shift @alleles;
		    $all_alleles{$ind->unique_id} = $allele;
		    my $AA = $table->translate($allele);
		    next if( $AA eq 'X' || $AA eq '*' || $allele =~ /N/i);

		    my $label = $t;
		    if( $t eq 'outgroup' ) {
			$label = $t.$outcount++;
		    }
		    $codonvals{$label}->{$allele}++;
		    $codonvals{all}->{$allele}++;
		}
	    }
	}
	my $total = sum ( values %{$codonvals{'ingroup'}} );
	next if( $total && $total < 2 ); # skip sites with < alleles
	# process all the seen alleles (codons) 
	# this is a vertical slide through the alignment
	if( keys %{$codonvals{all}} <= 1 ) {
	    # no changes or no VALID codons - monomorphic
	} else { 
	    # grab only the first outgroup codon (what to do with rest?)
	    my ($outcodon) = keys %{$codonvals{'outgroup1'}};
            if( ! $outcodon ) { 
		$status{"no outgroup codon $codon"}++;
		next;
	    }
	    my $out_AA = $table->translate($outcodon);
	    my ($outcodon2) = keys %{$codonvals{'outgroup2'}};
	    if( ($polarized && ($outcodon ne $outcodon2)) ||
		$out_AA eq 'X' || $out_AA eq '*' ) {
		# skip if outgroup codons are different 
		# (when polarized option is on)
		# or skip if the outcodon is STOP or 'NNN'
		if( $verbose > 0 ) {
		    $self->debug("skipping $out_AA and $outcodon $outcodon2\n");
		}
		$status{'outgroup codons different'}++;
		next;
	    }

	    # check if ingroup is actually different from outgroup -
	    # if there are the same number of alleles when considering
	    # ALL or just the ingroup, then there is nothing new seen
	    # in the outgroup so it must be a shared allele (codon)

	    # so we just count how many total alleles were seen
	    # if this is the same as the number of alleles seen for just 
	    # the ingroup then the outgroup presents no new information

	    my @ingroup_codons = keys %{$codonvals{'ingroup'}};
	    my $diff_from_out = ! exists $codonvals{'ingroup'}->{$outcodon};

	    if( $verbose > 0 ) {
		$self->debug("alleles are in: ", join(",", @ingroup_codons),
			     " out: ", join(",", keys %{$codonvals{outgroup1}}),
			     " diff_from_out=$diff_from_out\n");

		for my $ind ( sort keys %all_alleles ) {
		    $self->debug( "$ind\t$all_alleles{$ind}\n");
		}
	    }
	    # are all the ingroup alleles the same and diferent from outgroup?
	    # fixed differences between species
	    if( $diff_from_out ) {
		if( scalar @ingroup_codons == 1 ) { 
		    # fixed differences
		    if( $outcodon =~ /^$gapchar/ ) {
			$status{'outgroup codons with gaps'}++;
			next;
		    } elsif( $ingroup_codons[0] =~ /$gapchar/) {
			$status{'ingroup codons with gaps'}++;
			next;
		    }
		    my $path = $codon_path->{uc $ingroup_codons[0].$outcodon};
		    $two_by_two{fixed_N} += $path->[0];
		    $two_by_two{fixed_S} += $path->[1];
		    if( $verbose > 0 ) {
			$self->debug("ingroup is @ingroup_codons outcodon is $outcodon\n");
			$self->debug("path is ",join(",",@$path),"\n");
			$self->debug
			    (sprintf("%-15s fixeddiff - %s;%s(%s) %d,%d\tNfix=%d Sfix=%d Npoly=%d Spoly=%s\n",$codon,$ingroup_codons[0], $outcodon,$out_AA,
				     @$path, map { $two_by_two{$_} } 
				     qw(fixed_N fixed_S poly_N poly_S)));
		    }
		} else { 
		    # polymorphic and all are different from outgroup
		    # Here we find the minimum number of NS subst
		    my ($Ndiff,$Sdiff) = (3,0);	# most different path
		    for my $c ( @ingroup_codons ) {
			next if( $c =~ /$gapchar/ || $outcodon =~ /$gapchar/);
			my $path = $codon_path->{uc $c.$outcodon};
			my ($tNdiff,$tSdiff) = @$path;
			if( $path->[0] < $Ndiff ||
			    ($tNdiff == $Ndiff &&
			     $tSdiff  <= $Sdiff)) {
			    ($Ndiff,$Sdiff) = ($tNdiff,$tSdiff);
			}
		    }
		    $two_by_two{fixed_N} += $Ndiff;
		    $two_by_two{fixed_S} += $Sdiff;
	            if( @ingroup_codons > 2 ) { 
			$status{"more than 2 ingroup codons $codon"}++;
			warn("more than 2 ingroup codons (@ingroup_codons)\n");	
		    } else {
		    	my $path = $codon_path->{uc join('',@ingroup_codons)};

		    	$two_by_two{poly_N} += $path->[0];
		    	$two_by_two{poly_S} += $path->[1];
		    	if( $verbose > 0 ) {
			    $self->debug(sprintf("%-15s polysite_all - %s;%s(%s) %d,%d\tNfix=%d Sfix=%d Npoly=%d Spoly=%s\n",$codon,join(',',@ingroup_codons), $outcodon,$out_AA,@$path, map { $two_by_two{$_} } qw(fixed_N fixed_S poly_N poly_S)));
			}
		    } 
		} 
	    } else {
		my %unq = map { $_ => 1 } @ingroup_codons;
		delete $unq{$outcodon};
		my @unique_codons = keys %unq;

		# calc path for diff add to poly
		# Here we find the minimum number of subst bw
		# codons
		my ($Ndiff,$Sdiff) = (3,0); # most different path
		for my $c ( @unique_codons ) {
		    my $path = $codon_path->{uc $c.$outcodon };
		    if( ! defined $path ) {
			die " cannot get path for ", $c.$outcodon, "\n";
		    }
		    my ($tNdiff,$tSdiff) = @$path;
		    if( $path->[0] < $Ndiff ||
			($tNdiff == $Ndiff &&
			 $tSdiff  <= $Sdiff)) {
			($Ndiff,$Sdiff) = ($tNdiff,$tSdiff);
		    }
		}

		if( @unique_codons == 2 ) {
		    my $path = $codon_path->{uc join('',@unique_codons)};
		    if( ! defined $path ) {
			$self->throw("no path for @unique_codons\n");
		    }
		    $Ndiff += $path->[0];
		    $Sdiff += $path->[1];
		}
		$two_by_two{poly_N} += $Ndiff;
		$two_by_two{poly_S} += $Sdiff;
		if( $verbose > 0 ) {
		    $self->debug(sprintf("%-15s polysite - %s;%s(%s) %d,%d\tNfix=%d Sfix=%d Npoly=%d Spoly=%s\n",$codon,join(',',@ingroup_codons), $outcodon,$out_AA,
					 $Ndiff, $Sdiff, map { $two_by_two{$_} } 
					 qw(fixed_N fixed_S poly_N poly_S)));
		}
	    }
	}	    
    }
    return ( $two_by_two{'poly_N'},
	     $two_by_two{'fixed_N'},
	     $two_by_two{'poly_S'},
	     $two_by_two{'fixed_S'},
	     {%status});
    
}

*MK = \&mcdonald_kreitman;


=head2 mcdonald_kreitman_counts

 Title   : mcdonald_kreitman_counts
 Usage   : my $MK = $statistics->mcdonald_kreitman_counts(

             N_poly -> integer of count of non-syn polymorphism
             N_fix  -> integer of count of non-syn fixed substitutions
             S_poly -> integer of count of syn polymorphism
             S_fix  -> integer of count of syn fixed substitutions
							  );
 Function:
 Returns : decimal number
 Args    : 

=cut


sub mcdonald_kreitman_counts {
    my ($self,$Npoly,$Nfix,$Spoly,$Sfix) = @_;
    if( $has_twotailed ) {
	return &Text::NSP::Measures::2D::Fisher2::twotailed::calculateStatistic 
	    (n11=>$Npoly,
	     n1p=>$Npoly+$Spoly,
	     np1=>$Npoly+$Nfix,
	     npp=>$Npoly+$Nfix+$Spoly+$Sfix);
    } else {
	$self->warn("cannot call mcdonald_kreitman_counts because no Fisher's exact is available - install Text::NSP::Measures::2D::Fisher2::twotailed");
	return 0;
    }
}


1;
