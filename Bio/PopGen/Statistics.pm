# $Id$
#
# BioPerl module for Bio::PopGen::Statistics
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
  use Bio::PopGen::Simulation::Coalescent;
 
  my $sim = new Bio::PopGen::Simulation::Coalescent( -samplesample_size => 12);

  my $tree = $sim->next_tree;
  
  $factory->add_Mutations($tree,20);

  my $stats = new Bio::PopGen::Statistics();
  my $pi = $stats->pi($tree);
  my $D = $stats->tajima_d($tree);
  

=head1 DESCRIPTION

This object is intended to provide implementations some standard
population genetics statistics about alleles in populations.

This module was previously named Bio::Tree::Statistics.

This object is a place to accumulate routines for calculating various
statistics from the coalescent simulation, marker/allele, or from
aligned sequence data given that you can calculate alleles, number of
segregating sites.

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

=head1 AUTHOR - Jason Stajich, Matthew Hahn

Email jason-at-bioperl-dot-org
Matt Hahn E<lt>matthew.hahn-at-duke.dukeE<gt>

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::PopGen::Statistics;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;

@ISA = qw(Bio::Root::Root );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::PopGen::Statistics();
 Function: Builds a new Bio::PopGen::Statistics object 
 Returns : an instance of Bio::PopGen::Statistics
 Args    : none


=cut


=head2 fu_and_li_D

 Title   : fu_and_li_D
 Usage   : my $D = $statistics->fu_an_li_D(\@ingroup,$extmutations);
 Function: Fu and Li D statistic for a list of individuals
           given an outgroup and the number of external mutations
           (either provided or calculated from list of outgroup individuals)
 Returns : decimal
 Args    : $individuals - array refence which contains ingroup individuals 
           (L<Bio::PopGen::Individual> or derived classes)
           $extmutations - number of external mutations OR
           arrayref of outgroup individuals
=cut

sub fu_and_li_D{
    my ($self,$ingroup,$outgroup) = @_;
    if( ref($ingroup) !~ /ARRAY/i   ) {
	$self->throw("expected array references which is list of Bio::PopGen::IndividualI to fu_and_li_D");
    }
    my $ext_mutations;
    if( ! defined $outgroup ) {
	$self->warn("Need to provide either an array ref to the outgroup individuals or the number of external mutations");
	return 0;
    } elsif( ref($outgroup) =~ /ARRAY/i ) {
	$self->warn("Currently cannot calculate fu_and_li_F with a set of outgroup individuals");
        return 0;
    } else { 
	$ext_mutations = $outgroup;
    }
       
   my $seg_sites    = scalar $ingroup->[0]->get_marker_names;
   my $sample_size  = scalar(@$ingroup);

   my $tipmutcount;
   
   if( $seg_sites <= 0 ) { 
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
    my $D = ( $seg_sites - (  $a * $ext_mutations) ) /
            ( sqrt ( ($u * $seg_sites ) + ( $v * $seg_sites **2) ) );
 
    return $D;
}

=head2 fu_and_li_D_star

 Title   : fu_and_li_D_star
 Usage   : my $D = $statistics->fu_an_li_D_star(\@individuals);
 Function: Fu and Li's D star statistic for this set of samples
            Without an outgroup
 Returns : decimal number
 Args    : array ref of L<Bio::PopGen::IndividualI> objects

=cut

#'
# fu_and_li_D*

sub fu_and_li_D_star {
    my ($self,$individuals) = @_;
    if( ref($individuals) !~ /ARRAY/i ) {
	$self->throw("expected an array reference of a list of Bio::PopGen::IndividualI to fu_and_li_D_star");
    }
    my $sample_size = scalar @$individuals;
    my $pi          = $self->pi($individuals);
    my $seg_sites   = scalar ($individuals->[0]->get_marker_names());
    my ($total,$unique,$singleton) = $self->allele_count($individuals);

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

=head2 fu_and_li_F

 Title   : fu_and_li_F
 Usage   : my $D = Bio::PopGen::Statistics->fu_and_li_F(\@ingroup,$ext_muts);
 Function: Calculate Fu and Li's F on an ingroup with either the set of 
           outgroup individuals, or the number of external mutations
 Returns : decimal number
 Args    : array ref of L<Bio::PopGen::IndividualI> objects for the ingroup
           number of external mutations OR list of individuals for the outgroup

=cut
#'

sub fu_and_li_F {
    my ($self,$ingroup,$outgroup) = @_;

    if( ref($ingroup) !~ /ARRAY/i ) {
	$self->throw("expected an array reference of a list of Bio::PopGen::IndividualI to fu_and_li_F");
    }
    my $ext_mutations;
    if( ! defined $outgroup ) {
	$self->warn("Need to provide either an array ref to the outgroup individuals or the number of external mutations");
	return 0;
    } elsif( ref($outgroup) =~ /ARRAY/i ) {
	$self->warn("Currently cannot calculate fu_and_li_F with a set of outgroup individuals");
        return 0;
    } else { 
	$ext_mutations = $outgroup;
    }

    my $sample_size = scalar @$ingroup;
    my $pi          = $self->pi($ingroup);
    my $seg_sites   = scalar $ingroup->[0]->get_marker_names();
    my ($total,$unique,$singleton) = $self->allele_count($ingroup);
    
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

=head2 tajima_D

 Title   : tajima_D
 Usage   : my $D = Bio::PopGen::Statistics->tajima_D(\@samples);
 Function: Calculate Tajima's D on a set of samples 
 Returns : decimal number
 Args    : array ref of L<Bio::PopGen::IndividualI> objects 


=cut

#'

sub tajima_D {
    my ($self,$individuals) = @_;
    if( ref($individuals) !~ /ARRAY/i ) {
	$self->throw("expected an array reference of a list of Bio::PopGen::IndividualI to tajima_D");
    }
    my $sample_size = scalar @$individuals;
    # pi - all pairwise differences 
    my $pi          = $self->pi($individuals);  
    my $seg_sites   = scalar ($individuals->[0]->get_marker_names);
    my ($total,$unique,$singleton) = $self->allele_count($individuals);

    my $a1 = 0; 
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
    
    my $D = ( $pi - ( $seg_sites / $a1 ) ) / 
	sqrt ( ($e1 * $seg_sites) + (( $e2 * $seg_sites) * ( $seg_sites - 1)));

    return $D;
}

=head2 pi

 Title   : pi
 Usage   : my $pi = Bio::PopGen::Statistics->pi(\@inds)
 Function: Calculate pi (...explain here...) given a list of individuals 
           which have the same number of markers/sites/mutation as 
           available from the get_Genotypes() call in 
           L<Bio::PopGen::IndividualI>
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
    my (%data,@marker_names,$sample_size);

    if( ref($individuals) =~ /ARRAY/i ) {
	# one possible argument is an arrayref of Bio::PopGen::IndividualI objs
	@marker_names = $individuals->[0]->get_marker_names;
	$sample_size = scalar @$individuals;

	# Here we're calculating the allele frequencies
	my %marker_total;
	foreach my $ind ( @$individuals ) {
	    if( ! $ind->isa('Bio::PopGen::IndividualI') ) {
		$self->warn("Expected an arrayref of Bio::PopGen::IndividualI objects, this is a ".ref($ind)."\n");
	    }
	    foreach my $m ( @marker_names ) {
		foreach my $a (map { $_->get_Alleles} $ind->get_Genotypes($m) ) {
		    $data{$m}->{$a}++;
		    $marker_total{$m}++;
		}
	    }
	}
	while( my ($marker,$count) =  each %marker_total ) {
	    foreach my $c ( values %{$data{$marker}} ) {
		$c /= $count;
	    }
	}
	# %data will contain allele frequencies for each marker, allele
    } elsif( ref($individuals) && 
	     $individuals->isa('Bio::PopGen::PopulationI') ) {
	my $pop = $individuals;
	$sample_size = $pop->number_individuals;
	foreach my $marker( $pop->get_Markers ) {
	    push @marker_names, $marker->name;
	    $data{$marker->name} = [$marker->get_Allele_Frequencies];
	}
    } else { 
	$self->throw("expected an array reference of a list of Bio::PopGen::IndividualI to pi");
    }
    # doing all pairwise combinations

    # For now we assume that all individuals have the same markers
    my ($diffcount,$totalcompare) = (0,0);
    my $pi = 0;
    foreach my $markerdat ( values %data ) {
	my $totalalleles; # this will only be different among markers
	                  # when there is missing data
	my @alleles = keys %$markerdat;
	foreach my $al ( @alleles ) { $totalalleles += $markerdat->{$al} }
	for( my $i =0; $i < scalar @alleles -1; $i++ ) {
	    my ($a1,$a2) = ( $alleles[$i], $alleles[$i+1]);
	    $pi += $self->heterozygosity($sample_size, 
					 $markerdat->{$a1} / $totalalleles,
					 $markerdat->{$a2} / $totalalleles);
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
 Function: Calculates theta (...explanation here... ) from the sample size 
           and the number of segregating sites.
           Providing the third parameter, total number of sites will
           return theta per site          
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

sub theta {
    my $self = shift;    
    my ( $sample_size, $seg_sites,$totalsites) = @_;
    if( ref($sample_size) =~ /ARRAY/i ) {
	my $samps = $sample_size;
	$totalsites = $seg_sites; # only 2 arguments if one is an array
	my %data;
	my @marker_names = $samps->[0]->get_marker_names;
	# we need to calculate number of polymorphic sites
	foreach my $ind ( @$samps ) {
	    foreach my $m ( @marker_names ) {
		foreach my $a (map { $_->get_Alleles} 
			       $ind->get_Genotypes($m) ) {
		    $data{$m}->{$a}++;
		}
	    }
	}
	# if there is >1 allele then it is polymorphic
	$seg_sites = 0;
	foreach my $marker ( @marker_names ) {
	    $seg_sites++ if( keys %{$data{$marker}} > 1 );
	}
	$sample_size = scalar @$samps;
    } elsif(ref($sample_size) &&
	    $sample_size->isa('Bio::PopGen::PopulationI') ) {
	# This will handle the case when we pass in a PopulationI object
	my $pop = $sample_size;
	$totalsites = $seg_sites; # shift the arguments over by one
	$sample_size = $pop->number_individuals;
	$seg_sites = 0;
	foreach my $marker( $pop->get_Markers ) {
	    $seg_sites if ( scalar $marker->get_Alleles > 1 );
	}
    }
    my $a1 = 0; 
    for(my $k= 1; $k < $sample_size; $k++ ) {
	$a1 += ( 1 / $k );
    }    
    if( $totalsites ) { # 0 and undef are the same can't divide by them
	$seg_sites /= $totalsites;
    }
    return $seg_sites / $a1;
}


=head2 allele_count

 Title   : allele_count
 Usage   : my ($totalct,$unqiuect,
	       $singletons) = Bio::PopGen::Statistics->allele_count(\@inds)
 Function: Calculate the number of mutations/alleles which only occur once in
           a list of individuals for all sites/markers
 Returns : A triple, 
               first is total number of alleles (2* num markers for diploids)
                                             or just num of sites for haploids)
               second is total number of unique alleles (so if a single allele 
                    is only seen at a site (marker))
               third is number of alleles which only occur once (integer)

           Called in list context it will just return the total number 
              of alleles         
 Args    : arrayref of L<Bio::PopGen::IndividualI> objects


=cut

sub allele_count {
    my ($self,$individuals) = @_;
    my (%sites,$singleton_allele_ct,$total_allele_ct, $unique_allele_ct);
    # find number of sites where a particular allele is only seen once
    
    # first collect all the alleles into a hash structure
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
    foreach my $site ( values %sites ) { # don't really care what the namem is
	foreach my $allelect ( values %$site ) { # 
	    $total_allele_ct += $allelect;
            # find the sites which have an allele with only 1 copy
 	    $singleton_allele_ct++ if( $allelect == 1 );
	    $unique_allele_ct++;
	}
    }
    if( wantarray )  {
	return ($total_allele_ct, $unique_allele_ct,$singleton_allele_ct);
    } else {
	return $total_allele_ct;
    }
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

1;
