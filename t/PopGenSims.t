# -*-Perl-*- mode for emacs
# $Id$

# This will outline tests for the population genetics simulation
# in the Bio::PopGen::Simulation namespace
# Coalescent has its own tests though in t/Coalescent.t

my $error;

use vars qw($SKIPXML $LASTXMLTEST); 
use strict;
use lib '.';

BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use vars qw($NTESTS);
    $NTESTS = 22;
    $error = 0;

    use Test;
    plan tests => $NTESTS; 

}

if( $error == 1 ) {
    exit(0);
}


use Bio::PopGen::Simulation::GeneticDrift;

my $sim = new Bio::PopGen::Simulation::GeneticDrift(-popsize => 40,
						    -alleles => {A => 0.2,
								 B => 0.8});

my (@Afreqs,@Bfreqs);
for(my $i =0 ;$i < 10; $i++ ) {
    my %f = $sim->next_generation;
    push @Afreqs, $f{'A'};
    push @Bfreqs, $f{'B'};
    ok(($f{'A'}||0) + ($f{'B'}||0), 1, 'Allele freqs should sum to 1');
}

ok(@Afreqs, 10);
ok(($Afreqs[9]||0) <= 1, 1, 'All frequencies should be <= 1');

$sim = new Bio::PopGen::Simulation::GeneticDrift(-popsize => 50,
						 -alleles => {A => 0.2,
							      B => 0.3,
							      C => 0.5,
							  });

for(my $i =0 ;$i < 10; $i++ ) {
    my %f = $sim->next_generation;
    ok(($f{'A'}||0) + ($f{'B'}||0) + ($f{'C'}||0), 1, 'Allele freqs should sum to 1');
}
