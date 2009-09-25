# -*-Perl-*- Test Harness script for Bioperl
# $Id$

# This will outline tests for the population genetics simulation
# in the Bio::PopGen::Simulation namespace
# Coalescent has its own tests though in t/Coalescent.t

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 23);
	
    use_ok 'Bio::PopGen::Simulation::GeneticDrift';
}

my $sim = Bio::PopGen::Simulation::GeneticDrift->new(-popsize => 40,
						    -alleles => {A => 0.2,
								 B => 0.8});

my (@Afreqs,@Bfreqs);
for my $i (0..9) {
    my %f = $sim->next_generation;
    push @Afreqs, $f{'A'};
    push @Bfreqs, $f{'B'};
    is(($f{'A'}||0) + ($f{'B'}||0), 1, 'Allele freqs should sum to 1');
}

is(@Afreqs, 10);
cmp_ok (($Afreqs[9]||0), '<=', 1, 'All frequencies should be <= 1');

$sim = Bio::PopGen::Simulation::GeneticDrift->new(-popsize => 50,
						 -alleles => {A => 0.2,
							      B => 0.3,
							      C => 0.5,
							  });

for my $i (0..9) {
    my %f = $sim->next_generation;
    is(($f{'A'}||0) + ($f{'B'}||0) + ($f{'C'}||0), 1, 'Allele freqs should sum to 1');
}
