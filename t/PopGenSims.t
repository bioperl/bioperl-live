# -*-Perl-*- mode for emacs
# $Id$

# This will outline tests for the population genetics simulation
# in the Bio::PopGen::Simulation namespace
# Coalescent has its own tests though in t/Coalescent.t

my $error;

use strict;
use lib '.';

BEGIN {     
    eval { require Test::More; };
    if( $@ ) {
        use lib 't/lib';
    }
    use vars qw($NTESTS);
    $NTESTS = 23;
    $error = 0;

    use Test::More;
    plan tests => $NTESTS; 
    use_ok 'Bio::PopGen::Simulation::GeneticDrift';
}

my $sim = new Bio::PopGen::Simulation::GeneticDrift(-popsize => 40,
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

$sim = new Bio::PopGen::Simulation::GeneticDrift(-popsize => 50,
						 -alleles => {A => 0.2,
							      B => 0.3,
							      C => 0.5,
							  });

for my $i (0..9) {
    my %f = $sim->next_generation;
    is(($f{'A'}||0) + ($f{'B'}||0) + ($f{'C'}||0), 1, 'Allele freqs should sum to 1');
}
