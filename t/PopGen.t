# -*-Perl-*- mode for emacs
# $Id$

# This will outline many tests for the population genetics
# objects in the Bio::PopGen namespace

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
    $NTESTS = 19;
    $error = 0;

    use Test;
    plan tests => $NTESTS; 

}

if( $error == 1 ) {
    exit(0);
}


use Bio::PopGen::Individual;
use Bio::PopGen::Genotype;
use Bio::PopGen::Population;

my @individuals = ( new Bio::PopGen::Individual(-unique_id => '10a'));
ok($individuals[0]);

my @genotypes = ( new Bio::PopGen::Genotype(-marker_name    => 'Mkr1',
					    -individual_id  => '10a',
					    -alleles => [ qw(A a)]),
		  new Bio::PopGen::Genotype(-marker_name    => 'Mkr2',
					    -individual_id  => '10a',
					    -alleles => [ qw(B B)]),
		  new Bio::PopGen::Genotype(-marker_name    => 'Mkr3',
					    -individual_id  => '10a',
					    -alleles => [ qw(A a)]));
ok(($genotypes[1]->get_Alleles)[0], 'B');

$individuals[0]->add_Genotype(@genotypes);
ok($individuals[0]->get_Genotypes,3);
ok($individuals[0]->get_Genotypes(-marker => 'Mkr3')->get_Alleles(),2);
my @alleles = $individuals[0]->get_Genotypes(-marker => 'Mkr2')->get_Alleles();
ok($alleles[0], 'B');

					     
my $population = new Bio::PopGen::Population(-name        => 'TestPop1',
					     -source      => 'testjasondata',
					     -description => 'throw away example',
					     -individuals => \@individuals);

ok(scalar ($population->get_Individuals()), 1);
ok($population->name, 'TestPop1');
ok($population->source, 'testjasondata');
ok($population->description, 'throw away example');

my @genotypes2 = ( new Bio::PopGen::Genotype(-marker_name   => 'Mkr1',
					     -individual_id => '11',
					     -alleles       => [ qw(A A)]),
		   new Bio::PopGen::Genotype(-marker_name   => 'Mkr2',
					     -individual_id => '11',
					     -alleles       => [ qw(B B)]),
		   new Bio::PopGen::Genotype(-marker_name   => 'Mkr3',
					     -individual_id => '11',
					     -alleles       => [ qw(a a)]),
		   new Bio::PopGen::Genotype(-marker_name   => 'Mkr4',
					     -individual_id => '11',
					     -alleles       => [ qw(C C)])
		   );
push @individuals, new Bio::PopGen::Individual(-genotypes   => \@genotypes2,
					       -unique_id   => '11');
$population->add_Individual($individuals[1]);

ok(scalar ($population->get_Individuals()), 2);
my ($found_ind) = $population->get_Individuals(-unique_id => '10a');
ok($found_ind->unique_id, '10a');
ok(scalar($population->get_Individuals(-marker => 'Mkr4')) , 1);
ok(scalar($population->get_Individuals(-marker => 'Mkr3')) , 2);

my @g = $population->get_Genotypes(-marker => 'Mkr4');

ok($g[0]->individual_id, '11');
ok(($g[0]->get_Alleles())[0], 'C');

my $marker = $population->get_Marker('Mkr3');
ok($marker);

ok($marker->get_Alleles,2);
my %af = $marker->get_Allele_Frequencies();
ok($af{'a'}, 0.75);
ok($af{'A'}, 0.25);

