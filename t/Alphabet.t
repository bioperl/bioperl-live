# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

use strict;

BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;

    plan tests => 6;
}

use Bio::Symbol::Alphabet;
use Bio::Symbol::Symbol;

my $a = new Bio::Symbol::Symbol(-token => 'A' );
my $u = new Bio::Symbol::Symbol(-token => 'U' );
my $g = new Bio::Symbol::Symbol(-token => 'G' );
my $t = new Bio::Symbol::Symbol(-token => 'T' );

my $rna = new Bio::Symbol::Alphabet( -symbols => [ $a, $u, $g, $t ] );
				     
ok($rna);
my @symbols = $rna->symbols;
ok(scalar @symbols, 4);

ok($rna->contains($a));
ok($rna->contains($t));
ok($rna->contains($u));
ok($rna->contains($g));

