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

    plan tests => 8;
}

use Bio::Symbol::Symbol;

my $thymine = new Bio::Symbol::Symbol(-name => 'Arg',
				      -token=> 'R');
my $a = new Bio::Symbol::Symbol(-token => 'A' );
my $u = new Bio::Symbol::Symbol(-token => 'U' );
my $g = new Bio::Symbol::Symbol(-token => 'G' );

ok($thymine);
ok($thymine->name, "Arg");
ok($thymine->token, "R");
my $M = new Bio::Symbol::Symbol(-name  => 'Met',
				-token => 'M',
				-symbols => [ $a, $u, $g ]);

ok($M->name, "Met");
ok($M->token, 'M');
my @symbols = $M->symbols;
my @expected = ($a, $u, $g);
foreach ( 0..2 ) {
    ok($expected[$_], $symbols[$_]);
}
