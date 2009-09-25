# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 9);
	
	use_ok('Bio::Symbol::Symbol');
}

ok my $thymine = Bio::Symbol::Symbol->new(-name => 'Arg',
				      -token=> 'R');
my $a = Bio::Symbol::Symbol->new(-token => 'A' );
my $u = Bio::Symbol::Symbol->new(-token => 'U' );
my $g = Bio::Symbol::Symbol->new(-token => 'G' );

is($thymine->name, "Arg");
is($thymine->token, "R");
my $M = Bio::Symbol::Symbol->new(-name  => 'Met',
				-token => 'M',
				-symbols => [ $a, $u, $g ]);

is($M->name, "Met");
is($M->token, 'M');
my @symbols = $M->symbols;
my @expected = ($a, $u, $g);
foreach ( 0..2 ) {
    is($expected[$_], $symbols[$_]);
}
