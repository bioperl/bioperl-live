# -*-Perl-*- Test Harness script for Bioperl
# $Id: SearchIO_wise.t 11733 2007-10-26 18:22:10Z jason $

use strict;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 20);
	
    use_ok('Bio::SearchIO');
}

my $parser = Bio::SearchIO->new(-file => test_input_file('genewise.out'),
			    -format   => 'wise',
			    -wisetype => 'genewise');
my $result = $parser->next_result;
my $hit = $result->next_hit;
is($result->query_name, 'SINFRUP00000067802');
is($hit->name, 'Scaffold_2042.1');

is($hit->score, 2054.68);
my $hsp = $hit->next_hsp;

is($hsp->query->start,22265);
is($hsp->query->end,22396);
is($hsp->query->strand,1);
is($hsp->query->score, 2054.68);

is($hsp->hit->start,1);
is($hsp->hit->end,44);
is($hsp->hit->strand,0);
is($hsp->hit->score, 2054.68);

$hsp = $hit->next_hsp;

is($hsp->query->start,24224);
is($hsp->query->end,24328);

is($hsp->hit->start,45);
is($hsp->hit->end,79);

$hsp = $hit->next_hsp;

is($hsp->query->start,24471);
is($hsp->query->end,24513);

is($hsp->hit->start,80);
is($hsp->hit->end,93);
