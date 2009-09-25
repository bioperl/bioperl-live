# -*-Perl-*- Test Harness script for Bioperl
# $Id: SearchIO_megablast.t 14995 2008-11-16 06:20:00Z cjfields $

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 31);
    
    use_ok('Bio::SearchIO');
}

my ($searchio, $result, $hit, $hsp);

# this is megablast output type 0 
my $in = Bio::SearchIO->new(-file          => test_input_file('503384.MEGABLAST.0'),
			-report_format => 0,
			-format        => 'megablast'); 
my $r = $in->next_result;
my @dcompare = ( 
	      ['Contig634', 7620, 7941, 1, 1, 321, -1],
	      ['Contig1853', 6406, 6620, 1, 1691, 1905, 1],  
	      ['Contig3700', 8723,9434, 1, 4083, 4794, -1],
	      ['Contig3997', 1282, 1704, 1, 1546, 1968,-1 ],
	      );

is($r->query_name, '503384');

while( my $hit = $r->next_hit ) {
    my $d = shift @dcompare;
    is($hit->name, shift @$d);
    my $hsp = $hit->next_hsp;
    is($hsp->query->start, shift @$d);
    is($hsp->query->end, shift @$d);
    is($hsp->query->strand, shift @$d);
    is($hsp->hit->start, shift @$d);
    is($hsp->hit->end, shift @$d);
    is($hsp->hit->strand, shift @$d);
}
is(@dcompare, 0);

