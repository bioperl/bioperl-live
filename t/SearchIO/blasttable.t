# -*-Perl-*- Test Harness script for Bioperl
# $Id: SearchIO.t 14995 2008-11-16 06:20:00Z cjfields $

use strict;

BEGIN {
	use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 154);
	
	use_ok('Bio::SearchIO');
    use_ok('Bio::Search::SearchUtils');
}

my ($searchio, $result,$iter,$hit,$hsp);

# test blasttable output
my @eqset = qw( c200-vs-yeast.BLASTN.m9);
$searchio = Bio::SearchIO->new(-file => test_input_file('c200-vs-yeast.BLASTN'),
			      -format => 'blast');
$result = $searchio->next_result;
isa_ok($result,'Bio::Search::Result::ResultI');
my %ref = &Bio::Search::SearchUtils::result2hash($result);
is( scalar keys %ref, 67);
$searchio = Bio::SearchIO->new(-file => test_input_file('c200-vs-yeast.BLASTN.m8'),
			      -program_name => 'BLASTN',
			      -format => 'blasttable');
$result = $searchio->next_result;
my %tester = &Bio::Search::SearchUtils::result2hash($result);
is( scalar keys %tester, 67);
foreach my $key ( sort keys %ref ) {
    is($tester{$key}, $ref{$key},$key);
}

# test WU-BLAST blasttable output
$searchio = Bio::SearchIO->new(-file => test_input_file('test1.wublastp'),
			      			   -format => 'blast');
$result = $searchio->next_result;
isa_ok($result,'Bio::Search::Result::ResultI');
my %wuref = &Bio::Search::SearchUtils::result2hash($result);
is( scalar keys %wuref, 31);
$searchio = Bio::SearchIO->new(-file => test_input_file('test1.blasttab3'),
			      -program_name => 'BLASTP',
			      -format => 'blasttable');
$result = $searchio->next_result;
my %wutester = &Bio::Search::SearchUtils::result2hash($result);
is( scalar keys %wutester, 31);
foreach my $key ( sort keys %ref ) {
    is($wutester{$key}, $wuref{$key},$key);
}

# BLAST 2.2.18+ tabular output (has 13 columns instead of 12)
$searchio = Bio::SearchIO->new(-format => 'blasttable',
							  -file   => test_input_file('2008.blasttable'));

while(my $res = $searchio->next_result) {
    is($res->query_name, 'gi|1786183|gb|AAC73113.1|');
    is($res->algorithm, 'BLASTP');
    is($res->algorithm_version, '2.2.18+');
    my $hit = $res->next_hit;
    is($hit->name, 'gi|34395933|sp|P00561.2|AK1H_ECOLI');
    $hit = $res->next_hit;
    my $hsp = $hit->next_hsp;
    is($hsp->bits, 331);
    is($hsp->evalue, '2e-91');
    is($hsp->start('hit'), 16);
    is($hsp->end('hit'), 805);
    is($hsp->start('query'), 5);
    is($hsp->end('query'), 812);
    is($hsp->length, 821);
    is($hsp->gaps, 14);
}
