# -*-Perl-*- Test Harness script for Bioperl
# $Id: SearchIO_erpin.t 14672 2008-04-22 21:42:50Z cjfields $

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 91);
    
    use_ok('Bio::SearchIO');
}

my ($searchio, $result, $iter, $hit, $hsp, $algorithm, $meta);

$searchio = Bio::SearchIO->new( -format => 'erpin',
                                -file   => test_input_file('testfile.erpin'),
                                -model => 'stem-loop',
                                -query_acc => 'test',
                                -version => 5.5,
                                -verbose => 2
                               );
$result = $searchio->next_result;
isa_ok($result, 'Bio::Search::Result::ResultI');
$algorithm = $result->algorithm;
is($result->algorithm, 'ERPIN', "Result $algorithm");
is($result->algorithm_reference, undef, "Result $algorithm reference");
is($result->algorithm_version, 5.5, "Result $algorithm version");
is($result->available_parameters, 2, "Result parameters");
is($result->available_statistics, 1, "Result statistics");
is($result->database_entries, '', "Result entries");
is($result->database_letters, '', "Result letters");
is($result->database_name, 'AE016879.fna', "Result database_name");
is($result->num_hits, 1, "Result num_hits");
is($result->program_reference, undef, "Result program_reference");
is($result->query_accession, 'test', "Result query_accession");
is($result->query_description, '40 sequences of length 43', "Result query_description");
is($result->query_name, '/home/Administrator/pyrR.epn', "Result query_name");

$hit = $result->next_hit;
isa_ok($hit, 'Bio::Search::Hit::HitI');
is($hit->accession, 'AE016879.1', "Hit accession");
is($hit->ncbi_gi, '30260185', "Hit GI");
is($hit->algorithm, 'ERPIN', "Hit algorithm");
is($hit->bits, 31.64, "Hit bits");
is($hit->description, 'Bacillus anthracis str. Ames, complete genome',
   "Hit description"); # no hit descs yet
is($hit->length, 0, "Hit length");
is($hit->locus, '', "Hit locus");
is($hit->n, 4, "Hit n");
is($hit->name, 'gi|30260185|gb|AE016879.1|', "Hit name");
is($hit->num_hsps, 4, "Hit num_hsps"); 
is($hit->overlap, 0, "Hit overlap");
is($hit->query_length, undef, "Hit query_length");
is($hit->rank, 1, "Hit rank");
is($hit->raw_score, 31.64, "Hit raw_score");
is($hit->score, 31.64, "Hit score");
float_is($hit->significance, 4.44e-06);

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'ERPIN', "HSP algorithm");
float_is($hsp->evalue, 1.68e-05);
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
is($hsp->frame('query'), 0, "HSP frame");
is($hsp->gaps, 2, "HSP gaps");
isa_ok($hsp->hit, 'Bio::SeqFeature::Similarity', "HSP hit");
is($hsp->hit_string,
   'CTTT.aacc--.CAACC.CCGTGA.GGTTG.a.GAAG',
   "HSP hit_string");
is($hsp->homology_string,
   '',
   "HSP homology_string");
is($hsp->hsp_group, undef, "HSP hsp_group");
is($hsp->hsp_length, 37, "HSP hsp_length");
is($hsp->length, 37, "HSP length");
is($hsp->links, undef, "HSP links");
isa_ok($hsp->query, 'Bio::SeqFeature::Similarity', "HSP query");
is($hsp->query_string,
   '',
   "HSP query_string");
is($hsp->range, 37, "HSP range");
is($hsp->rank, 1, "HSP rank");
float_is($hsp->significance, 1.68e-05);
float_is($hsp->expect, '1.68e-05', "HSP expect");
isa_ok($hsp->seq, 'Bio::LocatableSeq');
is($hsp->seq_str,
   '',
   "HSP seq_str");
is($hsp->start, 1, "HSP start");
is($hsp->custom_score, undef, "HSP custom_score");
is($hsp->meta, undef, "HSP meta");
is($hsp->meta, undef);
is($hsp->meta, undef);
is($hsp->strand('hit'), 1, "HSP strand");
is($hsp->meta, undef);
is($hsp->meta, undef);

# ERPIN lacks sequence for query, will spit back a warning..
eval{$hsp->get_aln};
like($@, qr'Missing query string, can\'t build alignment','ERPIN get_aln warning');

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'ERPIN', "HSP algorithm");
float_is($hsp->evalue, 5.61e-05);
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
is($hsp->frame('query'), 0, "HSP frame");
is($hsp->gaps, 1, "HSP gaps");
isa_ok($hsp->hit, 'Bio::SeqFeature::Similarity', "HSP hit");
is($hsp->hit_string,
   'CTTT.taatt-.CAGTC.CTGTGA.GACCG.g.AAAG',
   "HSP hit_string");
is($hsp->homology_string,
   '',
   "HSP homology_string");
is($hsp->query_string,
   '',
   "HSP query_string");
is($hsp->hsp_group, undef, "HSP hsp_group");
is($hsp->hsp_length, 37, "HSP hsp_length");
is($hsp->length, 37, "HSP length");
is($hsp->links, undef, "HSP links");
isa_ok($hsp->query, 'Bio::SeqFeature::Similarity');
is($hsp->range, 37, "HSP range");
is($hsp->rank, 2, "HSP rank");
float_is($hsp->significance, 5.61e-05);
is($hsp->end, 37, "HSP end");
float_is($hsp->expect, '5.61e-05', "HSP expect");
isa_ok($hsp->seq, 'Bio::LocatableSeq');
is($hsp->seq_str,   '',   "HSP seq_str");
is($hsp->start, 1, "HSP start");
is($hsp->custom_score, undef, "HSP custom_score");
is($hsp->meta, undef);
is($hsp->meta, undef);
is($hsp->strand('hit'), -1, "HSP strand");
