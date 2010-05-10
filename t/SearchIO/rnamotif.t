# -*-Perl-*- Test Harness script for Bioperl
# $Id: SearchIO_rnamotif.t 14672 2008-04-22 21:42:50Z cjfields $

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 60);
    
    use_ok('Bio::SearchIO');
}

my ($searchio, $result, $iter, $hit, $hsp, $algorithm, $meta);

#### RNAMotif ####

# regular data

$searchio = Bio::SearchIO->new( -format => 'rnamotif',
                                -file   => test_input_file('trna.strict.rnamotif'),
                                -model => 'trna.descr',
                                -query_acc => 'test',
                                -database => 'gbrna.fas',
                                -verbose => 2
                               );

$result = $searchio->next_result;
isa_ok($result, 'Bio::Search::Result::ResultI');
$algorithm = $result->algorithm;
is($result->algorithm, 'RNAMOTIF', "Result $algorithm");
is($result->algorithm_reference, undef, "Result $algorithm reference");
is($result->algorithm_version, '3.0.3', "Result $algorithm version");
is($result->database_entries, '', "Result entries");
is($result->database_letters, '', "Result letters");
is($result->database_name, 'gbrna.fas', "Result database_name");
is($result->num_hits, 28, "Result num_hits");
is($result->program_reference, undef, "Result program_reference");
is($result->query_accession, 'test', "Result query_accession");
is($result->query_description, 'h5 ss h5 ss h3 ss h5 ss h3 ss h5 ss h3 h3 ss', "Result query_description");
is($result->query_length, 0, "Result query_length");
is($result->query_name, 'trna.strict.descr', "Result query_name");

$hit = $result->next_hit;
isa_ok($hit, 'Bio::Search::Hit::HitI');
is($hit->accession, 'M10671', "Hit accession");
is($hit->ncbi_gi, '173683', "Hit GI");
is($hit->algorithm, 'RNAMOTIF', "Hit algorithm");
is($hit->description, 'Avian oncornavirus Trp-tRNA',
   "Hit description"); # no hit descs yet
is($hit->length, 0, "Hit length");
is($hit->locus, 'ACSTRW', "Hit locus");
is($hit->n, 8, "Hit n");
is($hit->name, 'gi|173683|gb|M10671|ACSTRW', "Hit name");
is($hit->num_hsps, 8, "Hit num_hsps"); 
is($hit->overlap, 0, "Hit overlap");
is($hit->rank, 1, "Hit rank");
is($hit->raw_score, 0, "Hit raw_score");
is($hit->score, 0, "Hit score");
float_is($hit->significance, undef);

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'RNAMOTIF', "HSP algorithm");
float_is($hsp->evalue, undef);
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
is($hsp->frame('query'), 0, "HSP frame");
is($hsp->gaps, 0, "HSP gaps");

# RNAMotif cannot build alignments
eval{$hsp->get_aln};
like($@, qr'Missing query string, can\'t build alignment','RNAMotif get_aln warning');
isa_ok($hsp->hit, 'Bio::SeqFeature::Similarity', "HSP hit");
is($hsp->hit_string,
   'gacctcgtggcgcaacggtagcgcgtctgactccagatcagaaggctgcgtgttcgaatcacgtcggggtcacca',
   "HSP hit_string");
is($hsp->homology_string,
   '',
   "HSP homology_string");
is($hsp->hsp_group, undef, "HSP hsp_group");
is($hsp->hsp_length, 75, "HSP hsp_length");
is($hsp->length, 75, "HSP length");
is($hsp->links, undef, "HSP links");
is($hsp->n, 1,"HSP n");
isa_ok($hsp->query, 'Bio::SeqFeature::Similarity', "HSP query");
is($hsp->query_string,
   '',
   "HSP query_string");
is($hsp->range, 75, "HSP range");
is($hsp->rank, 1, "HSP rank");
float_is($hsp->significance, undef);
is($hsp->end, 75, "HSP end");
float_is($hsp->expect, undef, "HSP expect");
isa_ok($hsp->seq, 'Bio::LocatableSeq');
is($hsp->seq_str,
   '',
   "HSP seq_str");
is($hsp->start, 1, "HSP start");
is($hsp->custom_score, undef, "HSP custom_score");
is($hsp->meta,
   '<<<<<<<..<<<.........>>>.<<<<<.......>>>>>.......<<<<.......>>>>>>>>>>>....',
   "HSP meta");
is($hsp->strand('hit'), 1, "HSP strand");
($meta) = $hsp->feature1->get_tag_values('meta');
is($meta, '<<<<<<<..<<<.........>>>.<<<<<.......>>>>>.......<<<<.......>>>>>>>>>>>....');
($meta) = $hsp->feature2->get_tag_values('meta');
is($meta, '<<<<<<<..<<<.........>>>.<<<<<.......>>>>>.......<<<<.......>>>>>>>>>>>....');

