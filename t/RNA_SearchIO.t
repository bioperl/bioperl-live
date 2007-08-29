# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 326);
    
    use_ok('Bio::SearchIO');
}

my ($searchio, $result,$iter,$hit,$hsp, $algorithm, $meta);

### Infernal ####

$searchio = Bio::SearchIO->new( -format => 'infernal',
                                -file   => test_input_file('test.infernal'),
                                # version is reset to the correct one by parser
                                -version => 0.7, 
                                -model => 'Purine',
                                -query_acc => 'RF00167',
                                -query_desc => 'Purine riboswitch',
                                -database => 'b_sub.fas',
                                -hsp_minscore => 40,
                                -convert_meta => 0,
                                -verbose => 2
                               );

$result = $searchio->next_result;
isa_ok($result, 'Bio::Search::Result::ResultI');
$algorithm = $result->algorithm;
is($result->algorithm, 'CMSEARCH', "Result $algorithm");
is($result->algorithm_reference, undef, "Result $algorithm reference");
is($result->algorithm_version, 0.7, "Result $algorithm version");
is($result->available_parameters, 0, "Result parameters");
is($result->available_statistics, 0, "Result statistics");
is($result->database_entries, '', "Result entries");
is($result->database_letters, '', "Result letters");
is($result->database_name, 'b_sub.fas', "Result database_name");
is($result->num_hits, 2, "Result num_hits");
is($result->program_reference, undef, "Result program_reference");
is($result->query_accession, 'RF00167', "Result query_accession");
is($result->query_description, 'Purine riboswitch', "Result query_description");
is($result->query_length, 102, "Result query_length");
is($result->query_name, 'Purine', "Result query_name");

$hit = $result->next_hit;
isa_ok($hit, 'Bio::Search::Hit::HitI');
is($hit->ncbi_gi, '2239287', "Hit GI");
is($hit->accession, 'U51115.1', "Hit accession");
is($hit->algorithm, 'CMSEARCH', "Hit algorithm");
is($hit->bits, 81.29, "Hit bits");
is($hit->description, '', "Hit description"); # no hit descs yet
is($hit->locus, 'BSU51115', "Hit locus");
is($hit->n, 2, "Hit n");
is($hit->name, 'gi|2239287|gb|U51115.1|BSU51115', "Hit name");
is($hit->num_hsps, 2, "Hit num_hsps");

# These Bio::Search::Hit::HitI methods are currently unimplemented in
# Bio::Search::Hit::ModelHit; they may be integrated over time but will require
# some reconfiguring for Model-based searches

eval { $hit->length_aln() };
like($@, qr'length_aln not implemented for Model-based searches',
     "Hit length_aln() not implemented");
eval {$hit->num_unaligned_hit};
like($@, qr'num_unaligned_hit/num_unaligned_sbjct not implemented for Model-based searches',
     "Hit num_unaligned_hit() not implemented");
eval {$hit->num_unaligned_query};
like($@, qr'num_unaligned_query not implemented for Model-based searches',
     "Hit num_unaligned_query() not implemented");
eval {$hit->num_unaligned_sbjct};
like($@, qr'num_unaligned_hit/num_unaligned_sbjct not implemented for Model-based searches',
     "Hit num_unaligned_sbjct() not implemented");
eval {$hit->start};
like($@, qr'start not implemented for Model-based searches','Hit start not implemented');
eval {$hit->end};
like($@, qr'end not implemented for Model-based searches','Hit end not implemented');
eval {$hit->strand};
like($@, qr'strand not implemented for Model-based searches','Hit strand not implemented');
eval {$hit->logical_length};
like($@, qr'logical_length not implemented for Model-based searches','Hit logical_length not implemented');
eval {$hit->frac_aligned_hit};
like($@, qr'frac_aligned_hit not implemented for Model-based searches','Hit frac_aligned_hit not implemented');
eval{$hit->frac_aligned_query};
like($@, qr'frac_aligned_query not implemented for Model-based searches','Hit frac_aligned_query not implemented');
eval {$hit->frac_conserved};
like($@, qr'frac_conserved not implemented for Model-based searches','Hit frac_conserved not implemented');
eval{$hit->frac_identical};
like($@, qr'frac_identical not implemented for Model-based searches','Hit frac_identical not implemented');
eval{$hit->matches};
like($@, qr'matches not implemented for Model-based searches','Hit matches not implemented');
eval{$hit->gaps};
like($@, qr'gaps not implemented for Model-based searches','Hit gaps not implemented');
eval{$hit->frame};
like($@, qr'frame not implemented for Model-based searches','Hit frame not implemented');
eval {$hit->range};
like($@, qr'range not implemented for Model-based searches','Hit range not implemented');
eval {$hit->seq_inds};
like($@, qr'seq_inds not implemented for Model-based searches','Hit seq_inds not implemented');

# p() works but there are no evalues yet for Infernal output, so catch and check...
eval {$hit->p};
like($@, qr'P-value not defined. Using expect\(\) instead',
     "No p values");

is($hit->length, 0, "Hit length");
is($hit->overlap, 0, "Hit overlap");
is($hit->query_length, 102, "Hit query_length");
is($hit->rank, 1, "Hit rank");
is($hit->raw_score, 81.29, "Hit raw_score");
is($hit->score, 81.29, "Hit score");
is($hit->significance, undef, "Hit significance");

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'CMSEARCH', "HSP algorithm");
is($hsp->evalue, undef, "HSP evalue");
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
($meta) = $hsp->feature1->get_tag_values('meta');
is($meta, ':::::::::::::::::((((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,)))).))))::::::::::::::');
($meta) = $hsp->feature2->get_tag_values('meta');
is($meta, ':::::::::::::::::((((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,)))).))))::::::::::::::');

is($hsp->frame, 0, "HSP frame");
is($hsp->gaps, 0, "HSP gaps");
is($hit->length, 0, "Hit length");
isa_ok($hsp->get_aln, 'Bio::Align::AlignI');
isa_ok($hsp->hit, 'Bio::SeqFeature::Similarity', "HSP hit");
is($hsp->hit_string,
   'CAUGAAAUCAAAACACGACCUCAUAUAAUCUUGGGAAUAUGGCCCAUAAGUUUCUACCCGGCAACCGUAAAUUGCCGGACUAUGcAGGGAAGUGAUCGAUAAA',
   "HSP hit_string");
is($hsp->homology_string,
   ' A+ A+A+ AAAA A   :CUC:UAUAAU: :GGGAAUAUGGCCC: :AGUUUCUACC:GGCAACCGUAAAUUGCC:GACUA:G AG: AA + ++  +++++',
   "HSP homology_string");
is($hsp->hsp_group, undef, "HSP hsp_group");
is($hsp->hsp_length, 103, "HSP hsp_length");
is($hsp->length, 103, "HSP length");
is($hsp->links, undef, "HSP links");
is($hsp->n, '', "HSP n");
is($hsp->pvalue, undef, "HSP pvalue");
isa_ok($hsp->query, 'Bio::SeqFeature::Similarity', "HSP query");
is($hsp->query_string,
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGACUAcG.aGuaAauauuaaauauuu',
   "HSP query_string");
is($hsp->range, 102, "HSP range");
is($hsp->rank, 1, "HSP rank");
is($hsp->seq_inds, 67, "HSP seq_inds");
is($hsp->significance, undef, "HSP significance");
is($hsp->end, 102, "HSP end");
is($hsp->expect, undef, "HSP expect");

# These Bio::Search::HSP::HSPI methods are currently unimplemented in
# Bio::Search::HSP::ModelHSP; they may be integrated over time but will require
# some reconfiguring for Model-based searches

eval {$hsp->matches};
like($@, qr'matches not implemented for Model-based searches','HSP matches not implemented');
eval {$hsp->frac_conserved};
like($@, qr'frac_conserved not implemented for Model-based searches','HSP frac_conserved not implemented');
eval {$hsp->frac_identical};
like($@, qr'frac_identical not implemented for Model-based searches','HSP frac_identical not implemented');
eval {$hsp->num_conserved};
like($@, qr'num_conserved not implemented for Model-based searches','HSP num_conserved not implemented');
eval {$hsp->num_identical};
like($@, qr'num_identical not implemented for Model-based searches','HSP num_identical not implemented');
eval {$hsp->percent_identity};
like($@, qr'percent_identity not implemented for Model-based searches','HSP percent_identity not implemented');
eval {$hsp->cigar_string};
like($@, qr'cigar_string not implemented for Model-based searches','HSP cigar_string not implemented');
eval {$hsp->generate_cigar_string};
like($@, qr'generate_cigar_string not implemented for Model-based searches','HSP cigar_string not implemented');

isa_ok($hsp->seq, 'Bio::LocatableSeq');
is($hsp->seq_str,
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGACUAcG.aGuaAauauuaaauauuu',
   "HSP seq_str");
is($hsp->start, 1, "HSP start");
is($hsp->custom_score, undef, "HSP custom_score");
is($hsp->meta,
   ':::::::::::::::::((((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,)))).))))::::::::::::::',
   "HSP meta");
is($hsp->strand, 1, "HSP strand");

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'CMSEARCH', "HSP algorithm");
is($hsp->evalue, undef, "HSP evalue");
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
is($hsp->frame, 0, "HSP frame");
is($hsp->gaps, 0, "HSP gaps");
# infernal can return alignment data
isa_ok($hsp->get_aln, 'Bio::Align::AlignI');
isa_ok($hsp->hit, 'Bio::SeqFeature::Similarity', "HSP hit");
is($hsp->hit_string,
   'AGAAAUCAAAUAAGAUGAAUUCGUAUAAUCGCGGGAAUAUGGCUCGCAAGUCUCUACCAAGCUACCGUAAAUGGCUUGACUACGUAAACAUUUCUUUCGUUU',
   "HSP hit_string");
is($hsp->homology_string,
   'A AAAU AAA+AA A+   : CGUAUAAU::CG:GAAUAUGGC:CG::AGU UCUACCA:GC ACCGUAAAU GC:UGACUACG :   AU+U +++  UUU',
   "HSP homology_string");
is($hsp->hsp_group, undef, "HSP hsp_group");
is($hsp->hsp_length, 103, "HSP hsp_length");
is($hsp->length, 103, "HSP length");
is($hsp->links, undef, "HSP links");
is($hsp->n, '', "HSP n");
is($hsp->pvalue, undef, "HSP pvalue");
isa_ok($hsp->query, 'Bio::SeqFeature::Similarity', "HSP query");
is($hsp->query_string,
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGACUAcGaGuaAauauuaaauauuu',
   "HSP query_string");
is($hsp->range, 102, "HSP range");
is($hsp->rank, 2, "HSP rank");
is($hsp->seq_inds, 69, "HSP seq_inds");
is($hsp->significance, undef, "HSP significance");
is($hsp->end, 102, "HSP end");
is($hsp->expect, undef, "HSP expect");
#is($hsp->matches, 2, "HSP matches");
isa_ok($hsp->seq, 'Bio::LocatableSeq');
# this should probably default to the hit string
is($hsp->seq_str,
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGACUAcGaGuaAauauuaaauauuu',
   "HSP seq_str");
is($hsp->start, 1, "HSP start");
is($hsp->custom_score, undef, "HSP custom_score");
is($hsp->meta,
   ':::::::::::::::::((((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,))))))))::::::::::::::',
   "HSP meta");
is($hsp->strand, 1, "HSP strand");

# one more hit...

$hit = $result->next_hit;
isa_ok($hit, 'Bio::Search::Hit::HitI');
is($hit->accession, 'X83878.1', "Hit accession");
is($hit->ncbi_gi, '633168', "Hit GI");
is($hit->algorithm, 'CMSEARCH', "Hit algorithm");
is($hit->bits, 79.36, "Hit bits");
is($hit->description, '', "Hit description"); # no hit descs yet
is($hit->length, 0, "Hit length");
is($hit->locus, '', "Hit locus");
is($hit->n, 1, "Hit n");
is($hit->name, 'gi|633168|emb|X83878.1|', "Hit name");
is($hit->num_hsps, 1, "Hit num_hsps"); 
is($hit->overlap, 0, "Hit overlap");
is($hit->query_length, 102, "Hit query_length");
is($hit->rank, 2, "Hit rank");
is($hit->raw_score, 79.36, "Hit raw_score");
is($hit->score, 79.36, "Hit score");
is($hit->significance, undef, "Hit significance");

# one more HSP...

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'CMSEARCH', "HSP algorithm");
is($hsp->evalue, undef, "HSP evalue");
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
is($hsp->frame, 0, "HSP frame");
is($hsp->gaps, 2, "HSP gaps");
isa_ok($hsp->get_aln, 'Bio::Align::AlignI');
isa_ok($hsp->hit, 'Bio::SeqFeature::Similarity', "HSP hit");
is($hsp->hit_string,
   'UUACAAUAUAAUAGGAACACUCAUAUAAUCGCGUGGAUAUGGCACGCAAGUUUCUACCGGGCA-CCGUAAA-UGUCCGACUAUGGGUGAGCAAUGGAACCGC',
   "HSP hit_string");
is($hsp->homology_string,
   '+ A A++A AA A  AA:AC+C:UAUAAU::CG:G AUAUGGC:CG::AGUUUCUACC:G CA CCGUAAA UG C:GACUA:G+GU:A  A+U  A+    ',
   "HSP homology_string");
is($hsp->hsp_group, undef, "HSP hsp_group");
is($hsp->hsp_length, 103, "HSP hsp_length");
is($hsp->length, 103, "HSP length");
is($hsp->links, undef, "HSP links");
is($hsp->n, '', "HSP n");
isa_ok($hsp->query, 'Bio::SeqFeature::Similarity', "HSP query");
is($hsp->query_string,
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGACUAcGaGuaAauauuaaauauuu',
   "HSP query_string");
is($hsp->range, 102, "HSP range");
is($hsp->rank, 1, "HSP rank");
is($hsp->seq_inds, 64, "HSP seq_inds");
is($hsp->significance, undef, "HSP significance");
is($hsp->end, 102, "HSP end");
is($hsp->expect, undef, "HSP expect");
isa_ok($hsp->seq, 'Bio::LocatableSeq');
is($hsp->seq_str,
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGACUAcGaGuaAauauuaaauauuu',
   "HSP seq_str");
is($hsp->start, 1, "HSP start");
is($hsp->custom_score, undef, "HSP custom_score");
is($hsp->meta,
   ':::::::::::::::::((((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,))))))))::::::::::::::',
   "HSP meta");
is($hsp->strand, 1, "HSP strand");

my $symbols = {
            '5-prime'        => '(',
            '3-prime'        => ')',
            'single-strand'  => ':',
            'unknown'        => '?',
            'gap'            => '-'
             # may add more for quartets, triplets
              };

$searchio = Bio::SearchIO->new( -format => 'infernal',
                                -file   => test_input_file('test.infernal'),
                                # version is reset to the correct one by parser
                                -version => 0.7, 
                                -model => 'Purine',
                                -query_acc => 'RF00167',
                                -query_desc => 'Purine riboswitch',
                                -database => 'b_sub.fas',
                                -hsp_minscore => 40,
                                -convert_meta => 1,
                                -symbols => $symbols
                               );

$result = $searchio->next_result;
$hit = $result->next_hit;
$hsp = $hit->next_hsp;
is($hsp->meta,
   ':::::::::::::::::((((((((:::(((((((:::::::)))))))::::::::(((((((:::::::)))))))::))))-))))::::::::::::::',
   "HSP meta gap bug");
$hsp = $hit->next_hsp;
is($hsp->meta,
   ':::::::::::::::::((((((((:::(((((((:::::::)))))))::::::::(((((((:::::::)))))))::))))))))::::::::::::::',
   "HSP meta");
$hit = $result->next_hit;
$hsp = $hit->next_hsp;
is($hsp->meta,
   ':::::::::::::::::((((((((:::(((((((:::::::)))))))::::::::(((((((:::::::)))))))::))))))))::::::::::::::',
   "HSP meta");
($meta) = $hsp->feature1->get_tag_values('meta');
is($meta, ':::::::::::::::::((((((((:::(((((((:::::::)))))))::::::::(((((((:::::::)))))))::))))))))::::::::::::::');
($meta) = $hsp->feature2->get_tag_values('meta');
is($meta, ':::::::::::::::::((((((((:::(((((((:::::::)))))))::::::::(((((((:::::::)))))))::))))))))::::::::::::::');

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
is($hit->raw_score, '', "Hit raw_score");
is($hit->score, '', "Hit score");
is($hit->significance, undef, "Hit significance");

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'RNAMOTIF', "HSP algorithm");
is($hsp->evalue, undef, "HSP evalue");
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
is($hsp->frame, 0, "HSP frame");
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
is($hsp->n, '',"HSP n");
isa_ok($hsp->query, 'Bio::SeqFeature::Similarity', "HSP query");
is($hsp->query_string,
   '',
   "HSP query_string");
is($hsp->range, 75, "HSP range");
is($hsp->rank, 1, "HSP rank");
is($hsp->significance, undef, "HSP significance");
is($hsp->end, 75, "HSP end");
is($hsp->expect, undef, "HSP expect");
isa_ok($hsp->seq, 'Bio::LocatableSeq');
is($hsp->seq_str,
   '',
   "HSP seq_str");
is($hsp->start, 1, "HSP start");
is($hsp->custom_score, undef, "HSP custom_score");
is($hsp->meta,
   '<<<<<<<..<<<.........>>>.<<<<<.......>>>>>.......<<<<.......>>>>>>>>>>>....',
   "HSP meta");
is($hsp->strand, 1, "HSP strand");
($meta) = $hsp->feature1->get_tag_values('meta');
is($meta, '<<<<<<<..<<<.........>>>.<<<<<.......>>>>>.......<<<<.......>>>>>>>>>>>....');
($meta) = $hsp->feature2->get_tag_values('meta');
is($meta, '<<<<<<<..<<<.........>>>.<<<<<.......>>>>>.......<<<<.......>>>>>>>>>>>....');

#### ERPIN ####

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
is($hit->significance, '4.44e-06', "Hit significance");

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'ERPIN', "HSP algorithm");
is($hsp->evalue, '1.68e-05', "HSP evalue");
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
is($hsp->frame, 0, "HSP frame");
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
is($hsp->significance, '1.68e-05', "HSP significance");
is($hsp->expect, '1.68e-05', "HSP expect");
isa_ok($hsp->seq, 'Bio::LocatableSeq');
is($hsp->seq_str,
   '',
   "HSP seq_str");
is($hsp->start, 1, "HSP start");
is($hsp->custom_score, undef, "HSP custom_score");
TODO: {
    local $TODO = 'Working on meta string building';
    isnt($hsp->meta, undef, "HSP meta");
    ($meta) = $hsp->feature1->has_tag('meta') ? $hsp->feature1->get_tag_values('meta') : undef;
    isnt($meta, undef);
    ($meta) = $hsp->feature2->has_tag('meta') ? $hsp->feature2->get_tag_values('meta') : undef;
    isnt($meta, undef);
}
is($hsp->strand, 1, "HSP strand");
($meta) = $hsp->feature1->has_tag('meta') ? $hsp->feature1->get_tag_values('meta') : undef;
is($meta, undef);
($meta) = $hsp->feature2->has_tag('meta') ? $hsp->feature2->get_tag_values('meta') : undef;
is($meta, undef);

# ERPIN lacks sequence for query, will spit back a warning..
eval{$hsp->get_aln};
like($@, qr'Missing query string, can\'t build alignment','ERPIN get_aln warning');

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'ERPIN', "HSP algorithm");
is($hsp->evalue, '5.61e-05', "HSP evalue");
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
is($hsp->frame, 0, "HSP frame");
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
is($hsp->significance, '5.61e-05', "HSP significance");
is($hsp->end, 37, "HSP end");
is($hsp->expect, '5.61e-05', "HSP expect");
isa_ok($hsp->seq, 'Bio::LocatableSeq');
is($hsp->seq_str,   '',   "HSP seq_str");
is($hsp->start, 1, "HSP start");
is($hsp->custom_score, undef, "HSP custom_score");
TODO: {
    local $TODO = 'Working on meta string building';
    isnt($hsp->meta, undef, "HSP meta");
    ($meta) = $hsp->feature1->has_tag('meta') ? $hsp->feature1->get_tag_values('meta') : undef;
    isnt($meta, undef);
    ($meta) = $hsp->feature2->has_tag('meta') ? $hsp->feature2->get_tag_values('meta') : undef;
    isnt($meta, undef);
}
is($hsp->strand, -1, "HSP strand");
