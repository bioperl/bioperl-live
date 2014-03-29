# -*-Perl-*- Test Harness script for Bioperl
# $Id: SearchIO_infernal.t 14672 2008-04-22 21:42:50Z cjfields $

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 412);
    
    use_ok('Bio::SearchIO');
}

my ($result, $iter, $hit, $hsp, $algorithm, $meta);

### Infernal v. 1.0 ####

my $searchio = Bio::SearchIO->new( -format => 'infernal',
                                -file   => test_input_file('test2.infernal'),
                                -model => 'tRNAtest',
                                -query_acc => 'RF01234',
                                -query_desc => 'tRNA',
                               );

$result = $searchio->next_result;
isa_ok($result, 'Bio::Search::Result::ResultI');
is($result->algorithm, 'CMSEARCH', "Result");
is($result->algorithm_reference, undef, "Result reference");
is($result->algorithm_version, '1.0', "Result version");
is($result->available_parameters, 0, "Result parameters");
is($result->available_statistics, 0, "Result statistics");
is($result->database_entries, '', "Result entries");
is($result->database_letters, 600000, "Result letters");
is($result->database_name, 'tosearch.300Kb.db',
   "Result database_name");
is($result->num_hits, 1, "Result num_hits");
is($result->program_reference, undef, "Result program_reference");
is($result->query_accession, 'RF01234', "Result query_accession");
is($result->query_description, 'tRNA', "Result query_description");
is($result->query_length, 72, "Result query_length");
is($result->query_name, 'trna.5-1', "Result query_name");

$hit = $result->next_hit;

isa_ok($hit, 'Bio::Search::Hit::HitI');
is($hit->ncbi_gi, '', "Hit GI");
is($hit->accession, 'example', "Hit accession");
is($hit->algorithm, 'CMSEARCH', "Hit algorithm");
is($hit->bits, '78.06', "Hit bits");
is($hit->description, '', "Hit description"); # no hit descs yet
is($hit->locus, '', "Hit locus");
is($hit->n, 3, "Hit n");
is($hit->name, 'example', "Hit name");
is($hit->num_hsps, 3, "Hit num_hsps");

# These Bio::Search::Hit::HitI methods are currently unimplemented in
# Bio::Search::Hit::ModelHit; they may be integrated over time but will require
# some reconfiguring for Model-based searches

# these need to be replaced by dies_ok() or warnings_like()
warning_like { $hit->length_aln() }
    qr'length_aln not implemented for Model-based searches',
    "Hit length_aln() not implemented";
warning_like {$hit->num_unaligned_hit}
    qr'num_unaligned_hit/num_unaligned_sbjct not implemented for Model-based searches',
    "Hit num_unaligned_hit() not implemented";
warning_like {$hit->num_unaligned_query}
    qr'num_unaligned_query not implemented for Model-based searches',
    "Hit num_unaligned_query() not implemented";
warning_like {$hit->num_unaligned_sbjct}
    qr'num_unaligned_hit/num_unaligned_sbjct not implemented for Model-based searches',
    "Hit num_unaligned_sbjct() not implemented";
warning_like {$hit->start}
    qr'start not implemented for Model-based searches',
    'Hit start not implemented';
warning_like {$hit->end}
    qr'end not implemented for Model-based searches',
    'Hit end not implemented';
warning_like {$hit->strand}
    qr'strand not implemented for Model-based searches',
    'Hit strand not implemented';
warning_like {$hit->logical_length}
    qr'logical_length not implemented for Model-based searches',
    'Hit logical_length not implemented';
warning_like {$hit->frac_aligned_hit}
    qr'frac_aligned_hit not implemented for Model-based searches',
    'Hit frac_aligned_hit not implemented';
warning_like {$hit->frac_aligned_query}
    qr'frac_aligned_query not implemented for Model-based searches',
    'Hit frac_aligned_query not implemented';
warning_like {$hit->frac_conserved}
    qr'frac_conserved not implemented for Model-based searches',
    'Hit frac_conserved not implemented';
warning_like {$hit->frac_identical}
    qr'frac_identical not implemented for Model-based searches',
    'Hit frac_identical not implemented';
warning_like {$hit->matches}
    qr'matches not implemented for Model-based searches',
    'Hit matches not implemented';
warning_like {$hit->gaps}
    qr'gaps not implemented for Model-based searches',
    'Hit gaps not implemented';
warning_like {$hit->frame}
    qr'frame not implemented for Model-based searches',
    'Hit frame not implemented';
warning_like {$hit->range}
    qr'range not implemented for Model-based searches',
    'Hit range not implemented';
warning_like {$hit->seq_inds} 
    qr'seq_inds not implemented for Model-based searches',
    'Hit seq_inds not implemented';

is($hit->length, 0, "Hit length");
is($hit->overlap, 0, "Hit overlap");
is($hit->query_length, 72, "Hit query_length");
is($hit->rank, 1, "Hit rank");
is($hit->raw_score, '78.06', "Hit raw_score");
is($hit->score, '78.06', "Hit score");
float_is($hit->p, '2.906e-26', "Hit p");
float_is($hit->significance, '3.133e-21');

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'CMSEARCH', "HSP algorithm");
float_is($hsp->evalue, '3.133e-21');
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
($meta) = $hsp->feature1->get_tag_values('meta');
is($meta, '(((((((,,<<<<___.____>>>>,<<<<<_______>>>>>,,,,,<<<<<_______>>>>>))))))):');
($meta) = $hsp->feature2->get_tag_values('meta');
is($meta, '(((((((,,<<<<___.____>>>>,<<<<<_______>>>>>,,,,,<<<<<_______>>>>>))))))):');

is($hsp->frame('query'), 0, "HSP frame");
is($hsp->gaps, 1, "HSP gaps");
is($hit->length, 0, "Hit length");
isa_ok($hsp->get_aln, 'Bio::Align::AlignI');
isa_ok($hsp->hit, 'Bio::SeqFeature::Similarity', "HSP hit");
is($hsp->hit_string,
   'GCGGAUUUAGCUCAGUuGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA',
   "HSP hit_string");
is($hsp->homology_string,
   'GC::A::UAGC:CAGU GG AG:GCGCCAG:CUG+++A:CUGGAGGUCC:G:GUUCGAU C:C:G::U::GCA',
   "HSP homology_string");
is($hsp->hsp_group, undef, "HSP hsp_group");
is($hsp->hsp_length, 73, "HSP hsp_length");
is($hsp->length, 73, "HSP length");
is($hsp->links, undef, "HSP links");
is($hsp->n, 1, "HSP n");
float_is($hsp->pvalue, 2.906e-26, "HSP pvalue");
isa_ok($hsp->query, 'Bio::SeqFeature::Similarity', "HSP query");
is($hsp->query_string,
   'gCcgacAUaGcgcAgU.GGuAgcgCgccagccUgucAagcuggAGgUCCgggGUUCGAUuCcccGUgucgGca',
   "HSP query_string");
is($hsp->range, 72, "HSP range");
is($hsp->rank, 1, "HSP rank");
float_is($hsp->significance, 3.133e-21);
is($hsp->end, 72, "HSP end");
float_is($hsp->expect, '3.133e-21', "HSP expect");

# These Bio::Search::HSP::HSPI methods are currently unimplemented in
# Bio::Search::HSP::ModelHSP; they may be integrated over time but will require
# some reconfiguring for Model-based searches

warning_like {$hsp->seq_inds}
    qr'seq_inds not implemented for Model-based searches',
    'HSP seq_inds not implemented';
warning_like {$hsp->matches}
    qr'matches not implemented for Model-based searches',
    'HSP matches not implemented';
warning_like {$hsp->frac_conserved}
    qr'frac_conserved not implemented for Model-based searches',
    'HSP frac_conserved not implemented';
warning_like {$hsp->frac_identical}
    qr'frac_identical not implemented for Model-based searches',
    'HSP frac_identical not implemented';
warning_like {$hsp->num_conserved}
    qr'num_conserved not implemented for Model-based searches',
    'HSP num_conserved not implemented';
warning_like {$hsp->num_identical}
    qr'num_identical not implemented for Model-based searches',
    'HSP num_identical not implemented';
warning_like {$hsp->percent_identity}
    qr'percent_identity not implemented for Model-based searches',
    'HSP percent_identity not implemented';
warning_like {$hsp->cigar_string}
    qr'cigar_string not implemented for Model-based searches',
    'HSP cigar_string not implemented';
warning_like {$hsp->generate_cigar_string}
    qr'generate_cigar_string not implemented for Model-based searches',
    'HSP cigar_string not implemented';

isa_ok($hsp->seq, 'Bio::LocatableSeq');
is($hsp->seq_str,
   'gCcgacAUaGcgcAgU.GGuAgcgCgccagccUgucAagcuggAGgUCCgggGUUCGAUuCcccGUgucgGca',
   "HSP seq_str");
is($hsp->start, 1, "HSP start");
is($hsp->custom_score, undef, "HSP custom_score");
is($hsp->meta,
   '(((((((,,<<<<___.____>>>>,<<<<<_______>>>>>,,,,,<<<<<_______>>>>>))))))):',
   "HSP meta");
is($hsp->strand('hit'), 1, "HSP strand");

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'CMSEARCH', "HSP algorithm");
float_is($hsp->evalue, 0.6752);
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
is($hsp->frame('query'), 0, "HSP frame");
is($hsp->gaps, 4, "HSP gaps");
# infernal can return alignment data
isa_ok($hsp->get_aln, 'Bio::Align::AlignI');
isa_ok($hsp->hit, 'Bio::SeqFeature::Similarity', "HSP hit");
is($hsp->hit_string,
   'UCUGCUAUGGCGUAAUGGCCACGCGC----CCAUCAACAAAGAUAUC*[19]*UAACAGGA',
   "HSP hit_string");
is($hsp->homology_string,
   ' C:G :AU+GCG:A+UGG  :CGCGC    C  UCAA +++GA +UC      U: C:G A',
   "HSP homology_string");
is($hsp->hsp_group, undef, "HSP hsp_group");
is($hsp->hsp_length, 73, "HSP hsp_length");
is($hsp->length, 73, "HSP length");
is($hsp->links, undef, "HSP links");
is($hsp->n, 1, "HSP n");
float_is($hsp->pvalue, 6.263e-06, "HSP pvalue");
isa_ok($hsp->query, 'Bio::SeqFeature::Similarity', "HSP query");
is($hsp->query_string,
   'gCcgacAUaGcgcAgUGGuAgcgCgccagccUgucAagcuggAGgUC*[17]*UgucgGca',
   "HSP query_string");
is($hsp->range, 72, "HSP range");
is($hsp->rank, 2, "HSP rank");
float_is($hsp->significance, 0.6752);
is($hsp->end, 72, "HSP end");
float_is($hsp->expect, 0.6752, "HSP expect");
isa_ok($hsp->seq, 'Bio::LocatableSeq');
# this should probably default to the hit string
is($hsp->seq_str,
   'gCcgacAUaGcgcAgUGGuAgcgCgccagccUgucAagcuggAGgUC*[17]*UgucgGca',
   "HSP seq_str");
is($hsp->start, 1, "HSP start");
is($hsp->custom_score, undef, "HSP custom_score");
is($hsp->meta,
   '(((((((,,<<<<_______>>>>,<<<<<_______>>>>>,,,,,~~~~~~))))))):',
   "HSP meta");
is($hsp->strand('hit'), 1, "HSP strand");

### Infernal pre-v. 1.0 ####

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

# p() works but there are no evalues yet for Infernal output, so catch and check...
warning_like {$hit->p}
    qr'P-value not defined. Using significance\(\) instead',
    "No p values";

is($hit->length, 0, "Hit length");
is($hit->overlap, 0, "Hit overlap");
is($hit->query_length, 102, "Hit query_length");
is($hit->rank, 1, "Hit rank");
is($hit->raw_score, 81.29, "Hit raw_score");
is($hit->score, 81.29, "Hit score");
float_is($hit->significance, undef);

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'CMSEARCH', "HSP algorithm");
float_is($hsp->evalue, undef);
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
($meta) = $hsp->feature1->get_tag_values('meta');
is($meta, ':::::::::::::::::((((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,)))).))))::::::::::::::');
($meta) = $hsp->feature2->get_tag_values('meta');
is($meta, ':::::::::::::::::((((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,)))).))))::::::::::::::');

is($hsp->frame('query'), 0, "HSP frame");
is($hsp->gaps, 1, "HSP gaps");
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
is($hsp->n, 1, "HSP n");
float_is($hsp->pvalue, undef, "HSP pvalue");
isa_ok($hsp->query, 'Bio::SeqFeature::Similarity', "HSP query");
is($hsp->query_string,
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGACUAcG.aGuaAauauuaaauauuu',
   "HSP query_string");
is($hsp->range, 102, "HSP range");
is($hsp->rank, 1, "HSP rank");
float_is($hsp->significance, undef);
is($hsp->end, 102, "HSP end");
float_is($hsp->expect, undef, "HSP expect");

isa_ok($hsp->seq, 'Bio::LocatableSeq');
is($hsp->seq_str,
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGACUAcG.aGuaAauauuaaauauuu',
   "HSP seq_str");
is($hsp->start, 1, "HSP start");
is($hsp->custom_score, undef, "HSP custom_score");
is($hsp->meta,
   ':::::::::::::::::((((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,)))).))))::::::::::::::',
   "HSP meta");
is($hsp->strand('hit'), 1, "HSP strand");

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'CMSEARCH', "HSP algorithm");
float_is($hsp->evalue, undef);
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
is($hsp->frame('query'), 0, "HSP frame");
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
is($hsp->n, 1, "HSP n");
float_is($hsp->pvalue, undef, "HSP pvalue");
isa_ok($hsp->query, 'Bio::SeqFeature::Similarity', "HSP query");
is($hsp->query_string,
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGACUAcGaGuaAauauuaaauauuu',
   "HSP query_string");
is($hsp->range, 102, "HSP range");
is($hsp->rank, 2, "HSP rank");
float_is($hsp->significance, undef);
is($hsp->end, 102, "HSP end");
float_is($hsp->expect, undef, "HSP expect");
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
is($hsp->strand('hit'), 1, "HSP strand");

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
float_is($hit->significance, undef);

# one more HSP...

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'CMSEARCH', "HSP algorithm");
float_is($hsp->evalue, undef);
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
is($hsp->frame('query'), 0, "HSP frame");
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
is($hsp->n, 1, "HSP n");
isa_ok($hsp->query, 'Bio::SeqFeature::Similarity', "HSP query");
is($hsp->query_string,
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGACUAcGaGuaAauauuaaauauuu',
   "HSP query_string");
is($hsp->range, 102, "HSP range");
is($hsp->rank, 1, "HSP rank");
float_is($hsp->significance, undef);
is($hsp->end, 102, "HSP end");
float_is($hsp->expect, undef, "HSP expect");
isa_ok($hsp->seq, 'Bio::LocatableSeq');
is($hsp->seq_str,
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGACUAcGaGuaAauauuaaauauuu',
   "HSP seq_str");
is($hsp->start, 1, "HSP start");
is($hsp->custom_score, undef, "HSP custom_score");
is($hsp->meta,
   ':::::::::::::::::((((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,))))))))::::::::::::::',
   "HSP meta");
is($hsp->strand('hit'), 1, "HSP strand");

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
                                -symbols => $symbols,
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

## Infernal 0.81 parsing ##

$searchio = Bio::SearchIO->new( -format => 'infernal',
                                -file   => test_input_file('purine_v081.infernal'),
                                # version is reset to the correct one by parser
                                -query_acc => 'RF00167',
                                -query_desc => 'Purine riboswitch',
                                -database => 'b_sub.fas',
                                -convert_meta => 0,
                               );

$result = $searchio->next_result;

isa_ok($result, 'Bio::Search::Result::ResultI');
$algorithm = $result->algorithm;
is($result->algorithm, 'CMSEARCH', "Result $algorithm");
is($result->algorithm_reference, undef, "Result $algorithm reference");
is($result->algorithm_version, 0.81, "Result $algorithm version");
is($result->available_parameters, 0, "Result parameters");
is($result->available_statistics, 0, "Result statistics");
is($result->database_entries, '', "Result entries");
is($result->database_letters, '', "Result letters");
is($result->database_name, 'b_sub.fas', "Result database_name");
is($result->num_hits, 3, "Result num_hits");
is($result->program_reference, undef, "Result program_reference");
is($result->query_accession, 'RF00167', "Result query_accession");
is($result->query_description, 'Purine riboswitch', "Result query_description");
is($result->query_length, 102, "Result query_length");
is($result->query_name, 'Purine', "Result query_name");

$hit = $result->next_hit;
isa_ok($hit, 'Bio::Search::Hit::HitI');
is($hit->ncbi_gi, '633168', "Hit GI");
is($hit->accession, 'X83878.1', "Hit accession");
is($hit->algorithm, 'CMSEARCH', "Hit algorithm");
is($hit->bits, 79.36, "Hit bits");
is($hit->description, '', "Hit description"); # no hit descs yet
is($hit->locus, '', "Hit locus");
is($hit->n, 2, "Hit n");
is($hit->name, 'gi|633168|emb|X83878.1|', "Hit name");
is($hit->num_hsps, 2, "Hit num_hsps");

# p() works but there are no evalues yet for Infernal output, so catch and check...
warnings_like {$hit->p} qr'P-value not defined. Using significance\(\) instead',
     "No p values";

is($hit->length, 0, "Hit length");
is($hit->overlap, 0, "Hit overlap");
is($hit->query_length, 102, "Hit query_length");
is($hit->rank, 1, "Hit rank");
is($hit->raw_score, 79.36, "Hit raw_score");
is($hit->score, 79.36, "Hit score");
float_is($hit->significance, 1.945e-07);

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'CMSEARCH', "HSP algorithm");
float_is($hsp->evalue, 1.945e-07);
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
($meta) = $hsp->feature1->get_tag_values('meta');
is($meta, ':::::::::::::::::((((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,))))))))::::::::::::::');
($meta) = $hsp->feature2->get_tag_values('meta');
is($meta, ':::::::::::::::::((((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,))))))))::::::::::::::');

is($hsp->frame('query'), 0, "HSP frame");
is($hsp->gaps, 2, "HSP gaps");
is($hit->length, 0, "Hit length");
isa_ok($hsp->get_aln, 'Bio::Align::AlignI');
isa_ok($hsp->hit, 'Bio::SeqFeature::Similarity', "HSP hit");
is($hsp->hit_string,
   'UUACAAUAUAAUAGGAACACUCAUAUAAUCGCGUGGAUAUGGCACGCAAGUUUCUACCGGGCA-CCGUAAA-UGUCCGACUAUGGGUGAGCAAUGGAACCGC',
   "HSP hit_string");
is($hsp->homology_string,
   '+ A A++A AA A  AA:AC+C:UAUAAU::CG:G AUAUGGC:CG::AGUUUCUACC:G CA CCGUAAA UG C:GACUA:G+GU:A  A+U  A+    ',
   "HSP homology_string");
is($hsp->hsp_group, undef, "HSP hsp_group");
is($hsp->hsp_length,102, "HSP hsp_length");
is($hsp->length, 102, "HSP length");
is($hsp->links, undef, "HSP links");
is($hsp->n, 1, "HSP n");
float_is($hsp->pvalue, 1.945e-07, "HSP pvalue");
isa_ok($hsp->query, 'Bio::SeqFeature::Similarity', "HSP query");
is($hsp->query_string,
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGACUAcGaGuaAauauuaaauauuu',
   "HSP query_string");
is($hsp->range, 102, "HSP range");
is($hsp->rank, 1, "HSP rank");
float_is($hsp->significance, 1.945e-07);
is($hsp->end, 102, "HSP end");
float_is($hsp->expect, 1.945e-07, "HSP expect");

isa_ok($hsp->seq, 'Bio::LocatableSeq');
is($hsp->seq_str,
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGACUAcGaGuaAauauuaaauauuu',
   "HSP seq_str");
is($hsp->start, 1, "HSP start");
is($hsp->custom_score, undef, "HSP custom_score");
is($hsp->meta,
   ':::::::::::::::::((((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,))))))))::::::::::::::',
   "HSP meta");
is($hsp->strand('hit'), 1, "HSP strand");

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'CMSEARCH', "HSP algorithm");
float_is($hsp->evalue, 6.802);
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
is($hsp->frame('query'), 0, "HSP frame");
is($hsp->gaps, 4, "HSP gaps");
# infernal can return alignment data
isa_ok($hsp->get_aln, 'Bio::Align::AlignI');
isa_ok($hsp->hit, 'Bio::SeqFeature::Similarity', "HSP hit");
is($hsp->hit_string,
   'CGUGCGGUUCCAUUGCUCACCCAUA-GUCGGACAU-UUACGG-UGCCCGGUAGAAACUUGCGUGCCAUAUCCACGCGAUUaUAUGAGUGUUCCUAUUAUAUUG',
   "HSP hit_string");
is($hsp->homology_string,
   '  +    +   A    +:AC C:UA  +::: ::   UA GG :: :::GU    AC: G::::CC UA  ::::C :   UA:G GU: +  U+++AUAUU ',
   "HSP homology_string");
is($hsp->hsp_group, undef, "HSP hsp_group");
is($hsp->hsp_length, 102, "HSP hsp_length");
is($hsp->length, 102, "HSP length");
is($hsp->links, undef, "HSP links");
is($hsp->n, 1, "HSP n");
float_is($hsp->pvalue, 0.9989, "HSP pvalue");
isa_ok($hsp->query, 'Bio::SeqFeature::Similarity', "HSP query");
is($hsp->query_string,
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGAC.UAcGaGuaAauauuaaauauuu',
   "HSP query_string");
is($hsp->range, 102, "HSP range");
is($hsp->rank, 2, "HSP rank");
float_is($hsp->significance, 6.802);
is($hsp->end, 102, "HSP end");
float_is($hsp->expect, 6.802, "HSP expect");
#is($hsp->matches, 2, "HSP matches");
isa_ok($hsp->seq, 'Bio::LocatableSeq');
# this should probably default to the hit string
is($hsp->seq_str,
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGAC.UAcGaGuaAauauuaaauauuu',
   "HSP seq_str");
is($hsp->start, 1, "HSP start");
is($hsp->custom_score, undef, "HSP custom_score");
is($hsp->meta,
   ':::::::::::::::::((((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,.))))))))::::::::::::::',
   "HSP meta");
is($hsp->strand('hit'), -1, "HSP strand");

# one more hit...

$hit = $result->next_hit;
isa_ok($hit, 'Bio::Search::Hit::HitI');
is($hit->accession, 'U51115.1', "Hit accession");
is($hit->ncbi_gi, '2239287', "Hit GI");
is($hit->algorithm, 'CMSEARCH', "Hit algorithm");
is($hit->bits, 81.29, "Hit bits");
is($hit->description, '', "Hit description"); # no hit descs yet
is($hit->length, 0, "Hit length");
is($hit->locus, 'BSU51115', "Hit locus");
is($hit->n, 11, "Hit n");
is($hit->name, 'gi|2239287|gb|U51115.1|BSU51115', "Hit name");
is($hit->num_hsps, 11, "Hit num_hsps"); 
is($hit->overlap, 0, "Hit overlap");
is($hit->query_length, 102, "Hit query_length");
is($hit->rank, 2, "Hit rank");
is($hit->raw_score, 81.29, "Hit raw_score");
is($hit->score, 81.29, "Hit score");
float_is($hit->significance, 1.259e-07);

# one more HSP...

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'CMSEARCH', "HSP algorithm");
float_is($hsp->evalue, 1.259e-07);
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
is($hsp->frame('query'), 0, "HSP frame");
is($hsp->gaps, 0, "HSP gaps");
isa_ok($hsp->get_aln, 'Bio::Align::AlignI');
isa_ok($hsp->hit, 'Bio::SeqFeature::Similarity', "HSP hit");
is($hsp->hit_string,
   'AGAAAUCAAAUAAGAUGAAUUCGUAUAAUCGCGGGAAUAUGGCUCGCAAGUCUCUACCAAGCUACCGUAAAUGGCUUGACUACGUAAACAUUUCUUUCGUUU',
   "HSP hit_string");
is($hsp->homology_string,
   'A AAAU AAA+AA A+   : CGUAUAAU::CG:GAAUAUGGC:CG::AGU UCUACCA:GC ACCGUAAAU GC:UGACUACG :   AU+U +++  UUU',
   "HSP homology_string");
is($hsp->hsp_group, undef, "HSP hsp_group");
is($hsp->hsp_length, 102, "HSP hsp_length");
is($hsp->length, 102, "HSP length");
is($hsp->links, undef, "HSP links");
is($hsp->n, 1, "HSP n");
isa_ok($hsp->query, 'Bio::SeqFeature::Similarity', "HSP query");
is($hsp->query_string,
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGACUAcGaGuaAauauuaaauauuu',
   "HSP query_string");
is($hsp->range, 102, "HSP range");
is($hsp->rank, 1, "HSP rank");
float_is($hsp->significance, 1.259e-07);
is($hsp->end, 102, "HSP end");
float_is($hsp->expect, 1.259e-07, "HSP expect");
isa_ok($hsp->seq, 'Bio::LocatableSeq');
is($hsp->seq_str,
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGACUAcGaGuaAauauuaaauauuu',
   "HSP seq_str");
is($hsp->start, 1, "HSP start");
is($hsp->custom_score, undef, "HSP custom_score");
is($hsp->meta,
   ':::::::::::::::::((((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,))))))))::::::::::::::',
   "HSP meta");
is($hsp->strand('hit'), 1, "HSP strand");
