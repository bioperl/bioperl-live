# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error;
use strict;
use Dumpvalue();
my $dumper = Dumpvalue->new();

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test::More; };
    if( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
    
    use vars qw($NTESTS);
    $NTESTS = 311;
    $error = 0;

    plan tests => $NTESTS; 
}

use_ok('Bio::SearchIO');
use_ok('Bio::Root::IO');

my ($searchio, $result,$iter,$hit,$hsp, $algorithm);

#### Infernal ####

$searchio = Bio::SearchIO->new( -format => 'infernal',
                                -file   => Bio::Root::IO->catfile
                                          ('t','data','test.infernal'),
                                # version is reset to the correct one by parser
                                -version => 0.7, 
                                -model => 'Purine',
                                -query_acc => 'RF00167',
                                -query_desc => 'Purine riboswitch',
                                -database => 'b_sub.fas',
                                -hsp_minscore => 40,
                                -convert_meta => 0,
                                -verbose => -1
                               );

$result = $searchio->next_result;
isa_ok($result, 'Bio::Search::Result::ResultI');
#$dumper->dumpValue($result);
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
#is($hit->end, '102', "Hit end");
#is($hit->frac_aligned_hit, '', "Hit frac_aligned_hit");
#is($hit->frac_aligned_query, '', "Hit frac_aligned_query");
#is($hit->frac_conserved, '', "Hit frac_conserved");
#is($hit->frac_identical, '', "Hit frac_identical");
#is($hit->frame, 0, "Hit frame");
#is($hit->gaps, 0, "Hit gaps");
#is($hit->length, '', "Hit length");
is($hit->length_aln, 102, "Hi length_aln");
is($hit->locus, 'BSU51115', "Hit locus");
is($hit->logical_length, 102, "Hit logical_length");
is($hit->matches, 2, "Hit matches");
is($hit->n, 2, "Hit n");
is($hit->name, 'gi|2239287|gb|U51115.1|BSU51115', "Hit name");
is($hit->num_hsps, 2, "Hit num_hsps"); 
#is($hit->num_unaligned_hit, '', "Hit num_unaligned_hit");
#is($hit->num_unaligned_query, '', "Hit num_unaligned_query");
#is($hit->num_unaligned_sbjct, '', "Hit num_unaligned_subject");
is($hit->overlap, 0, "Hit overlap");
#is($hit->p, '', "Hit p");
is($hit->query_length, 102, "Hit query_length");
is($hit->range, 102, "Hit range");
is($hit->rank, 1, "Hit rank");
is($hit->raw_score, 81.29, "Hit raw_score");
is($hit->score, 81.29, "Hit score");
is($hit->seq_inds, 84, "Hit seq_inds");
is($hit->significance, undef, "Hit significance");
is($hit->start, 1, "Hit start");
is($hit->strand, 0, "Hit strand");

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'CMSEARCH', "HSP algorithm");
is($hsp->cigar_string, '103M', "HSP cigar_string");
is($hsp->evalue, undef, "HSP evalue");
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
#is($hsp->frac_conserved, '', "HSP frac_conserved");
#is($hsp->frac_identical, '', "HSP frac_identical");
is($hsp->frame, 0, "HSP frame");
is($hsp->gaps, 0, "HSP gaps");
is($hsp->generate_cigar_string, '0M', "HSP generate_cigar_string");
#is($hsp->get_aln, '', "HSP aln");
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
#is($hsp->n, 0, "HSP n");
#is($hsp->num_conserved, '', "HSP num_conserved");
#is($hsp->num_identical, '', "HSP num_identical");
#is($hsp->percent_identity, '', "HSP percent_identity");
#is($hsp->pvalue, '', "HSP pvalue");
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
is($hsp->matches, 2, "HSP matches");
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
is($hsp->cigar_string, '102M', "HSP cigar_string");
is($hsp->evalue, undef, "HSP evalue");
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
#is($hsp->frac_conserved, '', "HSP frac_conserved");
#is($hsp->frac_identical, '', "HSP frac_identical");
is($hsp->frame, 0, "HSP frame");
is($hsp->gaps, 0, "HSP gaps");
is($hsp->generate_cigar_string, '0M', "HSP generate_cigar_string");
#is($hsp->get_aln, '', "HSP aln");
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
#is($hsp->n, 0, "HSP n");
#is($hsp->num_conserved, '', "HSP num_conserved");
#is($hsp->num_identical, '', "HSP num_identical");
#is($hsp->percent_identity, '', "HSP percent_identity");
#is($hsp->pvalue, '', "HSP pvalue");
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
is($hsp->matches, 2, "HSP matches");
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

# one more hit...

$hit = $result->next_hit;
isa_ok($hit, 'Bio::Search::Hit::HitI');
is($hit->accession, 'X83878.1', "Hit accession");
is($hit->ncbi_gi, '633168', "Hit GI");
is($hit->algorithm, 'CMSEARCH', "Hit algorithm");
is($hit->bits, 79.36, "Hit bits");
is($hit->description, '', "Hit description"); # no hit descs yet
is($hit->end, '102', "Hit end");
#is($hit->frac_aligned_hit, '', "Hit frac_aligned_hit");
#is($hit->frac_aligned_query, '', "Hit frac_aligned_query");
#is($hit->frac_conserved, '', "Hit frac_conserved");
#is($hit->frac_identical, '', "Hit frac_identical");
is($hit->frame, 0, "Hit frame");
is($hit->gaps, 2, "Hit gaps");
#is($hit->length, '', "Hit length");
is($hit->length_aln, 102, "Hi length_aln");
is($hit->locus, '', "Hit locus");
is($hit->logical_length, 102, "Hit logical_length");
is($hit->matches, 2, "Hit matches");
is($hit->n, 1, "Hit n");
is($hit->name, 'gi|633168|emb|X83878.1|', "Hit name");
is($hit->num_hsps, 1, "Hit num_hsps"); 
#is($hit->num_unaligned_hit, '', "Hit num_unaligned_hit");
#is($hit->num_unaligned_query, '', "Hit num_unaligned_query");
#is($hit->num_unaligned_sbjct, '', "Hit num_unaligned_subject");
is($hit->overlap, 0, "Hit overlap");
#is($hit->p, '', "Hit p");
is($hit->query_length, 102, "Hit query_length");
is($hit->range, 102, "Hit range");
is($hit->rank, 2, "Hit rank");
is($hit->raw_score, 79.36, "Hit raw_score");
is($hit->score, 79.36, "Hit score");
is($hit->seq_inds, 64, "Hit seq_inds");
is($hit->significance, undef, "Hit significance");
is($hit->start, 1, "Hit start");
is($hit->strand, 0, "Hit strand");

# one more HSP...

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'CMSEARCH', "HSP algorithm");
is($hsp->cigar_string, '63MI7MI30M', "HSP cigar_string");
is($hsp->evalue, undef, "HSP evalue");
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
#is($hsp->frac_conserved, '', "HSP frac_conserved");
#is($hsp->frac_identical, '', "HSP frac_identical");
is($hsp->frame, 0, "HSP frame");
is($hsp->gaps, 2, "HSP gaps");
is($hsp->generate_cigar_string, '0M', "HSP generate_cigar_string");
#is($hsp->get_aln, '', "HSP aln");
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
#is($hsp->n, 0, "HSP n");
#is($hsp->num_conserved, '', "HSP num_conserved");
#is($hsp->num_identical, '', "HSP num_identical");
#is($hsp->percent_identity, '', "HSP percent_identity");
#is($hsp->pvalue, '', "HSP pvalue");
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
is($hsp->matches, 2, "HSP matches");
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
                                -file   => Bio::Root::IO->catfile
                                          ('t','data','test.infernal'),
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

#### ERPIN ####

$searchio = Bio::SearchIO->new( -format => 'erpin',
                                -file   => Bio::Root::IO->catfile
                                          ('t','data','testfile.erpin'),
                                -model => 'stem-loop',
                                -query_acc => 'test',
                                -version => 5.5,
                                -verbose => 2
                               );
$result = $searchio->next_result;
isa_ok($result, 'Bio::Search::Result::ResultI');
#$dumper->dumpValue($result);
$algorithm = $result->algorithm;
is($result->algorithm, 'ERPIN', "Result $algorithm");
is($result->algorithm_reference, undef, " Result $algorithm reference");
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
#is($result->query_length, 102, "Result query_length");
is($result->query_name, '/home/Administrator/pyrR.epn', "Result query_name");

$hit = $result->next_hit;
isa_ok($hit, 'Bio::Search::Hit::HitI');
is($hit->accession, 'AE016879.1', "Hit accession");
is($hit->ncbi_gi, '30260185', "Hit GI");
is($hit->algorithm, 'ERPIN', "Hit algorithm");
is($hit->bits, 31.64, "Hit bits");
is($hit->description, 'Bacillus anthracis str. Ames, complete genome',
   "Hit description"); # no hit descs yet
#is($hit->end, '37', "Hit end");
#is($hit->frac_aligned_hit, '', "Hit frac_aligned_hit");
#is($hit->frac_aligned_query, '', "Hit frac_aligned_query");
#is($hit->frac_conserved, '', "Hit frac_conserved");
#is($hit->frac_identical, '', "Hit frac_identical");
#is($hit->frame, 0, "Hit frame");
#is($hit->gaps, 6, "Hit gaps");
#is($hit->length, '', "Hit length");
#is($hit->length_aln, 37, "Hi length_aln");
is($hit->locus, '', "Hit locus");
#is($hit->logical_length, 102, "Hit logical_length");
#is($hit->matches, 2, "Hit matches");
is($hit->n, 4, "Hit n");
is($hit->name, 'gi|30260185|gb|AE016879.1|', "Hit name");
is($hit->num_hsps, 4, "Hit num_hsps"); 
#is($hit->num_unaligned_hit, '', "Hit num_unaligned_hit");
#is($hit->num_unaligned_query, '', "Hit num_unaligned_query");
#is($hit->num_unaligned_sbjct, '', "Hit num_unaligned_subject");
is($hit->overlap, 0, "Hit overlap");
#is($hit->p, '', "Hit p");
#is($hit->query_length, 102, "Hit query_length");
#is($hit->range, 37, "Hit range");
is($hit->rank, 1, "Hit rank");
is($hit->raw_score, 31.64, "Hit raw_score");
is($hit->score, 31.64, "Hit score");
is($hit->seq_inds, 0, "Hit seq_inds");
is($hit->significance, '4.44e-06', "Hit significance");
#is($hit->start, 1, "Hit start");
#is($hit->strand, 0, "Hit strand");

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'ERPIN', "HSP algorithm");
#is($hsp->cigar_string, '103M', "HSP cigar_string");
is($hsp->evalue, '1.68e-05', "HSP evalue");
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
#is($hsp->frac_conserved, '', "HSP frac_conserved");
#is($hsp->frac_identical, '', "HSP frac_identical");
is($hsp->frame, 0, "HSP frame");
is($hsp->gaps, 2, "HSP gaps");
is($hsp->generate_cigar_string, '0M', "HSP generate_cigar_string");
#is($hsp->get_aln, '', "HSP aln");
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
#is($hsp->n, 0, "HSP n");
#is($hsp->num_conserved, '', "HSP num_conserved");
#is($hsp->num_identical, '', "HSP num_identical");
#is($hsp->percent_identity, '', "HSP percent_identity");
#is($hsp->pvalue, '', "HSP pvalue");
isa_ok($hsp->query, 'Bio::SeqFeature::Similarity', "HSP query");
is($hsp->query_string,
   '',
   "HSP query_string");
is($hsp->range, 37, "HSP range");
is($hsp->rank, 1, "HSP rank");
#is($hsp->seq_inds, 67, "HSP seq_inds");
is($hsp->significance, '1.68e-05', "HSP significance");
#is($hsp->end, 37, "HSP end");
is($hsp->expect, '1.68e-05', "HSP expect");
#is($hsp->matches, 2, "HSP matches");
isa_ok($hsp->seq, 'Bio::LocatableSeq');
is($hsp->seq_str,
   '',
   "HSP seq_str");
is($hsp->start, 1, "HSP start");
is($hsp->custom_score, undef, "HSP custom_score");
SKIP: {
    skip('Working on meta string building; TODO', 1);
    is($hsp->meta,
       '',
       "HSP meta");
}
is($hsp->strand, 1, "HSP strand");

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'ERPIN', "HSP algorithm");
#is($hsp->cigar_string, '102M', "HSP cigar_string");
is($hsp->evalue, '5.61e-05', "HSP evalue");
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
#is($hsp->frac_conserved, '', "HSP frac_conserved");
#is($hsp->frac_identical, '', "HSP frac_identical");
is($hsp->frame, 0, "HSP frame");
is($hsp->gaps, 1, "HSP gaps");
is($hsp->generate_cigar_string, '0M', "HSP generate_cigar_string");
#is($hsp->get_aln, '', "HSP aln");
isa_ok($hsp->hit, 'Bio::SeqFeature::Similarity', "HSP hit");
is($hsp->hit_string,
   'CTTT.taatt-.CAGTC.CTGTGA.GACCG.g.AAAG',
   "HSP hit_string");
is($hsp->homology_string,
   '',
   "HSP homology_string");
is($hsp->hsp_group, undef, "HSP hsp_group");
is($hsp->hsp_length, 37, "HSP hsp_length");
is($hsp->length, 37, "HSP length");
is($hsp->links, undef, "HSP links");
#is($hsp->n, 0, "HSP n");
#is($hsp->num_conserved, '', "HSP num_conserved");
#is($hsp->num_identical, '', "HSP num_identical");
#is($hsp->percent_identity, '', "HSP percent_identity");
#is($hsp->pvalue, '', "HSP pvalue");
isa_ok($hsp->query, 'Bio::SeqFeature::Similarity', "HSP query");
is($hsp->query_string,
   '',
   "HSP query_string");
is($hsp->range, 37, "HSP range");
is($hsp->rank, 2, "HSP rank");
#is($hsp->seq_inds, 69, "HSP seq_inds");
is($hsp->significance, '5.61e-05', "HSP significance");
is($hsp->end, 37, "HSP end");
is($hsp->expect, '5.61e-05', "HSP expect");
#is($hsp->matches, 2, "HSP matches");
isa_ok($hsp->seq, 'Bio::LocatableSeq');
is($hsp->seq_str,
   '',
   "HSP seq_str");
is($hsp->start, 1, "HSP start");
is($hsp->custom_score, undef, "HSP custom_score");
SKIP: {
    skip('Working on meta string building; TODO', 1);
    is($hsp->meta,
       '',
       "HSP meta");
}
is($hsp->strand, -1, "HSP strand");

#### RNAMotif ####

# regular data

$searchio = Bio::SearchIO->new( -format => 'rnamotif',
                                -file   => Bio::Root::IO->catfile
                                          ('t','data','trna.strict.rnamotif'),
                                -model => 'trna.descr',
                                -query_acc => 'test',
                                -database => 'gbrna.fas',    
                               );

$result = $searchio->next_result;
isa_ok($result, 'Bio::Search::Result::ResultI');
#$dumper->dumpValue($result);
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
#is($result->query_length, 102, "Result query_length");
is($result->query_name, 'trna.strict.descr', "Result query_name");

$hit = $result->next_hit;
isa_ok($hit, 'Bio::Search::Hit::HitI');
is($hit->accession, 'M10671', "Hit accession");
is($hit->ncbi_gi, '173683', "Hit GI");
is($hit->algorithm, 'RNAMOTIF', "Hit algorithm");
#is($hit->bits, 31.64, "Hit bits");
is($hit->description, 'Avian oncornavirus Trp-tRNA',
   "Hit description"); # no hit descs yet
is($hit->end, '75', "Hit end");
#is($hit->frac_aligned_hit, '', "Hit frac_aligned_hit");
#is($hit->frac_aligned_query, '', "Hit frac_aligned_query");
#is($hit->frac_conserved, '', "Hit frac_conserved");
#is($hit->frac_identical, '', "Hit frac_identical");
is($hit->frame, 0, "Hit frame");
is($hit->gaps, 0, "Hit gaps");
#is($hit->length, '', "Hit length");
is($hit->length_aln, 75, "Hi length_aln");
is($hit->locus, 'ACSTRW', "Hit locus");
#is($hit->logical_length, 102, "Hit logical_length");
is($hit->matches, 2, "Hit matches");
is($hit->n, 8, "Hit n");
is($hit->name, 'gi|173683|gb|M10671|ACSTRW', "Hit name");
is($hit->num_hsps, 8, "Hit num_hsps"); 
#is($hit->num_unaligned_hit, '', "Hit num_unaligned_hit");
#is($hit->num_unaligned_query, '', "Hit num_unaligned_query");
#is($hit->num_unaligned_sbjct, '', "Hit num_unaligned_subject");
is($hit->overlap, 0, "Hit overlap");
#is($hit->p, '', "Hit p");
#is($hit->query_length, 75, "Hit query_length");
is($hit->range, 75, "Hit range");
is($hit->rank, 1, "Hit rank");
is($hit->raw_score, '', "Hit raw_score");
is($hit->score, '', "Hit score");
is($hit->seq_inds, 0, "Hit seq_inds");
is($hit->significance, undef, "Hit significance");
is($hit->start, 1, "Hit start");
is($hit->strand, 0, "Hit strand");

$hsp = $hit->next_hsp;
isa_ok($hsp, 'Bio::Search::HSP::HSPI');
is($hsp->algorithm, 'RNAMOTIF', "HSP algorithm");
#is($hsp->cigar_string, '103M', "HSP cigar_string");
is($hsp->evalue, undef, "HSP evalue");
isa_ok($hsp->feature1, 'Bio::SeqFeature::Similarity');
isa_ok($hsp->feature2, 'Bio::SeqFeature::Similarity');
#is($hsp->frac_conserved, '', "HSP frac_conserved");
#is($hsp->frac_identical, '', "HSP frac_identical");
is($hsp->frame, 0, "HSP frame");
is($hsp->gaps, 0, "HSP gaps");
is($hsp->generate_cigar_string, '0M', "HSP generate_cigar_string");
#is($hsp->get_aln, '', "HSP aln");
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
#is($hsp->n, 0, "HSP n");
#is($hsp->num_conserved, '', "HSP num_conserved");
#is($hsp->num_identical, '', "HSP num_identical");
#is($hsp->percent_identity, '', "HSP percent_identity");
#is($hsp->pvalue, '', "HSP pvalue");
isa_ok($hsp->query, 'Bio::SeqFeature::Similarity', "HSP query");
is($hsp->query_string,
   '',
   "HSP query_string");
is($hsp->range, 75, "HSP range");
is($hsp->rank, 1, "HSP rank");
#is($hsp->seq_inds, 67, "HSP seq_inds");
is($hsp->significance, undef, "HSP significance");
is($hsp->end, 75, "HSP end");
is($hsp->expect, undef, "HSP expect");
is($hsp->matches, 2, "HSP match es");
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

