# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 66);
	
	use_ok('Bio::Tools::TandemRepeatsFinder');
}

# first, open empty output file
# make sure no results get returned
my $trf =
   Bio::Tools::TandemRepeatsFinder->new( -file => test_input_file('tandem_repeats_finder.noresults')
  );
ok $trf,"Parser created  successfully" ;
my $feat = $trf->next_result;
ok( !defined($feat), "empty results file correctly returns no results" );


# now check some results
$trf =
   Bio::Tools::TandemRepeatsFinder->new( -file => test_input_file('tandem_repeats_finder.dat')
  );
ok $trf,"Second parser created  successfully" ;
my $feat1 = $trf->next_result();

# these are the parameters that should be parsed from 
# the following line:
# Parameters: 2 7 7 80 10 50 12
my $expected_run_parameters = {
                match_weight    => 2,
                mismatch_weight => 7,
                indel_weight    => 7,
                match_prob      => 80,
                indel_prob      => 10,
                min_score       => 50,
                max_period_size => 12
};

# test feature properties
is ( $feat1->seq_id(),      "DDB0169550",             "seq_id for first result correctly parsed");
is ( $feat1->start(),       13936,                    "start for first result correctly parsed");
is ( $feat1->end(),         13960,                    "end for first result correctly parsed");
is ( $feat1->source_tag(),  'Tandem Repeats Finder',  "source tag for first result correctly parsed");
is ( $feat1->primary_tag(), 'tandem repeat',          "primary tag for first result correctly parsed");
is ( $feat1->score(),       50,                       "score for first result correctly parsed");

# test tag values
# all of the data other than start, end, score, and seq_id 
# is stored in tags
my ($seqence_description) = $feat1->get_tag_values( 'sequence_description' );
is ( $seqence_description, "|Masked Chromosomal Sequence| on chromosome: M", "sequence description correctly parsed.");
my ($parameters) = $feat1->get_tag_values( 'run_parameters' );
is_deeply ( $parameters, $expected_run_parameters ,"correctly parsed all run parameters");
my ($period_size) = $feat1->get_tag_values( 'period_size' );
is ( $period_size, 12 ,"correctly parsed period_size for first result");
my ($copy_number) = $feat1->get_tag_values( 'copy_number' );
is ( $copy_number, 2.1 ,"correctly parsed copy_number for first result");
my ($consensus_size) = $feat1->get_tag_values( 'consensus_size' );
is ( $consensus_size, 12 ,"correctly parsed consensus_size for first result");
my ($percent_matches) = $feat1->get_tag_values( 'percent_matches' );
is ( $percent_matches, 100 ,"correctly parsed percent_matches for first result");
my ($percent_indels) = $feat1->get_tag_values( 'percent_indels' );
is ( $percent_indels, 0 ,"correctly parsed percent_indels for first result");
my ($percent_a) = $feat1->get_tag_values( 'percent_a' );
is ( $percent_a, 16 ,"correctly parsed percent_a for first result");
my ($percent_c) = $feat1->get_tag_values( 'percent_c' );
is ( $percent_c, 8 ,"correctly parsed percent_c for first result");
my ($percent_g) = $feat1->get_tag_values( 'percent_g' );
is ( $percent_g, 52 ,"correctly parsed percent_g for first result");
my ($percent_t) = $feat1->get_tag_values( 'percent_t' );
is ( $percent_t, 24 ,"correctly parsed percent_t for first result");
my ($entropy) = $feat1->get_tag_values( 'entropy' );
is ( $entropy, "1.70", "correctly parsed entropy for first result");
my ($repeat_sequence) = $feat1->get_tag_values( 'repeat_sequence' );
is ( $repeat_sequence, "GGCGTAATGGGTGGCGTAATGGGTG", "correctly parsed repeat_sequence for first result");
my ($consensus_sequence) = $feat1->get_tag_values( 'consensus_sequence' );
is ( $consensus_sequence, "GGCGTAATGGGT", "correctly parsed consensus_sequence for first result");

# for the next result just check the sequence _id, start, and end
my $feat2 = $trf->next_result();
is ( $feat2->seq_id(),      "DDB0169550",        "seq_id for second result correctly parsed");
is ( $feat2->start(),       "16937",             "start for second result correctly parsed");
is ( $feat2->end(),         "16965",             "end for second result correctly parsed");
is ( $feat2->source_tag(),  'Tandem Repeats Finder',  "source tag for first result correctly parsed");
is ( $feat2->primary_tag(), 'tandem repeat',          "primary tag for first result correctly parsed");
is ( $feat2->score(),       58,                       "score for first result correctly parsed");

# test tag values
# all of the data other than start, end, score, and seq_id 
# is stored in tags
($seqence_description) = $feat2->get_tag_values( 'sequence_description' );
is ( $seqence_description, "|Masked Chromosomal Sequence| on chromosome: M", "sequence description correctly parsed.");

($parameters) = $feat2->get_tag_values( 'run_parameters' );
is_deeply ( $parameters, $expected_run_parameters ,"correctly reatained all run parameters for second feature");
($period_size) = $feat2->get_tag_values( 'period_size' );
is ( $period_size, 9 ,"correctly parsed period_size for second result");
($copy_number) = $feat2->get_tag_values( 'copy_number' );
is ( $copy_number, "3.2" ,"correctly parsed copy_number for second result");
($consensus_size) = $feat2->get_tag_values( 'consensus_size' );
is ( $consensus_size, 9 ,"correctly parsed consensus_size for second result");
($percent_matches) = $feat2->get_tag_values( 'percent_matches' );
is ( $percent_matches, 100 ,"correctly parsed percent_matches for second result");
($percent_indels) = $feat2->get_tag_values( 'percent_indels' );
is ( $percent_indels, 0 ,"correctly parsed percent_indels for second result");
($percent_a) = $feat2->get_tag_values( 'percent_a' );
is ( $percent_a, 44 ,"correctly parsed percent_a for second result");
($percent_c) = $feat2->get_tag_values( 'percent_c' );
is ( $percent_c, 0 ,"correctly parsed percent_c for second result");
($percent_g) = $feat2->get_tag_values( 'percent_g' );
is ( $percent_g, 10 ,"correctly parsed percent_g for second result");
($percent_t) = $feat2->get_tag_values( 'percent_t' );
is ( $percent_t, 44 ,"correctly parsed percent_t for second result");
($entropy) = $feat2->get_tag_values( 'entropy' );
is ( $entropy, "1.38", "correctly parsed entropy for second result");
($repeat_sequence) = $feat2->get_tag_values( 'repeat_sequence' );
is ( $repeat_sequence, "TATATAGTATATATAGTATATATAGTATA", "correctly parsed repeat_sequence for second result");
($consensus_sequence) = $feat2->get_tag_values( 'consensus_sequence' );
is ( $consensus_sequence, "TATATAGTA", "correctly parsed consensus_sequence for second result");

# now, check the full results again for last result (on a different sequence).
my $feat3 = $trf->next_result();
# test feature properties
is ( $feat3->seq_id(),      "DDB0215018",             "seq_id for first result correctly parsed");
is ( $feat3->start(),       1649,                    "start for first result correctly parsed");
is ( $feat3->end(),         1679,                    "end for first result correctly parsed");
is ( $feat3->source_tag(),  'Tandem Repeats Finder',  "source tag for first result correctly parsed");
is ( $feat3->primary_tag(), 'tandem repeat',          "primary tag for first result correctly parsed");
is ( $feat3->score(),       62,                       "score for first result correctly parsed");

# test tag values
# all of the data other than start, end, score, and seq_id 
# is stored in tags
($seqence_description) = $feat3->get_tag_values( 'sequence_description' );
is ( $seqence_description, "|Masked Chromosomal Sequence| on chromosome: 2F", "sequence description correctly parsed.");

($parameters) = $feat3->get_tag_values( 'run_parameters' );
is_deeply ( $parameters, $expected_run_parameters ,"correctly reatained all run parameters for third feature");
($period_size) = $feat3->get_tag_values( 'period_size' );
is ( $period_size, 1 ,"correctly parsed period_size for third result");
($copy_number) = $feat3->get_tag_values( 'copy_number' );
is ( $copy_number, "31.0" ,"correctly parsed copy_number for third result");
($consensus_size) = $feat3->get_tag_values( 'consensus_size' );
is ( $consensus_size, 1 ,"correctly parsed consensus_size for third result");
($percent_matches) = $feat3->get_tag_values( 'percent_matches' );
is ( $percent_matches, 100 ,"correctly parsed percent_matches for third result");
($percent_indels) = $feat3->get_tag_values( 'percent_indels' );
is ( $percent_indels, 0 ,"correctly parsed percent_indels for third result");
($percent_a) = $feat3->get_tag_values( 'percent_a' );
is ( $percent_a, 0 ,"correctly parsed percent_a for third result");
($percent_c) = $feat3->get_tag_values( 'percent_c' );
is ( $percent_c, 0 ,"correctly parsed percent_c for third result");
($percent_g) = $feat3->get_tag_values( 'percent_g' );
is ( $percent_g, 0 ,"correctly parsed percent_g for third result");
($percent_t) = $feat3->get_tag_values( 'percent_t' );
is ( $percent_t, 100 ,"correctly parsed percent_t for third result");
($entropy) = $feat3->get_tag_values( 'entropy' );
is ( $entropy, "0.00", "correctly parsed entropy for third result");
($repeat_sequence) = $feat3->get_tag_values( 'repeat_sequence' );
is ( $repeat_sequence, "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", "correctly parsed repeat_sequence for third result");
($consensus_sequence) = $feat3->get_tag_values( 'consensus_sequence' );
is ( $consensus_sequence, "T", "correctly parsed consensus_sequence for third result");

my $empty_feat = $trf->next_result();
ok( !defined($empty_feat), "correctly return undef when no features are left" );


 $trf =
   Bio::Tools::TandemRepeatsFinder->new( -file => test_input_file('tandem_repeats_finder_no_desc.dat')
  );

my $feat_with_seqid = $trf->next_result();
#ensuring that we can parse out seq_id in the file even in the absent of description
is( $feat_with_seqid->seq_id(),"DDB0169550", "Correctly parsed seq_id even if description does not exist" );
