use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    test_begin(-tests => 134);
    use_ok('Bio::Assembly::IO');
    use_ok('Bio::Assembly::Tools::ContigSpectrum');
}

my $in = Bio::Assembly::IO->new(
  -file => test_input_file("contigspectrumtest.asm"),
  -format => 'tigr'
);
isa_ok($in, 'Bio::Assembly::IO');
my $sc = $in->next_assembly;
isa_ok($sc, 'Bio::Assembly::Scaffold');

# Try all the get/set methods
ok(my $csp = Bio::Assembly::Tools::ContigSpectrum->new);
isa_ok($csp, 'Bio::Assembly::Tools::ContigSpectrum');
ok($csp->id('asdf'));
is($csp->id, 'asdf');
ok($csp->nof_seq(123));
is($csp->nof_seq, 123);
ok($csp->nof_rep(456));
is($csp->nof_rep, 456);
ok($csp->max_size(789));
is($csp->max_size, 789);
ok($csp->nof_overlaps(111));
is($csp->nof_overlaps, 111);
ok($csp->min_overlap(50));
is($csp->min_overlap, 50);
ok($csp->avg_overlap(54.3));
is($csp->avg_overlap, 54.3);
ok($csp->min_identity(89.1));
is($csp->min_identity, 89.1);
ok($csp->avg_identity(98.7));
is($csp->avg_identity, 98.7);
ok($csp->avg_seq_len(123.456));
is($csp->avg_seq_len, 123.456);
ok($csp->eff_asm_params(1));
is($csp->eff_asm_params, 1);

# contig spectrum based on simple spectrum
ok(my $spectrum_csp = Bio::Assembly::Tools::ContigSpectrum->new);
ok($spectrum_csp->spectrum({1=>1, 2=>2, 3=>3}));
is($spectrum_csp->eff_asm_params, 0);
is($spectrum_csp->nof_seq, 14);
is($spectrum_csp->max_size, 3);
is($spectrum_csp->nof_rep, 1);
is($spectrum_csp->nof_overlaps, 0);
is($spectrum_csp->min_overlap, undef);
is($spectrum_csp->avg_overlap, 0);
is($spectrum_csp->min_identity, undef);
is($spectrum_csp->avg_identity, 0);
is($spectrum_csp->avg_seq_len, 0);

ok(my $string = $spectrum_csp->to_string(1));
is($string, '1 2 3');
ok($string = $spectrum_csp->to_string(2));
is($string, "1\t2\t3");
ok($string = $spectrum_csp->to_string(3));
is($string, "1\n2\n3");

# mixed contig spectrum imported from assembly
ok(my $mixed_csp = Bio::Assembly::Tools::ContigSpectrum->new(
  -assembly       => $sc,
  -eff_asm_params => 1 ));
is(scalar @{$mixed_csp->assembly()}, 1);
is_deeply($mixed_csp->spectrum, {1=>0, 9=>1}); # [0 0 0 0 0 0 0 0 1]
is($mixed_csp->eff_asm_params, 1);
is($mixed_csp->nof_seq, 9);
is($mixed_csp->max_size, 9);
is($mixed_csp->nof_rep, 1);
is($mixed_csp->nof_overlaps, 8);
is($mixed_csp->min_overlap, 35);
is($mixed_csp->avg_overlap, 71.875);
float_is($mixed_csp->min_identity, 97.1153846153846);
float_is($mixed_csp->avg_identity, 99.6394230769231);
float_is($mixed_csp->avg_seq_len, 100.222222222222);

# dissolved contig spectrum
ok(my $dissolved_csp = Bio::Assembly::Tools::ContigSpectrum->new(
  -dissolve => [$mixed_csp, 'ZZZ'] ));
is_deeply($dissolved_csp->spectrum, {1=>1}); # [1]
is($dissolved_csp->eff_asm_params, 0);
is($dissolved_csp->nof_seq, 1);
is($dissolved_csp->max_size, 1);
is($dissolved_csp->nof_rep, 1);
is($dissolved_csp->nof_overlaps, 0);
is($dissolved_csp->min_overlap, undef);
is($dissolved_csp->avg_overlap, 0);
is($dissolved_csp->min_identity, undef);
is($dissolved_csp->avg_identity, 0);
is($dissolved_csp->avg_seq_len, 96);

ok($dissolved_csp = Bio::Assembly::Tools::ContigSpectrum->new(
  -dissolve => [$mixed_csp, 'ABC'] ));
is_deeply($dissolved_csp->spectrum, {1=>1, 7=>1}); # [1 0 0 0 0 0 1]
is($dissolved_csp->eff_asm_params, 0);
is($dissolved_csp->nof_seq, 8);
is($dissolved_csp->max_size, 7);
is($dissolved_csp->nof_rep, 1);
is($dissolved_csp->nof_overlaps, 0);
is($dissolved_csp->min_overlap, undef);
is($dissolved_csp->avg_overlap, 0);
is($dissolved_csp->min_identity, undef);
is($dissolved_csp->avg_identity, 0);
is($dissolved_csp->avg_seq_len, 100.75);

ok($dissolved_csp = Bio::Assembly::Tools::ContigSpectrum->new(
  -min_overlap  => 62,
  -min_identity => 1,
  -dissolve     => [$mixed_csp, 'ABC'] ));
is_deeply($dissolved_csp->spectrum, {1=>1, 7=>1}); # [1 0 0 0 0 0 1]

ok($dissolved_csp = Bio::Assembly::Tools::ContigSpectrum->new(
   -min_overlap  => 63,
   -min_identity => 1,
   -dissolve     => [$mixed_csp, 'ABC'] ));
is_deeply($dissolved_csp->spectrum, {1=>1, 2=>1, 5=>1}); # [1 1 0 0 1]

ok($dissolved_csp = Bio::Assembly::Tools::ContigSpectrum->new(
  -min_overlap  => 62,
  -min_identity => 97,
  -dissolve     => [$mixed_csp, 'ABC'] ));
is_deeply($dissolved_csp->spectrum, {1=>1, 7=>1}); # [1 0 0 0 0 0 1]

ok($dissolved_csp = Bio::Assembly::Tools::ContigSpectrum->new(
   -min_overlap  => 62,
   -min_identity => 98,
   -dissolve     => [$mixed_csp, 'ABC'] ));
is_deeply($dissolved_csp->spectrum, {1=>2, 6=>1}); # [2 0 0 0 0 1]

ok($dissolved_csp = Bio::Assembly::Tools::ContigSpectrum->new(
  -dissolve       => [$mixed_csp, 'ABC'],
  -eff_asm_params => 1 ));
is_deeply($dissolved_csp->spectrum, {1=>1, 7=>1}); # [1 0 0 0 0 0 1]
is($dissolved_csp->eff_asm_params, 1);
is($dissolved_csp->nof_seq, 8);
is($dissolved_csp->max_size, 7);
is($dissolved_csp->nof_rep, 1);
is($dissolved_csp->nof_overlaps, 6);
is($dissolved_csp->min_overlap, 62);
float_is($dissolved_csp->avg_overlap, 81.3333333333333);
float_is($dissolved_csp->min_identity, 97.1153846153846);
float_is($dissolved_csp->avg_identity, 99.5192307692308);
float_is($dissolved_csp->avg_seq_len, 100.75);

# cross contig spectrum
ok(my $cross_csp = Bio::Assembly::Tools::ContigSpectrum->new(
  -cross => $mixed_csp));
is_deeply($cross_csp->spectrum, {1=>2, 9=>1}); # [2 0 0 0 0 0 0 0 1]

# sum of contig spectra
ok(my $sum_csp = Bio::Assembly::Tools::ContigSpectrum->new(-eff_asm_params=>1));
ok($sum_csp->add($dissolved_csp));
ok($sum_csp->add($mixed_csp));
is_deeply($sum_csp->spectrum, {1=>1, 7=>1, 9=>1}); # [1 0 0 0 0 0 1 0 1]
is($sum_csp->eff_asm_params, 1);
is($sum_csp->nof_seq, 17);
is($sum_csp->max_size, 9);
is($sum_csp->nof_rep, 2);
is($sum_csp->nof_overlaps, 14);
is($sum_csp->min_overlap, 35);
float_is($sum_csp->avg_overlap, 75.9285714285714);
float_is($sum_csp->min_identity, 97.1153846153846);
float_is($sum_csp->avg_identity, 99.5879120879121);
float_is($sum_csp->avg_seq_len, 100.470588235294);

# average of contig spectra
ok(my $avg_csp = Bio::Assembly::Tools::ContigSpectrum->new);
ok($avg_csp = $avg_csp->average([$dissolved_csp, $mixed_csp]));
is_deeply($avg_csp->spectrum, {1=>0.5, 7=>0.5, 9=>0.5}); # [0.5 0 0 0 0 0 0.5 0 0.5]
is($avg_csp->eff_asm_params, 1);
is($avg_csp->nof_seq, 8.5);
is($avg_csp->max_size, 9);
is($avg_csp->nof_rep, 2);
is($avg_csp->nof_overlaps, 7);
is($avg_csp->min_overlap, 35);
float_is($avg_csp->avg_overlap, 75.9285714285714);
float_is($avg_csp->min_identity, 97.1153846153846);
float_is($avg_csp->avg_identity, 99.5879120879121);
float_is($avg_csp->avg_seq_len, 100.470588235294);

# drop assembly info from contig spectrum
ok($mixed_csp->drop_assembly());
is(scalar @{$mixed_csp->assembly()}, 0);
