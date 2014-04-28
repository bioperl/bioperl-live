use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    test_begin( -tests            => 239,
                -requires_modules => [ 'DB_File',
                                       qw(Bio::Assembly::Tools::ContigSpectrum)] );
    use_ok 'Bio::Assembly::IO';
    use_ok 'Bio::Assembly::Tools::ContigSpectrum';
}


my $in = Bio::Assembly::IO->new(
   -file => test_input_file('contigspectrumtest.tigr'),
   -format => 'tigr'
);
isa_ok $in, 'Bio::Assembly::IO';
my $sc = $in->next_assembly;
isa_ok $sc, 'Bio::Assembly::Scaffold';

# Try all the get/set methods
ok my $csp = Bio::Assembly::Tools::ContigSpectrum->new, 'get/set methods';
isa_ok $csp, 'Bio::Assembly::Tools::ContigSpectrum';
ok $csp->id('asdf');
is $csp->id, 'asdf';
ok $csp->nof_seq(123);
is $csp->nof_seq, 123;
ok $csp->nof_rep(456);
is $csp->nof_rep, 456;
ok $csp->max_size(789);
is $csp->max_size, 789;
ok $csp->nof_overlaps(111);
is $csp->nof_overlaps, 111;
ok $csp->min_overlap(50);
is $csp->min_overlap, 50;
ok $csp->avg_overlap(54.3);
is $csp->avg_overlap, 54.3;
ok $csp->min_identity(89.1);
is $csp->min_identity, 89.1;
ok $csp->avg_identity(98.7);
is $csp->avg_identity, 98.7;
ok $csp->avg_seq_len(123.456);
is $csp->avg_seq_len, 123.456;
ok $csp->eff_asm_params(1);
is $csp->eff_asm_params, 1;

# contig spectrum based on simple spectrum
ok my $spectrum_csp = Bio::Assembly::Tools::ContigSpectrum->new, 'simple spectrum';
ok $spectrum_csp->spectrum({1=>1, 2=>2, 3=>3});
is $spectrum_csp->eff_asm_params, 0;
is $spectrum_csp->nof_seq, 14;
is $spectrum_csp->max_size, 3;
is $spectrum_csp->nof_rep, 1;
is $spectrum_csp->nof_overlaps, 0;
is $spectrum_csp->min_overlap, undef;
is $spectrum_csp->avg_overlap, 0;
is $spectrum_csp->min_identity, undef;
is $spectrum_csp->avg_identity, 0;
is $spectrum_csp->avg_seq_len, 0;
is scalar $spectrum_csp->assembly, 0;

ok my $string = $spectrum_csp->to_string(1);
is $string, '1 2 3';
ok $string = $spectrum_csp->to_string(2);
is $string, "1\t2\t3";
ok $string = $spectrum_csp->to_string(3);
is $string, "1\n2\n3";

# score
my $test_csp;
my $spectrum;
ok $test_csp = Bio::Assembly::Tools::ContigSpectrum->new(-spectrum=>$spectrum), 'contig spectrum score';
is $test_csp->score, undef;
$spectrum = {1=>120};
ok $test_csp = Bio::Assembly::Tools::ContigSpectrum->new(-spectrum=>$spectrum);
is $test_csp->score, 0;
$spectrum = {120=>1};
ok $test_csp = Bio::Assembly::Tools::ContigSpectrum->new(-spectrum=>$spectrum);
is $test_csp->score, 1;
float_is $test_csp->score(240), 0.248953974895397;
$spectrum = {1=>120, 120=>1};
ok $test_csp = Bio::Assembly::Tools::ContigSpectrum->new(-spectrum=>$spectrum);
float_is $test_csp->score, 0.248953974895397;

SKIP: {
   test_skip( -tests => 183, -requires_module => 'Graph::Undirected' ); #####

   # mixed contig spectrum imported from assembly
   ok my $mixed_csp = Bio::Assembly::Tools::ContigSpectrum->new(
      -assembly       => $sc,
      -eff_asm_params => 1 ), 'mixed contig spectrum';
   is_deeply $mixed_csp->spectrum, {1=>0, 2=>3, 6=>1, 9=>1}; # [0 3 0 0 0 1 0 0 1]
   is $mixed_csp->eff_asm_params, 1;
   is $mixed_csp->max_size, 9;
   is $mixed_csp->nof_rep, 1;
   is $mixed_csp->nof_seq, 21;
   float_is $mixed_csp->avg_seq_len, 303.81;
   is $mixed_csp->nof_overlaps, 16;
   is $mixed_csp->min_overlap, 35;
   float_is $mixed_csp->avg_overlap, 155.875;
   float_is $mixed_csp->min_identity, 96.8421;
   float_is $mixed_csp->avg_identity, 98.8826;
   is scalar $mixed_csp->assembly, 1;

   # dissolved contig spectrum
   ok my $dissolved_csp = Bio::Assembly::Tools::ContigSpectrum->new(
      -dissolve => [$mixed_csp, 'ZZZ'] ), 'dissolved contig spectrum';
   is_deeply $dissolved_csp->spectrum, {1=>2, 2=>1}; # [2 1]
   is $dissolved_csp->eff_asm_params, 0;
   is $dissolved_csp->max_size, 2;
   is $dissolved_csp->nof_rep, 1;
   is $dissolved_csp->nof_seq, 4;
   float_is $dissolved_csp->avg_seq_len, 321;
   # eff_asm_params haven't been requested
   is $dissolved_csp->nof_overlaps, 0;
   is $dissolved_csp->min_overlap, undef;
   is $dissolved_csp->avg_overlap, 0;
   is $dissolved_csp->min_identity, undef;
   is $dissolved_csp->avg_identity, 0;

   ok $dissolved_csp = Bio::Assembly::Tools::ContigSpectrum->new(
      -dissolve => [$mixed_csp, 'sdsu'] );
   is_deeply $dissolved_csp->spectrum, {1=>3, 6=>1}; # [3 0 0 0 0 1]
   is $dissolved_csp->eff_asm_params, 0;
   is $dissolved_csp->max_size, 6;
   is $dissolved_csp->nof_rep, 1;
   is $dissolved_csp->nof_seq, 9;
   float_is $dissolved_csp->avg_seq_len, 441.222222222222;
   # eff_asm_params haven't been requested
   is $dissolved_csp->nof_overlaps, 0;
   is $dissolved_csp->min_overlap, undef;
   is $dissolved_csp->avg_overlap, 0;
   is $dissolved_csp->min_identity, undef;
   is $dissolved_csp->avg_identity, 0;

   ok $dissolved_csp = Bio::Assembly::Tools::ContigSpectrum->new(
      -dissolve => [$mixed_csp, 'ABC'] );
   is_deeply $dissolved_csp->spectrum, {1=>2, 6=>1}; # [2 0 0 0 0 1]
   is $dissolved_csp->eff_asm_params, 0;
   is $dissolved_csp->max_size, 6;
   is $dissolved_csp->nof_rep, 1;
   is $dissolved_csp->nof_seq, 8;
   float_is $dissolved_csp->avg_seq_len, 140.625;
   # eff_asm_params haven't been requested
   is $dissolved_csp->nof_overlaps, 0;
   is $dissolved_csp->min_overlap, undef;
   is $dissolved_csp->avg_overlap, 0;
   is $dissolved_csp->min_identity, undef;
   is $dissolved_csp->avg_identity, 0;

   ok $dissolved_csp = Bio::Assembly::Tools::ContigSpectrum->new(
      -min_overlap  => 62,
      -min_identity => 1,
      -dissolve     => [$mixed_csp, 'ABC'] );
   is_deeply $dissolved_csp->spectrum, {1=>2, 6=>1}; # [2 0 0 0 0 1]

   ok $dissolved_csp = Bio::Assembly::Tools::ContigSpectrum->new(
      -min_overlap  => 63,
      -min_identity => 1,
      -dissolve     => [$mixed_csp, 'ABC'] );
   is_deeply $dissolved_csp->spectrum, {1=>3, 5=>1}; # [3 0 0 0 1]

   # after dissolving, the remaining assembly objects should be 3 singlets and 1 6-contig
   my @contigs = ($dissolved_csp->assembly);
   is scalar @contigs, 4;
   my @contig_ids = sort qw( 144 652_1 652_2 652_3 );
   is_deeply [sort map($_->id, @contigs)], \@contig_ids;
   my @contig_sizes = sort qw( 1 1 1 5 );
   is_deeply [sort map($_->num_sequences, @contigs)], \@contig_sizes;
   my @contig_isas = sort qw( Bio::Assembly::Singlet Bio::Assembly::Singlet
   Bio::Assembly::Singlet Bio::Assembly::Contig );
   is_deeply [sort map(ref $_, @contigs)], \@contig_isas;
   my @reads = ($contigs[1])->each_seq;
   my @read_ids = sort qw(ABC|9980040 ABC|9937790 ABC|9956706 ABC|9960711 ABC|9976538);
   is_deeply [sort map($_->id, @reads)], \@read_ids;

   ok $dissolved_csp = Bio::Assembly::Tools::ContigSpectrum->new(
      -min_overlap  => 62,
      -min_identity => 97,
      -dissolve     => [$mixed_csp, 'ABC'] );
   is_deeply $dissolved_csp->spectrum, {1=>2, 6=>1}; # [2 0 0 0 0 1]

   ok $dissolved_csp = Bio::Assembly::Tools::ContigSpectrum->new(
      -min_overlap  => 62,
      -min_identity => 98,
      -dissolve     => [$mixed_csp, 'ABC'] );
   is_deeply $dissolved_csp->spectrum, {1=>2, 6=>1}; # [2 0 0 0 0 1]

   ok $dissolved_csp = Bio::Assembly::Tools::ContigSpectrum->new(
      -dissolve       => [$mixed_csp, 'ABC'],
      -eff_asm_params => 1 );
   is_deeply $dissolved_csp->spectrum, {1=>2, 6=>1}; # [2 0 0 0 0 1]
   is $dissolved_csp->eff_asm_params, 1;
   is $dissolved_csp->max_size, 6;
   is $dissolved_csp->nof_rep, 1;
   is $dissolved_csp->nof_seq, 8;
   float_is $dissolved_csp->avg_seq_len, 140.625;
   is $dissolved_csp->nof_overlaps, 5;
   float_is $dissolved_csp->avg_overlap, 76.8;
   float_is $dissolved_csp->avg_identity, 100.0;
   # min_overlap and min_identity not explicitely specified for the dissolved csp
   # min_overlap and min_identity are thus taken from the mixed csp
   is $dissolved_csp->min_overlap, 35;
   float_is $dissolved_csp->min_identity, 96.8421;

   # cross contig spectrum
   ok my $cross_csp = Bio::Assembly::Tools::ContigSpectrum->new(
      -cross => $mixed_csp), 'cross-contig spectrum';
   is_deeply $cross_csp->spectrum, {1=>7, 2=>2, 9=>1}; # [7 2 0 0 0 0 0 0 1]

   # assembly should have 2 singlets and 1 9-contig
   @contigs = $cross_csp->assembly;
   is scalar @contigs, 3;
   @contig_sizes = sort qw( 2 2 9 );
   is_deeply [sort map($_->num_sequences, @contigs)], \@contig_sizes;
   @contig_isas = sort qw( Bio::Assembly::Contig Bio::Assembly::Contig Bio::Assembly::Contig);
   is_deeply [sort map(ref $_, @contigs)], \@contig_isas;
   @read_ids = sort qw(sdsu|SDSU_RFPERU_006_E04.x01.phd.1 ZZZ|SDSU_RFPERU_010_B05.x01.phd.1);
   is_deeply [sort map($_->id, $contigs[0]->each_seq)], \@read_ids;
   @read_ids = sort qw(sdsu|SDSU_RFPERU_013_H05.x01.phd.1 ABC|SDSU_RFPERU_005_F02.x01.phd.1);
   is_deeply [sort map($_->id, $contigs[1]->each_seq)], \@read_ids;
   @read_ids = sort qw( ZZZ|9962187 ABC|9937790 ABC|9944760 ABC|9956706
           sdsu|9986984 ABC|9960711 ABC|9970175 ABC|9976538 ABC|9980040);
   is_deeply [sort map($_->id, $contigs[2]->each_seq)], \@read_ids;

   # effective assembly params
   ok $cross_csp = Bio::Assembly::Tools::ContigSpectrum->new(
      -cross          => $mixed_csp,
      -eff_asm_params => 1 ), 'cross-contig spectrum';
   is_deeply $cross_csp->spectrum, {1=>7, 2=>2, 9=>1}; # [7 2 0 0 0 0 0 0 1]
   is $cross_csp->nof_rep, 1;
   is $cross_csp->eff_asm_params, 1;
   is $cross_csp->max_size, 9;
   is $cross_csp->nof_seq, 13;
   float_is $cross_csp->avg_seq_len, 206.308;
   is $cross_csp->nof_overlaps, 10;
   float_is $cross_csp->avg_overlap, 76.9;
   float_is $cross_csp->avg_identity, 99.2357;
   # min_overlap and min_identity not explicitly specified for the cross csp
   # min_overlap and min_identity are thus taken from the mixed csp
   is $cross_csp->min_overlap, 35;
   float_is $cross_csp->min_identity, 96.8421;

   # with a specified minimum overlap and identity
   ok $cross_csp = Bio::Assembly::Tools::ContigSpectrum->new(
      -cross          => $mixed_csp,
      -min_overlap    => 50,
      -min_identity   => 98 ), 'cross-contig spectrum';
   is_deeply $cross_csp->spectrum, {1=>3, 2=>1, 7=>1}; # [3 1 0 0 0 0 1]
   is $cross_csp->nof_rep, 1;
   is $cross_csp->eff_asm_params, 0;
   is $cross_csp->max_size, 7;
   is $cross_csp->nof_seq, 9;
   float_is $cross_csp->avg_seq_len, 191.222;
   is $cross_csp->min_overlap, 50;
   float_is $cross_csp->min_identity, 98;

   # sum of contig spectra
   ok my $sum_csp = Bio::Assembly::Tools::ContigSpectrum->new(-eff_asm_params=>1), 'contig spectrum sum';
   ok $sum_csp->add($dissolved_csp);
   ok $sum_csp->add($mixed_csp);
   is_deeply $sum_csp->spectrum, {1=>2, 2=>3, 6=>2, 9=>1}; # [2 3 0 0 0 2 0 0 1]
   is $sum_csp->eff_asm_params, 1;
   is $sum_csp->max_size, 9;
   is $sum_csp->nof_rep, 2;
   is $sum_csp->nof_seq, 29;
   float_is $sum_csp->avg_seq_len, 258.7934;
   is $sum_csp->nof_overlaps, 21;
   is $sum_csp->min_overlap, 35;
   float_is $sum_csp->avg_overlap, 137.0476;
   float_is $sum_csp->min_identity, 96.8421;
   float_is $sum_csp->avg_identity, 99.1487;
   is scalar $sum_csp->assembly, 4;

   # average of contig spectra
   ok my $avg_csp = Bio::Assembly::Tools::ContigSpectrum->new(-eff_asm_params=>1), 'average contig spectrum';
   ok $avg_csp = $avg_csp->average([$dissolved_csp, $mixed_csp]);
   is_deeply $avg_csp->spectrum, {1=>1, 2=>1.5, 6=>1, 9=>0.5}; # [1 1 0 0 0 1 0 0 0.5]
   is $avg_csp->eff_asm_params, 1;
   is $avg_csp->max_size, 9;
   is $avg_csp->nof_rep, 2;
   is $avg_csp->nof_seq, 14.5;
   float_is $avg_csp->avg_seq_len, 258.7934;
   is $avg_csp->nof_overlaps, 10.5;
   is $avg_csp->min_overlap, 35;
   float_is $avg_csp->avg_overlap, 137.0476;
   float_is $avg_csp->min_identity, 96.8421;
   float_is $avg_csp->avg_identity, 99.1487;
   is scalar $avg_csp->assembly, 4;

   # drop assembly info from contig spectrum
   ok $mixed_csp->drop_assembly(), 'drop assembly';
   is scalar $mixed_csp->assembly(), 0;

   # large contig (27 reads)
   $in = Bio::Assembly::IO->new(
      -file   => test_input_file('27-contig_Newbler.ace'),
      -format => 'ace-454'
   );
   isa_ok $in, 'Bio::Assembly::IO';
   $sc = $in->next_assembly;
   isa_ok $sc, 'Bio::Assembly::Scaffold';
   ok my $large_csp = Bio::Assembly::Tools::ContigSpectrum->new(
      -assembly       => $sc,
      -eff_asm_params => 1 ), 'large contig spectrum';
   is scalar $large_csp->assembly(), 1;
   is_deeply $large_csp->spectrum, {1=>0, 27=>1};
   is $large_csp->eff_asm_params, 1;
   is $large_csp->max_size, 27;
   is $large_csp->nof_rep, 1;
   is $large_csp->nof_seq, 27;
   float_is $large_csp->avg_seq_len, 100;
   is $large_csp->nof_overlaps, 26;
   is $large_csp->min_overlap, 54;
   # operation returns sometimes 88.76923... and sometimes 88.80769...
   ok $large_csp->avg_overlap >= 88.7692;
   ok $large_csp->avg_overlap <= 88.8077;
   float_is $large_csp->min_identity, 33.3333;
   cmp_ok($large_csp->avg_identity, '>=', 70, $large_csp->avg_identity);
   cmp_ok($large_csp->avg_identity, '<=', 80, $large_csp->avg_identity);

   ok my $large_xcsp = Bio::Assembly::Tools::ContigSpectrum->new(
      -cross          => $large_csp,
      -eff_asm_params => 1           ), 'large cross-contig spectrum';
   is $large_xcsp->nof_overlaps, 26;
   cmp_ok($large_xcsp->avg_overlap, '>=', 88.7, $large_xcsp->avg_overlap);
   cmp_ok($large_xcsp->avg_overlap, '<=', 88.9, $large_xcsp->avg_overlap);
   is_deeply $large_xcsp->spectrum, {1=>21, 27=>1};

   ok $large_xcsp = Bio::Assembly::Tools::ContigSpectrum->new(
      -cross          => $large_csp,
      -min_overlap    => 100);
   is_deeply $large_xcsp->spectrum, {1=>18, 2=>5, 3=>1, 7=>1};
   my @xcontigs = sort {$a->id cmp $b->id} $large_xcsp->assembly;
   is scalar @xcontigs, 7; # the cross-1-contigs are not included
   my @xcontig_ids = sort qw( contig00001_1 contig00001_2 contig00001_3 contig00001_4
      contig00001_5 contig00001_6 contig00001_7 );
   is_deeply [map($_->id, @xcontigs)], \@xcontig_ids;
   my @xcontig_sizes = sort qw( 2 2 2 2 2 3 7 );
   is_deeply [sort map($_->num_sequences, @xcontigs)], \@xcontig_sizes;

   # Examine largest cross-contig
   my $xcontig = (sort {$b->num_sequences <=> $a->num_sequences} $large_xcsp->assembly)[0];
   is $xcontig->num_sequences, 7;
   is $xcontig->get_seq_coord($xcontig->get_seq_by_name('species1635|5973'))->start, 1;
   is $xcontig->get_seq_coord($xcontig->get_seq_by_name('species158|7890'))->start, 1;
   is $xcontig->get_seq_coord($xcontig->get_seq_by_name('species2742|48'))->end, 140;

   # one contig at a time
   $in = Bio::Assembly::IO->new(
      -file => test_input_file('contigspectrumtest.tigr'),
      -format => 'tigr'
   );
   $sc = $in->next_assembly;
   ok $csp = Bio::Assembly::Tools::ContigSpectrum->new(
      -eff_asm_params => 1 ), 'one contig at a time';
   for my $contig ($sc->all_contigs) {
     ok $csp->assembly($contig);
   }

   is scalar $csp->assembly(), 5;
   is_deeply $csp->spectrum, {1=>0, 2=>3, 6=>1, 9=>1}; # [0 3 0 0 0 1 0 0 1]
   is $csp->eff_asm_params, 1;
   is $csp->max_size, 9;
   is $csp->nof_rep, 5;
   is $csp->nof_seq, 21;
   float_is $csp->avg_seq_len, 303.81;
   is $csp->nof_overlaps, 16;
   is $csp->min_overlap, 35;
   float_is $csp->avg_overlap, 155.875;
   float_is $csp->min_identity, 96.8421;
   float_is $csp->avg_identity, 98.8826;
}
