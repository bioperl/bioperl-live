# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 45,
	       -requires_module => 'Graph');

    use_ok('Bio::FeatureIO');
}

my $io;
my $f;
my $s;
my $fcount;
my $scount;

################################################################################
#
# use FeatureIO::gff to read a FASTA file.
#
$fcount = 0;
$scount = 0;

ok( $io = Bio::FeatureIO->new( -file => test_input_file('dna1.fa') ) );

#read features
while($f = $io->next_feature()){
    warn $f;
    $fcount++;
}
is($fcount, 0);

#then try to read sequences again.  should get seqs now
while($s = $io->next_seq()){
    $scount++;
    TODO: {
	local $TODO = 'How did this ever work?!?';
	if ($scount == 1) {
	    is($s->seq, 'Test1');
	}
    }
}
is($scount,  1);

################################################################################
#
# use FeatureIO::gff to read a GFF3 file.
#
$fcount = 0;
$scount = 0;

ok( $io = Bio::FeatureIO->new( -file => test_input_file('knownGene.gff3') ) );

#try to read sequences first.  should be undef
while($s = $io->next_seq()){
    $scount++;
}
is($scount,0);

#then read features
while($f = $io->next_feature()){
    $fcount++;
}
is($fcount, 15);

#then try to read sequences again.  should still be undef
while($s = $io->next_seq()){
    $scount++;
}
is($scount,0);

################################################################################
#
# use FeatureIO::gff to read a GFF3 file w/ directivized FASTA tail
#
$fcount = 0;
$scount = 0;

ok( $io = Bio::FeatureIO->new( -file => test_input_file('hybrid1.gff3') ) );

#try to read sequences first.  should be undef
while($s = $io->next_seq()){
    $scount++;
}
is($scount , 0);

#then read features
while($f = $io->next_feature()){
    $fcount++;
}
is($fcount , 6);

#then try to read sequences again.
while($s = $io->next_seq()){
    $scount++;
    TODO: {
	local $TODO = 'How did this ever work?!?';
	if ($scount == 1) {
	    is($s->seq, 'Test1');
	}
    }
}
is($scount , 1);

################################################################################
#
# use FeatureIO::gff to read a GFF3 file w/ non-directivized FASTA tail
#
$fcount = 0;
$scount = 0;

ok( $io = Bio::FeatureIO->new( -file => test_input_file('hybrid2.gff3') ) );

#try to read sequences first.  should be undef
while($s = $io->next_seq()){
    $scount++;
}
is($scount , 0);

#then read features
while($f = $io->next_feature()){
    $fcount++;
}
is($fcount , 6);

################################################################################
#
# use FeatureIO::gff to read a GFF3 file of directives
#
$fcount = 0;
$scount = 0;

ok( $io = Bio::FeatureIO->new(-file => test_input_file('directives.gff3'),
			      -verbose => test_debug() ? test_debug() : -1));

#read features
while($f = $io->next_feature()){
    $fcount++;
}
is($fcount , 1); #sequence-region

################################################################################
#
# use FeatureIO::gff to read a GFF3 file as aggregated feature groups
#
$fcount = 0;
$scount = 0;

ok( $io = Bio::FeatureIO->new( -file => test_input_file('hybrid1.gff3') ) );

#try to read sequences first.  should be undef
while($s = $io->next_seq()){
    $scount++;
}
is($scount , 0);

#read feature groups
$f = $io->next_feature_group();
is($f , 1);
$f = $io->next_feature_group();
is($f , 0);

#then try to read sequences again.
while($s = $io->next_seq()){
    $scount++;
    TODO: {
	local $TODO = 'How did this ever work?!?';
	if ($scount == 1) {
	    is($s->seq, 'Test1');
	}
    }
}
is($scount , 1);


################################################################################
#
# use FeatureIO::bed to read a bed file
#
ok($io = Bio::FeatureIO->new(-file => test_input_file('1.bed')));

ok($f = $io->next_feature);
# Check correct conversion of [0, feature-end+1) bed-coordinates into [1, feature-end]
# bioperl coordinates.  (here: bed [0, 10))
is($f->start, 1);
is($f->end, 10);

# Check field values.
my @tags = $f->get_tag_values("Name");
is(scalar(@tags), 1);
is($tags[0], "test-coordinates-1");

is($f->seq_id, "chr1");





################################################################################
#
# use FeatureIO::gff to read a PTT file.
#
$fcount = 0;

my $ptt_in = Bio::FeatureIO->new(
    -file => test_input_file('test.ptt'),
    -format => 'ptt',
    );
ok($ptt_in);

while (my $f = $ptt_in->next_feature) {
    $fcount++;
    if ($fcount==2) {
	# 2491..3423  + 310 24217063  metF  LB002 - COG0685E  5,10-methylenetetrahydrofolate reductase
	is( $f->start , 2491 );
	is( $f->end , 3423 );
	cmp_ok( $f->strand, '>', 0 );
	is( ($f->get_tag_values('PID'))[0],'24217063' );
	is( ($f->get_tag_values('Gene'))[0], 'metF' );
	is( ($f->get_tag_values('Synonym'))[0], 'LB002' );
	ok( ! $f->has_tag('Code') );
	is( ($f->get_tag_values('COG'))[0],'COG0685E' );
	is( ($f->get_tag_values('Product'))[0], '5,10-methylenetetrahydrofolate reductase' );
    }
}
is($fcount , 367);

################################################################################
#
# use FeatureIO::vecscreen_simple to read a vecscreen file
#


{
    my @expected_features =
	(
	 {
	     'seq_id' => 'C02HBa0072A04.1',
	     'primary_tag' => 'moderate_match',
	     'end' => '60548',
	     'start' => '60522'
	 },
	 {
	     'seq_id' => 'SL_FOS91h17_SP6_0',
	     'primary_tag' => 'strong_match',
	     'end' => '122',
	     'start' => '60'
	 },
	 {
	     'seq_id' => 'SL_FOS91h18_T7_0',
	     'primary_tag' => 'strong_match',
	     'end' => '102',
	     'start' => '35'
	 },
	 {
	     'seq_id' => 'SL_FOS91h18_T7_0',
	     'primary_tag' => 'moderate_match',
	     'end' => '103',
	     'start' => '76'
	 },
	 {
	     'seq_id' => 'SL_FOS91h18_T7_0',
	     'primary_tag' => 'weak_match',
	     'end' => '104',
	     'start' => '82'
	 },
	 {
	     'seq_id' => 'SL_FOS91h18_T7_0',
	     'primary_tag' => 'suspect_origin',
	     'end' => '34',
	     'start' => '1'
	 },
	 {
	     'seq_id' => 'SL_FOS91i01_SP6_0',
	     'primary_tag' => 'strong_match',
	     'end' => '110',
	     'start' => '46'
	 },
	 {
	     'seq_id' => 'SL_FOS91i01_SP6_0',
	     'primary_tag' => 'suspect_origin',
	     'end' => '45',
	     'start' => '1'
	 },
	 {
	     'seq_id' => 'SL_FOS92b12_T7_0',
	     'primary_tag' => 'strong_match',
	     'end' => '108',
	     'start' => '41'
	 },
	 {
	     'seq_id' => 'SL_FOS92b12_T7_0',
	     'primary_tag' => 'moderate_match',
	     'end' => '109',
	     'start' => '82'
	 },
	 {
	     'seq_id' => 'SL_FOS92b12_T7_0',
	     'primary_tag' => 'weak_match',
	     'end' => '110',
	     'start' => '88'
	 },
	 {
	     'seq_id' => 'SL_FOS92b12_T7_0',
	     'primary_tag' => 'weak_match',
	     'end' => '1329',
	     'start' => '1313'
	 },
	 {
	     'seq_id' => 'SL_FOS92b12_T7_0',
	     'primary_tag' => 'suspect_origin',
	     'end' => '40',
	     'start' => '1'
	 },
	 {
	     'seq_id' => 'SL_FOS92b12_T7_0',
	     'primary_tag' => 'suspect_origin',
	     'end' => '1334',
	     'start' => '1330'
	 }
	);
    my @vs_features;
    my $vs_in = Bio::FeatureIO->new( -file => test_input_file('vecscreen_simple.test_output'),
				     -format => 'vecscreen_simple',
	);
    ok( $vs_in );
    while(my $feat = $vs_in->next_feature) {
	push @vs_features,$feat;
    }

    #convert the array of feature objects to something that can more easily be checked with is_deeply
    @vs_features = map {
	my $f = $_;
	my $rec = { map {$_ => $f->$_()} qw/start end primary_tag seq_id/ };
    } @vs_features;

    is_deeply(\@vs_features,\@expected_features,'vecscreen_simple gets the correct features');
}
