# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 37,
			   -requires_modules => [qw(IO::String LWP::UserAgent)]);
	
	use_ok('Bio::DB::CUTG');
	use_ok('Bio::CodonUsage::Table');
	use_ok('Bio::CodonUsage::IO');
    use_ok('Bio::SeqIO');
    use_ok('Bio::Tools::SeqStats');
}

my $outfile = test_output_file();
my $verbose = test_debug();

# try reading from file
ok my $io = Bio::CodonUsage::IO->new
  (-file=> test_input_file('MmCT'));
ok  my $cut2 = $io->next_data();
is int($cut2->aa_frequency('LEU')), 10;

# write
ok $io = Bio::CodonUsage::IO->new(-file => ">$outfile");
$io->write_data($cut2);
ok -e $outfile;

# can we read what we've written?
ok $io = Bio::CodonUsage::IO->new(-file => "$outfile");
ok $cut2 = $io->next_data();
is int($cut2->aa_frequency('LEU')), 10;

# now try making a user defined CUT from a sequence
ok my $seqobj = Bio::SeqIO->new (-file =>test_input_file('HUMBETGLOA.fa'),
				                        -format => 'fasta')->next_seq;
is $seqobj->subseq(10,20), 'TTGACACCACT';
ok my $codcont_Ref = Bio::Tools::SeqStats->count_codons($seqobj);
is $codcont_Ref->{'TGA'}, 16;
ok my $cut = Bio::CodonUsage::Table->new(-data=>$codcont_Ref);
is $cut->codon_rel_frequency('CTG'), 0.18;
is $cut->codon_abs_frequency('CTG'), 2.6;
is $cut->codon_count('CTG'), 26;
is $cut->get_coding_gc(1), "39.70";
ok my $ref = $cut->probable_codons(20);

# Examples:
# http://www.kazusa.or.jp/codon/cgi-bin/spsearch.cgi?species=Pan+troglodytes&c=s
# http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=37011&aa=1&style=GCG
# requiring Internet access, set env BIOPERLDEBUG to 1 to run
SKIP: {
	test_skip(-tests => 14, -requires_networking => 1);
	ok my $tool = Bio::WebAgent->new(-verbose => $verbose);
	ok $tool->sleep;
	is $tool->delay(1), 1;
	ok $tool->sleep;

	# get CUT from web
	ok my $db = Bio::DB::CUTG->new();
	$db->verbose($verbose ? $verbose : -1);
	my $cdtable;
	eval {$cdtable = $db->get_request(-sp =>'Pan troglodytes');};
	skip "Server/network problems? Skipping those tests\n$@", 9 if $@;
	
	# tests for Table.pm, the answers seem to change with time, so not specific
	cmp_ok($cdtable->cds_count(), '>', 10);
	cmp_ok(int($cdtable->aa_frequency('LEU')), '>', 1);
	ok $cdtable->get_coding_gc('all');
	cmp_ok($cdtable->codon_rel_frequency('ttc'), '<', 1); 
    
	## now lets enter a non-existent species ans check handling..
	## should default to human...
	my $db2 = Bio::DB::CUTG->new();
	$db2->verbose($verbose ? $verbose : -1);
	eval {$cut2 = $db2->get_request(-sp =>'Wookie magnus');};
	skip "Server/network problems? Skipping those tests\n$@", 5 if $@;
	is $cut2->species(), 'Homo sapiens';
	
	$db = Bio::DB::CUTG->new();
	$db->verbose($verbose ? $verbose : -1);
	eval {$cdtable = $db->get_request(-sp =>'Homo sapiens');};
	skip "Server/network problems? Skipping those tests\n$@", 4 if $@;
	
	# tests for Table.pm, the answers seem to change with time, so not specific
	cmp_ok($cdtable->cds_count(), '>', 10);
	cmp_ok(int($cdtable->aa_frequency('LEU')), '>', 1);
	ok $cdtable->get_coding_gc('all');
	cmp_ok($cdtable->codon_rel_frequency('ttc'), '<', 1); 

	my $db3 = Bio::DB::CUTG->new(-sp =>'Bacillus subtilis', -gc => 1);
	$db3->verbose($verbose ? $verbose : -1);
	my $cut3;
	eval {$cut3 = $db3->get_request();};
	skip "Server/network problems? Skipping those tests\n$@", 5 if $@;
	print $cut3->codon_rel_frequency('ATG');
}
