# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$ 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $DEBUG);

BEGIN {
	$NUMTESTS = 34;
	$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
	
	eval {require Test::More;};
	if ($@) {
		use lib 't';
	}
	use Test::More;
	
	eval {
		require IO::String; 
		require LWP::UserAgent;
	};
	if ($@) {
		plan skip_all => 'IO::String or LWP::UserAgent not installed. This means that the module is not usable. Skipping tests';
	}
	else {
		plan tests => $NUMTESTS;
	}
	
	use_ok('Bio::DB::CUTG');
	use_ok('Bio::CodonUsage::Table');
	use_ok('Bio::CodonUsage::IO');
    use_ok('Bio::SeqIO');
    use_ok('Bio::Tools::SeqStats');
}

my $outfile = "t/data/cutg.out";
my $verbose = 1 if $DEBUG;

# try reading from file
ok my $io = Bio::CodonUsage::IO->new
  (-file=> Bio::Root::IO->catfile("t", "data", "MmCT"));
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
ok my $seqobj = Bio::SeqIO->new (-file =>
			 Bio::Root::IO->catfile("t", "data", "HUMBETGLOA.fa"),
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

# requiring Internet access, set env BIOPERLDEBUG to 1 to run
SKIP: {
	skip "Skipping tests which require remote servers, set BIOPERLDEBUG=1 to test", 11 unless $DEBUG;
	ok my $tool = Bio::WebAgent->new(-verbose =>$verbose);
	ok $tool->sleep;
	is $tool->delay(1), 1;
	ok $tool->sleep;

	# get CUT from web
	ok my $db = Bio::DB::CUTG->new();
	ok $db->verbose(1);
	my $cdtable;
	eval {$cdtable = $db->get_request(-sp =>'Pan troglodytes');};
	skip "Could not connect to server, server/network problems? Skipping those tests", 5 if $@;
	
	# tests for Table.pm
	is $cdtable->cds_count(), 617; # new value at CUD
	is int($cdtable->aa_frequency('LEU')), 10;
	ok $cdtable->get_coding_gc('all');
	is $cdtable->codon_rel_frequency('ttc'), "0.61"; 
    
	## now lets enter a non-existent species ans check handling..
	## should default to human...
	my $db2 = Bio::DB::CUTG->new();
	eval {$cut2 = $db2->get_request(-sp =>'Wookie magnus');};
	skip "Could not connect to server, server/network problems? Skipping those tests", 1 if $@;
	is $cut2->species(), 'Homo sapiens';
}
