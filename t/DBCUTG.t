# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$ 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $DEBUG $ERROR $verbose $outfile);

$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
BEGIN {
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test; };
	$ERROR = 0;
	if( $@ ) {
		use lib 't';
	}
	use Test;

	$NUMTESTS = 32;
	plan tests => $NUMTESTS;

	eval {
		require IO::String; 
		require LWP::UserAgent;
	}; 
	if( $@ ) {
		warn("IO::String or LWP::UserAgent not installed. This means that the module is not usable. Skipping tests");
		$ERROR = 1;
	}
}

END {
	unlink $outfile if (-e $outfile);
	foreach ( $Test::ntest..$NUMTESTS) {
		skip('unable to complete DBCUTG tests',1);
	}
}

exit 0 if $ERROR ==  1;

require Bio::DB::CUTG;
require Bio::CodonUsage::Table;
require Bio::CodonUsage::IO;
require Bio::SeqIO;
require Bio::Tools::SeqStats;
ok 1;

$outfile = "t/data/cutg.out";
$verbose = 1 if $DEBUG;

# try reading from file
ok my $io = Bio::CodonUsage::IO->new
  (-file=> Bio::Root::IO->catfile("t", "data", "MmCT"));
ok  my $cut2 = $io->next_data();
ok int($cut2->aa_frequency('LEU')), 10;

# write
ok $io = Bio::CodonUsage::IO->new(-file => ">$outfile");
$io->write_data($cut2);
ok -e $outfile;

# can we read what we've written?
ok $io = Bio::CodonUsage::IO->new(-file => "$outfile");
ok $cut2 = $io->next_data();
ok int($cut2->aa_frequency('LEU')), 10;

# now try making a user defined CUT from a sequence
ok my $seqobj = Bio::SeqIO->new (-file =>
			 Bio::Root::IO->catfile("t", "data", "HUMBETGLOA.fa"),
				                        -format => 'fasta')->next_seq;
ok $seqobj->subseq(10,20), 'TTGACACCACT';
ok my $codcont_Ref = Bio::Tools::SeqStats->count_codons($seqobj);
ok $codcont_Ref->{'TGA'}, 16;
ok my $cut = Bio::CodonUsage::Table->new(-data=>$codcont_Ref);
ok $cut->codon_rel_frequency('CTG'), 0.18;
ok $cut->codon_abs_frequency('CTG'), 2.6;
ok $cut->codon_count('CTG'), 26;
ok $cut->get_coding_gc(1), "39.70";
ok my $ref = $cut->probable_codons(20);

# requiring Internet access, set env BIOPERLDEBUG to 1 to run
if( $DEBUG ) {
	ok my $tool = Bio::WebAgent->new(-verbose =>$verbose);
	ok $tool->sleep;
	ok $tool->delay(1), 1;
	ok $tool->sleep;

	# get CUT from web
	ok my $db = Bio::DB::CUTG->new();
	ok $db->verbose(1);
	my $cdtable;
	eval {
		$cdtable = $db->get_request(-sp =>'Pan troglodytes');
	};
	if ($@) {
		foreach ( $Test::ntest..$NUMTESTS) { 
			skip('Could not connect to server, skipping tests requiring remote servers',1);
		}
		exit(0);
    }
	
	# tests for Table.pm
	ok $cdtable->cds_count(), 617; # new value at CUD
	ok int($cdtable->aa_frequency('LEU')), 10;
	ok $cdtable->get_coding_gc('all');
	ok $cdtable->codon_rel_frequency('ttc'), "0.61"; 
    
	## now lets enter a non-existent species ans check handling..
	## should default to human...
	my $db2 =  Bio::DB::CUTG->new();
	ok $cut2 = $db2->get_request(-sp =>'Wookie magnus');
	ok $cut2->species(), 'Homo sapiens';
	ok 1 ;
} else { 
   for ( $Test::ntest..$NUMTESTS) {
		skip("Skipping tests which require remote servers - set env variable BIOPERLDEBUG to test",1);
	}
}
