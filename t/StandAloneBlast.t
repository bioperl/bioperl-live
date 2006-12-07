# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#

use strict;
use constant NUMTESTS => 33;
BEGIN { 
	eval { require Test; };
	if ( $@ ) {
		use lib 't';
	}
	use Test;
	plan tests => NUMTESTS; 
}

END { 
	foreach( $Test::ntest..NUMTESTS) {
		skip('Blast or env variables not installed correctly', 1);
	}
	unlink('blastreport.out') if -e 'blastreport.out';
}

use Bio::Tools::BPlite;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Seq;
use Bio::Root::IO;
use Bio::SearchIO;

# Note: the swissprot and ecoli.nt data sets may be downloaded from
# ftp://ftp.ncbi.nih.gov/blast/db/FASTA
my $verbose = -1;
my $nt_database = 'ecoli.nt';
my $amino_database = 'swissprot';
my $evalue = 0.001;
my ($seq1,$seq2,$seq3,$seq4);

# Tests to check that "-attr" and "attr" and "a" all do the same thing
# http://bugzilla.open-bio.org/show_bug.cgi?id=1912

for my $p (qw(database db -d -database d)) {
  my $f = Bio::Tools::Run::StandAloneBlast->new($p => $nt_database);
  ok $f->d() eq $nt_database;
}
for my $p (qw(expect evalue -e -expect e)) {
  my $f = Bio::Tools::Run::StandAloneBlast->new($p => $evalue);
  ok $f->e() eq $evalue;
}

# Let's continue...

my @params = ('program'     => 'blastn',
				   'database'    => $nt_database , 
				   '_READMETHOD' => 'SearchIO', 
				   'output'      => 'blastreport.out',
				   'verbose'     => 0 );
my  $factory = Bio::Tools::Run::StandAloneBlast->new('-verbose' => $verbose,
						     @params);
ok $factory;

my $inputfilename = Bio::Root::IO->catfile("t","data","test.txt");

$factory->quiet(0);
$factory->q(-3);

ok($factory->q, -3);
ok($factory->quiet, 0);

if( ! $factory->executable('blastall') ) {
    skip('blastall not installed',1);
    exit;
}

if( ! defined $Bio::Tools::Run::StandAloneBlast::DATADIR ) {
    print STDERR "must have BLASTDIR and BLASTDB or BLASTDATADIR env variable set\n";
    exit(0);
}

my $nt_database_file = 
    Bio::Root::IO->catfile($Bio::Tools::Run::StandAloneBlast::DATADIR, 
			   $nt_database);
ok($nt_database_file, qr/$nt_database/);

my @testresults = qw(37 182 182  253 167 2);
my $testcount = 0;

# use ecoli.nt
if (-e $nt_database_file) {
	my $parser = $factory->blastall($inputfilename);
	my $blast_report = $parser->next_result;
	ok($blast_report->num_hits,$testresults[$testcount++]);

	$factory->_READMETHOD('BPlite');  # Note required leading underscore in _READMETHOD
	my $str = Bio::SeqIO->new('-file'  => Bio::Root::IO->catfile(qw(t data dna2.fa)),
									  '-format' => 'fasta');
	$seq1 = $str->next_seq();
	$seq2 = $str->next_seq();

	my $BPlite_report = $factory->blastall($seq1);
	my $sbjct = $BPlite_report->nextSbjct;
	my $hsp = $sbjct->nextHSP;
	ok($hsp->score, $testresults[$testcount]);

	$factory->_READMETHOD('Blast');
	my $searchio_report = $factory->blastall($seq1);
	$sbjct = $searchio_report->next_result->next_hit;
	$hsp = $sbjct->next_hsp;
	ok($hsp->score, $testresults[$testcount++]);

	my @seq_array =($seq1,$seq2);
	my $seq_array_ref = \@seq_array;
	$factory->_READMETHOD('BPlite');
	my $BPlite_report2 = $factory->blastall($seq_array_ref);
	$sbjct = $BPlite_report2->nextSbjct;
	$hsp = $sbjct->nextHSP;
	ok($hsp->score, $testresults[$testcount]);

	$factory->_READMETHOD('Blast');
	$searchio_report = $factory->blastall($seq_array_ref);
	$sbjct = $searchio_report->next_result->next_hit;
	$hsp = $sbjct->next_hsp;
	ok($hsp->score, $testresults[$testcount++]);

	$sbjct = $searchio_report->next_result->next_hit;
	ok($sbjct);
	$hsp = $sbjct->next_hsp;
	ok($hsp->score, $testresults[$testcount++]);

	@params = ('-verbose' => $verbose,
				  'program'  => 'blastp'); 
	$factory = Bio::Tools::Run::StandAloneBlast->new(@params);

	$str = Bio::SeqIO->new(-file => Bio::Root::IO->catfile(qw(t data amino.fa)),
								  -format => 'Fasta' );
	$seq3 = $str->next_seq();
	$seq4 = $str->next_seq();
	$factory->_READMETHOD('BPlite');
	my $bl2seq_report = $factory->bl2seq($seq3, $seq4);
	$hsp = $bl2seq_report->next_feature;
	ok ($hsp->hit->start, $testresults[$testcount], 
		 " failed creating or parsing BPlite bl2seq report object");

	$factory->_READMETHOD('Blast');
	$bl2seq_report = $factory->bl2seq($seq3, $seq4);
	$hsp = $bl2seq_report->next_result->next_hit->next_hsp;
	ok( $hsp->hit->start, $testresults[$testcount++], 
		 " failed creating or parsing SearchIO bl2seq report object");
} else {
	for (1..14) {
		skip("Database $nt_database not found, skipping",1);
	}	
}

# use nr
my $amino_database_file = 
    Bio::Root::IO->catfile($Bio::Tools::Run::StandAloneBlast::DATADIR, 
			   $amino_database);

if (-e $amino_database) {
	@params = ('database' => $amino_database,
				  '-verbose' => $verbose);
	$factory = Bio::Tools::Run::StandAloneBlast->new(@params);

	my $iter = 2;
	$factory->j($iter);    # 'j' is blast parameter for # of iterations
	my $new_iter = $factory->j();
	
	ok $new_iter, 2, " failed setting blast parameter";
	my $blast_report = $factory->blastpgp($seq3)->next_result;
	ok($blast_report->number_of_iterations, $testresults[$testcount]);

	$factory->_READMETHOD('BPlite');
	$iter = 2;
	$factory->j($iter);    # 'j' is blast parameter for # of iterations
	$new_iter = $factory->j();

	ok($new_iter, $iter, " failed setting blast parameter");

	$blast_report = $factory->blastpgp($seq3);
	ok($blast_report->number_of_iterations, $testresults[$testcount]);
} else {
	for (1..4) {
		skip("Database $amino_database not found, skipping",1);
	}
}

# "dashed parameters"
$factory = Bio::Tools::Run::StandAloneBlast->new(
												  -verbose     => $verbose,
												  -program     => 'blastn',
												  -database    => $nt_database , 
												  -_READMETHOD => 'SearchIO', 
												  -output      => 'blastreport.out',
												  -verbose     => 0	 );
ok(defined $factory);
