# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#

use strict;
use constant NUMTESTS => 16;
BEGIN { 
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    plan tests => NUMTESTS; 
}

END { 
    foreach( $Test::ntest..NUMTESTS) {
	skip('Blast or env variables not installed correctly',1);
    }
    unlink('blastreport.out');
}

use Bio::Tools::BPlite;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Seq;
use Bio::Root::IO;
use Bio::SearchIO;

my $verbose = -1;

my $nt_database = 'ecoli.nt';
my $amino_database = 'swissprot';

my @params = ('program' => 'blastn', 'database' => $nt_database , 
	      '_READMETHOD' => 'SearchIO', 
	      'output' => 'blastreport.out');
my  $factory = Bio::Tools::Run::StandAloneBlast->new('-verbose' => $verbose,
						     @params);

ok $factory;

my $inputfilename = Bio::Root::IO->catfile("t","data","test.txt");
my $program = 'blastn';

my $blast_present = $factory->executable('blastall');
if( ! $blast_present ) {
    skip('Blast not installed',1);
    exit;
} else { 
    ok($blast_present);
}
if( ! defined $Bio::Tools::Run::StandAloneBlast::DATADIR ) {
    print STDERR "must have BLASTDIR and BLASTDB or BLASTDATADIR env variable set\n";
    exit(0);
}
my $nt_database_file = 
    Bio::Root::IO->catfile($Bio::Tools::Run::StandAloneBlast::DATADIR, 
			   $nt_database);
ok($nt_database_file, qr/$nt_database/);
my $amino_database_file = 
    Bio::Root::IO->catfile($Bio::Tools::Run::StandAloneBlast::DATADIR, 
			   $amino_database);
my $file_present = -e $nt_database_file;

my $exit;
unless ($file_present) {
   warn "Blast Database $nt_database not found";
   $exit = 1;
}
my $file_present2 = -e $amino_database_file;

unless ($file_present2) {
    warn "Blast Database $amino_database not found";
    $exit=1;
}

if ($exit) {
   warn"Blast databases(s) not found, skipping remaining  tests";
   exit(0);
}

my @testresults = qw(37 182 182  253 167 2);

my $testcount = 0;
my $parser = $factory->blastall($inputfilename);
my $blast_report = $parser->next_result;
ok($blast_report->num_hits,$testresults[$testcount++]);

$factory->_READMETHOD('BPlite');  # Note required leading underscore in _READMETHOD
my $str = Bio::SeqIO->new('-file'  => Bio::Root::IO->catfile(qw(t data
								dna2.fa)),
			  '-format' => 'fasta');
my $seq1 = $str->next_seq();
my $seq2 = $str->next_seq();

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
	   'program' => 'blastp'); # This used to be blastp but atleast on my implementation it should be T
$factory = Bio::Tools::Run::StandAloneBlast->new(@params);

$str = Bio::SeqIO->new(-file=>Bio::Root::IO->catfile(qw(t data amino.fa)),
		       '-format' => 'Fasta' );
my $seq3 = $str->next_seq();
my $seq4 = $str->next_seq();
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

@params = ('database' => $amino_database,
	   '-verbose' => $verbose);
$factory = Bio::Tools::Run::StandAloneBlast->new(@params);

my $iter = 2;
$factory->j($iter);    # 'j' is blast parameter for # of iterations
my $new_iter = $factory->j();

ok $new_iter, 2, " failed setting blast parameter";
$blast_report = $factory->blastpgp($seq3)->next_result;
ok($blast_report->number_of_iterations, $testresults[$testcount]);

$factory->_READMETHOD('BPlite');
$iter = 2;
$factory->j($iter);    # 'j' is blast parameter for # of iterations
$new_iter = $factory->j();

ok($new_iter, $iter, " failed setting blast parameter");

$blast_report = $factory->blastpgp($seq3);
ok($blast_report->number_of_iterations, $testresults[$testcount]);
