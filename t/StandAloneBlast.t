# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#

use strict;
use constant NUMTESTS => 10;
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

ok(1);
my $verbose = -1;
my ($blast_report, $hsp, @testresults);

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
    exit();
}
my $nt_database_file = Bio::Root::IO->catfile($Bio::Tools::Run::StandAloneBlast::DATADIR, $nt_database);
ok($nt_database_file, qr/$nt_database/);
my $amino_database_file = Bio::Root::IO->catfile($Bio::Tools::Run::StandAloneBlast::DATADIR, $amino_database);
my $file_present = -e $nt_database_file;
my $exit;
unless ($file_present) {
   warn "Blast Database $nt_database not found";
   $exit=1;
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
if ($nt_database eq 'ecoli.nt') {	
	$testresults[3] = '$blast_report->num_hits == 1' ;
	$testresults[4] = '$hsp->score == 182';
	$testresults[5] = '$hsp->score == 182';
} else {
	$testresults[3] =  '$blast_report->num_hits';
	$testresults[4]  =  '$hsp->score';
	$testresults[5]  =  '$hsp->score';
}
if ($amino_database eq 'swissprot') {	
	$testresults[8]  =  '$blast_report->number_of_iterations == 2';
} else {
	$testresults[8] =  '$blast_report->number_of_iterations';

}
 $blast_report = $factory->blastall($inputfilename);
ok $testresults[3];

$factory->_READMETHOD('BPlite');    # Note required leading underscore in _READMETHOD

my $str = Bio::SeqIO->new(-file=>Bio::Root::IO->catfile("t","data","dna2.fa") , '-format' => 'Fasta', );
my $seq1 = $str->next_seq();
my $seq2 = $str->next_seq();

my $BPlite_report = $factory->blastall($seq1);
my $sbjct = $BPlite_report->nextSbjct;
 $hsp = $sbjct->nextHSP;
ok $testresults[4];


my @seq_array =($seq1,$seq2);
my $seq_array_ref = \@seq_array;

my $BPlite_report2 = $factory->blastall(\@seq_array);
 $sbjct = $BPlite_report2->nextSbjct;
 $hsp = $sbjct->nextHSP;
ok $testresults[5];

@params = ('-verbose' => $verbose,
	   'program' => 'blastp'); # This used to be blastp but atleast on my implementation it should be T
$factory = Bio::Tools::Run::StandAloneBlast->new(@params);

$str = Bio::SeqIO->new(-file=>Bio::Root::IO->catfile("t","data","amino.fa") , '-format' => 'Fasta', );
my $seq3 = $str->next_seq();
my $seq4 = $str->next_seq();
my $bl2seq_report = $factory->bl2seq($seq3, $seq4);
$hsp = $bl2seq_report->next_feature;
ok $hsp->hit()->start, 167, " failed creating or parsing bl2seq report object";


@params = ('database' => $amino_database);
$factory = Bio::Tools::Run::StandAloneBlast->new(@params);


my $iter = 2;
$factory->j($iter);    # 'j' is blast parameter for # of iterations
my $new_iter = $factory->j();

ok $new_iter, 2, " failed setting blast parameter";

$blast_report = $factory->blastpgp($seq3);
ok $testresults[8];
