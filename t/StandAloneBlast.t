## Bioperl Test Harness Script for Modules
##


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------

## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..8\n";
	use vars qw($loaded); }
# Edit the following line to point to the location of your local blast files...
# BEGIN {$ENV{BLASTDIR} = '/home/peter/blast/'; }
END {print "not ok 1\n" unless $loaded;}

use Bio::Tools::Blast;
use Bio::Tools::BPlite;
use Bio::Tools::StandAloneBlast;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Seq;
use strict;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 


sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

my ($blast_report, $hsp, @testresults);

# If you have different BLAST databases installed the following line will need to
# be changed to refer to those databases
my $nt_database = 'ecoli.nt';
my $amino_database = 'swissprot';

## Create blast factory
my @params = ('program' => 'blastn', 'database' => $nt_database , '_READMETHOD' => 'Blast');
my  $factory = Bio::Tools::StandAloneBlast->new(@params);

test 2, $factory, " couldn't create blast factory";

my $inputfilename = 't/test.txt';
my $program = 'blastn';

# If the blast program isn't found and executable at the expected location,
# there is no point in executing the remaining tests...

my $blast_present = Bio::Tools::StandAloneBlast->exists_blast();
unless ($blast_present) {
	warn "blast program not found. Skipping tests 3 to 8\n";
   	print "ok 3\n";  print "ok 4\n"; print "ok 5\n";
   	print "ok 6\n";  print "ok 7\n"; print "ok 8\n";
	exit 0;
}

# If the ecoli.nt and swissprot databases are not installed (but other ones presumably are),
# the tests only check that blast reports are being generated without checking
# their contents...


if ($nt_database eq 'ecoli.nt') {	
	$testresults[3] = '$blast_report->num_hits == 1' ;
	$testresults[4] = '$hsp->score == 182';
	$testresults[5] = '$hsp->score == 182';
} else {
	$testresults[3] =  '$blast_report->num_hits';
	$testresults[4]  =  '$hsp->score';
	$testresults[5]  =  '$hsp->score';
}
if ($nt_database eq 'swissprot') {	
	$testresults[8]  =  '$blast_report->number_of_iterations == 2';
} else {
	$testresults[8] =  '$blast_report->number_of_iterations';

}
 $blast_report = $factory->blastall($inputfilename);
test 3,$testresults[3] , " failed blastall blastn test";

# This time use a BioSeq object as the query and BPlite as _READMETHOD
$factory->_READMETHOD('BPlite');    # Note required leading underscore in _READMETHOD
$factory->outfile('BP1test.out');

# Get dna sequences from file
my $str = Bio::SeqIO->new(-file=>'t/dna2.fa' , '-format' => 'Fasta', );
my $seq1 = $str->next_seq();
my $seq2 = $str->next_seq();

my $BPlite_report = $factory->blastall($seq1);
my $sbjct = $BPlite_report->nextSbjct;
 $hsp = $sbjct->nextHSP;
test 4, $testresults[4] , " failed BPlite nt test";

# This time use a BioSeq object array ref as the query

my @seq_array =($seq1,$seq2);
my $seq_array_ref = \@seq_array;

my $BPlite_report2 = $factory->blastall(\@seq_array);
 $sbjct = $BPlite_report2->nextSbjct;
 $hsp = $sbjct->nextHSP;
test 5,$testresults[5], " failed Seq array input test";

# Bl2seq testing
# first create factory for bl2seq
@params = ('program' => 'blastp', 'outfile' => 'bl2seqtest.out');
$factory = Bio::Tools::StandAloneBlast->new(@params);

# Get protein sequences from file
$str = Bio::SeqIO->new(-file=>'t/amino.fa' , '-format' => 'Fasta', );
my $seq3 = $str->next_seq();
my $seq4 = $str->next_seq();


my $bl2seq_report = $factory->bl2seq($seq3, $seq4);
test 6, $bl2seq_report->subject->start == 167, " failed creating or parsing bl2seq report object";


# Psiblast testing
## Create factory for psiblast
@params = ('database' => $amino_database, 'outfile' => 'psiblast.out');
$factory = Bio::Tools::StandAloneBlast->new(@params);

# set psiblast iteration parameter

my $iter = 2;
$factory->j($iter);    # 'j' is blast parameter for # of iterations
my $new_iter = $factory->j();

test 7, $new_iter == 2, " failed setting blast parameter";


$blast_report = $factory->blastpgp($seq3);
test 8,$testresults[8] , " failed creating or parsing psiblast report object ";
