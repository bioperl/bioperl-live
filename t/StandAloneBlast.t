# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS);
BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;

    $NUMTESTS = 8; 
    plan tests => 8; 
}

# Edit the following line to point to the location of your local blast files...
# BEGIN {$ENV{BLASTDIR} = '/home/peter/blast/'; }
END { unlink('blastreport.out') }

use Bio::Tools::Blast;
use Bio::Tools::BPlite;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Seq;

ok(1);

my ($blast_report, $hsp, @testresults);

# If you have different BLAST databases installed the following line
# will need to be changed to refer to those databases

my $nt_database = 'ecoli.nt';
my $amino_database = 'swissprot';

## Create blast factory
my @params = ('program' => 'blastn', 'database' => $nt_database , 
	      '_READMETHOD' => 'Blast', 'output' => 'blastreport.out');
my  $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

ok $factory;

my $inputfilename = 't/test.txt';
my $program = 'blastn';

# If the blast program isn't found and executable at the expected location,
# there is no point in executing the remaining tests...

my $blast_present = Bio::Tools::Run::StandAloneBlast->exists_blast();
unless ($blast_present) {
    warn "blast program not found. Skipping tests $Test::ntest to $NUMTESTS\n";
    foreach ($Test::ntest..$NUMTESTS) {
	skip(1,1);
    }
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
ok $testresults[3];

# This time use a BioSeq object as the query and BPlite as _READMETHOD
$factory->_READMETHOD('BPlite');    # Note required leading underscore in _READMETHOD
#$factory->outfile('BP1test.out');

# Get dna sequences from file
my $str = Bio::SeqIO->new(-file=>'t/dna2.fa' , '-format' => 'Fasta', );
my $seq1 = $str->next_seq();
my $seq2 = $str->next_seq();

my $BPlite_report = $factory->blastall($seq1);
my $sbjct = $BPlite_report->nextSbjct;
 $hsp = $sbjct->nextHSP;
ok $testresults[4];

# This time use a BioSeq object array ref as the query

my @seq_array =($seq1,$seq2);
my $seq_array_ref = \@seq_array;

my $BPlite_report2 = $factory->blastall(\@seq_array);
 $sbjct = $BPlite_report2->nextSbjct;
 $hsp = $sbjct->nextHSP;
ok $testresults[5];

# Bl2seq testing
# first create factory for bl2seq
@params = ('program' => 'blastp');
#@params = ('program' => 'blastp', 'outfile' => 'bl2seqtest.out');
$factory = Bio::Tools::Run::StandAloneBlast->new(@params);

# Get protein sequences from file
$str = Bio::SeqIO->new(-file=>'t/amino.fa' , '-format' => 'Fasta', );
my $seq3 = $str->next_seq();
my $seq4 = $str->next_seq();

my $bl2seq_report = $factory->bl2seq($seq3, $seq4);
ok $bl2seq_report->subject->start, 167, " failed creating or parsing bl2seq report object";


# Psiblast testing
## Create factory for psiblast
@params = ('database' => $amino_database);
#@params = ('database' => $amino_database, 'outfile' => 'psiblast.out');
$factory = Bio::Tools::Run::StandAloneBlast->new(@params);

# set psiblast iteration parameter

my $iter = 2;
$factory->j($iter);    # 'j' is blast parameter for # of iterations
my $new_iter = $factory->j();

ok $new_iter, 2, " failed setting blast parameter";

$blast_report = $factory->blastpgp($seq3);
ok $testresults[8];
