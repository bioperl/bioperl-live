

use strict;
use vars qw($NUMTESTS);
BEGIN { 
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;

    $NUMTESTS = 8; 
    plan tests => 8; 
}

END { unlink('blastreport.out') }

use Bio::Tools::Blast;
use Bio::Tools::BPlite;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Seq;
use Bio::Root::IO;

ok(1);

my ($blast_report, $hsp, @testresults);


my $nt_database = 'ecoli.nt';
my $amino_database = 'swissprot';

my @params = ('program' => 'blastn', 'database' => $nt_database , 
	      '_READMETHOD' => 'Blast', 'output' => 'blastreport.out');
my  $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

ok $factory;

my $inputfilename = Bio::Root::IO->catfile("t","test.txt");
my $program = 'blastn';


my $blast_present = Bio::Tools::Run::StandAloneBlast->exists_blast();
unless ($blast_present) {
    warn "blast program not found. Skipping tests $Test::ntest to $NUMTESTS\n";
    foreach ($Test::ntest..$NUMTESTS) {
	skip(1,1);
    }
    exit 0;
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
if ($nt_database eq 'swissprot') {	
	$testresults[8]  =  '$blast_report->number_of_iterations == 2';
} else {
	$testresults[8] =  '$blast_report->number_of_iterations';

}
 $blast_report = $factory->blastall($inputfilename);
ok $testresults[3];

$factory->_READMETHOD('BPlite');    # Note required leading underscore in _READMETHOD

my $str = Bio::SeqIO->new(-file=>Bio::Root::IO->catfile("t","dna2.fa") , '-format' => 'Fasta', );
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

@params = ('program' => 'blastp');
$factory = Bio::Tools::Run::StandAloneBlast->new(@params);

$str = Bio::SeqIO->new(-file=>Bio::Root::IO->catfile("t","amino.fa") , '-format' => 'Fasta', );
my $seq3 = $str->next_seq();
my $seq4 = $str->next_seq();

my $bl2seq_report = $factory->bl2seq($seq3, $seq4);
ok $bl2seq_report->subject->start, 167, " failed creating or parsing bl2seq report object";


@params = ('database' => $amino_database);
$factory = Bio::Tools::Run::StandAloneBlast->new(@params);


my $iter = 2;
$factory->j($iter);    # 'j' is blast parameter for # of iterations
my $new_iter = $factory->j();

ok $new_iter, 2, " failed setting blast parameter";

$blast_report = $factory->blastpgp($seq3);
ok $testresults[8];
