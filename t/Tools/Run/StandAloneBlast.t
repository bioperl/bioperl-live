# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;
use warnings;
use File::Spec;

BEGIN { 
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 45);
	
	use_ok('Bio::Tools::Run::StandAloneBlast');
	use_ok('Bio::SeqIO');
}

# Note: the swissprot and ecoli.nt data sets may be downloaded from
# ftp://ftp.ncbi.nih.gov/blast/db/FASTA
my $verbose = test_debug() || -1;
my $nt_database = 'ecoli.nt';
my $amino_database = 'swissprot';
my $evalue = 0.001;
my ($seq1,$seq2,$seq3,$seq4);

# Tests to check that "-attr" and "attr" and "a" all do the same thing
# http://bugzilla.open-bio.org/show_bug.cgi?id=1912
for my $p (qw(database db -d -database d)) {
  my $f = Bio::Tools::Run::StandAloneBlast->new($p => $nt_database);
  is $f->d(), $nt_database;
}
for my $p (qw(expect evalue -e -expect e)) {
  my $f = Bio::Tools::Run::StandAloneBlast->new($p => $evalue);
  is $f->e(), $evalue;
}

# NCBI blast params are case-sensitive, wublast aren't
my $ncbi_factory = Bio::Tools::Run::StandAloneBlast->new(-program => 'blastn');
isa_ok $ncbi_factory, 'Bio::Tools::Run::StandAloneBlast';
isa_ok $ncbi_factory, 'Bio::Tools::Run::StandAloneNCBIBlast';
my $wu_factory = Bio::Tools::Run::StandAloneBlast->new(-program => 'wublastn');
isa_ok $wu_factory, 'Bio::Tools::Run::StandAloneBlast';
isa_ok $wu_factory, 'Bio::Tools::Run::StandAloneWUBlast';
for my $p (qw(e E)) {
  $ncbi_factory->$p($p);
  $wu_factory->$p($p);
}
is $ncbi_factory->e, 'e';
is $ncbi_factory->E, 'E';
is $wu_factory->e, 'E';
is $wu_factory->E, 'E';

# blastall switches like -I should take boolean but return 'T' or 'F' once set
is $ncbi_factory->I, undef;
is $ncbi_factory->I(1), 'T';
is $ncbi_factory->I(0), 'F';
is $ncbi_factory->I('T'), 'T';
is $ncbi_factory->I('F'), 'F';

# We should be able to set -F "m D" in an intuitive way, and also by manually
# quoting the value ourselves
$ncbi_factory->F('m D');
my $param_string = $ncbi_factory->_setparams('blastall');
like $param_string, qr/-F ['"]m D['"]/;
$ncbi_factory->F('"m S"');
$param_string = $ncbi_factory->_setparams('blastall');
like $param_string, qr/-F ["']m S['"]/;
$ncbi_factory->F("'m D'");
$param_string = $ncbi_factory->_setparams('blastall');
like $param_string, qr/-F ['"]m D["']/;

# dashed parameters should work
my $outfile = test_output_file();
ok my $factory = Bio::Tools::Run::StandAloneBlast->new(-verbose     => $verbose,
		-program     => 'blastn',
		-database    => $nt_database , 
		-_READMETHOD => 'SearchIO', 
		-output      => $outfile,
		-verbose     => 0);
is $factory->database, $nt_database;

# Setup and then do tests that actually run blast

my @params = ('program'     => 'blastn',
			  'database'    => $nt_database , 
			  '_READMETHOD' => 'SearchIO', 
			  'output'      => $outfile,
			  'verbose'     => 0 );
ok $factory = Bio::Tools::Run::StandAloneBlast->new('-verbose' => $verbose, @params);

my $inputfilename = test_input_file('test.txt');

is $factory->quiet(0), 0;
is $factory->q(-3), -3;

SKIP: {
	skip 'blastall not installed, skipping tests', 12 unless $factory->executable('blastall');
	skip 'must have BLASTDIR, BLASTDB or BLASTDATADIR env variable set, skipping tests', 12 unless defined $Bio::Tools::Run::StandAloneBlast::DATADIR;
	
	my @testresults = qw(37 182 182  253 167 2);
	my $testcount = 0;
	
	# use ecoli.nt
	my $nt_database_file = File::Spec->catfile($Bio::Tools::Run::StandAloneBlast::DATADIR, $nt_database);
	like $nt_database_file, qr/$nt_database/;
	SKIP: {
		skip "Database $nt_database not found, skipping tests on it", 8 unless -e $nt_database_file;
		
		my $parser = $factory->blastall($inputfilename);
		my $blast_report = $parser->next_result;
		is $blast_report->num_hits, $testresults[$testcount++];
		
		$factory->_READMETHOD('blast_pull');  # Note required leading underscore in _READMETHOD
		my $str = Bio::SeqIO->new('-file' => test_input_file('dna2.fa'),
								  '-format' => 'fasta');
		$seq1 = $str->next_seq();
		$seq2 = $str->next_seq();
		
		my $pull_report = $factory->blastall($seq1);
		my $sbjct = $pull_report->next_result->next_hit;
		my $hsp = $sbjct->next_hsp;
		is $hsp->score, $testresults[$testcount];
		
		$factory->_READMETHOD('Blast');
		my $searchio_report = $factory->blastall($seq1);
		$sbjct = $searchio_report->next_result->next_hit;
		$hsp = $sbjct->next_hsp;
		is $hsp->score, $testresults[$testcount++];
		
		my @seq_array =($seq1,$seq2);
		my $seq_array_ref = \@seq_array;
		$factory->_READMETHOD('blast_pull');
		$pull_report = $factory->blastall($seq_array_ref);
		$sbjct = $pull_report->next_result->next_hit;
		$hsp = $sbjct->next_hsp;
		is $hsp->score, $testresults[$testcount];
		
		$factory->_READMETHOD('Blast');
		$searchio_report = $factory->blastall($seq_array_ref);
		$sbjct = $searchio_report->next_result->next_hit;
		$hsp = $sbjct->next_hsp;
		is $hsp->score, $testresults[$testcount++];
		
		ok $sbjct = $searchio_report->next_result->next_hit;
		$hsp = $sbjct->next_hsp;
		is $hsp->score, $testresults[$testcount++];
		
		@params = ('-verbose' => $verbose, 'program'  => 'blastp'); 
		$factory = Bio::Tools::Run::StandAloneBlast->new(@params);
		
		$str = Bio::SeqIO->new(-file => test_input_file('amino.fa'),
							   -format => 'Fasta' );
		$seq3 = $str->next_seq();
		$seq4 = $str->next_seq();
		
		$factory->_READMETHOD('Blast');
		my $bl2seq_report = $factory->bl2seq($seq3, $seq4);
		$hsp = $bl2seq_report->next_result->next_hit->next_hsp;
		is $hsp->hit->start, $testresults[$testcount++], "creating/parsing SearchIO bl2seq report object";
	}
	
	# use swissprot
	my $amino_database_file = File::Spec->catfile($Bio::Tools::Run::StandAloneBlast::DATADIR, $amino_database);
	SKIP: {
		skip "Database $amino_database not found, skipping tests on it", 3 unless -e $amino_database_file;
		
		@params = ('database' => $amino_database, '-verbose' => $verbose);
		$factory = Bio::Tools::Run::StandAloneBlast->new(@params);
		
		my $iter = 2;
		$factory->j($iter);    # 'j' is blast parameter for # of iterations
		my $new_iter = $factory->j();
		is $new_iter, 2, "set blast parameter";
		
		my $blast_report = $factory->blastpgp($seq3)->next_result;
		is $blast_report->number_of_iterations, $testresults[$testcount];
		
		$factory->_READMETHOD('blast_pull');
		$iter = 2;
		$factory->j($iter);    # 'j' is blast parameter for # of iterations
		$new_iter = $factory->j();
		is $new_iter, $iter, "set blast parameter";
		
	}
}
