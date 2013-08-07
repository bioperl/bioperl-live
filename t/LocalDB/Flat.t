# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
	use Bio::Root::Test;

	test_begin(-tests => 25,
                -requires_module => 'DB_File');

	use_ok('Bio::DB::Flat');
}

my $verbose = test_debug();

# First of all we need to create an flat db

my $tmpdir = test_output_dir();

my $db = Bio::DB::Flat->new(-directory  => $tmpdir,
                            -index      => 'bdb',
									 -dbname     => 'mydb',
									 -format     => 'fasta',
									 -verbose    => $verbose,
									 -write_flag => 1 );
ok($db);
my $dir = test_input_file('AAC12660.fa');
my $result = $db->build_index(glob($dir));
ok($result);

# Now let's get the sequence out again
my $seq = $db->get_Seq_by_id('AAC12660');
ok($seq);
is($seq->length,504);
undef $db;

$db = Bio::DB::Flat->new(-directory  => $tmpdir,
                         -index      => 'bdb',
                         -format     => 'embl',
						 -dbname     => 'myembl',
                         -verbose    => $verbose,
                         -write_flag => 1 );

$dir= test_input_file('cds_sample.embl');

$result = $db->build_index(glob($dir));

is ($db->get_all_primary_ids, 1);
#is ($db->get_all_accs, 1);
ok($result);
$seq = $db->get_Seq_by_id('EAL24309');
ok($seq);
is($seq->length,192);

# deal with wantarray conditions
$seq = $db->get_Seq_by_acc('CH236947.1');
ok($seq && ref($seq));
is($seq->length,192);


undef $db;

$db = Bio::DB::Flat->new(-directory  => $tmpdir,
			 -index      => 'binarysearch',
			 -format     => 'fasta',
			 -dbname     => 'mybinfa',
			 -verbose    => $verbose,
			 -write_flag => 1
			 );

$dir= test_input_file('dbfa', '1.fa');
$result = $db->build_index($dir);
ok($result);
$seq = $db->get_Seq_by_id('AW057119');
ok($seq);
is($seq->length,808);
undef $db;

$db = Bio::DB::Flat->new(-directory  => $tmpdir,
			 -index      => 'binarysearch',
			 -format     => 'swiss',
			 -dbname     => 'mybinswiss',
			 -verbose    => $verbose,
			 -write_flag => 1
			 );
$dir= test_input_file('swiss.dat');
$result = $db->build_index($dir);

ok($result);
$seq = $db->get_Seq_by_id('ACON_CAEEL');
ok($seq);
is($seq->length,788);

$seq = $db->get_Seq_by_id('ACON_CAEEL');
ok($seq && ref($seq));

undef $db;

$db = Bio::DB::Flat->new(-directory  => $tmpdir,
                         -index      => 'binarysearch',
                         -format     => 'fasta',
								 -dbname     => 'myfasta',
                         -verbose    => $verbose,
                         -write_flag => 1 );

$dir = test_input_file('tmp.fst');
$result = $db->build_index(glob($dir));
ok($result);
$seq = $db->get_Seq_by_id('TEST00004');
is($seq->length,98);

undef $db;

$db = Bio::DB::Flat->new(-directory  => $tmpdir,
                         -index      => 'bdb',
                         -format     => 'fasta',
								 -dbname     => 'mybfasta',
                         -verbose    => $verbose,
                         -write_flag => 1 );

$dir = test_input_file('tmp.fst');
$result = $db->build_index(glob($dir));
ok($result);
for my $id ( qw(TEST00001 TEST00002 TEST00003 TEST00004) ) {
	$seq = $db->get_Seq_by_id($id);
	is($seq->length,98);
}
