#-*-Perl-*-
# $Id$
use strict;

BEGIN {     
    use lib 't/lib';
	use BioperlTest;
	
	test_begin(-tests => 14,
			   -requires_modules => [qw(Bio::DB::Fasta)]);
	
	use_ok('Bio::Root::IO');
	use_ok('File::Copy');
}

my $DEBUG = test_debug();

# this obfuscation is to deal with lockfiles by GDBM_File which can
# only be created on local filesystems apparently so will cause test
# to block and then fail when the testdir is on an NFS mounted system

my $io = Bio::Root::IO->new(-verbose => $DEBUG);
my $tempdir = $io->tempdir('CLEANUP' => 1);
my $test_dbdir = $io->catfile($tempdir, 'dbfa');
mkdir($test_dbdir); # make the directory
my $indir = $io->catfile(qw(. t data dbfa)); 
opendir(INDIR,$indir) || die("cannot open dir $indir");
# effectively do a cp -r but only copy the files that are in there, no subdirs
for my $file ( map { $io->catfile($indir,$_) } readdir(INDIR) ) {
	next unless (-f $file );
	copy($file, $test_dbdir);
}
closedir(INDIR);

# now use this temporary dir for the db file
my $db = Bio::DB::Fasta->new($test_dbdir, -reindex => 1);
ok($db);
cmp_ok($db->length('CEESC13F'), '>', 0);
is(length $db->seq('CEESC13F:1,10'), 10);
is(length $db->seq('AW057119',1,10), 10);
my $primary_seq = $db->get_Seq_by_id('AW057119');
ok($primary_seq);
cmp_ok(length($primary_seq->seq), '>', 0);
ok(!defined $db->get_Seq_by_id('foobarbaz'));
undef $db;
undef $primary_seq;

my (%h,$dna1,$dna2);
ok(tie(%h,'Bio::DB::Fasta',$test_dbdir));
ok($h{'AW057146'});
ok($dna1 = $h{'AW057146:1,10'});
ok($dna2 = $h{'AW057146:10,1'});

my $revcom = reverse $dna1;
$revcom =~ tr/gatcGATC/ctagCTAG/;
is($dna2, $revcom);
