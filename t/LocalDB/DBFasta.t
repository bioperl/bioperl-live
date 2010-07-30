# -*-Perl-*- Test Harness script for Bioperl
# $Id$


BEGIN {     
    use lib '.';
	use Bio::Root::Test;
	
    test_begin(-tests => 18,
	       -requires_modules => [qw(Bio::DB::Fasta Bio::SeqIO)]);
	use_ok('Bio::Root::IO');
	use_ok('File::Copy');
}

my $DEBUG = test_debug();

# this obfuscation is to deal with lockfiles by GDBM_File which can
# only be created on local filesystems apparently so will cause test
# to block and then fail when the testdir is on an NFS mounted system

my $io = Bio::Root::IO->new(-verbose => $DEBUG);
my $tempdir = test_output_dir();
my $test_dbdir = $io->catfile($tempdir, 'dbfa');
mkdir($test_dbdir); # make the directory
my $indir = test_input_file('dbfa');
opendir(my $INDIR,$indir) || die("cannot open dir $indir");
# effectively do a cp -r but only copy the files that are in there, no subdirs
for my $file ( map { $io->catfile($indir,$_) } readdir($INDIR) ) {
	next unless (-f $file );
	copy($file, $test_dbdir);
}
closedir($INDIR);

# now use this temporary dir for the db file
my $db = Bio::DB::Fasta->new($test_dbdir, -reindex => 1);
ok($db);
cmp_ok($db->length('CEESC13F'), '>', 0);
is(length $db->seq('CEESC13F:1,10'), 10);
is(length $db->seq('AW057119',1,10), 10);
my $primary_seq = $db->get_Seq_by_id('AW057119');
ok($primary_seq);
cmp_ok(length($primary_seq->seq), '>', 0);
is($primary_seq->trunc(1,10)->length, 10);
is($primary_seq->description, 'test description', 'bug 3126');
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

# test out writing the Bio::PrimarySeq::Fasta objects with SeqIO

$db = Bio::DB::Fasta->new($test_dbdir, -reindex => 1);
my $out = Bio::SeqIO->new(-format => 'genbank',
			  -file  => '>'.test_output_file());
$primary_seq = Bio::Seq->new(-primary_seq => $db->get_Seq_by_acc('AW057119'));
eval {
    #warn(ref($primary_seq),"\n");
    $out->write_seq($primary_seq) 
};
ok(!$@);

$out = Bio::SeqIO->new(-format => 'embl', -file  => '>'.test_output_file());

eval {
    $out->write_seq($primary_seq) 
};
ok(!$@);

