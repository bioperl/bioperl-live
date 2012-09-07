BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin( -tests => 43,
                -requires_modules => [qw(Bio::DB::Fasta Bio::SeqIO)]);
}
use strict;
use warnings;
use Bio::Root::Root;
use File::Copy;
my $DEBUG = test_debug();

{

my $test_dbdir = setup_temp_dir('dbfa');

# now use this temporary dir for the db file
my ($db, $primary_seq);
ok $db = Bio::DB::Fasta->new($test_dbdir, -reindex => 1), 'Index a directory';
isa_ok $db, 'Bio::DB::Fasta';
cmp_ok $db->length('CEESC13F'), '>', 0;
is length($db->seq('CEESC13F:1,10')), 10;
is length($db->seq('AW057119',1,10)), 10;
ok $primary_seq = $db->get_Seq_by_id('AW057119');
isa_ok $primary_seq, 'Bio::PrimarySeqI';
cmp_ok length($primary_seq->seq), '>', 0;
is $primary_seq->trunc(1,10)->length, 10;
is $primary_seq->description, 'test description', 'bug 3126';
is $db->get_Seq_by_id('foobarbaz'), undef;
undef $db;
undef $primary_seq;

my (%h, $dna1, $dna2);
ok tie(%h, 'Bio::DB::Fasta', $test_dbdir);
ok $h{'AW057146'};
ok $dna1 = $h{'AW057146:1,10'};
ok $dna2 = $h{'AW057146:10,1'};

my $revcom = reverse $dna1;
$revcom =~ tr/gatcGATC/ctagCTAG/;
is $dna2, $revcom;

# test out writing the Bio::PrimarySeq::Fasta objects with SeqIO

ok $db = Bio::DB::Fasta->new($test_dbdir, -reindex => 1);
my $out = Bio::SeqIO->new(
    -format => 'genbank',
    -file   => '>'.test_output_file()
);
$primary_seq = Bio::Seq->new(-primary_seq => $db->get_Seq_by_acc('AW057119'));
eval {
    $out->write_seq($primary_seq)
};
ok !$@;

$out = Bio::SeqIO->new(-format => 'embl', -file  => '>'.test_output_file());

eval {
    $out->write_seq($primary_seq)
};
is $@, '';


# Test alphabet
my $test_file = test_input_file('dbfa', 'mixed_alphabet.fasta');

ok $db = Bio::DB::Fasta->new( $test_file, -reindex => 1), 'Index a single file';
is $db->alphabet('gi|352962132|ref|NG_030353.1|'), 'dna';
is $db->alphabet('gi|352962148|ref|NM_001251825.1|'), 'rna';
is $db->alphabet('gi|194473622|ref|NP_001123975.1|'), 'protein';
is $db->alphabet('gi|61679760|pdb|1Y4P|B'), 'protein';
is $db->alphabet('123'), '';


# Test stream
ok my $stream = $db->get_PrimarySeq_stream;
ok $stream = $db->get_Seq_stream;
isa_ok $stream, 'Bio::DB::Indexed::Stream';
my $count = 0;
while (my $seq = $stream->next_seq) {
    $count++;
}
is $count, 5;
unlink "$test_file.index";


# Test an arbitrary index filename and cleaning
my $name = 'arbitrary.idx';
ok $db = Bio::DB::Fasta->new( $test_file,
    -reindex => 1, -index_name => $name, -clean => 1 );
is $db->index_name, $name;
ok -f $name;
unlink $name;
undef $db;
ok ! -f $name;


# Test opening set of files and test IDs
my $test_files = [
    test_input_file('dbfa', 'mixed_alphabet.fasta'),
    test_input_file('dbfa', '6.fa')
];
ok $db = Bio::DB::Fasta->new( $test_files, -reindex => 1), 'Index a set of files';
my @ids = sort $db->get_all_ids();
is_deeply \@ids, [ qw(
    123
    CEESC12R
    CEESC13F
    CEESC13R
    CEESC14F
    CEESC14R
    CEESC15F
    CEESC15R
    CEESC15RB
    CEESC16F
    CEESC17F
    CEESC17RB
    CEESC18F
    CEESC18R
    CEESC19F
    CEESC19R
    CEESC20F
    CEESC21F
    CEESC21R
    CEESC22F
    CEESC23F
    CEESC24F
    CEESC25F
    CEESC26F
    CEESC27F
    CEESC28F
    CEESC29F
    CEESC30F
    CEESC32F
    CEESC33F
    CEESC33R
    CEESC34F
    CEESC35R
    CEESC36F
    CEESC37F
    CEESC39F
    CEESC40R
    CEESC41F
    gi|194473622|ref|NP_001123975.1|
    gi|352962132|ref|NG_030353.1|
    gi|352962148|ref|NM_001251825.1|
    gi|61679760|pdb|1Y4P|B
)];
like $db->index_name, qr/^fileset_.+\.index$/;
unlink $db->index_name;


{
    # squash warnings locally
    local $SIG{__WARN__} = sub {};

    # Issue 3172
    $test_dbdir = setup_temp_dir('bad_dbfa');
    throws_ok {$db = Bio::DB::Fasta->new($test_dbdir, -reindex => 1)}
        qr/FASTA header doesn't match/;

    # Issue 3237

    # Empty lines within a sequence is bad...
    throws_ok {$db = Bio::DB::Fasta->new(test_input_file('badfasta.fa'), -reindex => 1)}
        qr/Blank lines can only precede header lines/;
}

# again, Issue 3237

# but empty lines preceding headers are okay, but let's check the seqs just in case
lives_ok {$db = Bio::DB::Fasta->new(test_input_file('spaced_fasta.fa'), -reindex => 1)};
is length($db->seq('CEESC39F')), 375, 'length is correct in sequences past spaces';
is length($db->seq('CEESC13F')), 389;

is $db->subseq('CEESC39F', 51, 60), 'acatatganc', 'subseq is correct';
is $db->subseq('CEESC13F', 146, 155), 'ggctctccct', 'subseq is correct';

# Remove temporary test file
my $outfile = test_input_file('spaced_fasta.fa').'.index';
unlink $outfile;

exit;

}


sub setup_temp_dir {
    # this obfuscation is to deal with lockfiles by GDBM_File which can
    # only be created on local filesystems apparently so will cause test
    # to block and then fail when the testdir is on an NFS mounted system

    my $data_dir = shift;

    my $io = Bio::Root::IO->new();
    my $tempdir = test_output_dir();
    my $test_dbdir = $io->catfile($tempdir, $data_dir);
    mkdir($test_dbdir); # make the directory
    my $indir = test_input_file($data_dir);
    opendir(my $INDIR,$indir) || die("cannot open dir $indir");
    # effectively do a cp -r but only copy the files that are in there, no subdirs
    for my $file ( map { $io->catfile($indir,$_) } readdir($INDIR) ) {
        next unless (-f $file );
        copy($file, $test_dbdir);
    }
    closedir($INDIR);
    return $test_dbdir
}
