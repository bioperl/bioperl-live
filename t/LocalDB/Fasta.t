BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin( -tests => 46,
                -requires_modules => [qw(Bio::DB::Fasta Bio::SeqIO)]);
}
use strict;
use warnings;
use Bio::Root::Root;
use File::Copy;
my $DEBUG = test_debug();

my $test_dir  = setup_temp_dir('dbfa');
my $test_file = test_input_file('dbfa', 'mixed_alphabet.fasta');
my $test_files = [
    test_input_file('dbfa', 'mixed_alphabet.fasta'),
    test_input_file('dbfa', '6.fa')
];


{
    # Test basic functionalities
    ok my $db = Bio::DB::Fasta->new($test_dir, -reindex => 1), 'Index a directory';
    isa_ok $db, 'Bio::DB::Fasta';
    is $db->length('CEESC13F'), 389;
    is $db->seq('CEESC13F:1,10'), 'cttgcttgaa';
    is $db->seq('AW057119', 1, 10), 'tcatgttggc';
    is $db->header('AW057119'), 'AW057119 test description';
    is $db->seq('foobarbaz'), undef;
    is $db->get_Seq_by_id('foobarbaz'), undef;
    ok my $primary_seq = $db->get_Seq_by_id('AW057119');
    isa_ok $primary_seq, 'Bio::PrimarySeqI';
    is $primary_seq->trunc(11, 20)->length, 10;
    is $primary_seq->trunc(11, 20)->seq, 'ttctcggggt';
    is $primary_seq->description, 'test description', 'bug 3126';
    is $primary_seq->seq, 'tcatgttggcttctcggggtttttatggattaatacattttccaaacgattctttgcgccttctgtggtgccgccttctccgaaggaactgacgaaaaatgacgtggatttgctgacaaatccaggcgaggaatatttggacggattgatgaaatggcacggcgacgagcgacccgtgttcaaaagagaggacatttatcgttggtcggatagttttccagaatatcggctaagaatgatttgtctgaaagacacgacaagggtcattgcagtcggtcaatattgttactttgatgctctgaaagaaaggagagcagccattgttcttcttaggattgggatggacggatcctgaatatcgtaatcgggcagttatggagcttcaagcttcgatggcgctggaggagagggatcggtatccgactgccaacgcggcatcgcatccaaataagttcatgaaacgattttggcacatattcaacggcctcaaagagcacgaggacaaaggtcacaaggctgccgctgtttcatacaagagcttctacgacctcanagacatgatcattcctgaaaatctggatgtcagtggtattactgtaaatgatgcacgaaaggtgccacaaagagatataatcaactacgatcaaacatttcatccatatcatcgagaaatggttataatttctcacatgtatgacaatgatgggtttggaaaagtgcgtatgatgaggatggaaatgtacttggaattgtctagcgatgtctttanaccaacaagactgcacattagtcaattatgcagatagcc';

}


{
    # Test tied hash access
    my %h;
    ok tie(%h, 'Bio::DB::Fasta', $test_dir), 'Tied hash access';
    ok exists $h{'AW057146'};
    is $h{'AW057146:1,10'}, 'aatgtgtaca';
    is $h{'AW057146:10,1'}, 'tgtacacatt'; # reverse complement
}


{
    # Test writing the Bio::PrimarySeq::Fasta objects with SeqIO
    ok my $db = Bio::DB::Fasta->new($test_dir, -reindex => 1);
    my $out = Bio::SeqIO->new(
        -format => 'genbank',
        -file   => '>'.test_output_file()
    );
    my $primary_seq = Bio::Seq->new(-primary_seq => $db->get_Seq_by_acc('AW057119'));
    eval {
        $out->write_seq($primary_seq)
    };
    is $@, '';

    $out = Bio::SeqIO->new(-format => 'embl', -file  => '>'.test_output_file());
    eval {
        $out->write_seq($primary_seq)
    };
    is $@, '';
}


{
    # Test alphabet
    ok my $db = Bio::DB::Fasta->new( $test_file, -reindex => 1), 'Index a single file';
    is $db->alphabet('gi|352962132|ref|NG_030353.1|'), 'dna';
    is $db->alphabet('gi|352962148|ref|NM_001251825.1|'), 'rna';
    is $db->alphabet('gi|194473622|ref|NP_001123975.1|'), 'protein';
    is $db->alphabet('gi|61679760|pdb|1Y4P|B'), 'protein';
    is $db->alphabet('123'), '';
}


{
    # Test stream
    ok my $db = Bio::DB::Fasta->new( $test_file, -reindex => 1);
    ok my $stream = $db->get_PrimarySeq_stream;
    ok $stream = $db->get_Seq_stream;
    isa_ok $stream, 'Bio::DB::Indexed::Stream';
    my $count = 0;
    while (my $seq = $stream->next_seq) {
        $count++;
    }
    is $count, 5;
    unlink "$test_file.index";
}


{
    # Test an arbitrary index filename and cleaning
    my $name = 'arbitrary.idx';
    ok my $db = Bio::DB::Fasta->new( $test_file,
        -reindex => 1, -index_name => $name, -clean => 1 );
    is $db->index_name, $name;
    ok -f $name;
    unlink $name;
    undef $db;
    ok ! -f $name;
}


{
    # Test opening set of files and test IDs
    ok my $db = Bio::DB::Fasta->new( $test_files, -reindex => 1), 'Index a set of files';
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
}


{
    # Squash warnings locally
    local $SIG{__WARN__} = sub {};

    # Issue 3172
    my $test_dir = setup_temp_dir('bad_dbfa');
    throws_ok {my $db = Bio::DB::Fasta->new($test_dir, -reindex => 1)}
        qr/FASTA header doesn't match/;

    # Issue 3237
    # Empty lines within a sequence is bad...
    throws_ok {my $db = Bio::DB::Fasta->new(test_input_file('badfasta.fa'), -reindex => 1)}
        qr/Blank lines can only precede header lines/;
}


{
    # again, Issue 3237
    # but empty lines preceding headers are okay, but let's check the seqs just in case
    my $db;
    lives_ok {$db = Bio::DB::Fasta->new(test_input_file('spaced_fasta.fa'), -reindex => 1)};
    is length($db->seq('CEESC39F')), 375, 'length is correct in sequences past spaces';
    is length($db->seq('CEESC13F')), 389;

    is $db->subseq('CEESC39F', 51, 60)  , 'acatatganc', 'subseq is correct';
    is $db->subseq('CEESC13F', 146, 155), 'ggctctccct', 'subseq is correct';

    # Remove temporary test file
    my $outfile = test_input_file('spaced_fasta.fa').'.index';
    unlink $outfile;
}

exit;



sub setup_temp_dir {
    # this obfuscation is to deal with lockfiles by GDBM_File which can
    # only be created on local filesystems apparently so will cause test
    # to block and then fail when the testdir is on an NFS mounted system

    my $data_dir = shift;

    my $io = Bio::Root::IO->new();
    my $tempdir = test_output_dir();
    my $test_dir = $io->catfile($tempdir, $data_dir);
    mkdir($test_dir); # make the directory
    my $indir = test_input_file($data_dir);
    opendir(my $INDIR,$indir) || die("cannot open dir $indir");
    # effectively do a cp -r but only copy the files that are in there, no subdirs
    for my $file ( map { $io->catfile($indir,$_) } readdir($INDIR) ) {
        next unless (-f $file );
        copy($file, $test_dir);
    }
    closedir($INDIR);
    return $test_dir
}
