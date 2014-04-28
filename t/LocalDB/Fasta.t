BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin( -tests => 107,
                -requires_modules => [qw(Bio::DB::Fasta Bio::SeqIO)] );
}
use strict;
use warnings;
use Bio::Root::Root;
use File::Copy;
my $DEBUG = test_debug();


# Test Bio::DB::Fasta, but also the underlying module, Bio::DB::IndexedBase

my $test_dir  = setup_temp_dir('dbfa');
my $test_file = test_input_file('dbfa', 'mixed_alphabet.fasta');
my $test_files = [
    test_input_file('dbfa', 'mixed_alphabet.fasta'),
    test_input_file('dbfa', '6.fa')
];


{
    # Test basic functionalities
    ok my $db = Bio::DB::Fasta->new($test_dir, -reindex => 1), 'Index a directory';
    is $db->glob, '*.{fa,FA,fasta,FASTA,fast,FAST,dna,DNA,fna,FNA,faa,FAA,fsa,FSA}';
    isa_ok $db, 'Bio::DB::Fasta';
    is $db->length('CEESC13F'), 389;
    is $db->seq('CEESC13F:1,10'), 'cttgcttgaa';
    is $db->seq('CEESC13F:1-10'), 'cttgcttgaa';
    is $db->seq('CEESC13F:1..10'), 'cttgcttgaa';
    is $db->seq('CEESC13F:1..10/1'), 'cttgcttgaa';
    is $db->seq('CEESC13F:1..10/+1'), 'cttgcttgaa';
    is $db->seq('CEESC13F:1..10/-1'), 'ttcaagcaag';
    is $db->seq('CEESC13F/1'), 'cttgcttgaaaaatttatataaatatttaagagaagaaaaataaataatcgcatctaatgacgtctgtccttgtatccctggtttccattgactggtgcactttcctgtctttgaggacatggacaatattcggcatcagttcctggctctccctcctctcctggtgctccagcagaaccgttctctccattatctcccttgtctccacgtggtccacgctctcctggtgctcctggaataccttgagctccctcgtgccgaattcctgcagcccgggggatccactagttctagagcggccgccaccgcggtgggagctccagcttttgttncctttagtgagggttaatttcgagcttggcgtaatcatggtcatagctgtttcctg';
    is $db->seq('CEESC13F/-1'), 'caggaaacagctatgaccatgattacgccaagctcgaaattaaccctcactaaaggnaacaaaagctggagctcccaccgcggtggcggccgctctagaactagtggatcccccgggctgcaggaattcggcacgagggagctcaaggtattccaggagcaccaggagagcgtggaccacgtggagacaagggagataatggagagaacggttctgctggagcaccaggagaggagggagagccaggaactgatgccgaatattgtccatgtcctcaaagacaggaaagtgcaccagtcaatggaaaccagggatacaaggacagacgtcattagatgcgattatttatttttcttctcttaaatatttatataaatttttcaagcaag';
    is $db->seq('AW057119', 1, 10), 'tcatgttggc';
    is $db->seq('AW057119', 1, 10, 1), 'tcatgttggc';
    is $db->seq('AW057119', 1, 10, -1), 'gccaacatga';
    is $db->seq('AW057119', 10, 1), 'gccaacatga';
    is $db->seq('AW057119', 10, 1, -1), 'tcatgttggc';
    is $db->header('AW057119'), 'AW057119 test description';
    is $db->seq('foobarbaz'), undef;
    is $db->get_Seq_by_id('foobarbaz'), undef;
    is $db->file('AW057119'), '1.fa';
    is $db->file('AW057410'), '3.fa';
    is $db->file('CEESC13F'), '6.fa';

    # Bio::DB::RandomAccessI and Bio::DB::SeqI methods
    ok my $primary_seq = $db->get_Seq_by_id('AW057119');
    ok $primary_seq = $db->get_Seq_by_acc('AW057119');
    ok $primary_seq = $db->get_Seq_by_version('AW057119');
    ok $primary_seq = $db->get_Seq_by_primary_id('AW057119');
    isa_ok $primary_seq, 'Bio::PrimarySeq::Fasta';
    isa_ok $primary_seq, 'Bio::PrimarySeqI';

    # Bio::PrimarySeqI methods
    is $primary_seq->id, 'AW057119';
    is $primary_seq->display_id, 'AW057119';
    like $primary_seq->primary_id, qr/^Bio::PrimarySeq::Fasta=HASH/;
    is $primary_seq->alphabet, 'dna';
    is $primary_seq->accession_number, 'unknown';
    is $primary_seq->is_circular, undef;
    is $primary_seq->subseq(11, 20), 'ttctcggggt';
    is $primary_seq->description, 'test description', 'bug 3126';
    is $primary_seq->seq, 'tcatgttggcttctcggggtttttatggattaatacattttccaaacgattctttgcgccttctgtggtgccgccttctccgaaggaactgacgaaaaatgacgtggatttgctgacaaatccaggcgaggaatatttggacggattgatgaaatggcacggcgacgagcgacccgtgttcaaaagagaggacatttatcgttggtcggatagttttccagaatatcggctaagaatgatttgtctgaaagacacgacaagggtcattgcagtcggtcaatattgttactttgatgctctgaaagaaaggagagcagccattgttcttcttaggattgggatggacggatcctgaatatcgtaatcgggcagttatggagcttcaagcttcgatggcgctggaggagagggatcggtatccgactgccaacgcggcatcgcatccaaataagttcatgaaacgattttggcacatattcaacggcctcaaagagcacgaggacaaaggtcacaaggctgccgctgtttcatacaagagcttctacgacctcanagacatgatcattcctgaaaatctggatgtcagtggtattactgtaaatgatgcacgaaaggtgccacaaagagatataatcaactacgatcaaacatttcatccatatcatcgagaaatggttataatttctcacatgtatgacaatgatgggtttggaaaagtgcgtatgatgaggatggaaatgtacttggaattgtctagcgatgtctttanaccaacaagactgcacattagtcaattatgcagatagcc';
    ok my $trunc = $primary_seq->trunc(11,20);
    isa_ok $trunc, 'Bio::PrimarySeq::Fasta';
    isa_ok $trunc, 'Bio::PrimarySeqI';
    is $trunc->length, 10;
    is $trunc->seq, 'ttctcggggt';
    ok my $rev = $trunc->revcom;
    isa_ok $rev, 'Bio::PrimarySeq::Fasta';
    isa_ok $rev, 'Bio::PrimarySeqI';
    is $rev->seq, 'accccgagaa';
    is $rev->length, 10;
}


{
    # Re-open an existing index.
    # Doing this test properly involves unloading and reloading Bio::DB::Fasta.

    SKIP: {
        test_skip(-tests => 1, -requires_modules => [qw(Class::Unload)]);
        use_ok('Class::Unload');
        Class::Unload->unload( 'Bio::DB::Fasta' );
        Class::Unload->unload( 'Bio::DB::IndexedBase' );
        require Bio::DB::Fasta;
    }

    ok my $db = Bio::DB::Fasta->new($test_dir), 'Re-open an existing index';
    is $db->seq('AW057119', 1, 10), 'tcatgttggc';
}


{
    # Test tied hash access
    my %h;
    ok tie(%h, 'Bio::DB::Fasta', $test_dir), 'Tied hash access';
    ok exists $h{'AW057146'};
    is $h{'AW057146:1,10'} , 'aatgtgtaca'; # in file 1.fa
    is $h{'AW057146:10,1'} , 'tgtacacatt'; # reverse complement
    is $h{'AW057443:11,20'}, 'gaaccgtcag'; # in file 4.fa
}


{
    # Test writing the Bio::PrimarySeq::Fasta objects with SeqIO
    ok my $db = Bio::DB::Fasta->new($test_dir, -reindex => 1), 'Writing with SeqIO';
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
    # Test alphabet and reverse-complement RNA
    ok my $db = Bio::DB::Fasta->new( $test_file, -reindex => 1), 'Index a single file';
    is $db->alphabet('gi|352962132|ref|NG_030353.1|'), 'dna';
    is $db->alphabet('gi|352962148|ref|NM_001251825.1|'), 'rna';
    is $db->alphabet('gi|194473622|ref|NP_001123975.1|'), 'protein';
    is $db->alphabet('gi|61679760|pdb|1Y4P|B'), 'protein';
    is $db->alphabet('123'), '';
    is $db->seq('gi|352962148|ref|NM_001251825.1|', 20, 29,  1), 'GUCAGCGUCC';
    is $db->seq('gi|352962148|ref|NM_001251825.1|', 20, 29, -1), 'GGACGCUGAC';

    # Test empty sequence
    is $db->seq('123'), '';

    is $db->file('gi|352962132|ref|NG_030353.1|'), 'mixed_alphabet.fasta';
}


{
    # Test stream
    ok my $db = Bio::DB::Fasta->new( $test_file, -reindex => 1);
    ok my $stream = $db->get_PrimarySeq_stream;
    isa_ok $stream, 'Bio::DB::Indexed::Stream';
    my $count = 0;
    while (my $seq = $stream->next_seq) {
        $count++;
    }
    is $count, 5;

    # ActivePerl will not allow deletion if the tie-hash is still active
    $db->DESTROY;
    # Strawberry Perl temporary file
    unlink "$test_file.index" if -e "$test_file.index";
    # ActivePerl temporary files
    unlink "$test_file.index.dir" if -e "$test_file.index.dir";
    unlink "$test_file.index.pag" if -e "$test_file.index.pag";
}


{
    # Concurrent databases (bug #3390)
    ok my $db1 = Bio::DB::Fasta->new( test_input_file('dbfa', '1.fa') );
    ok my $db3 = Bio::DB::Fasta->new( test_input_file('dbfa', '3.fa') );
    ok my $db4 = Bio::DB::Fasta->new( $test_dir );
    ok my $db2 = Bio::DB::Fasta->new( test_input_file('dbfa', '2.fa') );
    is $db4->file('AW057231'), '1.fa';
    is $db2->file('AW057302'), '2.fa';
    is $db4->file('AW057119'), '1.fa';
    is $db3->file('AW057336'), '3.fa';
    is $db1->file('AW057231'), '1.fa';
    is $db4->file('AW057410'), '3.fa';

    # ActivePerl will not allow deletion if the tie-hash is still active
    $db1->DESTROY;
    $db2->DESTROY;
    $db3->DESTROY;
    # Strawberry Perl temporary file
    unlink $db1->index_name if -e $db1->index_name;
    unlink $db2->index_name if -e $db2->index_name;
    unlink $db3->index_name if -e $db3->index_name;
    # ActivePerl temporary files
    unlink $db1->index_name().'.dir' if -e $db1->index_name().'.dir';
    unlink $db2->index_name().'.dir' if -e $db2->index_name().'.dir';
    unlink $db3->index_name().'.dir' if -e $db3->index_name().'.dir';
    unlink $db1->index_name().'.pag' if -e $db1->index_name().'.pag';
    unlink $db2->index_name().'.pag' if -e $db2->index_name().'.pag';
    unlink $db3->index_name().'.pag' if -e $db3->index_name().'.pag';
}


{
    # Test an arbitrary index filename and cleaning
    my $name = 'arbitrary.idx';
    ok my $db = Bio::DB::Fasta->new( $test_file,
        -reindex => 1, -index_name => $name, -clean => 1,
    );
    is $db->index_name, $name;

    # Tied-hash in Strawberry Perl produce $name,
    # while in ActivePerl produce "$name.dir" and "$name.pag"
    if (-e "$name.pag") {
        ok -f "$name.pag";
        # ActivePerl will not allow deletion if the tie-hash is still active
        $db->DESTROY;
        unlink "$name.dir" if -e "$name.dir";
        unlink "$name.pag" if -e "$name.pag";
        ok ! -f "$name.pag";
    }
    else {
        ok -f $name;
        # ActivePerl will not allow deletion if the tie-hash is still active
        $db->DESTROY;
        unlink $name if -e $name;
        ok ! -f $name;
    }
    undef $db;
}


{
    # Test makeid
    ok my $db = Bio::DB::Fasta->new( $test_file,
        -reindex => 1, -clean => 1, -makeid => \&extract_gi,
    ), 'Make single ID';
    is_deeply [sort $db->get_all_primary_ids], ['', 194473622, 352962132, 352962148, 61679760];
    is $db->get_Seq_by_id('gi|352962148|ref|NM_001251825.1|'), undef;
    isa_ok $db->get_Seq_by_id(194473622), 'Bio::PrimarySeqI';
}


{
    # Test makeid that generates several IDs, bug #3389
    ok my $db = Bio::DB::Fasta->new( $test_file,
        -reindex => 1, -clean => 1, -makeid => \&extract_gi_and_ref,
    ), 'Make multiple IDs, bug #3389';
    is_deeply [sort $db->get_all_primary_ids], ['', 194473622, 352962132, 352962148, 61679760, 'NG_030353.1',  'NM_001251825.1', 'NP_001123975.1'];
    is $db->get_Seq_by_id('gi|352962148|ref|NM_001251825.1|'), undef;
    isa_ok $db->get_Seq_by_id('NG_030353.1'), 'Bio::PrimarySeqI';
}


{
    # Test opening set of files and test IDs
    ok my $db = Bio::DB::Fasta->new( $test_files, -reindex => 1), 'Index a set of files';
    ok $db->ids;
    ok $db->get_all_ids;
    my @ids = sort $db->get_all_primary_ids();
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

    my $index = $db->index_name;
    # ActivePerl will not allow deletion if the tie-hash is still active
    $db->DESTROY;
    # Strawberry Perl temporary file
    unlink $index if -e $index;
    # ActivePerl temporary files
    unlink "$index.dir" if -e "$index.dir";
    unlink "$index.pag" if -e "$index.pag";
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
    # Issue 3237 again
    # but empty lines preceding headers are okay, but let's check the seqs just in case
    my $db;
    lives_ok {$db = Bio::DB::Fasta->new(test_input_file('spaced_fasta.fa'), -reindex => 1)};
    is length($db->seq('CEESC39F')), 375, 'length is correct in sequences past spaces';
    is length($db->seq('CEESC13F')), 389;

    is $db->subseq('CEESC39F', 51, 60)  , 'acatatganc', 'subseq is correct';
    is $db->subseq('CEESC13F', 146, 155), 'ggctctccct', 'subseq is correct';

    # Remove temporary test file
    my $outfile = test_input_file('spaced_fasta.fa').'.index';

    # ActivePerl will not allow deletion if the tie-hash is still active
    $db->DESTROY;
    # Strawberry Perl temporary file
    unlink $outfile if -e $outfile;
    # ActivePerl temporary files
    unlink "$outfile.dir" if -e "$outfile.dir";
    unlink "$outfile.pag" if -e "$outfile.pag";
}

exit;


sub extract_gi {
    # Extract GI from RefSeq
    my $header = shift;
    my ($id) = ($header =~ /gi\|(\d+)/m);
    return $id || '';
}


sub extract_gi_and_ref {
    # Extract GI and from RefSeq
    my $header = shift;
    my ($gi)  = ($header =~ /gi\|(\d+)/m);
    $gi ||= '';
    my ($ref) = ($header =~ /ref\|([^|]+)/m);
    $ref ||= '';
    return $gi, $ref;
}


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
