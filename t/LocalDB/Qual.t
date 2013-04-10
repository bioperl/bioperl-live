BEGIN {     
    use lib '.';
    use Bio::Root::Test;

    test_begin( -tests => 56,
                -requires_module => 'Bio::DB::Qual');

    use_ok('Bio::Root::IO');
    use_ok('File::Copy');
}

my $DEBUG = test_debug();

{

my $test_dbdir = setup_temp_dir('dbqual');

# now use this temporary dir for the db file
ok my $db = Bio::DB::Qual->new($test_dbdir, -reindex => 1);
is $db->glob, '*.{qual,QUAL,qa,QA}';
isa_ok $db, 'Bio::DB::Qual';
ok my @ids = $db->ids;
is scalar(@ids), 15;
@ids = sort {$a <=> $b} @ids;
is $ids[0] , '17601976';
is $ids[14], '17601991';
my $seqid = '17601979';

# direct indexed qual file database access
is ref($db->qual($seqid)), 'ARRAY';
is_deeply $db->qual($seqid), [23, 32, 24, 27, 26, 27, 27, 27, 28, 23, 28, 31, 23, 27];
is $db->length($seqid), 14;
is $db->length($seqid, -1000, 1000), 14; # length() ignores start and stop
is $db->header($seqid), '17601979';
is_deeply $db->qual($seqid, 2, 11), [32, 24, 27, 26, 27, 27, 27, 28, 23, 28];
is_deeply $db->qual($seqid, 2, 11, 1), [32, 24, 27, 26, 27, 27, 27, 28, 23, 28];
is_deeply $db->qual($seqid, 11, 2), [28, 23, 28, 27, 27, 27, 26, 27, 24, 32];
is_deeply $db->qual($seqid, 2, 11, -1), [28, 23, 28, 27, 27, 27, 26, 27, 24, 32];
is_deeply $db->qual($seqid, 11, 2, -1), [32, 24, 27, 26, 27, 27, 27, 28, 23, 28];

# the bioperl  way
is $db->get_Qual_by_id('foobarbaz'), undef;
ok my $obj = $db->get_Qual_by_id($seqid);
isa_ok $obj, 'Bio::Seq::PrimaryQual::Qual';
isa_ok $obj, 'Bio::Seq::QualI';
is ref($obj->qual($seqid)), 'ARRAY';
is $obj->length, 14;
is $obj->id, '17601979';
is $obj->display_id, '17601979';
is $obj->accession_number, 'unknown';
like $obj->primary_id, qr/^Bio::Seq::PrimaryQual::Qual=HASH/;
is $obj->validate_qual( join(' ', @{$obj->qual($seqid)}) ), 1;
is $obj->translate, 0;
is $obj->qualat(12), 31;
is_deeply $obj->subqual(2, 11), [32, 24, 27, 26, 27, 27, 27, 28, 23, 28];
is $obj->header, undef;
is $obj->desc, undef;
ok my $truncobj = $obj->trunc(1,3);
isa_ok $truncobj, 'Bio::Seq::PrimaryQual::Qual';
isa_ok $obj, 'Bio::Seq::QualI';
is ref($truncobj->qual($seqid)), 'ARRAY';
is $truncobj->length, 3;
ok my $revobj = $obj->revcom;
isa_ok $revobj, 'Bio::Seq::PrimaryQual::Qual';
isa_ok $revobj, 'Bio::Seq::PrimaryQual';
is ref($revobj->qual), 'ARRAY';
is $revobj->length, 14;
undef $obj;
undef $truncobj;
undef $revobj;

# using get_PrimarySeq_stream streaming
ok my $stream = $db->get_PrimaryQual_stream;
ok $stream = $db->get_PrimarySeq_stream;
isa_ok $stream, 'Bio::DB::Indexed::Stream';
ok my $streamqual = $stream->next_seq;
isa_ok $streamqual, 'Bio::Seq::PrimaryQual';

# using newFh streaming
ok my $fh = Bio::DB::Qual->newFh($test_dbdir);
my $fhqual = <$fh>;
isa_ok $fhqual, 'Bio::Seq::PrimaryQual';
undef $fh;

# tied-hash access
my (%h,$dna1,$dna2);
ok tie(%h,'Bio::DB::Qual',$test_dbdir);
ok $h{$seqid};
ok $dna1 = $h{"$seqid:1,10"};
ok $dna2 = $h{"$seqid:10,1"};

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
