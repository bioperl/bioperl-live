# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use lib '../..';
    use Bio::Root::Test;

    test_begin(-tests => 96);

    use_ok('Bio::SeqIO::embl');
}

my $verbose = test_debug();

my $ast = Bio::SeqIO->new( -format => 'embl',
                           -verbose => $verbose,
                           -file => test_input_file('roa1.dat'));
$ast->verbose($verbose);
my $as = $ast->next_seq();
ok defined $as->seq;
is($as->display_id, 'HSHNCPA1');
is($as->accession_number, 'X79536');
is($as->seq_version, 1);
is($as->version, 1);
is($as->desc, 'H.sapiens mRNA for hnRNPcore protein A1');
is($as->molecule, 'RNA');
is($as->alphabet, 'rna');
is(scalar $as->all_SeqFeatures(), 4);
is($as->length, 1198);
is($as->species->binomial(), 'Homo sapiens');

# EMBL Release 87 changes (8-17-06)

$ast = Bio::SeqIO->new( -format => 'embl',
                        -verbose => $verbose,
                        -file => test_input_file('roa1_v2.dat'));
$ast->verbose($verbose);
$as = $ast->next_seq();
ok defined $as->seq;
# accession # same as display name now
is($as->display_id, 'X79536'); 
is($as->accession_number, 'X79536');
is($as->seq_version, 1);
is($as->version, 1);
is($as->desc, 'H.sapiens mRNA for hnRNPcore protein A1');
# mRNA instead of RNA
is($as->molecule, 'mRNA');
is($as->alphabet, 'rna');
is(scalar $as->all_SeqFeatures(), 4);
is($as->length, 1198);
is($as->species->binomial(), 'Homo sapiens');

my $ent = Bio::SeqIO->new( -file => test_input_file('test.embl'),
                           -format => 'embl');
my $seq = $ent->next_seq();

is(defined $seq->seq(), 1,
   'success reading Embl with ^ location and badly split double quotes');
is(scalar $seq->annotation->get_Annotations('reference'), 3);

my $out_file = test_output_file();
my $out = Bio::SeqIO->new(-file=> ">$out_file",
                          -format => 'embl');
is($out->write_seq($seq),1,
   'success writing Embl format with ^ < and > locations');

# embl with no FT
$ent = Bio::SeqIO->new( -file => test_input_file('test.embl'),
                        -format => 'embl');
$seq = $ent->next_seq();

ok($seq);
is(lc($seq->subseq(1,10)),'gatcagtaga');
is($seq->length, 4870);

# embl with no FH
my $noFH = Bio::SeqIO->new( -file => test_input_file('no_FH.embl'),
                            -format => 'embl');
is(scalar($noFH->next_seq->get_SeqFeatures), 4);

# bug 1571
$ent = Bio::SeqIO->new( -format => 'embl',
                        -file   => test_input_file('test.embl2sq'));
is($ent->next_seq->length,4877);

# embl repbase
$ent = Bio::SeqIO->new(-file => test_input_file('BEL16-LTR_AG.embl'), -format => 'embl');
$seq = $ent->next_seq;
is($seq->display_id,'BEL16-LTR_AG');

# test secondary accessions in EMBL (bug #1332)
my $seqio = Bio::SeqIO->new( -format => 'embl',
                             -file => test_input_file('ECAPAH02.embl'));
$seq = $seqio->next_seq;
is($seq->accession_number, 'D10483');
is($seq->seq_version, 2);
my @accs = $seq->get_secondary_accessions();
is($accs[0], 'J01597');
is($accs[-1], 'X56742');

### TPA TESTS - Thanks to Richard Adams ###
# test Third Party Annotation entries in EMBL/Gb format 
# to ensure compatability with parsers.
my $str = Bio::SeqIO->new( -format =>'embl',
                           -file => test_input_file('BN000066-tpa.embl'));
$seq = $str->next_seq;
ok(defined $seq);
is($seq->accession_number, 'BN000066');
is($seq->alphabet, 'dna');
is($seq->display_id, 'AGA000066');
is($seq->length, 5195);
is($seq->division, 'INV');
is($seq->get_dates, 2);
is($seq->keywords, 'acetylcholinesterase; achE1 gene; Third Party Annotation; TPA.');
is($seq->seq_version, 1);
is($seq->feature_count, 15);

my $spec_obj = $seq->species;
is ($spec_obj->common_name, 'African malaria mosquito');
is ($spec_obj->species, 'gambiae');
is ($spec_obj->genus, 'Anopheles');
is ($spec_obj->binomial, 'Anopheles gambiae');

my $ac = $seq->annotation;
my $reference =  ($ac->get_Annotations('reference') )[1];
is ($reference->title,'"A novel acetylcholinesterase gene in mosquitoes codes for the insecticide target and is non-homologous to the ace gene in Drosophila"');
is ($reference->authors,'Weill M., Fort P., Berthomi eu A., Dubois M.P., Pasteur N., Raymond M.');
my $cmmnt =  ($ac->get_Annotations('comment') )[0];
is($cmmnt->text, 'see also AJ488492 for achE-1 from Kisumu strain Third Party Annotation Database: This TPA record uses Anopheles gambiae trace archive data (http://trace.ensembl.org) ');


$ent = Bio::SeqIO->new( -file => test_input_file('test.embl'),
                        -format => 'embl');
$ent->verbose($verbose);
$seq = $ent->next_seq();
my $species = $seq->species();
my @cl = $species->classification();
is( $cl[3] ne $species->genus(), 1, 'genus duplication test');
$ent->close();

#
## read-write - test embl writing of a PrimarySeq
#
my $primaryseq = Bio::PrimarySeq->new( -seq => 'AGAGAGAGATA',
                                      -id  => 'myid',
                                      -desc => 'mydescr',
                                      -alphabet => 'DNA',
                                      -accession_number => 'myaccession');



$verbose = -1 unless $ENV{'BIOPERLDEBUG'};  # silence warnings unless we are debuggin

my $embl = Bio::SeqIO->new(-format => 'embl',
                          -verbose => $verbose,
                          -file => ">$out_file");

ok($embl->write_seq($primaryseq));

# this should generate a warning
my $scalar = "test";
eval {
    $embl->write_seq($scalar);
};
ok ($@);

# CDS records
# (which have nonstandard 'PA' and 'OX' tags)
# see http://bioperl.org/pipermail/bioperl-l/2009-February/029252.html
# and the rest of that thread
my $cds_file = Bio::SeqIO->new(-format =>'embl',
                               -file => test_input_file('cds_sample.embl'));
my $cds_seq = $cds_file->next_seq;
ok(defined $cds_seq);
is($cds_seq->display_id, 'EAL24309');
is($cds_seq->accession_number, 'CH236947.1', 'CDS - accession on PA line');
is($cds_seq->alphabet, 'dna');
is($cds_seq->length, 192);
is($cds_seq->species->binomial(), 'Homo sapiens');
is($cds_seq->seq_version, 1);
is($cds_seq->feature_count, 2);
my $cds_annot = $cds_seq->annotation;
ok(defined $cds_annot);
my $cds_dblink = ($cds_annot->get_Annotations('dblink'))[0];
ok(defined $cds_dblink);
is($cds_dblink->tagname, 'dblink', 'CDS - OX tagname');
is($cds_dblink->database, 'NCBI_TaxID', 'CDS - OX database');
is($cds_dblink->primary_id, '9606', 'CDS - OX primary_id');

#bug 2982 - parsing contig descriptions sans sequence data

ok( $embl = Bio::SeqIO->new( -file => test_input_file('bug2982.embl'),
                             -format => 'embl') );
my $i;
for ($i=0; my $seq = $embl->next_seq; $i++) {
    ok !$seq->seq;
    ok ( my $ann = ($seq->annotation->get_Annotations('contig'))[0] );
    like $ann->value, qr/join\(/;
}
is $i, 4;


# bug 3086 - parsing long lines correctly

ok( $embl = Bio::SeqIO->new(-file => test_input_file('bug3086.embl'),
                            -format => 'embl',
                            -verbose => '$verbose') );
$seq = $embl->next_seq;
foreach my $feature ($seq->top_SeqFeatures) {
    if ($feature->has_tag('product')) {
        my ($product) = $feature->get_tag_values('product');
        is($product,
           'bifunctional phosphoribosylaminoimidazolecarboxamide formyltransferase/IMP cyclohydrolase',
           'Check if product was parsed correctly');
    }
}

# long labels handled

{
    # Create sequence with feature with a long label qualifier
    my $seq=Bio::Seq->new(-seq=>'actg');
    my $feature=Bio::SeqFeature::Generic->new(-primary=>'CDS', -start=>1, -end=>4);
    $feature->add_tag_value(label=>'1 2 3 4 5 6 7 8 9 a b c d e f g h i j k l m n o p q r');
    $seq->add_SeqFeature($feature);

    # Write EMBL
    my $string;
    open my $str_fh, '>', \$string or skip("Could not write string, skipping", 2);
    
    my $out=Bio::SeqIO->new(-format=>'embl', -fh => $str_fh);
    $out->write_seq($seq);

    # Read EMBL
    my $in=Bio::SeqIO->new(-format=>'embl', -string => $string);
    my $ret=eval { my $embl=$in->next_seq };
    my $error;
    $error=$@ if (!$ret);
    ok($ret, 'Parse long qualifier');
    is($error, undef);
}

# NCBI TaxIDs should roundtrip
{
    my $seq=Bio::Seq->new(-seq=>'actg');
    my $species = Bio::Species->new(-ncbi_taxid => 7165, -classification=>
                                    [ 'Anopheles gambiae', 'Anopheles', 'Culicoidea',
                                      'Nematocera', 'Diptera', 'Endopterygota',
                                      'Neoptera', 'Pterygota', 'Insecta', 'Hexapoda',
                                      'Arthropoda', 'Metazoa', 'Eukaryota' ]);

    $seq->species($species);
    is($seq->species->ncbi_taxid, 7165, 'TaxID set correctly');

    # Write EMBL
    my $string;
    open my $str_fh, '>', \$string or skip("Could not write string, skipping", 2);

    my $out=Bio::SeqIO->new(-format=>'embl', -fh => $str_fh);
    $out->write_seq($seq);

    # Read EMBL
    my $in=Bio::SeqIO->new(-format=>'embl', -string => $string);
    my $embl_seq;
    my $ret=eval { $embl_seq=$in->next_seq };
    my $error;
    $error=$@ if (!$ret);

    # Check that TaxID has roundtripped
    my $embl_species = $embl_seq->species;
    ok(defined $embl_species, "The read sequence has a species object");
    is($embl_species->ncbi_taxid, 7165, "NCBI TaxID has roundtripped");
    is($embl_species->binomial(), 'Anopheles gambiae', "Name has roundtripped");
}

# a taxon db_xref on a source feature should override an OX line
{
    my $seq=Bio::Seq->new(-seq=>'actg');
    my $species = Bio::Species->new(-ncbi_taxid => 7165, -classification=>
                                    [ 'Anopheles gambiae', 'Anopheles', 'Culicoidea',
                                      'Nematocera', 'Diptera', 'Endopterygota',
                                      'Neoptera', 'Pterygota', 'Insecta', 'Hexapoda',
                                      'Arthropoda', 'Metazoa', 'Eukaryota' ]);

    $seq->species($species);
    is($seq->species->ncbi_taxid, 7165, 'TaxID set correctly');

    my $seq_feature = Bio::SeqFeature::Generic->new(-primary=>'source',
                                                    -start => 1,
                                                    -end=> length($seq->seq));

    $seq_feature->add_tag_value('db_xref', 'taxon:9606');
    $seq->add_SeqFeature($seq_feature);

    # Write EMBL
    my $string;
    open my $str_fh, '>', \$string or skip("Could not write string, skipping", 2);

    my $out=Bio::SeqIO->new(-format=>'embl', -fh => $str_fh);
    $out->write_seq($seq);

    # Read EMBL
    my $in=Bio::SeqIO->new(-format=>'embl', -string => $string);
    my $embl_seq;
    my $ret=eval { $embl_seq=$in->next_seq };
    my $error;
    $error=$@ if (!$ret);

    # Check that TaxID has roundtripped
    my $embl_species = $embl_seq->species;
    ok(defined $embl_species, "The read sequence has a species object");
    is($embl_species->ncbi_taxid, 9606, "The taxid of the source feature overrides that of the OX line");
    is($embl_species->binomial(), 'Anopheles gambiae', "Name has roundtripped");
}

# Handle Seq objects that only define an ID, not an accession number
{
    my $seq = Bio::Seq->new(-seq=>'actg', -id=>'test_id');

    my $string;
    open my $str_fh, '>', \$string or skip("Could not write string, skipping", 1);

    my $out = Bio::SeqIO->new(-format=>'embl', -fh=>$str_fh);
    $out->write_seq($seq);

    ok($string =~ m/ID   test_id;/, "The ID field was written correctly");
}
