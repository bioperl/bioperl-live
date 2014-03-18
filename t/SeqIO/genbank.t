# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    test_begin(-tests => 287);
    use_ok('Bio::SeqIO::genbank');
}

my $verbose = test_debug;

my $ast = Bio::SeqIO->new(-format  => 'genbank' ,
                          -verbose => $verbose,
                          -file    => test_input_file('roa1.genbank'));
isa_ok($ast, 'Bio::SeqIO');
$ast->verbose($verbose);
my $as = $ast->next_seq;
is $as->molecule, 'mRNA',$as->accession_number;
is $as->alphabet, 'dna';
is $as->division, 'EST';
is join(',',$as->get_dates), '27-OCT-1998';
is($as->primary_id, 3598416);
my @class = $as->species->classification;
is $class[$#class],'Eukaryota';

$ast = Bio::SeqIO->new(-format  => 'genbank',
                       -verbose => $verbose,
                       -file    => test_input_file('NT_021877.gbk'));
$ast->verbose($verbose);
$as = $ast->next_seq;
is $as->molecule, 'DNA',$as->accession_number;
is $as->alphabet, 'dna';
is $as->division, 'CON';
is join(',',$as->get_dates), '17-OCT-2003';
is($as->primary_id, 37539616);
is($as->accession_number, 'NT_021877');

my ($cds) = grep { $_->primary_tag eq 'CDS' } $as->get_SeqFeatures;
is(($cds->get_tag_values('transl_except'))[1],
   '(pos:complement(4224..4226),aa:OTHER)');

# test for a DBSOURCE line
$ast = Bio::SeqIO->new(-format  => 'genbank',
                       -verbose => $verbose,
                       -file    => test_input_file('BAB68554.gb'));
$ast->verbose($verbose);
$as = $ast->next_seq;
is $as->molecule, 'PRT',$as->accession_number;
is $as->alphabet, 'protein';
is $as->division, 'VRT';
is join(',',$as->get_dates), '11-APR-2002';
# Though older GenBank releases indicate SOURCE contains only the common name,
# this is no longer true.  In general, this line will contain an abbreviated
# form of the full organism name (but may contain the full length name),
# as well as the optional common name and organelle.  There is no get/set
# for the abbreviated name but it is accessible via name()
ok defined($as->species->name('abbreviated')->[0]);
is $as->species->name('abbreviated')->[0], 'Aldabra giant tortoise';
is($as->primary_id, 15824047);
my $ac = $as->annotation;
ok defined $ac;
my @dblinks = $ac->get_Annotations('dblink');
is(scalar @dblinks,1);
is($dblinks[0]->database, 'GenBank');
is($dblinks[0]->primary_id, 'AB072353');
is($dblinks[0]->version, '1');
is($dblinks[0]->display_text, 'GenBank:AB072353.1','operator overloading in AnnotationI is deprecated');

# test for multi-line SOURCE
$ast = Bio::SeqIO->new(-format  => 'genbank',
                       -verbose => $verbose,
                       -file    => test_input_file('NC_006346.gb'));
$as = $ast->next_seq;
is $as->species->binomial('FULL'), 'Bolitoglossa n. sp. RLM-2004',$as->accession_number;;
@class = $as->species->classification;
is($class[$#class],'Eukaryota');
is($as->species->common_name,'mushroomtongue salamander');

$ast = Bio::SeqIO->new(-format  => 'genbank',
                       -verbose => $verbose,
                       -file    => test_input_file('U71225.gb'));
$as = $ast->next_seq;
@class = $as->species->classification;
is($class[$#class],'Eukaryota',$as->accession_number);
is $as->species->common_name,'black-bellied salamander';

# test for unusual common name
$ast = Bio::SeqIO->new(-format  => 'genbank',
                       -verbose => $verbose,
                       -file    => test_input_file('AB077698.gb'));
$as = $ast->next_seq;
# again, this is not a common name but is in name('abbreviated')
ok defined($as->species->name('abbreviated')->[0]),$as->accession_number;
is $as->species->name('abbreviated')->[0],'Homo sapiens cDNA to mRNA';

# test for common name with parentheses
$ast = Bio::SeqIO->new(-format  => 'genbank',
                       -verbose => $verbose,
                       -file    => test_input_file('DQ018368.gb'));
$as = $ast->next_seq;
is $as->species->scientific_name,'(Populus tomentosa x P. bolleana) x P. tomentosa var. truncata',
$as->accession_number;;

# test secondary accessions
my $seqio = Bio::SeqIO->new(-format  => 'genbank',
                            -verbose => $verbose,
                            -file    => test_input_file('D10483.gbk'));
my $seq = $seqio->next_seq;
my @kw =  $seq->get_keywords;
is(scalar @kw, 118, $seq->accession_number);
is($kw[-1], 'yabO');
my @sec_acc = $seq->get_secondary_accessions;
is(scalar @sec_acc,14);
is($sec_acc[-1], 'X56742');

# bug #1487
my $str = Bio::SeqIO->new(-verbose => $verbose,
                          -file    => test_input_file('D12555.gbk'));
eval {
    $seq = $str->next_seq;
};

ok(! $@, 'bug 1487');

# bug 1647 rpt_unit sub-feature with multiple parens
$str = Bio::SeqIO->new(-format  => 'genbank',
                       -verbose => $verbose,
                       -file    => test_input_file('mini-AE001405.gb'));
ok($seq = $str->next_seq);
my @rpts = grep { $_->primary_tag eq 'repeat_region' }
  $seq->get_SeqFeatures;
is $#rpts, 2, 'bug 1647';
my @rpt_units = grep {$_->has_tag('rpt_unit')} @rpts;
is $#rpt_units, 0;
is(($rpt_units[0]->get_tag_values('rpt_unit'))[0],'(TG)10;A;(TG)7');

# test bug #1673 , RDB-II genbank files
$str = Bio::SeqIO->new(-format  => 'genbank',
                       -verbose => $verbose,
                       -file    => test_input_file('Mcjanrna_rdbII.gbk')
              );
ok($seq = $str->next_seq, 'bug 1673');
my @refs = $seq->annotation->get_Annotations('reference');
is(@refs, 1);
is($seq->display_id,'Mc.janrrnA');
is($seq->molecule ,'RNA');
is $as->division, 'PLN';
is join(',',$as->get_dates), '23-MAY-2005';

$str  = Bio::SeqIO->new(-format  => 'genbank',
                        -file    => test_input_file('AF165282.gb'),
                        -verbose => $verbose);
$seq = $str->next_seq;
my @features = $seq->all_SeqFeatures;
is(@features, 5, $seq->accession_number);
is($features[0]->start, 1);
is($features[0]->end, 226);
my $location = $features[1]->location;
ok($location->isa('Bio::Location::SplitLocationI'));
my @sublocs = $location->sub_Location;
is(@sublocs, 29);

# version and primary ID - believe it or not, this wasn't working
is ($seq->version, 1);
is ($seq->seq_version, 1);
is ($seq->primary_id, "5734104");

# streaming and Bio::RichSeq creation
my $stream = Bio::SeqIO->new(-file    => test_input_file('test.genbank'),
                             -verbose => $verbose,
                             -format  => 'genbank');
$stream->verbose($verbose);
my $seqnum = 0;
my $species;
my @cl;
my $lasts;
my @ids = qw(DDU63596 DDU63595 HUMBDNF);
my @tids = (44689, 44689, 9606);
my @tnames = ("Dictyostelium discoideum",
              "Dictyostelium discoideum",
              "Homo sapiens");
while($seq = $stream->next_seq) {
    if($seqnum < 3) {
        is $seq->display_id, $ids[$seqnum];
        $species = $seq->species;
        @cl = $species->classification;
        is( $species->binomial, $tnames[$seqnum],
            'species parsing incorrect for genbank');
        is( $cl[3] ne $species->genus, 1,
            'genus duplicated in genbank parsing');
        is( $species->ncbi_taxid, $tids[$seqnum] );
    }
    $seqnum++;
    $lasts = $seq;
}
is($seqnum, 5,'streaming');
is $lasts->display_id, "HUMBETGLOA";
my ($ref) = $lasts->annotation->get_Annotations('reference');
is($ref->medline, 94173918);
$stream->close;

$stream = Bio::SeqIO->new(-file    => test_input_file('test.genbank.noseq'),
                          -verbose => $verbose,
                          -format  => 'genbank' );
$seqnum = 0;
while($seq = $stream->next_seq) {
    if($seqnum < 3) {
        is $seq->display_id, $ids[$seqnum];
    }
    elsif( $seq->display_id eq 'M37762') {
        is( ($seq->get_keywords)[0], 'neurotrophic factor');
    }
    $seqnum++;
}
is $seqnum, 5, "Total number of sequences in test file";

# fuzzy
$seq = Bio::SeqIO->new( -file    => test_input_file('testfuzzy.genbank'),
                        -format  => 'genbank',
                        -verbose => $verbose );
ok(defined($as = $seq->next_seq));

@features = $as->all_SeqFeatures;
is(@features,21,'Fuzzy in');
my $lastfeature = pop @features;
# this is a split location; the root doesn't have strand
is($lastfeature->strand, undef);
$location = $lastfeature->location;
#$location->verbose(-1); # silence the warning of undef seq_id()
# see above; splitlocs roots do not have a strand really
is($location->strand, undef);
is($location->start, 83202);
is($location->end, 84996);

@sublocs = $location->sub_Location;

is(@sublocs, 2);
my $loc = shift @sublocs;
is($loc->start, 83202);
is($loc->end, 83329);
is($loc->strand, -1);

$loc = shift @sublocs;
is($loc->start, 84248);
is($loc->end, 84996);
is($loc->strand,1);

$seq = Bio::SeqIO->new(-format  => 'genbank',
                       -verbose => $verbose,
                       -file    => ">" .test_output_file);
$seq->verbose($verbose);
ok($seq->write_seq($as),'Fuzzy out');

## now genbank ##
$str = Bio::SeqIO->new(-format  => 'genbank',
                       -verbose => $verbose,
                       -file    => test_input_file('BK000016-tpa.gbk'));
$seq = $str->next_seq;
ok(defined $seq, $seq->accession_number);
ok(defined $seq->seq);
is($seq->accession_number, 'BK000016',$seq->accession_number);
is($seq->alphabet, 'dna');
is($seq->display_id, 'BK000016');
is($seq->length, 1162);
is($seq->division, 'ROD');
is($seq->get_dates, 1);
is($seq->keywords, 'Third Party Annotation; TPA');
is($seq->desc, 'TPA: Mus musculus pantothenate kinase 4 mRNA, partial cds.');
is($seq->seq_version, 1);
is($seq->feature_count, 2);
my $spec_obj = $seq->species;
is ($spec_obj->common_name, 'house mouse');
is ($spec_obj->species, 'musculus');
is ($spec_obj->genus, 'Mus');
is ($spec_obj->binomial, 'Mus musculus');
$ac = $seq->annotation;
my $reference =  ($ac->get_Annotations('reference') )[0];
is ($reference->pubmed, '11479594');
is ($reference->medline, '21372465',$seq->accession_number);

# validate that what is written is what is read
my $testfile = test_output_file;
my $out = Bio::SeqIO->new(-file   => ">$testfile",
                          -format => 'genbank');
$out->write_seq($seq);
$out->close;

$str = Bio::SeqIO->new(-format => 'genbank',
                       -file   => $testfile);
$seq = $str->next_seq;
ok(defined $seq,'roundtrip');
ok(defined $seq->seq);
is($seq->accession_number, 'BK000016');
is($seq->alphabet, 'dna');
is($seq->display_id, 'BK000016');
is($seq->length, 1162);
is($seq->division, 'ROD');
is($seq->get_dates, 1);
is($seq->keywords, 'Third Party Annotation; TPA');
is($seq->desc, 'TPA: Mus musculus pantothenate kinase 4 mRNA, partial cds.');
is($seq->seq_version, 1);
is($seq->feature_count, 2);
$spec_obj = $seq->species;
is ($spec_obj->common_name, 'house mouse');
is ($spec_obj->species, 'musculus');
is ($spec_obj->genus, 'Mus');
is ($spec_obj->binomial, 'Mus musculus');
$ac = $seq->annotation;
$reference =  ($ac->get_Annotations('reference') )[0];
is ($reference->pubmed, '11479594');
is ($reference->medline, '21372465');

# write revcomp split location
my $gb = Bio::SeqIO->new(-format => 'genbank',
                         # This sequence has an odd LOCUS line which sets off a warning, setting
                         # verbose to -1.
                         # The newest Ensembl seq lacks this.  Maybe update?  cjfields 6-5-07
                         -verbose => $verbose ? $verbose : -1,
                         -file   => test_input_file('revcomp_mrna.gb'));
$seq = $gb->next_seq;

$gb = Bio::SeqIO->new(-format => 'genbank',
                      -file   => ">$testfile");

$gb->write_seq($seq);
undef $gb;
ok(! -z $testfile, 'revcomp split location');

# bug 1925, continuation of long ORGANISM line ends up in @classification:
# ORGANISM  Salmonella enterica subsp. enterica serovar Paratyphi A str. ATCC
#           9150
#           Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales;
#           Enterobacteriaceae; Salmonella.
$gb = Bio::SeqIO->new(-format  => 'genbank',
                      -verbose => $verbose,
                      -file    => test_input_file('NC_006511-short.gbk'));
$seq = $gb->next_seq;
is $seq->species->common_name, undef, "Bug 1925";
is $seq->species->scientific_name, "Salmonella enterica subsp. enterica serovar Paratyphi A str. ATCC 9150";
@class = $seq->species->classification;
is $class[$#class], "Bacteria";

# WGS   tests
$gb = Bio::SeqIO->new(-format  => 'genbank',
                      -verbose => $verbose,
                      -file    => test_input_file('O_sat.wgs'));
$seq = $gb->next_seq;

my @tests = ('wgs'        => 'AAAA02000001-AAAA02050231',
             'wgs_scafld' => 'CM000126-CM000137',
             'wgs_scafld' => 'CH398081-CH401163');

my @wgs = map {$seq->annotation->get_Annotations(lc($_))} qw(WGS WGS_SCAFLD);

my $ct=0;

for my $wgs (@wgs) {
    my ($tagname, $value) = (shift @tests, shift @tests);
    is($wgs->tagname, $tagname, $tagname);
    is($wgs->value, $value);
    $ct++;
}

is ($ct, 3);

# make sure we can retrieve a feature with a primary tag of 'misc_difference'
$gb = Bio::SeqIO->new(-format  => 'genbank',
                      -verbose => $verbose,
                      -file    => test_input_file('BC000007.gbk'));
$seq = $gb->next_seq;
($cds) = grep { $_->primary_tag eq 'misc_difference' } $seq->get_SeqFeatures;
my @vals = $cds->get_tag_values('gene');
is $vals[0], 'PX19', $seq->accession_number;

# Check that the source,organism section is identical between input and output.
# - test an easy one where organism is species, then two different formats of
# subspecies, then a species with a format that used to be mistaken for
# subspecies, then a bacteria with no genus, and finally a virus with a genus.

# These tests are now somewhat out-of-date since we are moving to a Bio::Taxon-
# based system for verifying taxonomic information.  Right now they just verify
# changes so are really useless; I will change them to verify common name,
# organelle, scientific name, etc.

my $outfile = test_output_file;

# output always adds a period (GenBank std), but two of these files do not use them.

foreach my $in ('BK000016-tpa.gbk', 'ay116458.gb', 'ay149291.gb', 'NC_006346.gb', 'ay007676.gb', 'dq519393.gb') {
    my $infile =  test_input_file($in);

    $str = Bio::SeqIO->new(-format  => 'genbank',
                           -verbose => $verbose,
                           -file    => $infile);
    $seq = $str->next_seq;

    $out = Bio::SeqIO->new(-file   => ">$outfile",
                           -format => 'genbank');
    $out->write_seq($seq);
    $out->close;

    open my $IN, '<', $infile or die "Could not read file '$infile': $!\n";
    my @in = <$IN>;
    close $IN;

    open my $RESULT, '<', $outfile or die "Could not read file '$outfile': $!\n";
    my $line = 0;
    my $check = 0;
    my $is = 1;

  FILECHECK:
    while (my $result = <$RESULT>) {
        if ($result =~ /^KEYWORDS/) {
            $check = 1;
            next;
        }

        if ($result =~ /^REFERENCE/) {
            last FILECHECK;
        }

        if ($check) {
            # end periods don't count (not all input files have them)
            $result =~ s{\.$}{};
            $in[$line] =~ s{\.$}{};

            if ($result ne $in[$line]) {
                $is = 0;
                last;
            }
        }
    }
    continue {
        $line++;
    }
    close $RESULT;

    ok $is, $in;
}

# NB: there should probably be full testing on all lines to ensure that output
# matches input.

# 20061117: problem with *double* colon in some annotation-dblink values
$ct = 0;

foreach my $in ('P35527.gb') {
    my $infile =  test_input_file($in);
    $str = Bio::SeqIO->new(-format  => 'genbank',
                           -verbose => $verbose,
                           -file    => $infile);
    $seq = $str->next_seq;

    my $ac = $seq->annotation;      # Bio::AnnotationCollection
    foreach my $key ($ac->get_all_annotation_keys ) {
        my @values = $ac->get_Annotations($key);
        foreach my $ann (@values) {
            my $value = $ann->display_text;
            $ct++;
            if ($key eq 'dblink') {
                ok (index($value,'::') < 0);   # this should never be true
                ok ($value, $value);           # check value is not empty

                #  print "  ann/", sprintf('%12s  ',$key), '>>>', $value , '<<<', "\n";
                #  print "        index double colon: ",index($value   ,'::'), "\n";

                #  check db name:
                my @parts = split(/:/,$value);
                if ( $parts[0] =~ /^(?:
                        #  not an exhaustive list of databases;
                        #  just the db's referenced in P35527.gb:
                          swissprot | GenBank | GenPept      | HSSP | IntAct   | Ensembl | KEGG
                        | HGNC      | MIM     | ArrayExpress | GO   | InterPro | Pfam    | PRINTS
                        | PROSITE    )$/x
                ) {
                    ok 1;
                }
                else {
                    ok 0;
                }
                    ok ( $parts[1], "$parts[0]" );
            }
        }
    }
}

is($ct, 46);

# bug 2195

$str = Bio::SeqIO->new(-format  => 'genbank',
                       -verbose => $verbose,
                       -file    => test_input_file('AF305198.gb')
                     );

$species = $str->next_seq->species;

is($species->scientific_name, 'Virginia creeper phytoplasma', 'Bug 2195');
is(join(', ',$species->classification),
     'Virginia creeper phytoplasma, 16SrV (Elm yellows group), '
   . 'Candidatus Phytoplasma, Acholeplasmataceae, Acholeplasmatales, '
   . 'Mollicutes, Firmicutes, Bacteria',
   'Bug 2195');

# bug 2569, PROJECT line support, read and write, round-tripping

$str = Bio::SeqIO->new(-format  => 'genbank',
                       -verbose => $verbose,
                       -file    => test_input_file('NC_008536.gb'));
$seq = $str->next_seq;

my $project = ($seq->annotation->get_Annotations('project'))[0];
isa_ok($project, 'Bio::Annotation::SimpleValue');

if ($project) {
    is($project->value, 'GenomeProject:12638');
}
else {
    ok(0, "PROJECT not parsed");
}

$outfile = test_output_file;
$gb = Bio::SeqIO->new(-format  => 'genbank',
                      -verbose => $verbose,
                      -file    => ">$outfile");
$gb->write_seq($seq);

$str = Bio::SeqIO->new(-format  => 'genbank',
                       -verbose => $verbose,
                       -file    => $outfile);
$seq = $str->next_seq;

$project = ($seq->annotation->get_Annotations('project'))[0];
isa_ok($project, 'Bio::Annotation::SimpleValue');

if ($project) {
    is($project->value, 'GenomeProject:12638');
}
else {
    ok(0, "Roundtrip test failed");
}

# test for swissprot/UniProt/UniProtKB DBSOURCE line (Bug : RT 44536)
$ast = Bio::SeqIO->new(-format  => 'genbank',
                       -verbose => $verbose,
                       -file    => test_input_file('P39765.gb'));
$as = $ast->next_seq;
is $as->molecule, 'PRT',$as->accession_number;
is $as->division, 'BCT',$as->accession_number;
is join(',',$as->get_dates), '03-MAR-2009',$as->accession_number;
is $as->alphabet, 'protein';
# Though older GenBank releases indicate SOURCE contains only the common name,
# this is no longer true.  In general, this line will contain an abbreviated
# form of the full organism name (but may contain the full length name),
# as well as the optional common name and organelle.  There is no get/set
# for the abbreviated name but it is accessible via name()
ok defined($as->species->name('abbreviated')->[0]);
is $as->species->name('abbreviated')->[0], 'Bacillus subtilis';
is($as->primary_id, 20141743);
$ac = $as->annotation;
ok defined $ac;
@dblinks = $ac->get_Annotations('dblink');
is(scalar @dblinks,31);
is($dblinks[0]->database, 'UniProtKB');
is($dblinks[0]->primary_id, 'PYRR_BACSU');
is($dblinks[0]->version, undef);
is($dblinks[0]->display_text, 'UniProtKB:PYRR_BACSU','operator overloading in AnnotationI is deprecated');

#bug 2982 embl/genbank contig handling

$ast = Bio::SeqIO->new( -file   => test_input_file('bug2982.gb'),
                        -format => 'genbank' );
$seq = $ast->next_seq;

ok my @ctg = $seq->annotation->get_Annotations('contig');
like $ctg[0]->value, qr/join\(.*?gap.*?complement/;

# write_seq() and FTHelper duplicate specific tags, need to check a round-trip
$ast = Bio::SeqIO->new(-format  => 'genbank' ,
                       -verbose => $verbose,
                       -file    => test_input_file('singlescore.gbk'));
$as = $ast->next_seq;
($cds) = grep { $_->primary_tag eq 'CDS' } $as->get_SeqFeatures;
my @notes = $cds->get_tag_values('note');
is(scalar @notes, 2);
$testfile = test_output_file;
$out = Bio::SeqIO->new(-file   => ">$testfile",
                       -format => 'genbank');
$out->write_seq($as);
$out->close;
$ast = Bio::SeqIO->new(-format  => 'genbank' ,
                       -verbose => $verbose,
                       -file    => $testfile );
$as = $ast->next_seq;
($cds) = grep { $_->primary_tag eq 'CDS' } $as->get_SeqFeatures;
@notes = $cds->get_tag_values('note');
is(scalar @notes, 2);


#bug 3375
my $in = Bio::SeqIO->new(-format => 'genbank',
                         -file   => test_input_file('NC_002058_multDBLINK_bug3375.gb'));
$seq = $in->next_seq;     # should not throw a warning now
@dblinks = $seq->annotation->get_Annotations('dblink');    # contains 5 dblink references
# testing DBLINK      BioProject: PRJNA15288
is($dblinks[0]->database, 'BioProject', 'bug3375 database is BioProject');
is($dblinks[0]->primary_id, 'PRJNA15288', 'bug3375 primary_id is PRJNA15288');
# testing DBLINK      Project:100,200,300
is($dblinks[3]->database, 'Project');
is($dblinks[3]->primary_id, '300');
# testing DBLINK      NC_002058.3
is($dblinks[4]->database, 'GenBank');
is($dblinks[4]->primary_id, 'NC_002058');
is($dblinks[4]->version, '3');

# long labels handled
{
    # Create sequence with feature with a long label qualifier
    my $seq=Bio::Seq->new(-seq => 'actg',
                          -id  => 'abacab');
    my $feature=Bio::SeqFeature::Generic->new(-primary=>'CDS', -start=>1, -end=>4);
    my $label='1 2 3 4 5 6 7 8 9 a b c d e f g h i j k l m n o p q r';
    $feature->add_tag_value(label => $label);
    $seq->add_SeqFeature($feature);

    # Write genbank
    my $string;
    open my $str_fh, '>', \$string or skip("Could not write string, skipping", 2);
    my $out = Bio::SeqIO->new(-format => 'genbank',
                              -fh     => $str_fh);
    $out->write_seq($seq);

    # Read genbank
    my $in = Bio::SeqIO->new(-format => 'genbank',
                             -string => $string);
    my $genbank = $in->next_seq;
    my ($read_feature) = $genbank->get_SeqFeatures;
    my ($read_label) = $read_feature->get_tag_values('label');
    is($read_label, $label, 'Label is the same');
}

# bug 3448
$in = Bio::SeqIO->new(-format  => 'genbank',
                      -file    => test_input_file('YP_007988852.gp'),
                      -verbose => $verbose);
$seq = $in->next_seq;     # should not throw a warning now
is($seq->length, 205);

my @anns = $seq->annotation->get_Annotations('contig');
is(@anns, 1);
isa_ok($anns[0], 'Bio::Annotation::SimpleValue');
is($anns[0]->value, 'join(WP_015639704.1:1..205)');

is($seq->seq, 'MENRKFGYIRVSSKDQNEGRQLEAMRKIGITERDIYLDKQSGKNFERANYQLLKRIIRKGDI'
            . 'LYIHSLDRFGRNKEEILQEWNDLTKNIEADIVVLDMPLLDTTQYKDSMGTFIADLVLQILSWMAEEERERIRK'
            . 'RQREGIDLALQNGIQFGRSPVVVSDEFKEVYRKWKAKELTAVEAMQEAGVKKTSFYKLVKAHENSIKVNS');

