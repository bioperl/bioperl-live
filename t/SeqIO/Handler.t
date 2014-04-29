# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests           => 561,
               -requires_module => 'Data::Stag');
	
    use_ok('Bio::SeqIO');
}


my $verbose = test_debug();

################################## GenBank ##################################

my $ast = Bio::SeqIO->new(-format => 'gbdriver' ,
                        -verbose => $verbose,
                        -file => test_input_file("roa1.genbank"));
$ast->verbose($verbose);
my $as = $ast->next_seq();
is $as->molecule, 'mRNA',$as->accession_number;
is $as->alphabet, 'dna';
is($as->primary_id, 3598416);
my @class = $as->species->classification;
is $class[$#class],'Eukaryota';

$ast = Bio::SeqIO->new(-format => 'gbdriver',
                              -verbose => $verbose,
                       -file => test_input_file("NT_021877.gbk"));
$ast->verbose($verbose);
$as = $ast->next_seq();
is $as->molecule, 'DNA',$as->accession_number;
is $as->alphabet, 'dna';
is($as->primary_id, 37539616);
is($as->accession_number, 'NT_021877');

my ($cds) = grep { $_->primary_tag eq 'CDS' } $as->get_SeqFeatures();
is(($cds->get_tag_values('transl_except'))[1],
   '(pos:complement(4224..4226),aa:OTHER)');

# test for a DBSOURCE line
$ast = Bio::SeqIO->new(-format => 'gbdriver',
                              -verbose => $verbose,
                       -file => test_input_file("BAB68554.gb"));
$ast->verbose($verbose);
$as = $ast->next_seq();
is $as->molecule, 'linear',$as->accession_number;;
is $as->alphabet, 'protein';
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
$ast = Bio::SeqIO->new(-format => 'gbdriver',
                              -verbose => $verbose,
                       -file => test_input_file("NC_006346.gb"));
$as = $ast->next_seq;
is $as->species->binomial('FULL'), 'Bolitoglossa n. sp. RLM-2004',$as->accession_number;;
@class = $as->species->classification;
is($class[$#class],'Eukaryota');
is($as->species->common_name,'mushroomtongue salamander');

$ast = Bio::SeqIO->new(-format => 'gbdriver',
                              -verbose => $verbose,
                       -file => test_input_file("U71225.gb"));
$as = $ast->next_seq;
@class = $as->species->classification;
is($class[$#class],'Eukaryota',$as->accession_number);
is $as->species->common_name,'black-bellied salamander';

# test for unusual common name
$ast = Bio::SeqIO->new(-format => 'gbdriver',
                              -verbose => $verbose,
                       -file => test_input_file("AB077698.gb"));
$as = $ast->next_seq;
# again, this is not a common name but is in name('abbreviated')
ok defined($as->species->name('abbreviated')->[0]),$as->accession_number;
is $as->species->name('abbreviated')->[0],'Homo sapiens cDNA to mRNA';

# test for common name with parentheses
$ast = Bio::SeqIO->new(-format => 'gbdriver',
                              -verbose => $verbose,
                       -file => test_input_file("DQ018368.gb"));
$as = $ast->next_seq;
is $as->species->scientific_name,'(Populus tomentosa x P. bolleana) x P. tomentosa var. truncata',
$as->accession_number;;

# test secondary accessions
my $seqio = Bio::SeqIO->new(-format => 'gbdriver',
                                    -verbose => $verbose,
                                    -file => test_input_file('D10483.gbk'));
my $seq = $seqio->next_seq;
my @kw =  $seq->get_keywords;
is(scalar @kw, 118, $seq->accession_number);
is($kw[-1], 'yabO');
my @sec_acc = $seq->get_secondary_accessions();
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
$str = Bio::SeqIO->new(-format => 'gbdriver',
                              -verbose => $verbose,
                       -file => test_input_file('mini-AE001405.gb'));
ok($seq = $str->next_seq);
my @rpts = grep { $_->primary_tag eq 'repeat_region' }
  $seq->get_SeqFeatures;
is $#rpts, 2, 'bug 1647';
my @rpt_units = grep {$_->has_tag('rpt_unit')} @rpts;
is $#rpt_units, 0;
is(($rpt_units[0]->get_tag_values('rpt_unit'))[0],'(TG)10;A;(TG)7');

# test bug #1673 , RDB-II genbank files
$str = Bio::SeqIO->new(-format => 'gbdriver',
                              -verbose => $verbose,
                       -file => test_input_file('Mcjanrna_rdbII.gbk')
              );
ok($seq = $str->next_seq, 'bug 1673');
my @refs = $seq->annotation->get_Annotations('reference');
is(@refs, 1);
is($seq->display_id,'Mc.janrrnA');
is($seq->molecule ,'RNA');

$str  = Bio::SeqIO->new(-format => 'gbdriver',
                              -file   => test_input_file("AF165282.gb"),
                              -verbose => $verbose);
$seq = $str->next_seq;
my @features = $seq->all_SeqFeatures();
is(@features, 5, $seq->accession_number);
is($features[0]->start, 1);
is($features[0]->end, 226);
my $location = $features[1]->location;
ok($location->isa('Bio::Location::SplitLocationI'));
my @sublocs = $location->sub_Location();
is(@sublocs, 29);

# version and primary ID - believe it or not, this wasn't working
is ($seq->version, 1);
is ($seq->seq_version, 1);
is ($seq->primary_id, "5734104");

# streaming and Bio::RichSeq creation
my $stream = Bio::SeqIO->new(-file => test_input_file("test.genbank"),
                                      -verbose => $verbose,
                             -format => 'gbdriver');
$stream->verbose($verbose);
my $seqnum = 0;
my $species;
my @cl;
my $lasts;
my @ids = qw(DDU63596 DDU63595 HUMBDNF);
my @tids = (44689, 44689, 9606);
my @tnames = ("Dictyostelium discoideum","Dictyostelium discoideum",
                  "Homo sapiens");
while($seq = $stream->next_seq()) {
    if($seqnum < 3) {
        is $seq->display_id(), $ids[$seqnum];
        $species = $seq->species();
        @cl = $species->classification();
        is( $species->binomial(), $tnames[$seqnum],
             'species parsing incorrect for genbank');
        is( $cl[3] ne $species->genus(), 1,
             'genus duplicated in genbank parsing');
        is( $species->ncbi_taxid, $tids[$seqnum] );
    }
    $seqnum++;
    $lasts = $seq;
}
is($seqnum, 5,'streaming');
is $lasts->display_id(), "HUMBETGLOA";
my ($ref) = $lasts->annotation->get_Annotations('reference');
is($ref->medline, 94173918);
$stream->close();

$stream = Bio::SeqIO->new(-file => test_input_file("test.genbank.noseq"),
                                  -verbose => $verbose,
                                  -format => 'gbdriver' );
$seqnum = 0;
while($seq = $stream->next_seq()) {
    if($seqnum < 3) {
        is $seq->display_id(), $ids[$seqnum];
    } elsif( $seq->display_id eq 'M37762') {
        is( ($seq->get_keywords())[0], 'neurotrophic factor');
    }
    $seqnum++;
}
is $seqnum, 5, "Total number of sequences in test file";

# fuzzy
$seq = Bio::SeqIO->new( -format => 'gbdriver',
                                -verbose => $verbose,
                        -file =>test_input_file("testfuzzy.genbank"));
$seq->verbose($verbose);
ok(defined($as = $seq->next_seq()));

@features = $as->all_SeqFeatures();
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

@sublocs = $location->sub_Location();

is(@sublocs, 2);
my $loc = shift @sublocs;
is($loc->start, 83202);
is($loc->end, 83329);
is($loc->strand, -1);

$loc = shift @sublocs;
is($loc->start, 84248);
is($loc->end, 84996);
is($loc->strand,1);

my $outfile = test_output_file();
$seq = Bio::SeqIO->new(-format => 'genbank',
                              -verbose => $verbose,
                       -file=> ">$outfile");
$seq->verbose($verbose);
ok($seq->write_seq($as),'Fuzzy out');

## now genbank ##
$str = Bio::SeqIO->new(-format =>'gbdriver',
                             -verbose => $verbose,
                             -file => test_input_file('BK000016-tpa.gbk'));
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
my $testfile = test_output_file();
my $out = Bio::SeqIO->new(-file => ">$testfile",
                             -format => 'genbank');
$out->write_seq($seq);
$out->close();

$str = Bio::SeqIO->new(-format =>'gbdriver',
                             -file => $testfile);
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
my $gb = Bio::SeqIO->new(-format => 'gbdriver',
                        -verbose => $verbose,
                        -file   => test_input_file('revcomp_mrna.gb'));
$seq = $gb->next_seq();

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
$gb = Bio::SeqIO->new(-format => 'gbdriver',
                     -verbose => $verbose,
                        -file   => test_input_file('NC_006511-short.gbk'));
$seq = $gb->next_seq;
is $seq->species->common_name, undef, "Bug 1925";
is $seq->species->scientific_name, "Salmonella enterica subsp. enterica serovar Paratyphi A str. ATCC 9150";
@class = $seq->species->classification;
is $class[$#class], "Bacteria";

# WGS tests
$gb = Bio::SeqIO->new(-format => 'gbdriver',
                      -verbose => $verbose,
                    -file   => test_input_file('O_sat.wgs'));
$seq = $gb->next_seq;

my @tests = ('wgs'        => 'AAAA02000001-AAAA02050231',
            'wgs_scafld' => 'CM000126-CM000137',
            'wgs_scafld' => 'CH398081-CH401163');

my @wgs = map {$seq->annotation->get_Annotations(lc($_))} (qw(WGS WGS_SCAFLD));

my $ct=0;

for my $wgs (@wgs) {
    my ($tagname, $value) = (shift @tests, shift @tests);
    is($wgs->tagname, $tagname, $tagname);
    is($wgs->value, $value);
    $ct++;
}

is ($ct, 3);

# make sure we can retrieve a feature with a primary tag of 'misc_difference'
$gb = Bio::SeqIO->new(-format => 'gbdriver',
                     -verbose => $verbose,
                    -file   => test_input_file('BC000007.gbk'));
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


# output always adds a period (GenBank std), but two of these files do not use them.

foreach my $in ('BK000016-tpa.gbk', 'ay116458.gb', 'ay149291.gb', 'NC_006346.gb', 'ay007676.gb', 'dq519393.gb') {
    my $infile =  test_input_file($in);
	$outfile = test_output_file();
    
    $str = Bio::SeqIO->new(-format =>'genbank',
                          -verbose => $verbose,
                          -file => $infile);
    $seq = $str->next_seq;
    
    $out = Bio::SeqIO->new(-file => $outfile, -format => 'genbank');
    $out->write_seq($seq);
    $out->close();
    
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
    } continue { $line++ }
    close $RESULT;
    
    ok $is, $in;
}

# NB: there should probably be full testing on all lines to ensure that output
# matches input.

# 20061117: problem with *double* colon in some annotation-dblink values
$ct = 0;

foreach my $in ('P35527.gb') {
    my $infile =  test_input_file($in);
    $str = Bio::SeqIO->new(-format =>'genbank',
                         -verbose => $verbose,
                         -file => $infile);
    $seq = $str->next_seq;
    my $ac      = $seq->annotation();      # Bio::AnnotationCollection
    foreach my $key ($ac->get_all_annotation_keys() ) {
        my @values = $ac->get_Annotations($key);
        foreach my $ann (@values) {
            my $value = $ann->display_text;
            $ct++;
            if ($key eq 'dblink') {

                ok (index($value,'::') < 0);   # this should never be true

                ok ($value, $value);   # check value is not empty

                #  print "  ann/", sprintf('%12s  ',$key), '>>>', $value , '<<<', "\n";
                #  print "        index double colon: ",index($value   ,'::'), "\n";

                #  check db name:
                my @parts = split(/:/,$value);
                if ( $parts[0] =~ /^(?:
                        #  not an exhaustive list of databases;
                        #  just the db's referenced in P35527.gb:
                        swissprot | GenBank | GenPept  | HSSP| IntAct | Ensembl | KEGG | HGNC | MIM | ArrayExpress
                                  | GO      | InterPro | Pfam| PRINTS | PROSITE
                                     )$/x )
                {
                    ok 1;
                }
                else {
                    ok 0;
                }
                    ok ( $parts[1], "$parts[0]" );
            }
                # elsif ($key eq 'reference') { }
        }
    }
}

is($ct, 46);

# bug 2195
    
$str = Bio::SeqIO->new(-format =>'gbdriver',
                      -verbose => $verbose,
                      -file => test_input_file('AF305198.gb')
                     );

$species = $str->next_seq->species;

is($species->scientific_name, 'Virginia creeper phytoplasma', 'Bug 2195');
is(join(', ',$species->classification), 'Virginia creeper phytoplasma, '.
   '16SrV (Elm yellows group), Candidatus Phytoplasma, '.
   'Acholeplasmataceae, Acholeplasmatales, Mollicutes, '.
   'Firmicutes, Bacteria', 'Bug 2195');

# bug 2569, PROJECT line support, read and write, round-tripping
    
$str = Bio::SeqIO->new(-format =>'gbdriver',
                      -verbose => $verbose,
                      -file => test_input_file('NC_008536.gb'));

$seq = $str->next_seq;

my $project = ($seq->annotation->get_Annotations('project'))[0];

isa_ok($project, 'Bio::Annotation::SimpleValue');

if ($project) {
	is($project->value, 'GenomeProject:12638');
} else {
	ok(0, "PROJECT not parsed");
}

$outfile = test_output_file();

$gb = Bio::SeqIO->new(-format => 'genbank',
                              -verbose => $verbose,
                       -file=> ">$outfile");

$gb->write_seq($seq);

$str = Bio::SeqIO->new(-format =>'gbdriver',
                      -verbose => $verbose,
                      -file => $outfile);

$seq = $str->next_seq;

$project = ($seq->annotation->get_Annotations('project'))[0];

isa_ok($project, 'Bio::Annotation::SimpleValue');

if ($project) {
	is($project->value, 'GenomeProject:12638');
} else {
	ok(0, "Roundtrip test failed");
}

################################## EMBL ##################################

# Set to -1 for release version, so warnings aren't printed

$ast = Bio::SeqIO->new( -format => 'embldriver',
			   -verbose => $verbose,
			   -file => test_input_file("roa1.dat"));
$ast->verbose($verbose);
$as = $ast->next_seq();
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
is($as->get_dates, 2);

# EMBL Release 87 changes (8-17-06)

$ast = Bio::SeqIO->new( -format => 'embldriver',
			   -verbose => $verbose,
			   -file => test_input_file("roa1_v2.dat"));
$ast->verbose($verbose);
$as = $ast->next_seq();
ok defined $as->seq;
# accession # same as display name now
is($as->display_id, 'X79536');
is($as->get_dates, 2);
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

my $ent = Bio::SeqIO->new( -file => test_input_file("test.embl"),
			   -format => 'embldriver');
$seq = $ent->next_seq();

is(defined $seq->seq(), 1,
   'success reading Embl with ^ location and badly split double quotes');
is(scalar $seq->annotation->get_Annotations('reference'), 3);
is($seq->get_dates, 0);

$out = Bio::SeqIO->new(-file=> ">$outfile",
			  -format => 'embl');
is($out->write_seq($seq),1,
   'success writing Embl format with ^ < and > locations');

# embl with no FT
$ent = Bio::SeqIO->new( -file => test_input_file("test.embl"),
			-format => 'embldriver');
$seq = $ent->next_seq();

ok($seq);
is(lc($seq->subseq(1,10)),'gatcagtaga');
is($seq->length, 4870);

# embl with no FH
my $noFH = Bio::SeqIO->new(-file => test_input_file("no_FH.embl"),
			-format => 'embldriver');
$seq = $noFH->next_seq;
is(scalar($seq->get_SeqFeatures), 4);
is($seq->display_id, 'AE000001');
is($seq->get_dates, 0);

# bug 1571
$ent = Bio::SeqIO->new(-format => 'embldriver',
		       -file   => test_input_file('test.embl2sq'));
$seq = $ent->next_seq;
is($seq->length,4870);
is($seq->get_dates, 0);

# embl repbase
$ent = Bio::SeqIO->new(-file => test_input_file("BEL16-LTR_AG.embl"), -format => 'embldriver');
$seq = $ent->next_seq;
is($seq->display_id,'BEL16-LTR_AG');
is($seq->get_dates, 2);

# test secondary accessions in EMBL (bug #1332)
$seqio = Bio::SeqIO->new(-format => 'embldriver',
			   -file => test_input_file('ECAPAH02.embl'));
$seq = $seqio->next_seq;
is($seq->accession_number, 'D10483');
is($seq->seq_version, 2);
my @accs = $seq->get_secondary_accessions();
is($accs[0], 'J01597');
is($accs[-1], 'X56742');
is($seq->get_dates, 2);

### TPA TESTS - Thanks to Richard Adams ###
# test Third Party Annotation entries in EMBL/Gb format 
# to ensure compatability with parsers.
$str = Bio::SeqIO->new(-verbose => $verbose,
                         -format =>'embldriver',
			 -file => test_input_file('BN000066-tpa.embl'));
$seq = $str->next_seq;
ok(defined $seq);
is($seq->accession_number, 'BN000066');
is($seq->alphabet, 'dna');
is($seq->display_id, 'AGA000066');
is($seq->length, 5195);
is($seq->division, 'INV');
is($seq->keywords, 'acetylcholinesterase; achE1 gene; Third Party Annotation; TPA');
is($seq->seq_version, 1);
is($seq->feature_count, 15);
is($seq->get_dates, 2);

$spec_obj = $seq->species;
is ($spec_obj->common_name, 'African malaria mosquito');
is ($spec_obj->species, 'gambiae');
is ($spec_obj->genus, 'Anopheles');
is ($spec_obj->binomial, 'Anopheles gambiae');

$ac = $seq->annotation;
$reference =  ($ac->get_Annotations('reference') )[1];
is ($reference->title,'"A novel acetylcholinesterase gene in mosquitoes codes for the insecticide target and is non-homologous to the ace gene in Drosophila"');
is ($reference->authors,'Weill M., Fort P., Berthomi eu A., Dubois M.P., Pasteur N., Raymond M.');
my $cmmnt =  ($ac->get_Annotations('comment') )[0];
is($cmmnt->text, 'see also AJ488492 for achE-1 from Kisumu strain Third Party Annotation Database: This TPA record uses Anopheles gambiae trace archive data (http://trace.ensembl.org)');


$ent = Bio::SeqIO->new( -file => test_input_file("test.embl"),
                        -format => 'embldriver');
$ent->verbose($verbose);
$seq = $ent->next_seq();
$species = $seq->species();
@cl = $species->classification();
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
                          -file => ">$outfile");

ok($embl->write_seq($primaryseq));

# this should generate a warning
my $scalar = "test";
eval {
	$embl->write_seq($scalar);
};
ok ($@);

############################## Swiss/UniProt ##############################

$seqio = Bio::SeqIO->new( -verbose => $verbose,
                                     -format => 'swissdriver',
                                     -file   => test_input_file('test.swiss'));

isa_ok($seqio, 'Bio::SeqIO');
$seq = $seqio->next_seq;
my @gns = $seq->annotation->get_Annotations('gene_name');

$outfile = test_output_file();
$seqio = Bio::SeqIO->new( -verbose => $verbose,
                                 -format => 'swiss',
                                 -file   => ">$outfile");

$seqio->write_seq($seq);

# reads it in once again
$seqio = Bio::SeqIO->new( -verbose => $verbose,
                                 -format => 'swissdriver',
                                 -file => $outfile);
    
$seq = $seqio->next_seq;
isa_ok($seq->species, 'Bio::Species');
is($seq->species->ncbi_taxid, 6239);

# version, seq_update, dates (5 tests)
is($seq->version, 40);
my ($ann) = $seq->annotation->get_Annotations('seq_update');
eval {is($ann->display_text, 35,'operator overloading in AnnotationI is deprecated')};
ok(!$@);

my @dates = $seq->get_dates;
my @date_check = qw(01-NOV-1997 01-NOV-1997 16-OCT-2001);

for my $date (@dates) {
    my $expdate = shift @date_check;
    is($date, $expdate,'dates');
}

my @gns2 = $seq->annotation->get_Annotations('gene_name');
# check gene name is preserved (was losing suffix in worm gene names)
ok($#gns2 == 0 && $gns[0]->value eq $gns2[0]->value);

# test swissprot multiple RP lines
$str = Bio::SeqIO->new(-file => test_input_file('P33897'));
$seq = $str->next_seq;
isa_ok($seq, 'Bio::Seq::RichSeqI');
@refs = $seq->annotation->get_Annotations('reference');
is( @refs, 23);
is($refs[20]->rp, 'VARIANTS X-ALD LEU-98; ASP-99; GLU-217; GLN-518; ASP-608; ILE-633 AND PRO-660, AND VARIANT THR-13.');

# version, seq_update, dates (5 tests)
is($seq->version, 44);
($ann) = $seq->annotation->get_Annotations('seq_update');
is($ann->display_text, 28,'operator overloading in AnnotationI is deprecated');
@dates = $seq->get_dates;
@date_check = qw(01-FEB-1994 01-FEB-1994 15-JUN-2004);
for my $date (@dates) {
    is($date, shift @date_check);
}

$ast = Bio::SeqIO->new(-verbose => $verbose,
                                  -format => 'swissdriver' ,
                                  -file => test_input_file("roa1.swiss"));
$as = $ast->next_seq();

ok defined $as->seq;
is($as->id, 'ROA1_HUMAN', "id is ".$as->id);
like($as->primary_id, qr(Bio::PrimarySeq));
is($as->length, 371);
is($as->alphabet, 'protein');
is($as->division, 'HUMAN');
is(scalar $as->all_SeqFeatures(), 16);
is(scalar $as->annotation->get_Annotations('reference'), 11);

# version, seq_update, dates (5 tests)
is($as->version, 35);
($ann) = $as->annotation->get_Annotations('seq_update');
is($ann->display_text, 15,'operator overloading in AnnotationI is deprecated');
@dates = $as->get_dates;
@date_check = qw(01-MAR-1989 01-AUG-1990 01-NOV-1997);
for my $date (@dates) {
    is($date, shift @date_check);
}

($ent,$out) = undef;
($as,$seq) = undef;

$seqio = Bio::SeqIO->new(-format => 'swissdriver' ,
                                 -verbose => $verbose,
                                 -file => test_input_file("swiss.dat"));
$seq = $seqio->next_seq;
isa_ok($seq, 'Bio::Seq::RichSeqI');

# more tests to verify we are actually parsing correctly
like($seq->primary_id, qr(Bio::PrimarySeq));
is($seq->display_id, 'MA32_HUMAN');
is($seq->length, 282);
is($seq->division, 'HUMAN');
is($seq->alphabet, 'protein');
my @f = $seq->all_SeqFeatures();
is(@f, 2);
is($f[1]->primary_tag, 'CHAIN');
is(($f[1]->get_tag_values('description'))[0], 'COMPLEMENT COMPONENT 1, Q SUBCOMPONENT BINDING PROTEIN');

# version, seq_update, dates (5 tests)
is($seq->version, 40);
($ann) = $seq->annotation->get_Annotations('seq_update');
is($ann->display_text, 31,'operator overloading in AnnotationI is deprecated');
@dates = $seq->get_dates;
@date_check = qw(01-FEB-1995 01-FEB-1995 01-OCT-2000);
for my $date (@dates) {
    is($date, shift @date_check);
}

my @genenames = qw(GC1QBP HABP1 SF2P32 C1QBP);
($ann) = $seq->annotation->get_Annotations('gene_name');
# use Data::Stag findval and element name to get values/nodes
foreach my $gn ( $ann->findval('Name') ) {
    ok ($gn, shift(@genenames));
}
foreach my $gn ( $ann->findval('Synonyms') ) {
    ok ($gn, shift(@genenames));
}
like($ann->value, qr/Name: GC1QBP/);

# test for feature locations like ?..N
$seq = $seqio->next_seq();
isa_ok($seq, 'Bio::Seq::RichSeqI');
like($seq->primary_id, qr(Bio::PrimarySeq));
is($seq->display_id, 'ACON_CAEEL');
is($seq->length, 788);
is($seq->division, 'CAEEL');
is($seq->alphabet, 'protein');
is(scalar $seq->all_SeqFeatures(), 5);

foreach my $gn ( $seq->annotation->get_Annotations('gene_name') ) {
    ok ($gn->value, 'F54H12.1');
}

# test species in swissprot -- this can be a n:n nightmare
$seq = $seqio->next_seq();
isa_ok($seq, 'Bio::Seq::RichSeqI');
like($seq->primary_id, qr(Bio::PrimarySeq));
@sec_acc = $seq->get_secondary_accessions();
is($sec_acc[0], 'P29360');
is($sec_acc[1], 'Q63631');
is($seq->accession_number, 'P42655');
@kw = $seq->get_keywords;
is( $kw[0], 'Brain');
is( $kw[1], 'Neurone');
is($kw[3], 'Multigene family');
is($seq->display_id, '143E_HUMAN');

# hybrid names from old sequences are no longer valid, these are chopped
# off at the first organism
is($seq->species->binomial, "Homo sapiens");
is($seq->species->common_name, "Human");
is($seq->species->ncbi_taxid, 9606);

$seq = $seqio->next_seq();
isa_ok($seq, 'Bio::Seq::RichSeqI');
like($seq->primary_id, qr(Bio::PrimarySeq));
is($seq->species->binomial, "Bos taurus");
is($seq->species->common_name, "Bovine");
is($seq->species->ncbi_taxid, 9913);

# multiple genes in swissprot
$seq = $seqio->next_seq();
isa_ok($seq, 'Bio::Seq::RichSeqI');
like($seq->primary_id, qr(Bio::PrimarySeq));

($ann) = $seq->annotation->get_Annotations("gene_name");
@genenames = qw(CALM1 CAM1 CALM CAM CALM2 CAM2 CAMB CALM3 CAM3 CAMC);
my $flatnames = "(CALM1 OR CAM1 OR CALM OR CAM) AND (CALM2 OR CAM2 OR CAMB) AND (CALM3 OR CAM3 OR CAMC)";

my @names = @genenames; # copy array

my @ann_names = $ann->get_all_values();
is(scalar(@ann_names), scalar(@names));

# do this in a layered way (nested tags)
for my $node ($ann->findnode('gene_name')) {
    for my $name ($node->findval('Name')) {
        is($name, shift(@names));
    }
    for my $name ($node->findval('Synonyms')) {
        is($name, shift(@names));
    }
}

is(scalar(@names),0);

# same entry as before, but with the new gene names format
$seqio = Bio::SeqIO->new(-format => 'swissdriver',
                                 -verbose => $verbose,
                         -file => test_input_file("calm.swiss"));
$seq = $seqio->next_seq();
isa_ok($seq, 'Bio::Seq::RichSeqI');
like($seq->primary_id, qr(Bio::PrimarySeq));

($ann) = $seq->annotation->get_Annotations("gene_name");
@names = @genenames; # copy array

my @ann_names2 = $ann->get_all_values(); #emulate StructuredValue's flattened array
is(scalar(@ann_names2), scalar(@names));

for my $node ($ann->findnode('gene_name')) {
    for my $name ($node->findval('Name')) {
        is($name, shift(@names));
    }
    for my $name ($node->findval('Synonyms')) {
        is($name, shift(@names));
    }
}

is(scalar(@names),0);

# test proper parsing of references
my @litrefs = $seq->annotation->get_Annotations('reference');
is(scalar(@litrefs), 17);

my @titles = (
    '"Complete amino acid sequence of human brain calmodulin."',
    '"Multiple divergent mRNAs code for a single human calmodulin."',
    '"Molecular analysis of human and rat calmodulin complementary DNA clones. Evidence for additional active genes in these species."',
    '"Isolation and nucleotide sequence of a cDNA encoding human calmodulin."',
    '"Structure of the human CALM1 calmodulin gene and identification of two CALM1-related pseudogenes CALM1P1 and CALM1P2."',
    undef,
    '"Characterization of the human CALM2 calmodulin gene and comparison of the transcriptional activity of CALM1, CALM2 and CALM3."',
    '"Cloning of human full-length CDSs in BD Creator(TM) system donor vector."',
    '"The DNA sequence and analysis of human chromosome 14."',
    '"Generation and initial analysis of more than 15,000 full-length human and mouse cDNA sequences."',
    '"Alpha-helix nucleation by a calcium-binding peptide loop."',
    '"Solution structure of Ca(2+)-calmodulin reveals flexible hand-like properties of its domains."',
    '"Calmodulin structure refined at 1.7 A resolution."',
    '"Drug binding by calmodulin: crystal structure of a calmodulin-trifluoperazine complex."',
    '"Structural basis for the activation of anthrax adenylyl cyclase exotoxin by calmodulin."',
    '"Physiological calcium concentrations regulate calmodulin binding and catalysis of adenylyl cyclase exotoxins."',
    '"Crystal structure of a MARCKS peptide containing the calmodulin-binding domain in complex with Ca2+-calmodulin."',
);

my @locs = (
    "Biochemistry 21:2565-2569(1982).",
    "J. Biol. Chem. 263:17055-17062(1988).",
    "J. Biol. Chem. 262:16663-16670(1987).",
    "Biochem. Int. 9:177-185(1984).",
    "Eur. J. Biochem. 225:71-82(1994).",
    "Submitted (FEB-1995) to the EMBL/GenBank/DDBJ databases.",
    "Cell Calcium 23:323-338(1998).",
    "Submitted (MAY-2003) to the EMBL/GenBank/DDBJ databases.",
    "Nature 421:601-607(2003).",
    "Proc. Natl. Acad. Sci. U.S.A. 99:16899-16903(2002).",
    "Proc. Natl. Acad. Sci. U.S.A. 96:903-908(1999).",
    "Nat. Struct. Biol. 8:990-997(2001).",
    "J. Mol. Biol. 228:1177-1192(1992).",
    "Biochemistry 33:15259-15265(1994).",
    "Nature 415:396-402(2002).",
    "EMBO J. 21:6721-6732(2002).",
    "Nat. Struct. Biol. 10:226-231(2003).",
);

my @positions = (
     undef, undef,
    undef, undef,
    undef, undef,
    undef, undef,
    undef, undef,
    undef, undef,
    undef, undef,
    undef, undef,
    undef, undef,
    undef, undef,
    94, 103,
    1, 76,
    undef, undef,
    undef, undef,
    5, 148,
    1, 148,
    undef, undef,
);

foreach my $litref (@litrefs) {
    is($litref->title, shift(@titles));
    is($litref->location, shift(@locs));
    is($litref->start, shift(@positions));
    is($litref->end, shift(@positions));
}

# format parsing changes (pre-rel 9.0)

$seqio = Bio::SeqIO->new( -verbose => $verbose,
                         -format => 'swissdriver',
                         -file   => test_input_file('pre_rel9.swiss'));

ok($seqio);
$seq = $seqio->next_seq;
isa_ok($seq->species, 'Bio::Species');
is($seq->species->ncbi_taxid, "6239");

# version, seq_update, dates (5 tests)
is($seq->version, 44);
($ann) = $seq->annotation->get_Annotations('seq_update');
is($ann->display_text, 1);
@dates = $seq->get_dates;
@date_check = qw(01-NOV-1997 01-NOV-1996 30-MAY-2006 );
for my $date (@dates) {
    is($date, shift @date_check);
}

my @idcheck = qw(Z66513 T22647 Cel.30446 Q06319 Q20772 F54D5.7 WBGene00010052
		 F54D5.7 GO:0005515 IPR006089 IPR006091 IPR006090
		 IPR006092 IPR009075 IPR009100 IPR013764 PF00441
		 PF02770 PF02771 PS00072 PS00073);

for my $dblink ( $seq->annotation->get_Annotations('dblink') ) {
    is($dblink->primary_id, shift @idcheck);
}

$seqio = Bio::SeqIO->new( -verbose => $verbose,
                         -format => 'swissdriver',
                         -file   => test_input_file('pre_rel9.swiss'));

my @namespaces = qw(Swiss-Prot TrEMBL TrEMBL);

while (my $seq = $seqio->next_seq) {
    is($seq->namespace, shift @namespaces);
}

# format parsing changes (rel 9.0, Oct 2006)

$seqio = Bio::SeqIO->new( -verbose => $verbose,
                         -format => 'swissdriver',
                         -file   => test_input_file('rel9.swiss'));

ok($seqio);
$seq = $seqio->next_seq;
isa_ok($seq->species, 'Bio::Species');
is($seq->species->ncbi_taxid, 6239);

is($seq->version, 47);
($ann) = $seq->annotation->get_Annotations('seq_update');
is($ann->display_text, 1,'operator overloading in AnnotationI is deprecated');
@dates = $seq->get_dates;
@date_check = qw(01-NOV-1997 01-NOV-1996 31-OCT-2006 );
for my $date (@dates) {
    is($date, shift @date_check);
}

@idcheck = qw(Z66513 T22647 Cel.30446 Q06319 Q20772 F54D5.7 cel:F54D5.7
         WBGene00010052 F54D5.7 GO:0005515 IPR006089 IPR006091 IPR006090
         IPR006092 IPR009075 IPR013786 IPR009100 IPR013764 PF00441 PF02770
         PF02771 PS00072 PS00073 );

for my $dblink ( $seq->annotation->get_Annotations('dblink') ) {
    is($dblink->primary_id, shift @idcheck);
}

$seqio = Bio::SeqIO->new( -verbose => $verbose,
                         -format => 'swissdriver',
                         -file   => test_input_file('rel9.swiss'));

@namespaces = qw(Swiss-Prot TrEMBL TrEMBL);

while (my $seq = $seqio->next_seq) {
    is($seq->namespace, shift @namespaces);
}

# bug 2288
# Q8GBD3.swiss
$seqio = Bio::SeqIO->new( -verbose => $verbose,
                         -format => 'swiss',
                         -file   => test_input_file('Q8GBD3.swiss'));

while (my $seq = $seqio->next_seq) {
    my $lineage = join(';', $seq->species->classification);
	is ($lineage, 'Acetobacter aceti;Acetobacter subgen. Acetobacter;'.
		'Acetobacter;Acetobacteraceae;Rhodospirillales;Alphaproteobacteria;'.
		'Proteobacteria;Bacteria');
}

# test for GenBank swissprot/UniProt/UniProtKB DBSOURCE line (Bug : RT 44536)
$ast = Bio::SeqIO->new(-format => 'genbank',
                              -verbose => $verbose,
                       -file => test_input_file('P39765.gb'));
$ast->verbose($verbose);
$as = $ast->next_seq();
is $as->molecule, 'PRT',$as->accession_number;;
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
