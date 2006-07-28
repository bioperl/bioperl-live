# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

use strict;

BEGIN {
	eval { require Test; };
	if ( $@ ) {
		use lib 't';
	}
	use Test;
	plan tests => 126;
}

use Bio::SeqIO;
use Bio::Root::IO;

ok(1);

my $verbose = $ENV{'BIOPERLDEBUG'} || 0;

END {
	unlink "tmp_revcomp_mrna.gb" if -e "tmp_revcomp_mrna.gb";
	unlink 'testsource.gb';
}

my $ast = Bio::SeqIO->new(-format => 'GenBank' ,
								  -verbose => $verbose,
								  -file => Bio::Root::IO->catfile
								  ("t","data","roa1.genbank"));
$ast->verbose($verbose);
my $as = $ast->next_seq();
ok $as->molecule, 'mRNA';
ok $as->alphabet, 'dna';
ok($as->primary_id, 3598416);
my @class = $as->species->classification;
ok $class[$#class],'Eukaryota';

$ast = Bio::SeqIO->new(-format => 'genbank',
							  -verbose => $verbose,
                       -file => Bio::Root::IO->catfile
							  ("t","data","NT_021877.gbk"));
$ast->verbose($verbose);
$as = $ast->next_seq();
ok $as->molecule, 'DNA';
ok $as->alphabet, 'dna';
ok($as->primary_id, 37539616);
ok($as->accession_number, 'NT_021877');

my ($cds) = grep { $_->primary_tag eq 'CDS' } $as->get_SeqFeatures();
ok(($cds->get_tag_values('transl_except'))[1],
   '(pos:complement(4224..4226),aa:OTHER)');

# test for a DBSOURCE line
$ast = Bio::SeqIO->new(-format => 'genbank',
							  -verbose => $verbose,
                       -file => Bio::Root::IO->catfile("t","data","BAB68554.gb"));
$ast->verbose($verbose);
$as = $ast->next_seq();
ok $as->molecule, 'linear';
ok $as->alphabet, 'protein';
ok $as->species->common_name, 'Aldabra giant tortoise';
ok($as->primary_id, 15824047);
my $ac = $as->annotation;
ok defined $ac;
my @dblinks = $ac->get_Annotations('dblink');
ok(scalar @dblinks,1);

# test for multi-line SOURCE
$ast = Bio::SeqIO->new(-format => 'genbank',
							  -verbose => $verbose,
                       -file => Bio::Root::IO->catfile("t","data",
                                                       "NC_006346.gb"));
$as = $ast->next_seq;
ok $as->species->binomial,'Bolitoglossa n. sp. RLM-2004';
@class = $as->species->classification;
ok($class[$#class],'Eukaryota');
ok($as->species->common_name,'mitochondrion Bolitoglossa n. sp. RLM-2004 (mushroomtongue salamander)');

$ast = Bio::SeqIO->new(-format => 'genbank',
							  -verbose => $verbose,
                       -file => Bio::Root::IO->catfile("t","data",
                                                       "U71225.gb"));
$as = $ast->next_seq;
@class = $as->species->classification;
ok($class[$#class],'Eukaryota');
ok $as->species->common_name,'mitochondrion Desmognathus quadramaculatus (black-bellied salamander)';

# test for unusual common name
$ast = Bio::SeqIO->new(-format => 'genbank',
							  -verbose => $verbose,
                       -file => Bio::Root::IO->catfile("t","data",
                                                       "AB077698.gb"));
$as = $ast->next_seq;
ok $as->species->common_name,'Homo sapiens cDNA to mRNA';

# test for common name with parentheses
$ast = Bio::SeqIO->new(-format => 'genbank',
							  -verbose => $verbose,
                       -file => Bio::Root::IO->catfile("t","data",
                                                       "DQ018368.gb"));
$as = $ast->next_seq;
ok $as->species->common_name,'(Populus tomentosa x P. bolleana) x P. tomentosa var. truncata';

# test secondary accessions
my $seqio = new Bio::SeqIO(-format => 'genbank',
									-verbose => $verbose,
									-file => Bio::Root::IO->catfile
									(qw(t data D10483.gbk)));
my $seq = $seqio->next_seq;
my @kw =  $seq->get_keywords;
ok(scalar @kw, 118);
ok($kw[-1], 'yabO');
my @sec_acc = $seq->get_secondary_accessions();
ok(scalar @sec_acc,14);
ok($sec_acc[-1], 'X56742');

# bug #1487
my $str = new Bio::SeqIO(-verbose => $verbose,
								 -file    => Bio::Root::IO->catfile
								 (qw(t data D12555.gbk)));
eval {
	$seq = $str->next_seq;
};

ok(! $@ );

# bug 1647 rpt_unit sub-feature with multiple parens
$str = Bio::SeqIO->new(-format => 'genbank',
							  -verbose => $verbose,
                       -file => Bio::Root::IO->catfile
							  (qw(t data mini-AE001405.gb) ));
ok($seq = $str->next_seq);
my @rpts = grep { $_->primary_tag eq 'repeat_region' }
  $seq->get_SeqFeatures;
ok $#rpts, 2;
my @rpt_units = map {$_->get_tag_values('rpt_unit')} @rpts;
ok $#rpt_units, 0;
ok $rpt_units[0],'(TG)10;A;(TG)7';

# test bug #1673 , RDB-II genbank files
$str = Bio::SeqIO->new(-format => 'genbank',
							  -verbose => $verbose,
                       -file => Bio::Root::IO->catfile
							  (qw(t data Mcjanrna_rdbII.gbk) )
		      );
ok($seq = $str->next_seq);
my @refs = $seq->annotation->get_Annotations('reference');
ok(@refs, 1);
ok($seq->display_id,'Mc.janrrnA');
ok($seq->molecule ,'RNA');

$str  = new Bio::SeqIO(-format => 'genbank',
							  -file   => Bio::Root::IO->catfile
							  ("t","data","AF165282.gb"),
							  -verbose => $verbose);
$seq = $str->next_seq;
my @features = $seq->all_SeqFeatures();
ok(@features, 5);
ok($features[0]->start, 1);
ok($features[0]->end, 226);
my $location = $features[1]->location;
ok($location->isa('Bio::Location::SplitLocationI'));
my @sublocs = $location->sub_Location();
ok(@sublocs, 29);

# version and primary ID - believe it or not, this wasn't working
ok ($seq->version, 1);
ok ($seq->seq_version, 1);
ok ($seq->primary_id, "5734104");

# streaming and Bio::RichSeq creation
my $stream = Bio::SeqIO->new(-file => Bio::Root::IO->catfile
									  ("t","data","test.genbank"),
									  -verbose => $verbose,
                             -format => 'GenBank');
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
		ok $seq->display_id(), $ids[$seqnum];
		$species = $seq->species();
		@cl = $species->classification();
		ok( $species->binomial(), $tnames[$seqnum],
			 'species parsing incorrect for genbank');
		ok( $cl[3] ne $species->genus(), 1,
			 'genus duplicated in genbank parsing');
		ok( $species->ncbi_taxid, $tids[$seqnum] );
	}
	$seqnum++;
	$lasts = $seq;
}
ok $lasts->display_id(), "HUMBETGLOA";
my ($ref) = $lasts->annotation->get_Annotations('reference');
ok($ref->medline, 94173918);
$stream->close();

$stream = Bio::SeqIO->new(-file => Bio::Root::IO->catfile
								  ("t","data","test.genbank.noseq"),
								  -verbose => $verbose,
								  -format => 'GenBank' );
$seqnum = 0;
while($seq = $stream->next_seq()) {
	if($seqnum < 3) {
		ok $seq->display_id(), $ids[$seqnum];
	} elsif( $seq->display_id eq 'M37762') {
		ok( ($seq->get_keywords())[0], 'neurotrophic factor');
	}
	$seqnum++;
}
ok $seqnum, 5, "Total number of sequences in test file";

# fuzzy
$seq = Bio::SeqIO->new( -format => 'GenBank',
								-verbose => $verbose,
                        -file =>Bio::Root::IO->catfile
								("t","data","testfuzzy.genbank"));
$seq->verbose($verbose);
ok(defined($as = $seq->next_seq()));

@features = $as->all_SeqFeatures();
ok(@features,21);
my $lastfeature = pop @features;
# this is a split location; the root doesn't have strand
ok($lastfeature->strand, undef);
$location = $lastfeature->location;
$location->verbose(-1); # silence the warning of undef seq_id()
# see above; splitlocs roots do not have a strand really
ok($location->strand, undef);
ok($location->start, 83202);
ok($location->end, 84996);

@sublocs = $location->sub_Location();

ok(@sublocs, 2);
my $loc = shift @sublocs;
ok($loc->start, 83202);
ok($loc->end, 83329);
ok($loc->strand, -1);

$loc = shift @sublocs;
ok($loc->start, 84248);
ok($loc->end, 84996);
ok($loc->strand,1);

$seq = Bio::SeqIO->new(-format => 'GenBank',
							  -verbose => $verbose,
                       -file=> ">" . Bio::Root::IO->catfile
							  ("t","data","genbank.fuzzyout"));
$seq->verbose($verbose);
ok($seq->write_seq($as));
unlink(Bio::Root::IO->catfile("t","data","genbank.fuzzyout"));

## now genbank ##
$str = new Bio::SeqIO(-format =>'genbank',
							 -verbose => $verbose,
							 -file => Bio::Root::IO->catfile
							 ( qw(t data BK000016-tpa.gbk)));
$seq = $str->next_seq;
ok(defined $seq);
ok(defined $seq->seq);
ok($seq->accession_number, 'BK000016');
ok($seq->alphabet, 'dna');
ok($seq->display_id, 'BK000016');
ok($seq->length, 1162);
ok($seq->division, 'ROD');
ok($seq->get_dates, 1);
ok($seq->keywords, 'Third Party Annotation; TPA');
ok($seq->desc, 'TPA: Mus musculus pantothenate kinase 4 mRNA, partial cds.');
ok($seq->seq_version, 1);
ok($seq->feature_count, 2);
my $spec_obj = $seq->species;
ok ($spec_obj->common_name, 'Mus musculus (house mouse)');
ok ($spec_obj->species, 'musculus');
ok ($spec_obj->genus, 'Mus');
ok ($spec_obj->binomial, 'Mus musculus');
$ac = $seq->annotation;
my $reference =  ($ac->get_Annotations('reference') )[0];
ok ($reference->pubmed, '11479594');
ok ($reference->medline, '21372465');

# validate that what is written is what is read
my $testfile = "testtpa.gbk";
my $out = new Bio::SeqIO(-file => ">$testfile",
							 -format => 'genbank');
$out->write_seq($seq);
$out->close();

$str = new Bio::SeqIO(-format =>'genbank',
							 -file => $testfile);
$seq = $str->next_seq;
ok(defined $seq);
ok(defined $seq->seq);
ok($seq->accession_number, 'BK000016');
ok($seq->alphabet, 'dna');
ok($seq->display_id, 'BK000016');
ok($seq->length, 1162);
ok($seq->division, 'ROD');
ok($seq->get_dates, 1);
ok($seq->keywords, 'Third Party Annotation; TPA');
ok($seq->desc, 'TPA: Mus musculus pantothenate kinase 4 mRNA, partial cds.');
ok($seq->seq_version, 1);
ok($seq->feature_count, 2);
$spec_obj = $seq->species;
ok ($spec_obj->common_name, 'Mus musculus (house mouse)');
ok ($spec_obj->species, 'musculus');
ok ($spec_obj->genus, 'Mus');
ok ($spec_obj->binomial, 'Mus musculus');
$ac = $seq->annotation;
$reference =  ($ac->get_Annotations('reference') )[0];
ok ($reference->pubmed, '11479594');
ok ($reference->medline, '21372465');

unlink($testfile);

# write revcomp split location
my $gb = new Bio::SeqIO(-format => 'genbank',
                        -file   => Bio::Root::IO->catfile
                        (qw(t data revcomp_mrna.gb)));
$seq = $gb->next_seq();

$gb = new Bio::SeqIO(-format => 'genbank',
                     -file   => ">tmp_revcomp_mrna.gb");

$gb->write_seq($seq);
undef $gb;
ok(! -z "tmp_revcomp_mrna.gb");
# INSERT DIFFING CODE HERE

# bug 1925, continuation of long ORGANISM line ends up in @classification:
# ORGANISM  Salmonella enterica subsp. enterica serovar Paratyphi A str. ATCC
#           9150
#           Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales;
#           Enterobacteriaceae; Salmonella.
$gb = new Bio::SeqIO(-format => 'genbank',
							-file   => Bio::Root::IO->catfile
							(qw(t data NC_006511-short.gbk)));
$seq = $gb->next_seq;
ok $seq->species->common_name, "Salmonella enterica subsp. enterica serovar Paratyphi A str. ATCC 9150";
@class = $seq->species->classification;
ok $class[$#class], "Bacteria";

# WGS tests
$gb = new Bio::SeqIO(-format => 'genbank',
							-file   => Bio::Root::IO->catfile
							(qw(t data O_sat.wgs)));
$seq = $gb->next_seq;

my @tests = ('WGS'        => 'AAAA02000001-AAAA02050231',
				  'WGS_SCAFLD' => 'CM000126-CM000137',
				  'WGS_SCAFLD' => 'CH398081-CH401163');

foreach my $wgs
(map {$seq->annotation->get_Annotations($_)} qw(WGS WGS_SCAFLD)) {
    my ($tagname, $value) = (shift @tests, shift @tests);
	ok($wgs->tagname, $tagname);
	ok($wgs->value, $value);
}

# make sure we can retrieve a feature with a primary tag of 'misc_difference'
$gb = new Bio::SeqIO(-format => 'genbank',
							-file   => Bio::Root::IO->catfile
							(qw(t data BC000007.gbk)));
$seq = $gb->next_seq;
($cds) = grep { $_->primary_tag eq 'misc_difference' } $seq->get_SeqFeatures;
my @vals = $cds->get_tag_values('gene');
ok $vals[0], 'PX19';

# Check that the source,organism section is identical between input and output.
# - test an easy one where organism is species, then two different formats of
# subspecies, and finally a species with a format that used to be mistaken for
# subspecies.
my $outfile = 'testsource.gb';
foreach my $infile ('BK000016-tpa.gbk', 'ay116458.gb', 'ay149291.gb', 'NC_006346.gb') {
	$str = new Bio::SeqIO(-format =>'genbank',
						  -verbose => $verbose,
						  -file => Bio::Root::IO->catfile("t/data/$infile"));
	$seq = $str->next_seq;
	
	$out = new Bio::SeqIO(-file => ">$outfile", -format => 'genbank');
	$out->write_seq($seq);
	$out->close();
	
	open (IN, "t/data/$infile");
	my @in = <IN>;
	close(IN);
	open (RESULT, "$outfile");
	my $line = 0;
	my $check = 0;
	my $ok = 1;
	while (<RESULT>) {
		if (/^KEYWORDS/) {
			$check = 1;
			next;
		}
		
		if ($check) {
			if ($_ ne $in[$line]) {
				$ok = 0;
				last;
			}
		}
		
		if (/^REFERENCE/) {
			last;
		}
	} continue { $line++ }
	close(RESULT);
	
	ok $ok;
	
	unlink($outfile);
}

# NB: there should probably be full testing on all lines to ensure that output
# matches input.
