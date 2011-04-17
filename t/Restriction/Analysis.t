# -*-Perl-*- Test Harness script for Bioperl
# $Id$


use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 182);
	
    use_ok('Bio::Restriction::Enzyme');
    use_ok('Bio::Restriction::Enzyme::MultiCut');
    use_ok('Bio::Restriction::Enzyme::MultiSite');
    use_ok('Bio::Restriction::EnzymeCollection');
    use_ok('Bio::Restriction::Analysis');
    use_ok('Bio::SeqIO');
}

#
# Bio::Restriction::Enzyme
#

my ($re, $seq, $iso, %meth, $microbe, $source, @vendors, @refs, $name);
ok $re=Bio::Restriction::Enzyme->new(-enzyme=>'EcoRI', -site=>'G^AATTC');
isa_ok($re, 'Bio::Restriction::EnzymeI');
is $re->cut, 1;
ok ! $re->cut(0);
is $re->complementary_cut, 6;
ok $re->cut(1);

is $re->complementary_cut,5;
is $re->site,'G^AATTC';
ok $seq = $re->seq;
isa_ok($seq, 'Bio::PrimarySeqI');
is $seq->seq, 'GAATTC';
is $re->string,'GAATTC';
is $re->revcom, 'GAATTC';
is $re->recognition_length, 6;
is $re->cutter, 6;
is $re->palindromic, 1;
is $re->overhang, "5'";
is $re->overhang_seq, 'AATT';
is $re->is_ambiguous, 0;

ok $re->compatible_ends($re);

ok $re->isoschizomers('BamHI', 'AvaI'); # not really true :)

is my @isos=$re->isoschizomers, 2;
is $isos[0],'BamHI';
ok $re->purge_isoschizomers;
is scalar($re->isoschizomers), 0;
ok $re->methylation_sites(2,5); # not really true :)
ok %meth = $re->methylation_sites;
is $meth{2}, 5;
ok $re->purge_methylation_sites;
is scalar($re->methylation_sites), 0;


ok $re->microbe('E. coli');
ok $microbe = $re->microbe;
is $microbe, "E. coli";
ok $re->source("Rob"); # not really true :)

ok $source = $re->source;
is $source, "Rob";

ok !$re->vendor;
ok $re->vendors('NEB'); # my favorite
ok @vendors = $re->vendors;
is $vendors[0], "NEB";
$re->purge_vendors;
is scalar($re->vendors), 0;

ok $re->references('Rob et al');
ok @refs = $re->references;
is $refs[0], "Rob et al";
$re->purge_references;
is scalar($re->references), 0;

ok $re->name('BamHI');
ok $name = $re->name;
is $name, "BamHI";

$re->verbose(2);

eval {$re->is_prototype};
ok($@);
like($@, qr/Can't unequivocally assign prototype based on input format alone/, 'bug 2179');
$re->verbose(2);

is $re->is_prototype(0), 0;
is $re->is_prototype, 0;
is $re->is_prototype(1), 1;
is $re->is_prototype, 1;

is $re->prototype_name, $re->name;
ok ! $re->is_prototype(0);
is $re->prototype_name('XxxI'), 'XxxI';
is $re->prototype_name, 'XxxI';


is $re->cutter, 6;
ok $re->seq->seq('RCATGY');
is $re->cutter, 5;

ok my $re2 = $re->clone;
isnt $re, $re2;
is $re->name, $re2->name;

ok $re = Bio::Restriction::Enzyme->new(-enzyme=>'AciI', 
										-site=>'C^CGC');
is $re->palindromic, 0;
is $re->is_palindromic, 0;

#
# Bio::Restriction::Enzyme::MultiSite
#

ok $re=Bio::Restriction::Enzyme::MultiSite->new(-enzyme=>'TaqII',
                                              -site=>'GACCGA',
                                              -cut=>17,
                                              -complementary_cut=>15
                                             );
ok $re2=Bio::Restriction::Enzyme::MultiSite->new(-enzyme=>'TaqII',
                                                -site=>'CACCCA',
                                                -cut=>17,
                                                -complementary_cut=>15
                                               );
isa_ok( $re, 'Bio::Restriction::EnzymeI');
isa_ok( $re2, 'Bio::Restriction::EnzymeI');
ok $re->others($re2);
ok $re2->others($re);

is $re->others, 1;
is $re2->others, 1;

ok my $re3 = $re->clone;
isnt $re, $re3;
is $re->name , $re3->name; # wouldn't this be a circular reference???
#print Dumper $re, $re3;exit;

#
# Bio::Restriction::Enzyme::MultiCut
#
#Hin4I has four cut sites [(8/13)GAYNNNNNVTC(13/8)],

ok $re = Bio::Restriction::Enzyme::MultiCut->new(-enzyme=>'Hin4I',
                                              -site=>'GAYNNNNNVTC',
                                              -cut=>-8,
                                              -complementary_cut=>-13
                                             );
ok $re2 = Bio::Restriction::Enzyme::MultiCut->new(-enzyme=>'Hin4I',
                                               -site=>'GAYNNNNNVTC',
                                               -cut=>13,
                                               -complementary_cut=>8
                                              );
isa_ok($re, 'Bio::Restriction::EnzymeI');
isa_ok($re2, 'Bio::Restriction::EnzymeI');
ok $re->others($re2);
ok $re2->others($re);

ok $re3 = $re->clone;
isnt $re, $re3;
is $re->name, $re3->name;
#print Dumper $re, $re3;exit;

#
# Bio::Restriction::EnzymeCollection
#

my ($collection, $enz, $new_set);

ok $collection = Bio::Restriction::EnzymeCollection->new(-empty=>1);
is $collection->each_enzyme, 0;
# default set
$collection = Bio::Restriction::EnzymeCollection->new;
is $collection->each_enzyme, 532;
is $collection->each_enzyme, 532;

ok $enz = $collection->get_enzyme('AclI');
isa_ok($enz, 'Bio::Restriction::Enzyme');
is my @enzymes=$collection->available_list, 532;

ok $new_set = $collection->blunt_enzymes;
isa_ok($enz, 'Bio::Restriction::Enzyme');
is $new_set->each_enzyme, 114;

#map {print $_->name, ": ", $_->cutter, "\n"; } $collection->each_enzyme;

ok $new_set = $collection->cutters(8);
is $new_set->each_enzyme, 17;

ok $new_set=$collection->cutters(-start => 8, -end => 8);
is $new_set->each_enzyme, 17;

ok $new_set=$collection->cutters(-start => 6, -end => 8);
is $new_set->each_enzyme, 293;

ok $new_set=$collection->cutters(-start => 6, -end => 8,  -exclusive => 1);
is $new_set->each_enzyme, 10;

ok $new_set = $collection->cutters([4,8]);
is $new_set->each_enzyme, 129;

# bug 2128; enhancement request to pass array ref of sizes

#
# Restriction::Analysis
#


ok my $seqio=Bio::SeqIO->new(-file=>test_input_file('dna1.fa'),
                         -format=>'fasta');
ok $seq=$seqio->next_seq;

ok my $ra = Bio::Restriction::Analysis->new(-seq=>$seq);
ok my $uniq = $ra->unique_cutters;

# test most objects
is $ra->unique_cutters->each_enzyme, 42, 'number of unique cutters';
is $ra->fragments('RsaI'), 2, 'number of RsaI fragments';
is $ra->max_cuts, 9, 'number of maximum cutters';
is $ra->zero_cutters->each_enzyme, 477, 'number of zero cutters';
is $ra->cutters->each_enzyme, 55, 'number of cutters';
is $ra->cutters(3)->each_enzyme, 8, 'number of 3x cutters';
is $ra->fragments('MseI'), 4, '4 MseI fragments';
is $ra->cuts_by_enzyme('MseI'), 3, '3 MseI cut sites';

#my $z = $ra->cutters(3);
#my $out=Bio::Restriction::IO->new;
#$out->write($z);


is $ra->fragments('PspEI'), 2, 'expected 2 PspEI fragments';
is $ra->cuts_by_enzyme('PspEI'), 1;
is $ra->cuts_by_enzyme('XxxI'), undef;


is my @ss = $ra->sizes('PspEI'), 2, 'expected 2 sizes for PspEI';
is $ss[0] + $ss[1], $seq->length;

# Issue 3157
$re = Bio::Restriction::Enzyme->new(-enzyme=>'PspEI', -site=>'G^GTNACC');
is @ss = $ra->sizes($re), 2, 'expected 2 sizes for PspEI';
is $ss[0] + $ss[1], $seq->length;

is $ra->fragments('MwoI'), 1, 'not circular expected 1 fragments for MwoI as it doesnt cut';

# circularise the sequence, regenerate the cuts and test again
# note that there is one less fragment now!
ok $seq->is_circular(1);

# we need to regenerate all the cuts
ok $ra->cut;

is $ra->fragments('RsaI'), 1, 'number of RsaI fragments';
is $ra->fragments('MseI'), 3, '3 circular MseI fragments';
is $ra->cuts_by_enzyme('MseI'), 3, '3 circular MseI cut sites';
is $ra->fragments('AciI'), 1, 'number for AciI a non-palindromic enzyme';

is $ra->fragments('MwoI'), 1, '1 fragment for MwoI as it cuts across the circ point';

ok my @rb=($collection->get_enzyme("AluI"), $collection->get_enzyme("MseI"), $collection->get_enzyme("MaeIII"));

# test multiple digests
ok my $rbc=Bio::Restriction::EnzymeCollection->new(-empty=>1);
ok $rbc->enzymes(@rb);
ok $ra->cut('multiple', $rbc);
is $ra->fragments('multiple_digest'),7, '7 fragments in the multiple digest';
is my @pos=$ra->positions('multiple_digest'),7, '7 positions in the multiple digest';
is my @ssm = $ra->sizes('multiple_digest'),7, '7 sizes in the multiple digest';
my $check_len;
map {$check_len+=$_}@ssm;
is $check_len, $seq->length;

# now test the non-palindromic part
# HindI is a non palindromic enzyme that cuts 9 times
is $ra->positions('HindI'), 9, ' expected 9 cuts for HindI';

# now we need to test the fragment maps
# lets do this for HindI
is my @fm=$ra->fragment_maps('HindI'), 9, 'expect 9 fragment maps for HindI';
foreach my $fm (@fm) {
 is exists $fm->{'seq'}, 1, "sequence for ".$fm->{'seq'};
 is exists $fm->{'start'}, 1, "start at ".$fm->{'start'};
 is exists $fm->{'end'}, 1, "end at ".$fm->{'end'};
}

# bug 2139

eval {$re = Bio::Restriction::Enzyme->new(
        -name    => 'Invalid',
        -site    => 'G^AATTE' );};

ok $@;
like($@, qr(Unrecognized characters in site), 'bug 2139');

# 0-pos bug (Elia Stupka)

$seq = Bio::Seq->new(
    -display_name   => 'foo',
    -alphabet       => 'dna',
    -seq            => 'GATCNNNNGATC'
);

$ra = Bio::Restriction::Analysis->new(-seq=>$seq);

is $ra->fragments('HindIII'), 1, 'number of HindIII fragments';
is $ra->fragments('BfuCI'), 2, 'number of EcoRI fragments';

# passing a bad enzyme name

is $ra->fragments('FooBarI'), 1, 'number of RsaI fragments';


