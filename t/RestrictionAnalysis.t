# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

# test for Bio::Restriction::Analysis.pm
# written by Rob Edwards & Heikki Lehvaslaiho

use strict;
my $NUMTESTS;
my $error;


BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't','..';
    }
    use Test;
    $NUMTESTS = 153;
    $error  = 0;

    plan tests => $NUMTESTS;

}

if( $error ==  1 ) {
    exit(0);
}

require Bio::Restriction::Enzyme;
require Bio::Restriction::Enzyme::MultiCut;
require Bio::Restriction::Enzyme::MultiSite;
require Bio::Restriction::EnzymeCollection;
require Bio::Restriction::Analysis;
use Bio::Root::IO;
use Bio::SeqIO;
use Data::Dumper;
ok(1);


#
# Bio::Restriction::Enzyme
#

my ($re, $seq, $iso, %meth, $microbe, $source, @vendors, @refs, $name);
ok $re=new Bio::Restriction::Enzyme(-enzyme=>'EcoRI', -site=>'G^AATTC');
ok $re->isa('Bio::Restriction::EnzymeI');
ok $re->cut, 1;
ok ! $re->cut(0);
ok $re->complementary_cut, 6;
$re->cut(1);

ok $re->complementary_cut,5;
ok $re->site,'G^AATTC';
ok $seq=$re->seq;
ok $seq->seq, 'GAATTC';
ok $re->string,'GAATTC';
ok $re->revcom, 'GAATTC';
ok $re->recognition_length, 6;
ok $re->cutter, 6;
ok $re->palindromic, 1;
ok $re->overhang, "5'";
ok $re->overhang_seq, 'AATT';
ok $re->is_ambiguous, 0;

ok $re->compatible_ends($re);

ok $re->isoschizomers('BamHI', 'AvaI'); # not really true :)
ok my @isos=$re->isoschizomers, 2;
ok $isos[0] eq 'BamHI';
ok $re->purge_isoschizomers;
ok $re->methylation_sites(2,5); # not really true :)
ok %meth=$re->methylation_sites;
ok $meth{2} == 5;
ok $re->microbe('E. coli');
ok $microbe=$re->microbe;
ok $microbe eq "E. coli";
ok $re->source("Rob"); # not really true :)
ok $source=$re->source;
ok $source eq "Rob";

ok !$re->vendor;
ok $re->vendors('NEB'); # my favorite
ok @vendors=$re->vendors;
ok $vendors[0] eq "NEB";
ok $re->purge_vendors;

ok $re->references('Rob et al');
ok @refs=$re->references;
ok $refs[0] eq "Rob et al";
ok $re->purge_references;

ok $re->name('BamHI');
ok $name=$re->name;
ok $name eq "BamHI";

ok !$re->is_prototype;
ok $re->is_prototype(1);

ok $re->prototype_name, $re->name;
ok ! $re->is_prototype(0);
ok $re->prototype_name('XxxI'), 'XxxI';
ok $re->prototype_name, 'XxxI';


ok $re->cutter, 6;
$re->seq->seq('RCATGY');
ok $re->cutter, 5;

ok my $re2 = $re->clone;
ok $re ne $re2;
ok $re->name eq $re2->name;

ok $re = new Bio::Restriction::Enzyme(-enzyme=>'AciI', 
												  -site=>'C^CGC');
ok $re->palindromic, 0;
ok $re->is_palindromic, 0;

#
# Bio::Restriction::Enzyme::MultiSite
#

ok $re=new Bio::Restriction::Enzyme::MultiSite(-enzyme=>'TaqII',
                                              -site=>'GACCGA',
                                              -cut=>17,
                                              -complementary_cut=>15
                                             );
ok $re2=new Bio::Restriction::Enzyme::MultiSite(-enzyme=>'TaqII',
                                                -site=>'CACCCA',
                                                -cut=>17,
                                                -complementary_cut=>15
                                               );
ok $re->isa('Bio::Restriction::EnzymeI');
ok $re->others($re2);
ok $re2->others($re);

ok $re->others, 1;

ok my $re3 = $re->clone;
ok $re ne $re3;
ok $re->name eq $re3->name;
#print Dumper $re, $re3;exit;

#
# Bio::Restriction::Enzyme::MultiCut
#
#Hin4I has four cut sites [(8/13)GAYNNNNNVTC(13/8)],

ok $re=new Bio::Restriction::Enzyme::MultiCut(-enzyme=>'Hin4I',
                                              -site=>'GAYNNNNNVTC',
                                              -cut=>-8,
                                              -complementary_cut=>-13
                                             );
ok $re2=new Bio::Restriction::Enzyme::MultiCut(-enzyme=>'Hin4I',
                                               -site=>'GAYNNNNNVTC',
                                               -cut=>13,
                                               -complementary_cut=>8
                                              );
ok $re->isa('Bio::Restriction::EnzymeI');
ok $re->others($re2);
ok $re2->others($re);


ok $re3 = $re->clone;
ok $re ne $re3;
ok $re->name eq $re3->name;
#print Dumper $re, $re3;exit;

#
# Bio::Restriction::EnzymeCollection
#

my ($collection, $enz, $new_set);

ok $collection = new Bio::Restriction::EnzymeCollection(-empty=>1);
ok $collection->each_enzyme, 0;
# default set
ok $collection = new Bio::Restriction::EnzymeCollection;
ok $collection->each_enzyme, 532;
ok $collection = new Bio::Restriction::EnzymeCollection;
ok $collection->each_enzyme, 532;

ok $enz=$collection->get_enzyme('AclI');
ok $enz->isa('Bio::Restriction::Enzyme');
ok my @enzymes=$collection->available_list, 532;

ok $new_set=$collection->blunt_enzymes;
ok $new_set->each_enzyme, 114;

#map {print $_->name, ": ", $_->cutter, "\n"; } $collection->each_enzyme;

ok $new_set=$collection->cutters(8);
ok $new_set->each_enzyme, 17;

ok $new_set=$collection->cutters(-start => 8, -end => 8);
ok $new_set->each_enzyme, 17;

ok $new_set=$collection->cutters(-start => 6, -end => 8);
ok $new_set->each_enzyme, 293;

ok $new_set=$collection->cutters(-start => 6, -end => 8,  -exclusive => 1);
ok $new_set->each_enzyme, 10;


#
# Restriction::Analysis
#


my $seqio=new Bio::SeqIO(-file=>Bio::Root::IO->catfile("t","data","dna1.fa"),
                         -format=>'fasta');
$seq=$seqio->next_seq;

ok my $ra=Bio::Restriction::Analysis->new(-seq=>$seq);
ok my $uniq = $ra->unique_cutters;

# test most objects
ok $ra->unique_cutters->each_enzyme, 42, 'wrong number of unique cutters';
ok $ra->fragments('RsaI'), 2, 'wrong number of RsaI fragments';
ok $ra->max_cuts, 9, 'wrong number of maximum cutters';
ok $ra->zero_cutters->each_enzyme, 477, 'wrong number of zero cutters';
ok $ra->cutters->each_enzyme, 55, 'wrong number of cutters';
ok $ra->cutters(3)->each_enzyme, 8, 'wrong number of 3x cutters';
ok $ra->fragments('MseI'), 4, 'expected 4 MseI fragments';
ok $ra->cuts_by_enzyme('MseI'), 3, 'expected 3 MseI cut sites';

#my $z = $ra->cutters(3);
#my $out=Bio::Restriction::IO->new;
#$out->write($z);


ok $ra->fragments('PspEI'), 2, 'expected 2 PspEI fragments';
ok $ra->cuts_by_enzyme('PspEI'), 1;
ok $ra->cuts_by_enzyme('XxxI'), undef;


ok my @ss = $ra->sizes('PspEI'), 2, 'expected 2 sizes for PspEI';
ok $ss[0] + $ss[1], $seq->length;

ok $ra->fragments('MwoI'), 1, 'not circular expected 1 fragments for MwoI as it doesnt cut';

# circularise the sequence, regenerate the cuts and test again
# note that there is one less fragment now!
$seq->is_circular(1);

# we need to regenerate all the cuts
ok $ra->cut;

ok $ra->fragments('RsaI'), 1, 'wrong number of RsaI fragments';
ok $ra->fragments('MseI'), 3, 'expected 3 circular MseI fragments';
ok $ra->cuts_by_enzyme('MseI'), 3, 'expected 3 circular MseI cut sites';
ok $ra->fragments('AciI'), 1, 'wrong number for AciI a non-palindromic enzyme';

ok $ra->fragments('MwoI'), 1, 'expected 1 fragments for MwoI as it cuts across the circ point';

ok my @rb=($collection->get_enzyme("AluI"), $collection->get_enzyme("MseI"), $collection->get_enzyme("MaeIII"));

# test multiple digests
ok my $rbc=Bio::Restriction::EnzymeCollection->new(-empty=>1);
ok $rbc->enzymes(@rb);
ok $ra->cut('multiple', $rbc);
ok $ra->fragments('multiple_digest'),7, 'expected 7 fragments in the multiple digest';
ok my @pos=$ra->positions('multiple_digest'),7, 'expected 7 positions in the multiple digest';
ok my @ssm = $ra->sizes('multiple_digest'),7, 'expected 7 sizes in the multiple digest';
my $check_len;
map {$check_len+=$_}@ssm;
ok $check_len == $seq->length;

# now test the non-palindromic part
# HindI is a non palindromic enzyme that cuts 9 times
ok $ra->positions('HindI'), 9, ' expected 9 cuts for HindI';

# now we need to test the fragment maps
# lets do this for HindI
ok my @fm=$ra->fragment_maps('HindI'), 9, 'expect 9 fragment maps for HindI';
foreach my $fm (@fm) {
 ok exists $fm->{'seq'}, '1', "no sequence for $fm";
 ok exists $fm->{'start'}, '1', "wrong start at ".$fm->{'start'};
 ok exists $fm->{'end'}, 1, "wrong end at ".$fm->{'end'};
}


