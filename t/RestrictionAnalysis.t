# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

# test for Bio::Restriction::Analysis.pm
# written by Rob Edwards & Heikki Lehvaslaiho

use strict;
use constant NUMTESTS => 84;

BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't','..';
    }
    use Test;

    plan tests => NUMTESTS;
}

use Bio::Restriction::Enzyme;
use Bio::Restriction::Enzyme::MultiCut;
use Bio::Restriction::Enzyme::MultiSite;

use Bio::Root::IO;
use Bio::Restriction::EnzymeCollection;


#use Bio::Restriction::Analysis;
use Bio::SeqIO;
use Data::Dumper;
ok(1);


#
# Bio::Restriction::Enzyme
#


my ($re, $seq, $iso, %meth, $microbe, $source, @vendors, @refs, $name);
ok $re=new Bio::Restriction::Enzyme(-enzyme=>'EcoRI', -site=>'G^AATTC');
ok $re->isa('Bio::Restriction::EnzymeI');
#exit;
ok $re->cut, 1;
ok ! $re->cut(0);
ok $re->complementary_cut, 0;
$re->cut(1);

#print Dumper $re;
#exit;
ok $re->complementary_cut,5;
ok $re->site,'G^AATTC';
ok $seq=$re->seq;
ok $seq->seq, 'GAATTC';
ok $re->string,'GAATTC';
ok $re->revcom, 'GAATTC';
ok $re->recognition_length, 6;
ok $re->non_ambiguous_length, 6;
ok $re->cutter, 6;
ok $re->palindromic, 1;
#ok $re->is_blunt, 0;
ok $re->overhang, "5'";
ok $re->overhang_seq, 'AATT';
ok $re->is_ambiguous, 0;

ok $re->compatible_ends($re);
#print Dumper $re;


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

#print Dumper $re, $re2;

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


#print Dumper $re, $re2;
#exit;


#
# Bio::Restriction::EnzymeCollection
#

my ($collection, $enz, $new_set);

ok $collection = new Bio::Restriction::EnzymeCollection(-empty=>1);
ok $collection->each_enzyme, 0;
# default set
ok $collection = new Bio::Restriction::EnzymeCollection;
ok $collection->each_enzyme, 530;
ok $enz=$collection->get_enzyme('AclI');
ok $enz->isa('Bio::Restriction::Enzyme');
ok my @enzymes=$collection->available_list, 530;

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

my $enzymes;
my $seqio=new Bio::SeqIO(-file=>Bio::Root::IO->catfile("t","data","dna1.fa"),
                         -format=>'fasta');
$seq=$seqio->next_seq;
my ($anal, $sizes, $cuts, $enzobj);

#ok $anal=Bio::Restriction::Analysis->new
#    (-seq=>$seq,
#     -collection => Bio::Restriction::EnzymeCollection->new);

#exit;

#ok $enzymes=$anal->unique_cutters;
#ok defined $enzymes, 1, 'could not find any unique cutters';
#$enz=$$enzymes[0];
#ok defined $anal->fragments($enz), 1, "Couldn't find a fragment for $enz";
#ok $enzymes=$anal->zero_cutters;
#ok defined $enzymes, 1, 'could not find any zero cutters';
#ok $sizes=$anal->sizes($enz);
#ok defined $sizes, 1, 'could not find any sizes';
#ok $cuts=$anal->cuts_by_enzyme($enz);
#ok defined $cuts, 1, 'could not find any cuts';
#ok $cuts=$anal->cuts_by_frequency(1);
#ok $$cuts[0] eq $enz; # these should be the same!
#ok $enzymes=$anal->all_cutters;
#ok defined $enzymes, 1, 'could not find all cutters';
#ok $enzobj=$anal->enzyme($enz);
#ok $enzobj->name eq $enz;
#ok ref($enzobj) eq "Bio::Restriction::Enzyme";
