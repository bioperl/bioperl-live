# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

# test for Bio::Restriction::Analysis.pm
# written by Rob Edwards

use strict;
use constant NUMTESTS => 14;

BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't','..';
    }
    use Test;

    plan tests => NUMTESTS;
}

use Bio::Restriction::IO;

use Bio::Root::IO;
use Bio::Restriction::EnzymeCollection;

use Data::Dumper;
ok(1);


#
# default enz set
#
ok my $in  = Bio::Restriction::IO->new();
#print Dumper $in;
ok my $renzs = $in->read;
ok $renzs->each_enzyme, 530;
#print Dumper $renzs;

ok my $e = $renzs->get_enzyme('AccI');
ok $e->name, 'AccI';


ok my $out  = Bio::Restriction::IO->new(-format=>'base',  -file=> ">/tmp/r");
#$out->write($renzs);

#map {print $_->name, "\t", $_->site, "\t", $_->overhang, "\n"} $renzs->each_enzyme;

#my @l = $renzs->available_list;
#ok @l, 530;
#ok $renzs->blunt_enzymes;
#ok $renzs->blunt_enzymes->each_enzyme, 1;



#
# withrefm, 31
#

ok $in  = Bio::Restriction::IO->new
    (-format=> 31, -verbose => 0,
     -file => Bio::Root::IO->catfile("t","data","rebase.withrefm"));
#     -file => Bio::Root::IO->catfile("t","data","rebase_link_withrefm.txt"));

#print Dumper $in;
ok $renzs = $in->read;
ok $renzs->each_enzyme, 17;
#print Dumper $renzs;



#
# itype2, 8
#
#enzyme name [tab] prototype [tab] recognition sequence with cleavage site
#  [tab] methylation site and type [tab] commercial source [tab] references

ok $in  = Bio::Restriction::IO->new
    (-format=> 8, -verbose => 0,
     -file => Bio::Root::IO->catfile("t","data","rebase.itype2"));
#     -file => Bio::Root::IO->catfile("t","data","rebase_link_withrefm.txt"));

#print Dumper $in;
ok $renzs = $in->read;
ok $renzs->each_enzyme, 19;
#print Dumper $renzs;
#print $renzs->get_enzyme('AacI')->cuts_after, "\n";

ok  $out  = Bio::Restriction::IO->new(-format=>'base');
#$out->write($renzs);
