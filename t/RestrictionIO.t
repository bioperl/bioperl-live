# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

# test for Bio::Restriction::Analysis.pm
# written by Rob Edwards

use strict;

my $NUMTESTS;

BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't','..';
    }
    use Test;
    $NUMTESTS = 14;
    plan tests => $NUMTESTS;
}

require Bio::Restriction::IO;
use Bio::Root::IO;
require Bio::Restriction::EnzymeCollection;

my $tmpdir = "t/tmp";
mkdir($tmpdir,0777);

ok(1);

#
# default enz set
#
ok my $in  = Bio::Restriction::IO->new();
ok my $renzs = $in->read;
ok $renzs->each_enzyme, 532;

ok my $e = $renzs->get_enzyme('AccI');
ok $e->name, 'AccI';

ok my $out = Bio::Restriction::IO->new(-format => 'base',
													-file   => ">$tmpdir/r");
#$out->write($renzs);
#map {print $_->name, "\t", $_->site, "\t", $_->overhang, "\n"} $renzs->each_enzyme;

#
# withrefm, 31
#

ok $in  = Bio::Restriction::IO->new
    (-format=> 31, -verbose => 0,
     -file => Bio::Root::IO->catfile("t","data","rebase.withrefm"));
ok $renzs = $in->read;
ok $renzs->each_enzyme, 11;

#
# itype2, 8
#
#enzyme name [tab] prototype [tab] recognition sequence with cleavage site
#  [tab] methylation site and type [tab] commercial source [tab] references

ok $in  = Bio::Restriction::IO->new
    (-format=> 8, -verbose => 0,
     -file => Bio::Root::IO->catfile("t","data","rebase.itype2"));

ok $renzs = $in->read;
ok $renzs->each_enzyme, 16;

ok  $out  = Bio::Restriction::IO->new(-format=>'base');

END { cleanup(); }

sub cleanup {
   eval {
      Bio::Root::IO->rmtree($tmpdir) if (-d $tmpdir);
   };
}
