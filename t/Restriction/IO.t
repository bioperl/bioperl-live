# -*-Perl-*- Test Harness script for Bioperl
# $Id$


use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 18);
	
    use_ok('Bio::Restriction::IO');
}

#
# default enz set
#
ok my $in  = Bio::Restriction::IO->new();
ok my $renzs = $in->read;
is $renzs->each_enzyme, 532;

ok my $e = $renzs->get_enzyme('AccI');
is $e->name, 'AccI';

my $outfile = test_output_file();
ok my $out = Bio::Restriction::IO->new(-format => 'base', -file => ">$outfile");
TODO: {
    local $TODO = "writing to a file doesn't seem to work? prints to STDOUT!";
    #$out->write($renzs);
    ok -s $outfile;
    #map {print $_->name, "\t", $_->site, "\t", $_->overhang, "\n"} $renzs->each_enzyme;
}

#
# withrefm, 31
#

ok $in  = Bio::Restriction::IO->new
  (-format=> 'withrefm',
	-verbose => 0,
	-file => test_input_file('rebase.withrefm'));
ok $renzs = $in->read;
is $renzs->each_enzyme, 11;

#
# itype2, 8
#
#enzyme name [tab] prototype [tab] recognition sequence with cleavage site
#  [tab] methylation site and type [tab] commercial source [tab] references

ok $in  = Bio::Restriction::IO->new
    (-format=> 'itype2', -verbose => 0,
     -file => test_input_file('rebase.itype2'));

ok $renzs = $in->read;
is $renzs->each_enzyme, 16;

ok  $out  = Bio::Restriction::IO->new(-format=>'base');

SKIP: {
   test_skip(-tests => 3, -requires_networking => 1);
   #test_skip(-tests => 2, -requires_module => 'LWP::UserAgent');
    
    ok  $in  = Bio::Restriction::IO->new(-format=>'prototype',
                                          -current => 1);
    
    ok my $coll = $in->read;
    cmp_ok $coll->each_enzyme, '>=', 307;
}

