# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 272);
	
    use_ok('Bio::Factory::FTLocationFactory');
}

my $simple_impl = "Bio::Location::Simple";
my $fuzzy_impl = "Bio::Location::Fuzzy";
my $split_impl = "Bio::Location::Split";

# Holds strings and results. The latter is an array of expected class name,
# min/max start position and position type, min/max end position and position
# type, location type, the number of locations, and the strand.
#
my %testcases = (
   # note: the following are directly taken from 
   # http://www.ncbi.nlm.nih.gov/collab/FT/#location
   "467" => [$simple_impl,
	    467, 467, "EXACT", 467, 467, "EXACT", "EXACT", 1, 1],
	"340..565" => [$simple_impl,
		 340, 340, "EXACT", 565, 565, "EXACT", "EXACT", 1, 1],
	"<345..500" => [$fuzzy_impl,
		 undef, 345, "BEFORE", 500, 500, "EXACT", "EXACT", 1, 1],
	"<1..888" => [$fuzzy_impl,
		 undef, 1, "BEFORE", 888, 888, "EXACT", "EXACT", 1, 1],
	"(102.110)" => [$fuzzy_impl,
		 102, 102, "EXACT", 110, 110, "EXACT", "WITHIN", 1, 1],
	"(23.45)..600" => [$fuzzy_impl,
		 23, 45, "WITHIN", 600, 600, "EXACT", "EXACT", 1, 1],
	"(122.133)..(204.221)" => [$fuzzy_impl,
		 122, 133, "WITHIN", 204, 221, "WITHIN", "EXACT", 1, 1],
	"123^124" => [$simple_impl,
		 123, 123, "EXACT", 124, 124, "EXACT", "IN-BETWEEN", 1, 1],
	"145^177" => [$fuzzy_impl,
		 145, 145, "EXACT", 177, 177, "EXACT", "IN-BETWEEN", 1, 1],
	"join(12..78,134..202)" => [$split_impl,
		 12, 12, "EXACT", 202, 202, "EXACT", "EXACT", 2, 1],
	"complement(join(4918..5163,2691..4571))" => [$split_impl,
		 2691, 2691, "EXACT", 5163, 5163, "EXACT", "EXACT", 2, -1],
	"complement(34..(122.126))" => [$fuzzy_impl,
		 34, 34, "EXACT", 122, 126, "WITHIN", "EXACT", 1, -1],
	"J00194:100..202" => [$simple_impl,
		 100, 100, "EXACT", 202, 202, "EXACT", "EXACT", 1, 1],
	# this variant is not really allowed by the FT definition
	# document but we want to be able to cope with it
	"J00194:(100..202)" => [$simple_impl,
		 100, 100, "EXACT", 202, 202, "EXACT", "EXACT", 1, 1],
	"((122.133)..(204.221))" => [$fuzzy_impl,
		 122, 133, "WITHIN", 204, 221, "WITHIN", "EXACT", 1, 1],
	"join(AY016290.1:108..185,AY016291.1:1546..1599)"=> [$split_impl,
		 108, 108, "EXACT", 185, 185, "EXACT", "EXACT", 2, undef],

	# UNCERTAIN locations and positions (Swissprot)
   "?2465..2774" => [$fuzzy_impl,
       2465, 2465, "UNCERTAIN", 2774, 2774, "EXACT", "EXACT", 1, 1],
   "22..?64" => [$fuzzy_impl,
       22, 22, "EXACT", 64, 64, "UNCERTAIN", "EXACT", 1, 1],
   "?22..?64" => [$fuzzy_impl,
       22, 22, "UNCERTAIN", 64, 64, "UNCERTAIN", "EXACT", 1, 1],
   "?..>393" => [$fuzzy_impl,
       undef, undef, "UNCERTAIN", 393, undef, "AFTER", "UNCERTAIN", 1, 1],
   "<1..?" => [$fuzzy_impl,
       undef, 1, "BEFORE", undef, undef, "UNCERTAIN", "UNCERTAIN", 1, 1],
   "?..536" => [$fuzzy_impl,
       undef, undef, "UNCERTAIN", 536, 536, "EXACT", "UNCERTAIN", 1, 1],
   "1..?" => [$fuzzy_impl,
       1, 1, "EXACT", undef, undef, "UNCERTAIN", "UNCERTAIN", 1, 1],
   "?..?" => [$fuzzy_impl,
       undef, undef, "UNCERTAIN", undef, undef, "UNCERTAIN", "UNCERTAIN", 1, 1],
   # Not working yet:
   #"12..?1" => [$fuzzy_impl,
   #    1, 1, "UNCERTAIN", 12, 12, "EXACT", "EXACT", 1, 1]
		 );

my $locfac = Bio::Factory::FTLocationFactory->new();
isa_ok($locfac,'Bio::Factory::LocationFactoryI');

# sorting is to keep the order constant from one run to the next
foreach my $locstr (keys %testcases) { 
	my $loc = $locfac->from_string($locstr);
	if($locstr eq "join(AY016290.1:108..185,AY016291.1:1546..1599)") {
		$loc->seq_id("AY016295.1");
	}
	my @res = @{$testcases{$locstr}};
	is(ref($loc), $res[0], $res[0]);
	is($loc->min_start(), $res[1]);
	is($loc->max_start(), $res[2]);
	is($loc->start_pos_type(), $res[3]);
	is($loc->min_end(), $res[4]);
	is($loc->max_end(), $res[5]);
	is($loc->end_pos_type(), $res[6]);
	is($loc->location_type(), $res[7]);
	my @locs = $loc->each_Location();
	is(@locs, $res[8]);
	my $ftstr = $loc->to_FTstring();
	# this is a somewhat ugly hack, but we want clean output from to_FTstring()
	# Umm, then these should really fail, correct?
	# Should we be engineering workarounds for tests?
	$locstr = "J00194:100..202" if $locstr eq "J00194:(100..202)";
	$locstr = "(122.133)..(204.221)" if $locstr eq "((122.133)..(204.221))";
	# now test
	is($ftstr, $locstr, "Location String: $locstr");
	# test strand production
	is($loc->strand(), $res[9]);
}

SKIP: {
    skip('nested matches in regex only supported in v5.6.1 and higher', 5) unless $^V gt v5.6.0;
    
	# bug #1674, #1765, 2101
	# EMBL-like 
	# join(20464..20694,21548..22763,join(complement(314652..314672),complement(232596..232990),complement(231520..231669)))
	# GenBank-like
	# join(20464..20694,21548..22763,complement(join(231520..231669,232596..232990,314652..314672)))
	# Note that
	# join(1000..2000,join(3000..4000,join(5000..6000,7000..8000)),9000..10000)
	# is the same as
	# join(1000..2000,3000..4000,5000..6000,7000..8000,9000..10000)
	# But I don't want to bother with it at this point
	my @expected = (# intentionally testing same expected string twice
					# as I am providing two different encodings
					# that should mean the same thing
	'join(11025..11049,complement(join(315036..315294,251354..251412,241499..241580,239890..240081)))',
	'join(11025..11049,complement(join(315036..315294,251354..251412,241499..241580,239890..240081)))',
	# ditto
	'join(20464..20694,21548..22763,complement(join(314652..314672,232596..232990,231520..231669)))',
	'join(20464..20694,21548..22763,complement(join(314652..314672,232596..232990,231520..231669)))',
	# this is just seen once
	'join(1000..2000,join(3000..4000,join(5000..6000,7000..8000)),9000..10000)',
	'order(S67862.1:72..75,join(S67863.1:1..788,1..19))'
   );

	for my $locstr (
		'join(11025..11049,join(complement(239890..240081),complement(241499..241580),complement(251354..251412),complement(315036..315294)))',
		'join(11025..11049,complement(join(315036..315294,251354..251412,241499..241580,239890..240081)))',
		'join(20464..20694,21548..22763,complement(join(314652..314672,232596..232990,231520..231669)))',
		'join(20464..20694,21548..22763,join(complement(231520..231669),complement(232596..232990),complement(314652..314672)))',
		'join(1000..2000,join(3000..4000,join(5000..6000,7000..8000)),9000..10000)',
		'order(S67862.1:72..75,join(S67863.1:1..788,1..19))'
	   ) {
		my $loc = $locfac->from_string($locstr);
		my $ftstr = $loc->to_FTstring();
		is($ftstr, shift @expected, $locstr);
	}
}
