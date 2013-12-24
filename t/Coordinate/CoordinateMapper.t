use strict;
use warnings;

BEGIN {
    use Bio::Root::Test;

    test_begin(-tests => 175);

    use_ok('Bio::Location::Simple');
    use_ok('Bio::Coordinate::Pair');
    use_ok('Bio::Coordinate::Result::Match');
    use_ok('Bio::Coordinate::Result::Gap');
    use_ok('Bio::Coordinate::Chain');
    use_ok('Bio::Coordinate::Collection');
}

my ($c, $value);

ok $c = Bio::Coordinate::Result::Match-> new;
ok $c = Bio::Coordinate::Result::Gap-> new;

# propepide
my $match1 = Bio::Location::Simple->new
    (-seq_id => 'propeptide', -start => 21, -end => 40, -strand=>1 );
# peptide
my $match2 = Bio::Location::Simple->new
    (-seq_id => 'peptide', -start => 1, -end => 20, -strand=>1 );

ok my $pair = Bio::Coordinate::Pair->new(-in => $match1,
                                         -out => $match2,
                                         -negative => 0, # false, default
                                        );

ok $pair->test;
is $pair->strand(), 1; #  = in->strand * out->strand
is $pair->in->seq_id(), 'propeptide';


my ($count, $pos, $pos2, $res, $match, $res2);


#
# match within
#
$pos = Bio::Location::Simple->new
    (-start => 25, -end => 25, -strand=> -1 );

# results are in Bio::Coordinate::Result
# they can be Matches and Gaps; are  Bio::LocationIs
ok $res = $pair->map($pos);
isa_ok $res, 'Bio::Coordinate::Result';
isa_ok $res, 'Bio::Location::SplitLocationI';
is $res->each_match, 1;
is $res->each_gap, 0;
is $res->each_Location, 1;

isa_ok $res->match, 'Bio::LocationI';
isa_ok $res->match, 'Bio::Coordinate::Result::Match';
is $res->match->start, 5;
is $res->match->end, 5;
is $res->match->strand, -1;
is $res->match->seq_id, 'peptide';
is $res->start, 5;
is $res->end, 5;
is $res->strand, -1;
#is $res->seq_id, 'peptide';

# lets do the reverse
$match = $res->match;
ok $pair->swap;
$res2 = $pair->map($match);
is $res2->match->start, $pos->start;
is $res2->match->end, $pos->end;
is $res2->match->strand, $pos->strand;
is $res2->match->seq_id, $pair->out->seq_id;
ok $pair->swap;

#
# match outside = Gap
#
$pos = Bio::Location::Simple->new (-start => 5, -end => 5 );

ok $res = $pair->map($pos);
#$res->verbose(2);
is $res->each_Location, 1;
is $res->each_gap, 1;

isa_ok $res->gap, 'Bio::Coordinate::Result::Gap';
isa_ok $res->gap, 'Bio::LocationI';
is $res->gap->strand, 1;
is $res->gap->start, 5;
is $res->gap->length, $pos->length;
is $res->gap->seq_id, 'propeptide';


#
# partial match = gap & match
#
$pos2 = Bio::Location::Simple->new
    (-start => 20, -end => 22, -strand=> -1 );

ok $res = $pair->map($pos2);

is $res->each_match, 1;
is $res->each_gap, 1;
is $res->each_Location, 2;
is $res->match->length + $res->gap->length, $pos2->length;

is $res->match->start, 1;
is $res->match->end, 2;
is $res->match->seq_id, 'peptide';
is $res->match->strand, -1;
is $res->gap->start, 20;
is $res->gap->end, 20;
is $res->gap->seq_id, 'propeptide';
is $res->gap->strand, -1;

#
# partial match =  match & gap
#
$pos2 = Bio::Location::Simple->new (-start => 40, -end => 41, -strand=> 1 );
ok $res = $pair->map($pos2);
is $res->match->length + $res->gap->length, $pos2->length;

#
#enveloping
#
$pos2 = Bio::Location::Simple->new (-start => 19, -end => 41, -strand=> 1 );
ok $res = $pair->map($pos2);
$count = 0; map {$count += $_->length} $res->each_Location;
is $count, $pos2->length;




#
# Testing insertions
#
#out
$pos = Bio::Location::Simple->new (-start => 5, -end => 6, -location_type=>'^');
$res = $pair->map($pos);
is $res->each_gap, 1;
is $res->each_Location, 1;

#in
$pos = Bio::Location::Simple->new (-start => 21, -end => 22, -location_type=>'^');
$res = $pair->map($pos);
is $res->each_match, 1;
is $res->each_Location, 1;

#just before
$pos = Bio::Location::Simple->new (-start => 20, -end => 21, -location_type=>'^');
$res = $pair->map($pos);
is $res->each_gap, 1;
is $res->each_Location, 1;

#just after
$pos = Bio::Location::Simple->new (-start => 40, -end => 41, -location_type=>'^');
$res = $pair->map($pos);
is $res->each_gap, 1;
is $res->each_Location, 1;

#
# strandness
#
#   11   6 4 2
#  -|--------|-
#  -|--------|-
#   2    7 9 11
#

# from
$match1 = Bio::Location::Simple->new
    (-seq_id => 'from', -start => 2, -end => 11, -strand=>1 );
# to
$match2 = Bio::Location::Simple->new
    (-seq_id => 'to', -start => 2, -end => 11, -strand=>-1 );
$pair = Bio::Coordinate::Pair->new(-in => $match1,
                                   -out => $match2
                                  );
#
# match within
#

ok $pair->test;
is $pair->strand(), -1;
$pos = Bio::Location::Simple->new
    (-seq_id => 'from', -start => 7, -end => 9, -strand=>1 );
$res = $pair->map($pos);
is $res->match->start, 4;
is $res->match->end, 6;
is $res->match->strand, -1;

$pos = Bio::Location::Simple->new
    (-seq_id => 'from', -start => 3, -end => 10, -strand=>-1 );
$res = $pair->map($pos);
is $res->match->start, 3;
is $res->match->end, 10;
is $res->match->strand, 1;

#
# match outside = Gap
#
$pos = Bio::Location::Simple->new
    (-seq_id => 'from', -start => 1, -end => 1, -strand=>1 );
$res = $pair->map($pos);
is $res->gap->start, 1;
is $res->gap->end, 1;
is $res->gap->strand, 1;
$pos = Bio::Location::Simple->new
    (-seq_id => 'from', -start => 12, -end => 12, -strand=>-1 );
$res = $pair->map($pos);
is $res->gap->start, 12;
is $res->gap->end, 12;
is $res->gap->strand, -1;


#
# partial match1 = gap & match
#
$pos = Bio::Location::Simple->new
    (-seq_id => 'from', -start => 1, -end => 7, -strand=>-1 );
$res = $pair->map($pos);
is $res->gap->start, 1;
is $res->gap->end, 1;
is $res->gap->strand, -1;
is $res->match->start, 6;
is $res->match->end, 11;
is $res->match->strand, 1;

#
# partial match2 =  match & gap
#

$pos = Bio::Location::Simple->new
    (-seq_id => 'from', -start => 9, -end => 12, -strand=>-1 );
$res = $pair->map($pos);
is $res->match->start, 2;
is $res->match->end, 4;
is $res->match->strand, 1;
is $res->gap->start, 12;
is $res->gap->end, 12;
is $res->gap->strand, -1;

#
#enveloping
#

$pos = Bio::Location::Simple->new
    (-seq_id => 'from', -start => 1, -end => 12, -strand=>-1 );
$res = $pair->map($pos);
is $res->match->start, 2;
is $res->match->end, 11;
is $res->match->strand, 1;

my ($gap1, $gap2) = $res->each_gap;
is $gap1->start, 1;
is $gap1->end, 1;
is $gap1->strand, -1;
is $gap2->start, 12;
is $gap2->end, 12;
is $gap2->strand, -1;

#
# Chain
#
# chain (two) mappers together
#

# propepide
$match1 = Bio::Location::Simple->new
    (-seq_id => 'propeptide', -start => 5, -end => 40, -strand=>1 );
# peptide
$match2 = Bio::Location::Simple->new
    (-seq_id => 'peptide', -start => 1, -end => 36, -strand=>1 );

ok $pair = Bio::Coordinate::Pair->new(-in => $match1,
                                      -out => $match2
                                      );


ok my $chain = Bio::Coordinate::Chain->new;
ok $chain->add_mapper($pair);
$chain->add_mapper($pair);


$pos = Bio::Location::Simple->new
    (-seq_id => 'from', -start => 6, -end => 21, -strand=> 1 );

#  6 ->  2 ->  1
# 21 -> 17 -> 13
$match = $chain->map($pos);
isa_ok $match, 'Bio::Coordinate::Result::Match';
is $match->start, 1;
is $match->end, 13;
is $match->strand, 1;



#
# Collection
#
#         1   5     6   10
#         |---|     |---|
#-----|-----------------------
#     1   5   9     15  19
#         pair1     pair2

# gene
$match1 = Bio::Location::Simple->new
    (-seq_id => 'gene', -start => 5, -end => 9, -strand=>1 );
# exon2
$match2 = Bio::Location::Simple->new
    (-seq_id => 'exon1', -start => 1, -end => 5, -strand=>1 );

ok my $pair1 = Bio::Coordinate::Pair->new(-in => $match1,
                                          -out => $match2,
                                        );
# gene
my $match3 = Bio::Location::Simple->new
    (-seq_id => 'gene', -start => 15, -end => 19, -strand=>1 );
# exon
my $match4 = Bio::Location::Simple->new
    (-seq_id => 'exon2', -start => 6, -end => 10, -strand=>1 );

ok my $pair2 = Bio::Coordinate::Pair->new(-in => $match3,
                                          -out => $match4,
                                         );

ok my $transcribe = Bio::Coordinate::Collection->new;
ok $transcribe->add_mapper($pair1);
ok $transcribe->add_mapper($pair2);


# simple match
$pos = Bio::Location::Simple->new (-start => 5, -end => 9 );
ok $res = $transcribe->map($pos);
is $res->match->start, 1;
is $res->match->end, 5;
is $res->match->seq_id, 'exon1';

# flank pre
$pos = Bio::Location::Simple->new (-start => 2, -end => 9 );
ok $res = $transcribe->map($pos);
is $res->each_gap, 1;
is $res->each_match, 1;
is $res->match->start, 1;
is $res->match->end, 5;

# flank post
$pos = Bio::Location::Simple->new (-start => 5, -end => 12 );
ok $res = $transcribe->map($pos);
is $res->each_gap, 1;
is $res->each_match, 1;
is $res->match->start, 1;
is $res->match->end, 5;

# match more than two
$pos = Bio::Location::Simple->new (-start => 5, -end => 19 );
ok $res = $transcribe->map($pos);
is $res->each_gap, 2;
is $res->each_match, 2;



# testing sorting
#
#         1   5     6   10    11  15
#         |---|     |---|     |---|
#-----|-----------------------|---|--
#     1   5   9     15  19    25  29
#         pair1     pair2     pair3
#
#
# create the third pair
# gene
my $match5 = Bio::Location::Simple->new
    (-seq_id => 'gene', -start => 25, -end => 29, -strand=>1 );
# exon
my $match6 = Bio::Location::Simple->new
    (-seq_id => 'exon3', -start => 11, -end => 15, -strand=>1 );

my $pair3 = Bio::Coordinate::Pair->new(-in => $match5,
                                       -out => $match6
                                      );

# create a new collection in wrong order
$transcribe = Bio::Coordinate::Collection->new;
$transcribe->add_mapper($pair3);
$transcribe->add_mapper($pair1);
$transcribe->add_mapper($pair2);
ok $transcribe->sort;
my @res;
map {push @res, $_->in->start } $transcribe->each_mapper;
ok compare_arrays ([5, 15, 25], \@res);


#
# Test using genomic data
#

my $mapper = Bio::Coordinate::Collection->new;

load_data($mapper, undef );

# transform a segment entirely within the first rawcontig
#test_transform ($mapper,
#               [627012, 2, 5, -1, "rawcontig"],
#               ["chr1", 2, 5, -1]);
$pos = Bio::Location::Simple->new (-start => 2, -end => 5, -strand => -1);
$res = $mapper->map($pos);
is $res->match->start, 2;
is $res->match->end, 5;
is $res->match->strand, -1;
is $res->match->seq_id, '627012';

## now a split coord
my @testres = (
             [314696, 31917, 31937, -1],
             [341, 126, 59773, -1],
             [315843, 5332, 5963, +1]
);
$pos = Bio::Location::Simple->new (-start => 383700, -end => 444000, -strand => 1);
$res = $mapper->map($pos);
 @res =  $res->each_match;
compare (shift @res, shift @testres);
compare (shift @res, shift @testres);
compare (shift @res, shift @testres);

## now a simple gap
@testres = (
            [627011, 7447, 7507, +1],
            ["chr1", 273762, 273781, 1]
           );
$pos = Bio::Location::Simple->new (-start => 273701, -end => 273781, -strand => 1);
$res = $mapper->map($pos);
is $res->each_match, 1;
is $res->each_gap, 1;
@res =  $res->each_Location;
compare (shift @res, shift @testres);
compare (shift @res, shift @testres);

ok $mapper->swap;
$pos = Bio::Location::Simple->new
    (-start => 2, -end => 5, -strand => -1, -seq_id => '627012');
$res = $mapper->map($pos);
is $res->match->start, 2;
is $res->match->end, 5;
is $res->match->strand, -1;
is $res->match->seq_id, 'chr1';

#
# tests for split locations
#

# testing a  simple pair
$match1 = Bio::Location::Simple->new
    (-seq_id => 'a', -start => 5, -end => 17, -strand=>1 );
$match2 = Bio::Location::Simple->new
    (-seq_id => 'b', -start => 1, -end => 13, -strand=>-1 );

$pair = Bio::Coordinate::Pair->new(-in => $match1,
                                   -out => $match2,
                                  );

# split location

ok my $split = Bio::Location::Split->new();
ok $split->add_sub_Location(Bio::Location::Simple->new(-start=>6,
                                                      -end=>8,
                                                      -strand=>1));
$split->add_sub_Location(Bio::Location::Simple->new(-start=>15,
                                                   -end=>16,
                                                   -strand=>1));

$res=$pair->map($split);
ok my @sublocs = $res->each_Location(1);
is @sublocs, 2;

#print Dumper \@sublocs;
is $sublocs[0]->start, 2;
is $sublocs[0]->end, 3;
is $sublocs[1]->start, 10;
is $sublocs[1]->end, 12;



#
# from Align
#

use Bio::Coordinate::Utils;
use Bio::LocatableSeq;
use Bio::SimpleAlign;

my $string;
#y $out = IO::String->new($string);

#AAA/3-10    --wtatgtng
#BBB/1-7     -aaaat-tt-

my $s1 = Bio::LocatableSeq->new(-id => 'AAA',
                                -seq => '--wtatgtng',
                                -start => 3,
                                -end => 10,
                                -alphabet => 'dna'
                               );
my $s2 = Bio::LocatableSeq->new(-id => 'BBB',
                                -seq => '-aaaat-tt-',
                                -start => 1,
                                -end => 7,
                                -alphabet => 'dna'
                               );
$a = Bio::SimpleAlign->new();
$a->add_seq($s1);
$a->add_seq($s2);
#use Data::Dumper;

ok my $uti = Bio::Coordinate::Utils->new;
$mapper = $uti->from_align($a);
#print Dumper $mapper;
is $mapper->return_match, 1;
is $mapper->return_match(1), 1;


$pos = Bio::Location::Simple->new
    (-start => 4, -end => 8, -strand => 1);
$res = $mapper->map($pos);
#print Dumper $res;

exit; # end of tests
#
# subroutines only after this
#

sub compare_arrays {
    my ($first, $second) = @_;

    return 0 unless @$first == @$second;
    for (my $i = 0; $i < @$first; $i++) {
        return 0 if $first->[$i] ne $second->[$i];
    }
    return 1;
}


sub compare {
    my ($match, $test) = @_;
    is $match->seq_id eq $test->[0], 1,
        "Match: |". $match->seq_id. "| Test: ". $test->[0]. "|";
    is $match->start, $test->[1];
    is $match->end, $test->[2];
    is $match->strand, $test->[3];
}


sub load_data {
    my ($map, $reverse) = @_;

#chr_name	raw_id	chr_start	chr_end	raw_start	raw_end	raw_ori
    my @sgp_dump = split ( /\n/, qq {
chr1	627012	1	31276	1	31276	1
chr1	627010	31377	42949	72250	83822	-1
chr1	2768	42950	180950	251	138251	1
chr1	10423	180951	266154	1	85204	-1
chr1	627011	266255	273761	1	7507	1
chr1	314698	273862	283122	1	9261	-1
chr1	627009	283223	331394	251	48422	-1
chr1	314695	331395	352162	1	20768	-1
chr1	314697	352263	359444	1	7182	-1
chr1	314696	359545	383720	31917	56092	-1
chr1	341	383721	443368	126	59773	-1
chr1	315843	443369	444727	5332	6690	1
chr1	315844	444828	453463	1	8636	-1
chr1	315834	453564	456692	1	3129	1
chr1	315831	456793	458919	1	2127	1
chr1	315827	459020	468965	251	10196	-1
chr1	544782	468966	469955	1	990	-1
chr1	315837	470056	473446	186	3576	-1
chr1	544807	473447	474456	1	1010	-1
chr1	315832	474557	477289	1	2733	1
chr1	544806	477390	477601	1086	1297	-1
chr1	315840	477602	482655	21	5074	1
chr1	544802	482656	483460	1	805	-1
chr1	544811	483561	484162	6599	7200	-1
chr1	315829	484163	498439	15	14291	-1
chr1	544813	498440	500980	1	2541	-1
chr1	544773	501081	502190	1217	2326	-1
chr1	315828	502191	513296	72	11177	1
chr1	544815	513297	517276	2179	6158	1
chr1	315836	517277	517662	2958	3343	1
chr1	544805	517663	520643	299	3279	1
chr1	315835	520744	521682	2462	3400	-1
chr1	544784	521683	526369	54	4740	1
chr1	544796	526470	527698	1	1229	1
chr1	315833	527799	528303	2530	3034	-1
chr1	544803	528304	531476	1	3173	-1
chr1	544821	531577	532691	1	1115	1
chr1	544810	532792	533843	1	1052	1
chr1	544800	533944	535249	1	1306	1
chr1	544786	535350	536652	1	1303	1
chr1	544814	536753	538358	1	1606	1
chr1	544812	538459	540004	1	1546	1
chr1	544818	540105	541505	1	1401	1
chr1	544816	541606	542693	1	1088	1
chr1	544778	542794	544023	1	1230	1
chr1	544779	544124	545709	1	1586	1
chr1	544804	545810	547660	1	1851	1
chr1	544774	547761	550105	1	2345	1
chr1	544817	550206	552105	1	1900	1
chr1	544781	552206	553640	1	1435	1
chr1	315830	553741	555769	1	2029	-1
chr1	544819	555870	558904	1	3035	-1
chr1	544777	559005	560670	1	1666	1
chr1	544795	560771	563092	1	2322	1
chr1	544809	563193	565523	1	2331	1
chr1	544808	565624	568113	1	2490	1
chr1	544798	568214	570324	1	2111	1
chr1	544783	570425	574640	1	4216	1
chr1	544824	574741	578101	1	3361	1
chr1	544775	578202	580180	1	1979	-1
chr1	544825	580281	581858	1	1578	-1
chr1	544772	581959	585312	1	3354	1
chr1	544793	585413	588740	1	3328	1
chr1	544785	588841	591656	1	2816	-1
chr1	544791	591757	594687	1	2931	1
chr1	544820	594788	597671	1	2884	1
chr1	544790	597772	601587	1	3816	1
chr1	544794	601688	603324	1	1637	-1
chr1	544823	603425	607433	1	4009	1
chr1	544789	607534	610856	1	3323	1
chr1	544799	610957	614618	1	3662	1
chr1	544776	614719	618674	1	3956	-1
chr1	544797	618775	624522	1	5748	-1
chr1	544787	624623	629720	1	5098	-1
chr1	544792	629821	637065	1	7245	1
chr1	622020	837066	851064	1	13999	-1
chr1	622021	851165	854101	1	2937	-1
chr1	622016	854202	856489	1	2288	-1
chr1	625275	856590	888524	420	32354	-1
chr1	622015	888525	891483	1	2959	-1
chr1	622024	891584	896208	8871	13495	-1
chr1	625537	896209	952170	1	55962	-1
chr1	625538	952271	1051812	251	99792	-1
chr1	625277	1051813	1055193	1	3381	-1
chr1	625266	1055294	1062471	1	7178	-1
chr1	598266	1062572	1086504	11	23943	-1
chr1	625271	1086505	1096571	3943	14009	1
chr1	625265	1096572	1100161	2436	6025	-1
chr1	173125	1100162	1106067	3329	9234	-1
chr1	598265	1106068	1112101	286	6319	1
chr1	625360	1112102	1172572	251	60721	1
chr1	173111	1172573	1172716	1	144	-1
chr1	173103	1172817	1173945	1	1129	1
chr1	170531	1174046	1174188	8791	8933	-1
chr1	625363	1174189	1183590	67	9468	1
chr1	173120	1183591	1183929	153	491	-1
chr1	170509	1183930	1184112	864	1046	1
chr1	173119	1184213	1189703	1	5491	-1
chr1	625357	1189804	1213915	1	24112	1
chr1	625359	1214016	1216330	1	2315	1
} );
    # test the auto-sorting feature
    #	@sgp_dump = reverse (@sgp_dump) if defined $reverse;

    my $first = 1;
    for my $line ( @sgp_dump ) {
        if( $first ) { $first = 0; next; }
        my ( $chr_name, $contig_id, $chr_start, $chr_end,
             $contig_start, $contig_end, $contig_strand ) =
                 split ( /\t/, $line );

        my $match1 = Bio::Location::Simple->new
            (-seq_id => $chr_name, -start => $chr_start,
             -end => $chr_end, -strand=>1 );
        my $match2 = Bio::Location::Simple->new
            (-seq_id => $contig_id, -start => $contig_start,
             -end => $contig_end, -strand=>$contig_strand );

        my $pair = Bio::Coordinate::Pair->new(-in => $match1,
                                              -out => $match2,
                                             );
        $map->add_mapper($pair);
    }
    return $map;
}
