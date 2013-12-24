use strict;
use warnings;

BEGIN {
    use Bio::Root::Test;

    test_begin(-tests => 116);

    use_ok('Bio::Location::Simple');
    use_ok('Bio::Coordinate::Pair');
    use_ok('Bio::Coordinate::ExtrapolatingPair');
    use_ok('Bio::Coordinate::GeneMapper');
}

#
# Extrapolating pairs
#
#    No gaps returned, matches extrapolated
#     returns always a match or undef
#     -strict
#


# the  reverse strand pair
my $inr = Bio::Location::Simple->new(-start=>2, -end=>5, -strand=>1);
my $outr = Bio::Location::Simple->new(-start=>10, -end=>13, -strand=>-1);
ok my $pairr = Bio::Coordinate::ExtrapolatingPair->
    new(-in => $inr,
        -out => $outr
       );

my $posr = Bio::Location::Simple->new
    (-start => 3, -end => 4, -strand=> 1 );
my $resr = $pairr->map($posr);
is $resr->start, 11;
is $resr->end, 12;
is $resr->strand, -1;



# propepide
my $match1 = Bio::Location::Simple->new
    (-seq_id => 'propeptide', -start => 21, -end => 40, -strand=>1 );
# peptide
my $match2 = Bio::Location::Simple->new
    (-seq_id => 'peptide', -start => 1, -end => 20, -strand=>1 );

ok my $pair = Bio::Coordinate::ExtrapolatingPair->
    new(-in => $match1,
        -out => $match2,
        -strict => 1
       );

ok $pair->test;
is $pair->strand(), 1; #  = in->strand * out->strand
is $pair->in->seq_id(), 'propeptide';
is $pair->strict(), 1;

my ($count, $pos, $pos2, $res, $match, $res2);

# match within
$pos = Bio::Location::Simple->new
    (-start => 25, -end => 25, -strand=> -1 );
$res = $pair->map($pos);

isa_ok $res, 'Bio::Location::Simple';
is $res->start, 5;
is $res->end, 5;
is $res->strand, -1;
is $res->seq_id, 'peptide';


# match outside = undef
$pos = Bio::Location::Simple->new (-start => 5, -end => 5 );
$res = $pair->map($pos);

is $res, undef;

#
# partial match = match
#
$pos2 = Bio::Location::Simple->new
    (-start => 20, -end => 22, -strand=> -1 );

ok $res = $pair->map($pos2);

is $res->start, 0;
is $res->end, 2;
is $res->seq_id, 'peptide';
is $res->strand, -1;


#
# partial match2 =  match & gap
#
$pos2 = Bio::Location::Simple->new (-start => 40, -end => 41, -strand=> 1 );
ok $res = $pair->map($pos2);
is $res->start, 20;
is $res->end, 20;

#
#enveloping
#
$pos2 = Bio::Location::Simple->new (-start => 19, -end => 41, -strand=> 1 );
ok $res = $pair->map($pos2);
is $res->start, 1;
is $res->end, 20;

#
# testing the changing the strand
#

# chr
$match1 = Bio::Location::Simple->new
    (-seq_id => 'chr', -start => 21, -end => 40, -strand=>1 );
# gene
$match2 = Bio::Location::Simple->new
    (-seq_id => 'gene', -start => 1, -end => 20, -strand=>-1 );

 $pair = Bio::Coordinate::ExtrapolatingPair->
#my $pair = Bio::Coordinate::Pair->
    new(-in => $match1,
        -out => $match2,
        -strict => 0
       );

$pos = Bio::Location::Simple->new
    (-start => 38, -end => 40, -strand=> 1 );
$res = $pair->map($pos);
is $res->start, 1;
is $res->end, 3;
is $res->strand, -1;

$pos = Bio::Location::Simple->new
    (-start => 1, -end => 3, -strand=> 1 );
$res = $pair->map($pos);
is $res->start, 38;
is $res->end, 40;
is $res->strand, -1;


#
#
# Gene Mapper
#
#

ok my $m = Bio::Coordinate::GeneMapper->new(-in => 'propeptide',
                                            -out => 'peptide');
#$m->verbose(2);

is $m->peptide_offset(5), 5;


# match within
$pos = Bio::Location::Simple->new
    (-start => 25, -end => 25, -strand=> 1 );
$res = $m->map($pos);

is $res->start, 20;
is $res->end, 20;
is $res->strand, 1;
is $res->seq_id, 'peptide';


#
# nozero
#

# match within
$pos = Bio::Location::Simple->new
    (-start => 4, -end => 5, -strand=> 1 );
$res = $m->map($pos);
is $res->start, -1;
is $res->end, 0;

is $m->nozero('in&out'), 'in&out';
$res = $m->map($pos);
is $res->start, -2;
is $res->end, -1;
is $m->nozero(0), 0;



ok $m->swap;
$pos = Bio::Location::Simple->new
    (-start => 5, -end => 5, -strand=> 1 );
$res = $m->map($pos);
is $res->start, 10;

# cds -> propeptide
is $m->in('cds'), 'cds';
is $m->out('propeptide'), 'propeptide';

$res = $m->map($pos);
is $res->start, 2;
ok $res = $m->_translate($pos);
is $res->start, 2;
ok $res = $m->_reverse_translate($pos);
is $res->start, 13;
is $res->end, 15;

$pos = Bio::Location::Simple->new
    (-start => 26, -end => 26, -strand=> 1 );
$m->out('peptide');
$res = $m->map($pos);
is $res->start, 4;


#
# frame
#

$pos = Bio::Location::Simple->new
    (-start => 1, -end => 3, -strand=> 1 );
$res = $m->_frame($pos);
is $res->start, 1;
is $res->end, 3;


# Collection representing exons
#
#  cds    1   5     6   10    11  15
#  exon   1   5     1   5     1   5
#  gene   1   5    11   15   21   25
#         |---|     |---|     |---|
#-----|-----------------------|---|--
# chr 1   5   9    15   19   25   29
#         pair1     pair2     pair3

# gene
my $e1 = Bio::Location::Simple->new
    (-seq_id => 'gene', -start => 5, -end => 9, -strand=>1 );
my $e2 = Bio::Location::Simple->new
    (-seq_id => 'gene', -start => 15, -end => 19, -strand=>1 );
my $e3 = Bio::Location::Simple->new
    (-seq_id => 'gene', -start => 25, -end => 29, -strand=>1 );
my @cexons = ($e1, $e2, $e3);

$m= Bio::Coordinate::GeneMapper->new();

$m->in('chr');
$m->out('gene');
my $off = $m->cds(5);
is $off->start, 5; # start of the coding region
is $m->exons(@cexons), 3;

$m->out('exon');
$pos = Bio::Location::Simple->new
    (-start => 6, -end => 7, -strand=> 1 );
$res = $m->map($pos);

is $res->start, 2;
is $res->end, 3;

$m->out('negative_intron');
$pos = Bio::Location::Simple->new
    (-start => 12, -end => 14, -strand=> 1 );
$res = $m->map($pos);
is $res->start, -3;
is $res->end, -1;
is $res->seq_id, 'intron1';


# cds
$m->out('cds');
$pos = Bio::Location::Simple->new
    (-start => 5, -end => 9, -strand=> 1 );
$res = $m->map($pos);
is $res->start, 1;
is $res->end, 5;

$pos = Bio::Location::Simple->new
    (-start => 15, -end => 25, -strand=> 1 );
$res = $m->map($pos);
is $res->start, 6;
is $res->end, 11;

$pos = Bio::Location::Simple->new
    (-start => 5, -end => 19, -strand=> 1 );
$res = $m->map($pos);
is $res->start, 1;
is $res->end, 10;


#
# chr to cds ; ranges into one
#
my $exons = Bio::Location::Split->new(-seq_id => 'gene');
$exons->add_sub_Location($e1);
$exons->add_sub_Location($e2);
$exons->add_sub_Location($e3);

$res = $m->map($exons);
isa_ok $res,'Bio::Location::Simple';
is $res->start, 1;
is $res->end, 15;

#
# cds to chr; single range into two
#
$m->in('cds');
$m->out('gene');

$pos = Bio::Location::Simple->new
    (-start => 4, -end => 7, -strand=> 1 );
$res = $m->map($pos);
is $res->start, 4;
is $res->end, 12;



# Collection representing exons
#
#  cds  -11  -7    -6  -2    -1   3  :27
#  cds   -6  -2    -1 1 3     4   8  :17
#  exon   1   5     1   5     1   5
#  gene -21  -17  -11  -7    -1 1 3  :27
#  gene -11  -7    -1 1 3     9   13 :17
#         |---|     |---|     |---|
#-----|-----------------------|---|--
# chr 1   5   9    15   19   25   29
#         pair1     pair2     pair3

$m= Bio::Coordinate::GeneMapper->new();

$m->in('chr');
$m->out('gene');
$off = $m->cds(17);
is $off->start, 17; # start of the coding region
is $m->exons(@cexons), 3;

# testing parameter handling in the constructor
ok $m = Bio::Coordinate::GeneMapper->new(-in => 'gene',
                                         -out => 'peptide',
                                         -cds => 3,
                                         -exons => @cexons,
                                         -utr => 7,
                                         -peptide_offset => 5
                                        );


#
# Real life data
# Mapping SNPs into  human serum protein MSE55 and
# human galecting LGALS2 from Ensembl:
#

#Ensembl Gene ID	Exon Start (Chr bp)	Exon End (Chr bp)	Exon Coding Start (Chr bp)
#	Exon Coding End (Chr bp)	Strand

my @gene1_dump = split ( /\n/, qq {
ENSG00000128283	34571058	34571126			1
ENSG00000128283	34576610	34577350	34576888	34577350	1
ENSG00000128283	34578646	34579858	34578646	34579355	1
});


my @gene2_dump = split ( /\n/, qq {
ENSG00000100079	34590438	34590464			-1
ENSG00000100079	34582387	34582469	34582387	34582469	-1
ENSG00000100079	34581114	34581273	34581114	34581273	-1
ENSG00000100079	34580784	34580950	34580804	34580950	-1
}); # exon start should be less than end or is this intentional?

#Chromosome Name	Location (bp)	Strand	Reference ID
my @snp_dump = split ( /\n/, qq {
22	34572694	1	2235335
22	34572799	1	2235336
22	34572843	1	2235337
22	34574896	1	2076087
22	34575256	1	2076088
22	34578830	1	2281098
22	34579111	1	2281099
22	34580411	1	2235338
22	34580591	1	2281097
22	34580845	1	2235339
22	34581963	1	2281100
22	34583722	1	140057
22	34585003	1	140058
22	34587726	1	968725
22	34588207	1	2284055
22	34591507	1	1969639
22	34591949	1	140059
});
shift @snp_dump;

my ($cdsr, @exons) = read_gene_data(@gene1_dump);

ok my $g1 = Bio::Coordinate::GeneMapper->new(-in=>'chr', -out=>'gene');
$g1->cds($cdsr);

#$pos = Bio::Location::Simple->new
#    (-start => 34576888, -end => 34579355, -strand=> 1 );
$res = $g1->map($cdsr);
is $res->start, 1;
is $res->end, 2468;

$g1->exons(@exons);
$g1->in('gene');
$g1->out('cds');
$res = $g1->map($res);
is $res->start, 1;
is $res->end, 1173;

#map_snps($g1, @snp_dump);


#gene 2 in reverse strand
($cdsr, @exons) = read_gene_data(@gene2_dump);
ok my $g2 = Bio::Coordinate::GeneMapper->new(-in=>'chr', -out=>'gene');
$g2->cds($cdsr);

$pos = Bio::Location::Simple->new
    (-start => $cdsr->end-2, -end => $cdsr->end, -strand=> 1 );
$res = $g2->map($pos);
is $res->start, 1;
is $res->end, 3;
is $res->strand, -1;


$g2->exons(@exons);

#map_snps($g2, @snp_dump);


$match1 = Bio::Location::Simple->new
    (-seq_id => 'a', -start => 5, -end => 17, -strand=>1 );
$match2 = Bio::Location::Simple->new
    (-seq_id => 'b', -start => 1, -end => 13, -strand=>-1 );
ok $pair = Bio::Coordinate::Pair->new(-in => $match1,
                                      -out => $match2,
                                     );

#
# split location
#

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

# testing  cds -> gene/chr which generates a split location from a simple one
# exons in reverse strand!
#
#  pept   33222     111
#  cds    8   4     3 1-1
#  exon   5   1     5   1
#  gene  13   9     3 1-2
#         |---|     |---|
#-----|-------------------
# chr 1   5   9    15   19
#           e1        e2

# gene
$e1 = Bio::Location::Simple->new
    (-seq_id => 'gene', -start => 5, -end => 9, -strand=>-1 );
$e2 = Bio::Location::Simple->new
    (-seq_id => 'gene', -start => 15, -end => 19, -strand=>-1 );
@cexons = ($e1, $e2);
my $cds= Bio::Location::Simple->new
    (-seq_id => 'gene', -start => 5, -end => 17, -strand=>-1 );

$m = Bio::Coordinate::GeneMapper->new(-in=>'cds', -out=>'chr');

$m->cds($cds); # this has to be set first!?
is $m->exons(@cexons), 2;


my $cds_f= Bio::Location::Simple->new
    (-start => 2, -end => 7, );
$res = $m->map($cds_f);

ok @sublocs = $res->each_Location(1);
is @sublocs, 2;

is $sublocs[0]->start, 6;
is $sublocs[0]->end, 9;
is $sublocs[1]->start, 15;
is $sublocs[1]->end, 16;


# test inex, exon & negative_intron

$m->in('gene');
$m->out('inex');

$pos = Bio::Location::Simple->new
    (-seq_id => 'gene', -start => 2, -end => 10, -strand=> 1 );

$res = $m->map($pos);
is $res->each_Location, 3;


$m->out('intron');
$res = $m->map($pos);
is $res->start, 1;
is $res->end, 5;
is $res->strand, 1;

$m->out('negative_intron');
$res = $m->map($pos);
is $res->start, -5;
is $res->end, -1;
is $res->strand, 1;

is $m->_mapper_code2string('1-2'), 'chr-gene';
is $m->_mapper_string2code('chr-gene'), '1-2';


#todo:
#  strict mapping mode
#  extrapolating pair code into Bio::Coordinate::Pair ?






sub read_gene_data {
    my ($self,@gene_dump) = @_;
    my ($cds_start, $cds_end, $strand, @exons);

    #one line per exon
    my ($first, $first_line);
    for my $line ( @gene_dump ) {

        my ($geneid, $exon_start, $exon_end, $exon_cstart,
            $exon_cend, $exon_strand) = split /\t/, $line;

        $strand = $exon_strand if $exon_strand;
        #print join (' ', $geneid, $exon_start, $exon_strand), "\n";

        # CDS location in chromosome coordinates
        $cds_start = $exon_cstart if !$cds_start and $exon_cstart;
        $cds_end = $exon_cend if $exon_cend;


        if ($exon_start > $exon_end) {
            ($exon_start, $exon_end) = ($exon_end, $exon_start);
        }

        my $exon = Bio::Location::Simple->new
            (-seq_id => 'gene', -start => $exon_start,
             -end => $exon_end, -strand=>$strand, -verbose=>2);
        push @exons, $exon;
        }

    if ($cds_start > $cds_end) {
        ($cds_start, $cds_end) = ($cds_end, $cds_start);
    }

    my $cdsr = Bio::Location::Simple->new (-start => $cds_start,
                                           -end => $cds_end,
                                           -strand=> $strand);

    return ($cdsr, @exons);
}


sub map_snps {
    my ($mapper, @snps) =@_;
    $mapper->in('chr');
    $mapper->out('cds');
    foreach my $line (@snps) {
        $mapper->out('cds');

        my ($chr, $start, $strand, $id) = split /\t/, $line;
        my $loc = Bio::Location::Simple->new
            ( -start => $start,
             -end => $start, -strand=>$strand );

        my $res = $mapper->map($loc);
        my $cds_start = 0;
        $cds_start = $res->start if defined $res;#defined $res->start;
        print $id, "\t", $cds_start, "\n";

        # coding
        if ($cds_start) {
            $mapper->out('propeptide');
            my $frame_obj = $mapper->_frame($res);
            my $res = $mapper->map($loc);
            my $cds_start = 0;
            $cds_start = $res->start if defined $res;#defined $res->start;
            print  "\t\t", $cds_start, " (", $frame_obj->start, ")\n";
        }
    }
}
