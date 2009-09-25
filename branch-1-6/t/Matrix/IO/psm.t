# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 63);
	
	use_ok ('Bio::Matrix::PSM::IO');
}

my $mmt= "chr04q	170164	170208	strong	-	0	Motif 3 occurrance in chr04q
chr04q	215755	215799	strong	+	0	Motif 4 occurrance in chr04q
chr04q	532530	532574	strong	+	2	Motif 2 occurrance in chr04q
chr04q	539492	539536	strong	-	1	Motif 1 occurrance in chr04q
chr04q	586113	586157	strong	+	2	Motif 2 occurrance in chr04q
chr04q	698245	698289	strong	-	0	Motif 4 occurrance in chr04q
chr04q	804412	804456	strong	-	0	Motif 3 occurrance in chr04q
chr04q	858870	858914	strong	-	2	Motif 3 occurrance in chr04q
chr04q	861561	861605	strong	-	2	Motif 3 occurrance in chr04q
chr04q	916898	916942	strong	-	1	Motif 1 occurrance in chr04q
chr04q	1146916	1146960	strong	-	0	Motif 1 occurrance in chr04q
chr04q	1315772	1315816	strong	+	1	Motif 1 occurrance in chr04q
chr04q	1636119	1636163	strong	+	2	Motif 3 occurrance in chr04q
chr04q	1636200	1636244	strong	+	2	Motif 1 occurrance in chr04q
chr04q	1636437	1636481	strong	+	2	Motif 4 occurrance in chr04q
chr04q	1637361	1637405	strong	+	2	Motif 2 occurrance in chr04q
chr04q	1652447	1652491	strong	+	1	Motif 4 occurrance in chr04q";
my @mmt=split(/\n/,$mmt);

#Let's try meme here
my $psmIO =  Bio::Matrix::PSM::IO->new(-format=>'meme', 
	     -file=>test_input_file('meme.dat'));
ok $psmIO;

my @inputfile=grep(/datafile/i,$psmIO->unstructured);
ok @inputfile;

my $release=$psmIO->release;
ok $release;

my @ids=$psmIO->hid;
is @ids,4;

my %weights=$psmIO->weight;
ok %weights;

my %seq = $psmIO->seq;
is %seq,'0';#Meme doesn't have seq

is $psmIO->version,'3.0';

my $psm = $psmIO->next_psm;
ok $psm;

#Lets try to compress and uncompress the log odds and the frequencies, see if there is no
#considerable loss of data.
my $fA=$psm->get_compressed_freq('A');
my @check=Bio::Matrix::PSM::SiteMatrix::_uncompress_string($fA,1,1);
my @A=$psm->get_array('A');
my ($var,$max) = (0,0);
for (my $i = 0; $i<@check;$i++) {
  my $diff=abs(abs($check[$i])-abs($A[$i]));
  $var += $diff;
  $max=$diff if ($diff>$max);
}
my $avg=$var/@check;
cmp_ok $avg,'<',0.01; #Loss of data under 1 percent
#print $avg,"\n";
is $psm->sequence_match_weight('CAGAAAAATAAAATGGCCACCACCC'),2015;

my $lA=$psm->get_compressed_logs('A');
@check=Bio::Matrix::PSM::SiteMatrix::_uncompress_string($lA,1000,2);
@A=$psm->get_logs_array('A');
($var,$max) = (0,0);
for (my $i = 0;$i<@check;$i++) {
  my $diff=abs(abs($check[$i])-abs($A[$i]));
  $var += $diff;
  $max=$diff if ($diff>$max);
}
$avg=$var/@check;
cmp_ok $avg,'<',10; #Loss of data under 1 percent

my $matrix=$psm->matrix;
ok $matrix;
my $psm2=$psm;
$psm2->matrix($matrix);
is $psm,$psm2;

my %psm_header=$psm->header;
is $psm_header{IC},38.1;
is $psm_header{sites},4;
is $psm_header{width},25;
is $psm_header{e_val},'1.2e-002';


#Quick check if returned object works
my $IUPAC=$psm->IUPAC;
is $IUPAC,'CMKWMAAAKWVAWTYCMCASCHCCM';
is $IUPAC,$psm2->IUPAC;
is $IUPAC,$matrix->IUPAC;

my $instances=$psm->instances;
ok $instances;

foreach my $instance (@{$instances}) {
  my $id=$instance->primary_id;
  is $instance->strand,1;
  last if (ok $id);
}

ok $psm->header('e_val');
#Meme parser should be OK if tests passed


#Now we are going to try transfac

$psmIO =  Bio::Matrix::PSM::IO->new(-format=>'transfac', 
	  -file=> test_input_file('transfac.dat'));
ok $psmIO;

my $version=$psmIO->version;
ok !$version;

is $psmIO->release, '6.4--2002-12-02';

$psm     = $psmIO->next_psm;
ok $psm;

# Lets try to compress and uncompress the the frequencies, see if
# there is no considerable loss of data.
$fA=$psm->get_compressed_freq('A');
@check=Bio::Matrix::PSM::SiteMatrix::_uncompress_string($fA,1,1);
@A=$psm->get_array('A');
($var,$max) = (0,0);
for (my $i = 0; $i<@check;$i++) {
  my $diff=abs(abs($check[$i])-abs($A[$i]));
  $var += $diff;
  $max=$diff if ($diff>$max);
}
$avg=$var/@check;
cmp_ok $avg,'<',0.01; #Loss of data under 1 percent

%weights = $psmIO->weight;
ok !$weights{''};

%seq     = $psmIO->seq;
is scalar keys %seq, 0;

#Quick check if returned object works
$IUPAC   = $psm->IUPAC;
is $IUPAC,'VVDCAKSTGBYD';

#Now we are going to try mast
$psmIO =  Bio::Matrix::PSM::IO->new(-format=>'mast', 
	  -file=>test_input_file('mast.dat'));
ok $psmIO;

@inputfile = grep(/datafile/i,$psmIO->unstructured);
ok !@inputfile;

is( $psmIO->release, '2002/04/02 0:11:59');

@ids     = $psmIO->hid;
is @ids,4;

%weights = $psmIO->weight;
ok !%weights; #Mast doesn't have weights

ok %seq    = $psmIO->seq;

foreach my $id ($psmIO->hid) {
    ok $seq{$id};
}
ok $psm=$psmIO->next_psm;

my %instances=$psmIO->instances;
ok %instances;

is $psmIO->version, '3.0';

my $mmastIO=Bio::Matrix::PSM::IO->new(-format=>'mast',-file=>test_input_file('mixedmast.dat'));

$psm = $mmastIO->next_psm; 
my $lastinstances = $psm->instances();
my $i=0;
foreach my $hit (@$lastinstances) {
    $hit -> end ( $hit-> start () + length ($hit->seq) - 1 ) ; # fix an old bug in InstanceSite.pm
    my $d=join("\t",$hit->{accession_number},$hit -> start () , $hit-> end (),$hit -> score (),
    $hit -> strand == 1 ? '+' : '-' , $hit -> frame,  $hit -> desc ( ));
    is $d,$mmt[$i];
    $i++;
    last if ($hit -> start == 1652447);
}
