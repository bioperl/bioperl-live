# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##

use Test;
use strict;

BEGIN { plan tests => 18 }

use Bio::SimpleAlign;
use Bio::AlignIO;

my ($str,$aln,$strout,$status);
## Now we test Bio::AlignIO::stockholm input
$str = Bio::AlignIO->new(-file=> 't/testaln.stockholm','-format' => 'stockholm');
ok defined($str) && ref($str) && $str->isa('Bio::AlignIO');
$aln = $str->next_aln();
ok $aln->{order}->{'0'}, '1433_LYCES-9-246';

## Now we test Bio::AlignIO::pfam

$str = Bio::AlignIO->new(-file=> 't/testaln.pfam');
ok defined($str) && ref($str) && $str->isa('Bio::AlignIO');
$aln = $str->next_aln();
ok $aln->{order}->{'0'}, '1433_LYCES-9-246', " failed pfam input test";

$strout = Bio::AlignIO->new(-file=> '>t/testout.pfam', '-format' => 'pfam');
$status = $strout->write_aln($aln);
ok $status, 1, " failed pfam output test";


## Now we test Bio::AlignIO::msf

$str = Bio::AlignIO->new(-file=> 't/testaln.msf');
$aln = $str->next_aln();
ok $aln->{order}->{'0'}, '1433_LYCES-9-246', " failed msf input test";



$strout = Bio::AlignIO->new(-file=> '>t/testout.msf', '-format' => 'msf');
$status = $strout->write_aln($aln);
ok $status, 1, "  failed msf output test";

## Now we test Bio::AlignIO::fasta

$str = Bio::AlignIO->new(-file=> 't/testaln.fasta', '-format' => 'fasta');
$aln = $str->next_aln();
ok $aln->{order}->{'0'}, 'AK1H_ECOLI-1-378', " failed fasta input test ";


$strout = Bio::AlignIO->new(-file=> '>t/testout.fasta', '-format' => 'fasta');
$status = $strout->write_aln($aln);
ok $status, 1,"  failed fasta output test";

## Now we test Bio::AlignIO::selex

$str = Bio::AlignIO->new(-file=> 't/testaln.selex','-format' => 'selex');
$aln = $str->next_aln();
ok $aln->{order}->{'0'}, 'AK1H_ECOLI-114-431', " failed selex format test ";

$strout = Bio::AlignIO->new(-file=> '>t/testout.selex', '-format' => 'selex');
$status = $strout->write_aln($aln);
ok $status, 1, "  failed selex output test";

## Now we test Bio::AlignIO::mase input
$str = Bio::AlignIO->new(-file=> 't/testaln.mase','-format' => 'mase');
$aln = $str->next_aln();
ok $aln->{order}->{'0'}, 'AK1H_ECOLI-1-378', " failed mase input test ";

## Now we test Bio::AlignIO::prodom input
$str = Bio::AlignIO->new(-file=> 't/testaln.prodom','-format' => 'prodom');
$aln = $str->next_aln();
ok $aln->{order}->{'0'}, 'P04777-1-33', " failed prodom input test ";

## Now we test Bio::AlignIO::clustalw output writing
$strout = Bio::AlignIO->new(-file=> '>t/testaln.clustal', '-format' => 'clustalw');
$status = $strout->write_aln($aln);
ok $status, 1, "  failed clustalw (.aln) output test";
undef $strout;
## Now we test Bio::AlignIO::clustalw input
$str = Bio::AlignIO->new('-file'=> 't/testaln.clustal', 
			 '-format' => 'clustalw');
$aln = $str->next_aln($aln);
ok $aln->{order}->{'0'}, 'Q40065-1-32', "  failed clustalw (.aln) output test - was " . $aln->{order}->{'0'};

# Testing filehandle manipulations

my $in  = Bio::AlignIO->newFh(-file => "t/testaln.fasta", '-format' => 'fasta');
my $out = Bio::AlignIO->newFh(-file => ">t/testout2.pfam", '-format' => 'pfam');
while ( $aln = <$in>) {
	ok $aln->{order}->{'0'}, 'AK1H_ECOLI-1-378', "  failed filehandle input test  ";
	$status = print $out $aln;
	last;
}
ok $status, 1, "  failed filehandle output test";

## Now we test Bio::AlignIO::bl2seq input
$str = Bio::AlignIO->new('-file'   => 't/bl2seq.out',
			 '-format' => 'bl2seq');
$aln = $str->next_aln();
ok $aln->{order}->{'1'}, 'ALEU_HORVU -60-360', "failed BLAST bl2seq format test";

