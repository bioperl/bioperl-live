
use strict;

BEGIN { 
    eval { require Test; };
    if( $@ ) { 
	use lib 't';
    }
    use Test;
    plan tests => 18 }

use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::Root::IO;

my ($str,$aln,$strout,$status);
$str = Bio::AlignIO->new(-file=> Bio::Root::IO->catfile("t","testaln.stockholm"),
			 '-format' => 'stockholm');
ok defined($str) && ref($str) && $str->isa('Bio::AlignIO');
$aln = $str->next_aln();
ok $aln->{order}->{'0'}, '1433_LYCES-9-246';


$str = Bio::AlignIO->new(-file=> Bio::Root::IO->catfile("t","testaln.pfam"));
ok defined($str) && ref($str) && $str->isa('Bio::AlignIO');
$aln = $str->next_aln();
ok $aln->{order}->{'0'}, '1433_LYCES-9-246', " failed pfam input test";

$strout = Bio::AlignIO->new(-file=> ">".Bio::Root::IO->catfile("t","testout.pfam"), '-format' => 'pfam');
$status = $strout->write_aln($aln);
ok $status, 1, " failed pfam output test";


$str = Bio::AlignIO->new(-file=> Bio::Root::IO->catfile("t","testaln.msf"));
$aln = $str->next_aln();
ok $aln->{order}->{'0'}, '1433_LYCES-9-246', " failed msf input test";



$strout = Bio::AlignIO->new(-file=> ">".Bio::Root::IO->catfile("t","testout.msf"), '-format' => 'msf');
$status = $strout->write_aln($aln);
ok $status, 1, "  failed msf output test";


$str = Bio::AlignIO->new(-file=> Bio::Root::IO->catfile("t","testaln.fasta"), '-format' => 'fasta');
$aln = $str->next_aln();
ok $aln->{order}->{'0'}, 'AK1H_ECOLI-114-431', " failed fasta input test ";


$strout = Bio::AlignIO->new(-file=> ">".Bio::Root::IO->catfile("t","testout.fasta"), '-format' => 'fasta');
$status = $strout->write_aln($aln);
ok $status, 1,"  failed fasta output test";


$str = Bio::AlignIO->new(-file=> Bio::Root::IO->catfile("t","testaln.selex"),'-format' => 'selex');
$aln = $str->next_aln();
ok $aln->{order}->{'0'}, 'AK1H_ECOLI-114-431', " failed selex format test ";

$strout = Bio::AlignIO->new(-file=> ">".Bio::Root::IO->catfile("t","testout.selex"), '-format' => 'selex');
$status = $strout->write_aln($aln);
ok $status, 1, "  failed selex output test";

$str = Bio::AlignIO->new(-file=> Bio::Root::IO->catfile("t","testaln.mase"),'-format' => 'mase');
$aln = $str->next_aln();
ok $aln->{order}->{'0'}, 'AK1H_ECOLI-1-318', " failed mase input test ";

$str = Bio::AlignIO->new(-file=> Bio::Root::IO->catfile("t","testaln.prodom"),'-format' => 'prodom');
$aln = $str->next_aln();
ok $aln->{order}->{'0'}, 'P04777-1-33', " failed prodom input test ";

$strout = Bio::AlignIO->new(-file=> ">".Bio::Root::IO->catfile("t","testaln.clustal"), 
			    '-format' => 'clustalw');
$status = $strout->write_aln($aln);
ok $status, 1, "  failed clustalw (.aln) output test";
undef $strout;
$str = Bio::AlignIO->new('-file'=> Bio::Root::IO->catfile("t","testaln.clustal"), 
			 '-format' => 'clustalw');
$aln = $str->next_aln($aln);
ok $aln->{order}->{'0'}, 'P04777-1-33', "  failed clustalw (.aln) output test - was " . $aln->{order}->{'0'};


my $in  = Bio::AlignIO->newFh(-file => Bio::Root::IO->catfile("t","testaln.fasta"), '-format' => 'fasta');
my $out = Bio::AlignIO->newFh(-file => ">".Bio::Root::IO->catfile("t","testout2.pfam"), '-format' => 'pfam');
while ( $aln = <$in>) {
	ok $aln->{order}->{'0'}, 'AK1H_ECOLI-114-431', "  failed filehandle input test  ";
	$status = print $out $aln;
	last;
}
ok $status, 1, "  failed filehandle output test";

$str = Bio::AlignIO->new('-file'   => Bio::Root::IO->catfile("t","bl2seq.out"),
			 '-format' => 'bl2seq');
$aln = $str->next_aln();
ok ($aln->{order}->{'1'}, 'ALEU_HORVU -60-360', 
    "failed BLAST bl2seq format test");

unlink(Bio::Root::IO->catfile("t","testout2.pfam"),Bio::Root::IO->catfile("t","testout.selex"),
       Bio::Root::IO->catfile("t","testout.pfam"),Bio::Root::IO->catfile("t","testout.msf"),
       Bio::Root::IO->catfile("t","testout.fasta"), Bio::Root::IO->catfile("t","testaln.clustal"));
