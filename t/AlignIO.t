# This is -*-Perl-*- code
# $Id$
use strict;

BEGIN { 
    eval { require Test; };
    if( $@ ) { 
	use lib 't';
    }
    use Test;
    plan tests => 25 }

use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::Root::IO;
ok 1;

my ($str,$aln,$strout,$status);

# STOCKHOLM
$str = Bio::AlignIO->new('-file' => Bio::Root::IO->catfile("t","data","testaln.stockholm"),
			   '-format' => 'stockholm');
ok defined($str) && ref($str) && $str->isa('Bio::AlignIO');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, '1433_LYCES/9-246';

# PFAM
$str = Bio::AlignIO->new('-file' => Bio::Root::IO->catfile("t","data","testaln.pfam"));
ok defined($str) && ref($str) && $str->isa('Bio::AlignIO');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, '1433_LYCES/9-246', " failed pfam input test";

$strout = Bio::AlignIO->new('-file' => ">".Bio::Root::IO->catfile("t","data","testout.pfam"), 
			      '-format' => 'pfam');
$status = $strout->write_aln($aln);
ok $status, 1, " failed pfam output test";

# MSF
$str = Bio::AlignIO->new('-file' => Bio::Root::IO->catfile("t","data","testaln.msf"));
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, '1433_LYCES/9-246', " failed msf input test";


$strout = Bio::AlignIO->new('-file' => ">".Bio::Root::IO->catfile("t","data","testout.msf"), 
			      '-format' => 'msf');
$status = $strout->write_aln($aln);
ok $status, 1, "  failed msf output test";

# FASTA
$str = Bio::AlignIO->new('-file' => Bio::Root::IO->catfile("t","data","testaln.fasta"), 
			   '-format' => 'fasta');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, 'AK1H_ECOLI/114-431', " failed fasta input test ";

$strout = Bio::AlignIO->new('-file' => ">".Bio::Root::IO->catfile("t","data","testout.fasta"), 
			      '-format' => 'fasta');
$status = $strout->write_aln($aln);
ok $status, 1,"  failed fasta output test";

# SELEX
$str = Bio::AlignIO->new('-file' => Bio::Root::IO->catfile("t","data","testaln.selex"),
			   '-format' => 'selex');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, 'AK1H_ECOLI/114-431', " failed selex format test ";

$strout = Bio::AlignIO->new('-file' => ">".Bio::Root::IO->catfile("t","data","testout.selex"), 
			      '-format' => 'selex');
$status = $strout->write_aln($aln);
ok $status, 1, "  failed selex output test";

# MASE
$str = Bio::AlignIO->new('-file' => Bio::Root::IO->catfile("t","data","testaln.mase"),
			   '-format' => 'mase');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, 'AK1H_ECOLI/1-318', " failed mase input test ";

# PRODOM
$str = Bio::AlignIO->new('-file' => Bio::Root::IO->catfile("t","data","testaln.prodom"),
			   '-format' => 'prodom');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, 'P04777/1-33', " failed prodom input test ";

# CLUSTAL
$strout = Bio::AlignIO->new('-file' => ">".Bio::Root::IO->catfile("t","data","testaln.clustal"), 
			      '-format' => 'clustalw');
$status = $strout->write_aln($aln);
ok $status, 1, "  failed clustalw (.aln) output test";
undef $strout;
$str = Bio::AlignIO->new('-file'=> Bio::Root::IO->catfile("t","data","testaln.clustal"), 
			   '-format' => 'clustalw');
$aln = $str->next_aln($aln);
ok $aln->get_seq_by_pos(1)->get_nse, 'P04777/1-33', "  failed clustalw (.aln) output test";

# FASTA
my $in = Bio::AlignIO->newFh('-file'  => Bio::Root::IO->catfile("t","data","testaln.fasta"), 
			       '-format' => 'fasta');
my $out = Bio::AlignIO->newFh('-file'  => ">".Bio::Root::IO->catfile("t","data","testout2.pfam"), 
				'-format' => 'pfam');
while ( $aln = <$in>) {
    ok $aln->get_seq_by_pos(1)->get_nse, 'AK1H_ECOLI/114-431', "  failed filehandle input test  ";
    $status = print $out $aln;
    last;
}
ok $status, 1, "  failed filehandle output test";

# BL2SEQ
$str = Bio::AlignIO->new('-file'   => Bio::Root::IO->catfile("t","data","bl2seq.out"),
			   '-format' => 'bl2seq');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(2)->get_nse, 'ALEU_HORVU/60-360', 
    "failed BLAST bl2seq format test";

# PHYLIP
$str = Bio::AlignIO->new('-file' => Bio::Root::IO->catfile("t","data","testaln.phylip"),
			  '-format' => 'phylip');
ok defined($str) && ref($str) && $str->isa('Bio::AlignIO');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, 'Homo_sapie/1-45';

$strout = Bio::AlignIO->new('-file'  => ">".Bio::Root::IO->catfile("t","data","testout.phylip"),
			      '-format' => 'phylip');
$status = $strout->write_aln($aln);
ok $status, 1, "  failed phylip output test";

# NEXUS
$str = Bio::AlignIO->new('-file' => Bio::Root::IO->catfile("t","data","testaln.nexus"),
			  '-format' => 'nexus');
ok defined($str) && ref($str) && $str->isa('Bio::AlignIO');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, 'Homo_sapiens/1-45';
$strout = Bio::AlignIO->new('-file'  => ">".Bio::Root::IO->catfile("testout.nexus"),
			    '-format' => 'nexus', );
$status = $strout->write_aln($aln);
ok $status, 1, "  failed nexus output test";


unlink(Bio::Root::IO->catfile("t","data","testout2.pfam"),
     Bio::Root::IO->catfile("t","data","testout.selex"),
     Bio::Root::IO->catfile("t","data","testout.pfam"),
     Bio::Root::IO->catfile("t","data","testout.msf"),
     Bio::Root::IO->catfile("t","data","testout.fasta"), 
     Bio::Root::IO->catfile("t","data","testaln.clustal"),
     Bio::Root::IO->catfile("t","data","testout.phylip"),
     Bio::Root::IO->catfile("t","data","testout.nexus")
       );

