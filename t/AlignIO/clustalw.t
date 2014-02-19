# -*-Perl-*- Test Harness script for Bioperl
# $Id: clustalw.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 6);
	
	use_ok('Bio::AlignIO::clustalw');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# CLUSTAL
my $io = Bio::AlignIO->new(
   -file => test_input_file("testaln.clustalw") );
$aln = $io->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->consensus_string, "MNEGEHQIKLDELFEKLLRARKIFKNKDVLRHSWEPKDLPHRHEQIEA".
"LAQILVPVLRGETMKIIFCGHHACELGEDRGTKGFVIDELKDVDEDRNGKVDVIEINCEHMDTHYRVLPNIAKLF".
"DDCTGIGVPMHGGPTDEVTAKLKQVIDMKERFVIIVLDEIDKLVKKSGDEVLYSLTRINTELKRAKVSVIGISND".
"LKFKEYLDPRVLSSLSEEEVVFPPYDANQLRDILTQRAEEAFYPGVLDEGVIPLCAALAAREHGDARKALDLLRV".
"AGEIAEREGASKVTEKHVWKAQEKIEQDMMEEVIKTLPLQSKVLLYAIVLLDENGDLPANTGDVYAVYRELCEYI".
"DLEPLTQRRISDLINELDMLGIINAKVVSKGRYGRTKEIRLMVTSYKIRNVLRYDYSIQPLLTISLKSEQRRLI",
"clustalw consensus_string test";

my $outfile = test_output_file();
$strout = Bio::AlignIO->new(
   '-file' => ">$outfile", 
			      '-format' => 'clustalw');
$status = $strout->write_aln($aln);
is $status, 1, "clustalw output test";
undef $strout;
$str = Bio::AlignIO->new(
   '-file'=> $outfile, 
			   '-format' => 'clustalw');
$aln = $str->next_aln($aln);
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'P84139/1-420', "clustalw input test";
