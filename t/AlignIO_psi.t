# -*-Perl-*- Test Harness script for Bioperl
# $Id: AlignIO.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 5);
	
	use_ok('Bio::AlignIO');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# PSI format  
$str  = Bio::AlignIO->new(
    '-file'	=> test_input_file("testaln.psi"),
    '-format'	=> 'psi');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse, 'QUERY/1-798');
is($aln->no_sequences, 56);