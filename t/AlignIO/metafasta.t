# -*-Perl-*- Test Harness script for Bioperl
# $Id: megafasta.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 4);
	
	use_ok('Bio::AlignIO::metafasta');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# METAFASTA
#print STDERR "Better Metafasta tests needed\n" if $DEBUG;
my $io = Bio::AlignIO->new(-verbose => $DEBUG ? $DEBUG : -1, 
   -file => test_input_file('testaln.metafasta'));
$aln = $io->next_aln;
isa_ok($aln,'Bio::Align::AlignI');
is $aln->consensus_string,'CDEFHIJKLMNOPQRSTUVWXYZ', "consensus_string on metafasta";
is $aln->symbol_chars,'23',"symbol_chars() using metafasta";

