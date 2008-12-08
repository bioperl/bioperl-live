# -*-Perl-*- Test Harness script for Bioperl
# $Id$


use strict;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 4);
    
	use_ok('Bio::SearchIO');
}

my $searchio = Bio::SearchIO->new(
	 -format => 'blast',
    -file => test_input_file('blast.report')
);

my @hsps = ();
while(my $result = $searchio->next_result){
    while(my $hit = $result->next_hit){
        while(my $hsp = $hit->next_hsp){
            push @hsps, $hsp;
        }
    }
}

my $first_hsp = shift @hsps;
my $first_hsp_cigar_string = '504M'; 
is $first_hsp->cigar_string, $first_hsp_cigar_string;
is $first_hsp->cigar_string, $first_hsp_cigar_string; # fetch from hash

my $second_hsp = $hsps[0];
my $second_hsp_cigar_string = '29M18I22M11I20MD33M4I22M3I25M5I21MI33MD14M';
is $second_hsp->cigar_string, $second_hsp_cigar_string;
