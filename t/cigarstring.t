# written by Juguang Xiao
use strict;

BEGIN {
    use lib 't';
    use Test;
    plan tests => 3;
}

END { }

use Bio::SearchIO;

my $searchio = new Bio::SearchIO(
    -format => 'blast',
    -file => 't/data/blast.report'
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
ok $first_hsp->cigar_string, $first_hsp_cigar_string;
ok $first_hsp->cigar_string, $first_hsp_cigar_string; # fetch from hash

my $second_hsp = $hsps[0];
my $second_hsp_cigar_string = '29M18I22M11I20MD33M4I22M3I25M5I21MI33MD14M';
ok $second_hsp->cigar_string, $second_hsp_cigar_string;

