# -*-Perl-*-
## Bioperl Test Harness Script for Modules

use strict;
BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) { 
        use lib 't';
    }
    use Test;

    plan tests => 2;
}

use Bio::PopGen::TagHaplotype;

my $hap = [
             [qw/0       0       0/],
             [qw/1       1       1/],
             [qw/0       0       1/],
             [qw/1       2       1/]
          ];

my $obj = Bio::PopGen::TagHaplotype -> new(-haplotype_block => $hap);


# check haplotype length 
ok( $obj->tag_length ,1); # Tag length for this set must be 1
 
# check the tag list
ok( (join ' ', @{($obj->tag_list)->[0]}) ,'1 2'); # the SNPs are 1 and 2 (1 2)
