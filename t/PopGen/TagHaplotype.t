# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 3);
	
    use_ok('Bio::PopGen::TagHaplotype');
}

my $hap = [
             [qw/0       0       0/],
             [qw/1       1       1/],
             [qw/0       0       1/],
             [qw/1       2       1/]
          ];

my $obj = Bio::PopGen::TagHaplotype->new(-haplotype_block => $hap);


# check haplotype length 
is( $obj->tag_length ,1); # Tag length for this set must be 1
 
# check the tag list
is( (join ' ', @{($obj->tag_list)->[0]}) ,'1 2'); # the SNPs are 1 and 2 (1 2)

# TODO? is this really enough tests?!
