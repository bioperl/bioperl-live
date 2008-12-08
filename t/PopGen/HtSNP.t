# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 8);
	
    use_ok('Bio::PopGen::HtSNP');
}

my $hap = [
     'acgt?cact',
     'acgt?ca-t',
     'cg?tag?gc',
     'cactcgtgc',
     'cgctcgtgc',
     'cggtag?gc',
     'ac?t?cact',
     ];

my $snp = [qw/s1 s2 s3 s4 s5 s6 s7 s8 s9/];

my $pop = [
     [qw/ uno    0.20/],
     [qw/ dos    0.20/],
     [qw/ tres   0.15/],
     [qw/ cuatro 0.15/],
     [qw/ cinco  0.10/],
     [qw/ seis   0.10/],
     [qw/ siete  0.10/],
       ];

my $obj = Bio::PopGen::HtSNP->new(-haplotype_block => $hap,
                                   -snp_ids         => $snp,
                                   -pattern_freq    => $pop,
);


# check lenght of the haplotype
is($obj->hap_length,9); # length of the haplotype must be 9 

# check silent SNPs
is( (join ' ', @{$obj->silent_snp}) ,'s4'); # the silent snp is in position 4 (counting from 1)

# check degenerated SNPs 
is( (join ' ', @{$obj->deg_snp}) ,'s7 s5 s3'); # degenerate SNPs 

# check useful SNP's
is( (join ' ', @{$obj->useful_snp}) ,'s1 s2 s6 s8 s9'); # degenerate SNPs 

# check the SNP code
is( (join ' ',@{$obj->snp_type_code}),'36 63 36 75 36'); # code for SNPs

# check the HtType 
is( (join ' ',@{$obj->ht_type}),'36 63 75'); # min snp_code 

my $tmp = $obj->deg_pattern();
my $err=0;

foreach my $family (keys %$tmp){
    if ($family eq '0'){
       unless ( (join ' ', @{$tmp->{$family}}) eq '0 6'){
           $err=1;
       }
    }
    if ($family eq '1'){
       unless ( (join ' ', @{$tmp->{$family}}) eq '1'){
           $err=1;
       }
    }
    if ($family eq '2'){
       unless ( (join ' ', @{$tmp->{$family}}) eq '2 4 5'){
           $err=1;
       }
    }
    if ($family eq '3'){
       unless ( (join ' ', @{$tmp->{$family}}) eq '3'){
           $err=1;
       }
    }
}

ok(! $err); # clustering degenerated haplotypes 
