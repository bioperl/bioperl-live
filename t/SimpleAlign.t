# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'
use Test;
use strict;

BEGIN { plan tests => 6 }
use Bio::SimpleAlign;
ok(1);

open(FH,"t/test.mase") || die "Could not open test.mase $!";
my $aln = Bio::SimpleAlign->new();
$aln->read_mase(\*FH);
close(FH);

ok( $aln );
open(OUT,">t/out.aln_fasta"); 
$aln->write_fasta(\*OUT);
close(OUT);
ok(1);

$aln = Bio::SimpleAlign->new();
open(FH,"t/test.pfam");
$aln->read_Pfam(\*FH);
close(FH);

ok ( $aln );

open(OUT,">t/out.pfam"); 
$aln->write_Pfam(\*OUT);
close(OUT);
ok(1);

$aln = Bio::SimpleAlign->new();
open(IN,"t/out.pfam");
$aln->read_Pfam(\*IN);
close(IN);

ok ( $aln );




