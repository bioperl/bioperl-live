# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

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

    plan tests => 6;
}
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

# make sure we can dogfood here
$aln = Bio::SimpleAlign->new();
open(IN,"t/out.pfam");
$aln->read_Pfam(\*IN);
close(IN);

ok ( $aln );

unlink('t/out.pfam', 't/out.aln_fasta');




