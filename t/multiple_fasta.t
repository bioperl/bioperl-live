# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##
# CVS Version
# $Id$


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use Test;
use strict;
BEGIN { plan tests => 5 }
use Bio::SeqIO;

my $in = Bio::SeqIO->new(-file => "<t/multifa.seq" , '-format' => 'Fasta');
ok $in;
my $c=0;
while ( my $seq = $in->next_seq() ) {
    ok($seq);
    $c++;
}
ok $c,3, " missing sequences in the file";
