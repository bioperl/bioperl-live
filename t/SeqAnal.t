#-*-Perl-*-
## Bioperl Test Harness Script for Modules

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use Test;
use strict;

BEGIN { plan test => 1 }

use Bio::Tools::SeqAnal;
use Bio::Tools::IUPAC;
use Bio::Seq;

# test IUPAC

my $ambiseq = new Bio::Seq (-seq => 'ARTCGTTGR', -type =>
			    'Dna'); 

my $stream  = new Bio::Tools::IUPAC($ambiseq);

my $b = 1; 
while (my $uniqueseq = $stream->next_seq()) {
    if( ! $uniqueseq->isa('Bio::Seq') ) {
	$b =0;
	last; # no point continuing if we get here
    }
}
ok $b;






