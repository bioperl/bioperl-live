#-*-Perl-*-
## Bioperl Test Harness Script for Modules

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
    plan test => 1;
}

#use Bio::Tools::SeqAnal; # deprecated, don't use any more
use Bio::Tools::IUPAC;
use Bio::Seq;

# test IUPAC

my $ambiseq = new Bio::Seq (-seq => 'ARTCGTTGR', -type =>
			    'Dna'); 

my $stream  = new Bio::Tools::IUPAC('-seq' => $ambiseq);

my $b = 1; 
while (my $uniqueseq = $stream->next_seq()) {
    if( ! $uniqueseq->isa('Bio::Seq') ) {
	$b =0;
	last; # no point continuing if we get here
    }
}
ok $b;
