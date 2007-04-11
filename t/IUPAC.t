#-*-Perl-*-
## Bioperl Test Harness Script for Modules

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;

BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;
    plan tests => 4;
	use_ok('Bio::Tools::IUPAC');
	use_ok('Bio::Seq');
}

# test IUPAC

my $ambiseq = new Bio::Seq (-seq => 'ARTCGTTGR',
			    -alphabet => 'dna'); 

my $stream  = new Bio::Tools::IUPAC('-seq' => $ambiseq);
is $stream->count(), 4;

my $b = 1; 
while (my $uniqueseq = $stream->next_seq()) {
    if( ! $uniqueseq->isa('Bio::Seq') ) {
	$b = 0;
	last; # no point continuing if we get here
    }
}
ok $b;

