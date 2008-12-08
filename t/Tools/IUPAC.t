# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 4);
	
	use_ok('Bio::Tools::IUPAC');
	use_ok('Bio::Seq');
}

# test IUPAC

my $ambiseq = Bio::Seq->new(-seq => 'ARTCGTTGR',
			    -alphabet => 'dna'); 

my $stream  = Bio::Tools::IUPAC->new('-seq' => $ambiseq);
is $stream->count(), 4;

my $b = 1; 
while (my $uniqueseq = $stream->next_seq()) {
    if( ! $uniqueseq->isa('Bio::Seq') ) {
	$b = 0;
	last; # no point continuing if we get here
    }
}
ok $b;
