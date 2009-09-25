# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 3);
		
	use_ok('Bio::SeqIO::MultiFile');
}

my $verbose = test_debug();

my $mf = Bio::SeqIO::MultiFile->new(-format => 'Fasta' ,
												-verbose => $verbose,
												-files =>
												[ test_input_file('multi_1.fa'),
												  test_input_file('multi_2.fa')]);
ok defined $mf;
my $count = 0;
eval {
	while (my $seq = $mf->next_seq() ) {
		$count++;
	}
};
is( $count,12 );
