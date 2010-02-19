# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 6);
	
	use_ok('Bio::SeqIO::metafasta');
}

my $verbose = test_debug();

my $io = Bio::SeqIO->new(-format => 'metafasta',
								 -verbose => $verbose,
								 -file => test_input_file('test.metafasta'));

isa_ok($io, 'Bio::SeqIO');
ok(my $seq = $io->next_seq);
isa_ok($seq, 'Bio::Seq::Meta');
is($seq->seq, "ABCDEFHIJKLMNOPQRSTUVWXYZ");
is($seq->display_id,'test');
