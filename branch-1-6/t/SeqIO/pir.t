# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 9);
	
	use_ok('Bio::SeqIO::pir');
}

my $verbose = test_debug();

my $str = Bio::SeqIO->new(-file => test_input_file('seqfile.pir'),
								  -verbose => $verbose,
								  -format => 'pir');

ok ( defined $str, 'new instance is defined ');
isa_ok ($str, 'Bio::SeqIO');

my $out = Bio::SeqIO->new(-format => 'pir',
								 -fh => \*STDOUT);

while (my $seq = $str->next_seq()) {
	ok( $seq->length > 1, 'checked length');
	$out->write_seq($seq) if $verbose > 0;
}
