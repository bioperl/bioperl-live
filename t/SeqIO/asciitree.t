# -*-Perl-*- Test Harness script for Bioperl
# $Id$

# `make test'. After `make install' it should work as `perl test.t'

use strict;

BEGIN {
	use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 2);
	
	use_ok('Bio::SeqIO');
}

my $verbose = test_debug();

# asciitree is a write-only format
my $in = Bio::SeqIO->new(-format => 'genbank',
						-verbose => $verbose,
						-file => test_input_file('AE003644_Adh-genomic.gb'));
my $seq = $in->next_seq;

my $out_file = test_output_file();
my $out = Bio::SeqIO->new(-file => ">$out_file",
								  -verbose => $verbose,
								  -format => 'asciitree');
$out->write_seq($seq);
ok (-s $out_file);
