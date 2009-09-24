# -*-Perl-*- Test Harness script for Bioperl
# $Id$

# `make test'. After `make install' it should work as `perl test.t'

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
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
my $out = Bio::SeqIO->new(-file => ">".$out_file,
						-verbose => $verbose,
						-format => 'asciitree');
$out->write_seq($seq);

# this is a bug and is failing on some systems like IRIX (not sure why, maybe
# File::Temp?)

if (-s $out_file) {
	ok(1, "File exists, has contents on ".$^O);
} else {
	TODO: {
        local $TODO = "Output doesn't exists on ".$^O;
		ok(-s $out_file);
	}
}
