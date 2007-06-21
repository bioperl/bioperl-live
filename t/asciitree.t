# -*-Perl-*-
# $Id$
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;

our $out_file;

BEGIN {
	use File::Spec;
	$out_file = File::Spec->catfile(qw(t data tmp-asciitree));
	
	use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 3);
	
	use_ok('Bio::SeqIO::asciitree');
	use_ok('Bio::Root::IO');
}

my $verbose = test_debug();

# asciitree is a write-only format
my $in = Bio::SeqIO->new(-format => 'genbank',
						-verbose => $verbose,
						-file => Bio::Root::IO->catfile
						qw(t data AE003644_Adh-genomic.gb) );
my $seq = $in->next_seq;

my $out = Bio::SeqIO->new(-file => ">$out_file",
								  -verbose => $verbose,
								  -format => 'asciitree');
$out->write_seq($seq);
ok (-e $out_file);

END {
	unlink $out_file if -e $out_file;
}
