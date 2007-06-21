# -*-Perl-*-
# $Id$

use strict;

our $out_file;

BEGIN {
	use File::Spec;
	$out_file = File::Spec->catfile(qw(t data tmp-chaosxml));
	
	use lib 't/lib';
	use BioperlTest;
	
	test_begin(-tests => 3,
			   -requires_modules => ['Data::Stag']);
	
	use_ok('Bio::SeqIO');
	use_ok('Bio::Root::IO');
}

END {
	unlink $out_file if -e $out_file;
}

my $verbose = test_debug();

# currently chaosxml is write-only
my $in = Bio::SeqIO->new(-format => 'genbank',
								 -verbose => $verbose,
								 -file => Bio::Root::IO->catfile
								 qw(t data AE003644_Adh-genomic.gb) );

my $seq = $in->next_seq;

my $out = Bio::SeqIO->new(-file => ">$out_file",
								  -verbose => $verbose,
								  -format => 'chaosxml');
$out->write_seq($seq);
ok (-e $out_file);
