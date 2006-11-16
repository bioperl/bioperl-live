# -*-Perl-*-
# $Id$
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($error $NUMTESTS);
BEGIN {
	$NUMTESTS = 5;
	$error = 0;
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test::More; };
	if ( $@ ) {
		use lib 't/lib';
	}
	# SeqIO modules abi.pm, ctf.pm, exp.pm, pln.pm, ztr.pm
	# all require Bio::SeqIO::staden::read, part of bioperl-ext
	use Test::More;
	plan tests => $NUMTESTS;
	use_ok('Bio::SeqIO::abi');
	use_ok('Bio::Root::IO');
}

exit(0) if ( $error == 1 );



my $verbose = $ENV{'BIOPERLDEBUG'};

SKIP: {
	eval {
		require Bio::SeqIO::staden::read;
	};
	skip("Bio::SeqIO::staden::read from bioperl-ext is not installed or".
		 " is installed incorrectly - skipping abi.t tests\n", 3) if $@;
	require_ok('Bio::SeqIO::staden::read');
	my $io = Bio::SeqIO->new(-format => 'abi',
									 -verbose => $verbose,
									 -file => Bio::Root::IO->catfile
									 (qw(t data readtest.abi) ));
	my $seq = $io->next_seq;
	isa_ok($seq, 'Bio::PrimarySeqI');
	is($seq->seq, "GCNTATGACGTGGATTNCGAATTCTNNNNNCGGTAGNNGAAAATCCCCGGNCAAGNTTNNCCCTGCAAANGGAANAANNTGGCCGAGCGCTACGGGCTGATCTGGGTGTGCCTGTTTCCCCCGGCCGGGGGGAGNGATGCAGGACATCCAAGTATCCCGCCNATGGNGGGCTGAGGACGAGGACGGCTTCCATCAGATCAGTGTGCCCGGNCTTCGACATCGGCGGCAGCGCCGCGCGCCAACTGGAAGGCTTCATCGACGTGNAGCATTTTGNCTTCNTGCGCACCGCTACCTTCACCCANCCGGACAAGCGCNAANTGCNGNCCTACACCACCACNGAAACACCGACCGGNTTNAATGCCGATTACCTGAGNNGCGTGGCAAATTATTCGGNGGACNTGCCGCTGNCGGACGTGGACCCGAACTTCCAATGGCTGCGTCATTNCTAGGTGAATCTGCCTTTCACCGCCACGCTCACCATCCACTTCCCGGTGCCGGGCAAGCGGTTGGTGATNATGAATGCCGCCAGACCGGTGTCCAAGCACACCANCCGCCTGNTGGTGCCGATCGNCCGCTAATTTCGACACCCATCTGCCNGNGGGAAGACGTACATGNGTTCAACCTTGCACNTNGTTCNAAAAAAACCNTGCCATGGTGGNAANCGCAAGCGGNCCGGAAATATCNGCCGGNTTGACCCGCNTGNTTGGAAAGTGCATATTCCCCNCCGATNCNCAATTTCGAT");
}
