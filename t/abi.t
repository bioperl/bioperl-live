# -*-Perl-*-
# $Id$
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($error $NUMTESTS);
BEGIN {
	$NUMTESTS = 3;
	$error = 0;
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test; };
	if ( $@ ) {
		use lib 't';
	}
	# SeqIO modules abi.pm, ctf.pm, exp.pm, pln.pm, ztr.pm
	# all require Bio::SeqIO::staden::read, part of bioperl-ext
	eval {
		require Bio::SeqIO::staden::read;
	};
	if ( $@ ) {
		$error = 1;
		warn "Bio::SeqIO::staden::read from bioperl-ext is not installed or is installed incorrectly - skipping abi.t tests\n";
   }
	use Test;
	plan tests => $NUMTESTS;
}

END { 
	foreach ( $Test::ntest..$NUMTESTS) {
		skip('Unable to run all of the abi tests',1);
   }
}

exit(0) if ( $error == 1 );

use Bio::SeqIO::abi;
use Bio::Root::IO;

my $verbose = $ENV{'BIOPERLDEBUG'};
ok(1);

my $io = Bio::SeqIO->new(-format => 'abi',
								 -verbose => $verbose,
								 -file => Bio::Root::IO->catfile
								 (qw(t data readtest.abi) ));
ok(my $seq = $io->next_seq);
ok($seq->seq, "GCNTATGACGTGGATTNCGAATTCTNNNNNCGGTAGNNGAAAATCCCCGGNCAAGNTTNNCCCTGCAAANGGAANAANNTGGCCGAGCGCTACGGGCTGATCTGGGTGTGCCTGTTTCCCCCGGCCGGGGGGAGNGATGCAGGACATCCAAGTATCCCGCCNATGGNGGGCTGAGGACGAGGACGGCTTCCATCAGATCAGTGTGCCCGGNCTTCGACATCGGCGGCAGCGCCGCGCGCCAACTGGAAGGCTTCATCGACGTGNAGCATTTTGNCTTCNTGCGCACCGCTACCTTCACCCANCCGGACAAGCGCNAANTGCNGNCCTACACCACCACNGAAACACCGACCGGNTTNAATGCCGATTACCTGAGNNGCGTGGCAAATTATTCGGNGGACNTGCCGCTGNCGGACGTGGACCCGAACTTCCAATGGCTGCGTCATTNCTAGGTGAATCTGCCTTTCACCGCCACGCTCACCATCCACTTCCCGGTGCCGGGCAAGCGGTTGGTGATNATGAATGCCGCCAGACCGGTGTCCAAGCACACCANCCGCCTGNTGGTGCCGATCGNCCGCTAATTTCGACACCCATCTGCCNGNGGGAAGACGTACATGNGTTCAACCTTGCACNTNGTTCNAAAAAAACCNTGCCATGGTGGNAANCGCAAGCGGNCCGGAAATATCNGCCGGNTTGACCCGCNTGNTTGGAAAGTGCATATTCCCCNCCGATNCNCAATTTCGAT");

