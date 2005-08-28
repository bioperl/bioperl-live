# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

use strict;

BEGIN {
	eval { require Test; };
	if ( $@ ) {
		use lib 't';
	}
	use Test;
	plan tests => 29;
}

use Bio::SeqIO;

ok(1);

# Set to -1 for release version, so warnings aren't printed
my $verbose = $ENV{'BIOPERLDEBUG'} || 0;

# Basic read and/or write tests for SeqIO. Specific tests for
# given module should go into their own file.

my @formats = qw(gcg fasta raw pir tab ace );
# The following files or formats are failing: swiss genbank interpro embl

foreach my $format (@formats) {
	print "======== $format ========\n" if $verbose;
	read_write($format);
}

sub read_write {
	my $format = shift;
	my $seq;
	my $str = Bio::SeqIO->new(-file=> Bio::Root::IO->catfile
									  ("t","data","test.$format"),
									  -format => $format);
	ok $seq = $str->next_seq();
	print "Sequence 1 of 2 from $format stream:\n", $seq->seq, "\n\n"
	  if  $verbose;
	unless ($format eq 'raw') {
		ok $seq->id, 'roa1_drome',"ID for format $format";
		ok $seq->length, 358;
	}

	unless ($format eq 'gcg') { # GCG file can contain only one sequence
		ok $seq = $str->next_seq();
		print "Sequence 2 of 2 from $format stream:\n", $seq->seq,
		  $seq->seq, "\n" if $verbose;
	}

	my $out = Bio::SeqIO->new(-file => ">". Bio::Root::IO->catfile
									  ("t","data","$format.out"),
									  -format => $format);
	ok $out->write_seq($seq);
	if ($format eq 'fasta') {
		my $id_type;
		ok($id_type = $out->preferred_id_type('accession.version'),
			'accession.version');
	}
}

END {
	map { unlink Bio::Root::IO->catfile("t","data","$_.out") } @formats
}
