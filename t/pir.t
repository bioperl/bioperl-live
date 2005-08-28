# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

use strict;

BEGIN {
	eval { require Test; };
	if ( $@ ) {
		use lib 't';
	}
	use Test;
	plan tests => 8;
}

use Bio::SeqIO;

ok(1);

my $verbose = $ENV{'BIOPERLDEBUG'} || 0;

my $str = Bio::SeqIO->new(-file => Bio::Root::IO->catfile
								  ("t","data","seqfile.pir"),
								  -verbose => $verbose,
								  -format => 'pir');
ok $str;
my $out = new Bio::SeqIO(-format => 'pir',
								 -fh => \*STDOUT);

while (my $seq = $str->next_seq()) {
	# ok($seq->id, qr /^[PF]1/ );
	ok($seq->length > 1);
	$out->write_seq($seq) if $verbose > 0;
}
