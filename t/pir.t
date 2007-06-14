# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

use strict;

BEGIN {
	eval { require Test::More; };
	if ( $@ ) {
		use lib 't/lib';
	}
	use Test::More;
	plan tests => 8;
}

use_ok('Bio::SeqIO', 'Bio::SeqIO can be used ' );

my $verbose = $ENV{'BIOPERLDEBUG'} || 0;

my $str = Bio::SeqIO->new(-file => Bio::Root::IO->catfile
								  ("t","data","seqfile.pir"),
								  -verbose => $verbose,
								  -format => 'pir');

ok ( defined $str, 'new instance is defined ');

my $out = Bio::SeqIO->new(-format => 'pir',
								 -fh => \*STDOUT);

while (my $seq = $str->next_seq()) {
	ok( $seq->length > 1, 'checked length');
	$out->write_seq($seq) if $verbose > 0;
}
