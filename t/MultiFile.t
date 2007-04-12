# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

use strict;

BEGIN {
	eval { require Test::More; };
	if ( $@ ) {
		use lib 't/lib';
	}
	use Test::More;
	plan tests => 4;	
	use_ok('Bio::SeqIO');
	use_ok('Bio::SeqIO::MultiFile');
}

my $verbose = $ENV{'BIOPERLDEBUG'} || 0;

my $mf = Bio::SeqIO::MultiFile->new(-format => 'Fasta' ,
												-verbose => $verbose,
												-files =>
												[ Bio::Root::IO->catfile
												("t","data","multi_1.fa"),
												Bio::Root::IO->catfile
												("t","data","multi_2.fa")]);
ok defined $mf;
my $count = 0;
eval {
	while (my $seq = $mf->next_seq() ) {
		$count++;
	}
};
is( $count,12 );
