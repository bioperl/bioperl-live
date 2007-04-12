# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

use strict;
BEGIN {     
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;
    plan tests => 10;
	use_ok('Bio::SeqIO');
	use_ok('Bio::Root::IO');
}

my $in = Bio::SeqIO->new(-file => "<".Bio::Root::IO->catfile("t","data", 
	"multifa.seq") , '-format' => 'Fasta');
ok $in;
my $c=0;
while ( my $seq = $in->next_seq() ) {
    ok($seq);
    $c++;
}
is $c,6, "all sequences in the file";
