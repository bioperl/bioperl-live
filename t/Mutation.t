# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    plan tests => 18;
}

use Bio::LiveSeq::Mutation;

my $a = Bio::LiveSeq::Mutation->new();
ok defined $a;

$a->seq('aaa');
ok $a->seq, 'aaa';

$a->seqori('ggg');
ok $a->seqori, 'ggg';

$a->pos(-4);
ok $a->pos, -4;

$a->pos(5);
ok $a->pos, 5;

ok ($a->len, 3);

$a->len(9);
ok ($a->len, 9);

$a->transpos(55);
ok $a->transpos, 55;

$a->issue(1);
ok $a->issue, 1;

$a->label(57);
ok $a->label, '57';

$a->prelabel(57);
ok $a->prelabel, '57';

$a->postlabel(57);
ok $a->postlabel, '57';

$a->lastlabel(57);
ok $a->lastlabel, '57';

#constuctor test
$b = Bio::LiveSeq::Mutation->new('-seq'=>'AC',
				 '-seqori' => 'GG',
				 '-pos' => 5,
				 '-len' => 2,
				 );
ok  defined $b;
ok $b->seqori, 'GG';
ok $b->len, 2;
ok $b->seq, 'AC';
ok $b->pos, 5;


