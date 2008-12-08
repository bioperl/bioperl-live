# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 19);
	
	use_ok('Bio::LiveSeq::Mutation');
}


my $a = Bio::LiveSeq::Mutation->new();
ok defined $a;

$a->seq('aaa');
is $a->seq, 'aaa';

$a->seqori('ggg');
is $a->seqori, 'ggg';

$a->pos(-4);
is $a->pos, -4;

$a->pos(5);
is $a->pos, 5;

is ($a->len, 3);

$a->len(9);
is ($a->len, 9);

$a->transpos(55);
is $a->transpos, 55;

$a->issue(1);
is $a->issue, 1;

$a->label(57);
is $a->label, '57';

$a->prelabel(57);
is $a->prelabel, '57';

$a->postlabel(57);
is $a->postlabel, '57';

$a->lastlabel(57);
is $a->lastlabel, '57';

#constuctor test
$b = Bio::LiveSeq::Mutation->new('-seq'=>'AC',
				 '-seqori' => 'GG',
				 '-pos' => 5,
				 '-len' => 2,
				 );
ok  defined $b;
is $b->seqori, 'GG';
is $b->len, 2;
is $b->seq, 'AC';
is $b->pos, 5;
