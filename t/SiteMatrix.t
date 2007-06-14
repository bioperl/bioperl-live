# -*-Perl-*-
# $Id$
#Some simple test, nothing fancy...

use strict;

CHECK {
  $ENV{PERL_HASH_SEED} = 0;
}

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;

    plan tests => 15;
}
use Bio::Matrix::PSM::SiteMatrix;

ok(1);

my $score;
my $A='a0501';
my $C='014a0';
my $G='01103';
my $T='08006';
my $eval=0.0001;
my %param=(pA=>$A,pC=>$C,pG=>$G,pT=>$T,e_val=>$eval, correction =>0);
my $matrix=Bio::Matrix::PSM::SiteMatrix->new(%param);
ok $matrix;

#Simple methods here
ok $matrix->IUPAC,'ABVCD';

ok $matrix->consensus,'ATACT';

ok $matrix->width,5;

ok $matrix->curpos,0;

ok $matrix->get_string('A'),$A;

my %x= (base=>'A',pA=>1,pC=>0,pG=>0,pT=>0,prob=>10,rel=>0, 
        lA=>undef,lC=>undef,lG=>undef,lT=>undef);
my %pos = $matrix->next_pos;
my ($all) = 1;
while(my ($k,$v) = each %x ) {
    my $r =$pos{$k};
    if( ! defined $v && ! defined $r) {
    } elsif($pos{$k} ne $v ) { 
	$all = 0;
	last;
    }
}
ok($all);

ok $matrix->curpos,1;

ok $matrix->e_val(0.0001);

ok $matrix->e_val,0.0001;

#Now some PSM specific methods like regexp and matrix info
ok $matrix->regexp,'[Aa][CcGgTtBb][AaCcGgVv][Cc][AaGgTtDd]';
my $regexp=$matrix->regexp;
ok 'ATCCT',"/$regexp/";

my @x=(1,0,0.5,0,0.1);
ok $matrix->get_array('A'),@x;

@x=qw(Aa Tt AaCc Cc GgTt);
ok $matrix->regexp_array,@x;

