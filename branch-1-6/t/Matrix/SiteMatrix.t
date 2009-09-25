# -*-Perl-*- Test Harness script for Bioperl
# $Id$


use strict;

BEGIN {
  use lib '.';
  use Bio::Root::Test;
  
  test_begin(-tests => 14);
  
  use_ok('Bio::Matrix::PSM::SiteMatrix');
}

my $score;
my $A='a0501';
my $C='014a0';
my $G='01103';
my $T='08006';
my $eval=0.0001;
my %param=(pA=>$A,pC=>$C,pG=>$G,pT=>$T,e_val=>$eval, correction =>0);
ok my $matrix=Bio::Matrix::PSM::SiteMatrix->new(%param);

#Simple methods here
is $matrix->IUPAC,'ABVCD';

is $matrix->consensus,'ATACT';

is $matrix->width,5;

is $matrix->curpos,0;

is $matrix->get_string('A'),$A;

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

is $matrix->curpos,1;

ok $matrix->e_val(0.0001);

is $matrix->e_val,0.0001;

#Now some PSM specific methods like regexp and matrix info
is $matrix->regexp,'[Aa][CcGgTtBb][AaCcGgVv][Cc][AaGgTtDd]';

my @x=(1,0,0.5,0,0.1);
is_deeply [$matrix->get_array('A')], \@x;

@x=qw([Aa] [CcGgTtBb] [AaCcGgVv] [Cc] [AaGgTtDd]);
is_deeply [$matrix->regexp_array], \@x;
