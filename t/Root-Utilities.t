#-*-Perl-*-

# Test file for Bio::Root::Utilities
# Author: Torsten Seemann

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use lib './';
BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;

    plan tests => 36;
}
use Bio::Root::Utilities;

# Test

ok(1);

# Object creation

my $u = Bio::Root::Utilities->new();
ok($u);

# month2num() and num2month()

my @month = qw(XXX Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
for my $i (1 .. 12) {
  ok $u->month2num($month[$i]), $i;
  ok $u->num2month($i), $month[$i];
}

# untaint()

# this is the behaviour the function has, but should it return '' ?
ok $u->untaint(''), undef;
ok $u->untaint('nice string'), 'nice string';
ok $u->untaint('bad *?&^$! string'), 'bad ';
ok $u->untaint( q{100% relaxed&;''\"|*?!~<>^()[]{}$}, 1 ), '100% relaxed';

# mean_stdev()

my($mu,$sd);

($mu,$sd) = $u->mean_stdev();
ok $mu, undef;
ok $sd, undef;

($mu,$sd) = $u->mean_stdev(42);
ok $mu, 42;
ok $sd, undef;

($mu,$sd) = $u->mean_stdev(-1,0,1);
ok $mu, 0;
ok $sd, 1;

# other functions still need testing
# but many depend on the filesystem and access to binaries

