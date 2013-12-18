# -*-Perl-*- Test Harness script for Bioperl
# $Id: RootIO.t 16840 2010-02-16 17:14:12Z cjfields $

use strict;
use warnings;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 29,
	       -requires_networking => 1);
	
    use_ok('Bio::Root::HTTPget');
}

my $TESTURL = 'http://www.google.com/index.html';

my $TEST_PROXY = 'http://myproxy';

my @TEST_AUTHENTICATION = qw(foo bar);

my ($fh, $proxy);

my @auth;

=head1 Bio::Root::HTTPget comments

This module is a bit schizophrenic in that it is called in three different
ways; as an instance method, a class method, or as an explicit subroutine.

These tests check for all call types.  They are by no means incomplete.

=cut 

# test object method calls
my $obj = Bio::Root::HTTPget->new();

ok defined($obj) && $obj->isa('Bio::Root::Root');

lives_ok {$obj->get($TESTURL)};
lives_ok {$fh = $obj->getFH($TESTURL)};
isa_ok($fh, 'IO::Socket::INET');

undef($fh);

is ($obj->proxy(), undef);
is_deeply([$obj->authentication], []);
$obj->proxy('http', $TEST_PROXY);
$obj->authentication(@TEST_AUTHENTICATION);
is ($obj->proxy(), $TEST_PROXY);
is_deeply([$obj->authentication], \@TEST_AUTHENTICATION);

# test class method calls; note that mixing class and sub calls pollutes the
# class attributes

lives_ok {Bio::Root::HTTPget->get($TESTURL)};
lives_ok {$fh = Bio::Root::HTTPget->getFH($TESTURL)};
isa_ok($fh, 'IO::Socket::INET');

undef($fh);

is (Bio::Root::HTTPget->proxy(), undef);
is_deeply([Bio::Root::HTTPget->authentication], []);
Bio::Root::HTTPget->proxy('http', $TEST_PROXY);
Bio::Root::HTTPget->authentication(@TEST_AUTHENTICATION);
is (Bio::Root::HTTPget->proxy('http'), $TEST_PROXY);
is_deeply([Bio::Root::HTTPget->authentication], \@TEST_AUTHENTICATION);

# test sub calls (not called as method)

lives_ok {Bio::Root::HTTPget::get($TESTURL)};
lives_ok {$fh = Bio::Root::HTTPget::getFH($TESTURL)};
isa_ok($fh, 'IO::Socket::INET');

undef($fh);

# note that mixing class and sub calls pollutes the class attributes, have to
# manually reset
Bio::Root::HTTPget->authentication(undef, undef);

my $old = Bio::Root::HTTPget->clear_proxy('http');
is (Bio::Root::HTTPget::proxy(), undef);
is ($old, $TEST_PROXY);

is_deeply([Bio::Root::HTTPget->authentication], [undef, undef]);
Bio::Root::HTTPget::proxy('http', $TEST_PROXY);
Bio::Root::HTTPget::authentication(@TEST_AUTHENTICATION);
is (Bio::Root::HTTPget::proxy('http'), $TEST_PROXY);
is_deeply([Bio::Root::HTTPget->authentication], \@TEST_AUTHENTICATION);

# check to make sure new instance attributes are not polluted by class attrbutes
# from previous tests

my $newobj = Bio::Root::HTTPget->new();

ok defined($newobj) && $obj->isa('Bio::Root::Root');

is ($newobj->proxy(), undef);
is_deeply([$newobj->authentication], []);
$newobj->proxy('http', $TEST_PROXY);
$newobj->authentication(@TEST_AUTHENTICATION);
is ($newobj->proxy(), $TEST_PROXY);
is_deeply([$newobj->authentication], \@TEST_AUTHENTICATION);
