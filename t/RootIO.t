# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;

BEGIN {
    # to handle systems with no installed Test::More module
    # we include the t dir (where a copy of Test/More.pm is located)
    # as a fallback
    eval { require Test::More; };
    if( $@ ) {
	use lib 't';
    }
    use Test::More;
    plan tests => 31;
}

my $DEBUG = $ENV{BIOPERLDEBUG} || 0;
$| = 1;

use_ok('Bio::Root::IO');

my $obj = Bio::Root::IO->new();
ok defined($obj) && $obj->isa('Bio::Root::IO');

#############################################
# tests for exceptions/debugging/verbosity
#############################################

eval { $obj->throw('Testing throw') };
like $@, qr/Testing throw/, 'throw()'; # 'throw failed';

$obj->verbose(-1);
eval { $obj->throw('Testing throw') };
like $@, qr/Testing throw/, 'throw() verbose(-1)'; # 'verbose(-1) throw did not work properly' . $@;

eval { $obj->warn('Testing warn') };
ok !$@, 'warn()';

$obj->verbose(1);
eval { $obj->throw('Testing throw') };
like $@, qr/Testing throw/, 'throw() verbose(1)'; # 'verbose(1) throw did not work properly' . $@;

my @stack = $obj->stack_trace();
ok scalar @stack == 2, 'stack_trace()';

my $verbobj = Bio::Root::IO->new(-verbose=>1,-strict=>1);
ok $verbobj->verbose() == 1, 'set verbosity to 1';

ok $obj->verbose(-1);

#############################################
# tests for handle read and write abilities
#############################################

ok my $TESTINFILE = Bio::Root::IO->catfile(qw(t data test.waba));

my($handle,$file) = $obj->tempfile;
ok $handle;
ok $file;

#test with files

ok my $rio = Bio::Root::IO->new(-file=>$TESTINFILE);
ok $rio->mode eq 'r', 'filename, read';

ok my $wio = Bio::Root::IO->new(-file=>">$file");
ok $wio->mode eq 'w', 'filename, write';

# test with handles

ok open(my $I, $TESTINFILE);
ok open(my $O, '>', $file);

ok $rio = Bio::Root::IO->new(-fh=>$I);
ok $rio->mode eq 'r', 'handle, read';

ok $wio = Bio::Root::IO->new(-fh=>$O);
ok $wio->mode eq 'w', 'handle, write';

##############################################
# tests _pushback for multi-line buffering
##############################################

my $line1 = $rio->_readline;
my $line2 = $rio->_readline;

ok $rio->_pushback($line1);
ok $rio->_pushback($line2);

my $line3 = $rio->_readline;
my $line4 = $rio->_readline;
my $line5 = $rio->_readline;

ok $line1 eq $line3;
ok $line2 eq $line4;
ok $line5 ne $line4;

ok close($I);
ok close($O);

##############################################
# tests http retrieval
##############################################

SKIP: {
  skip "Skipping tests which require network access, set BIOPERLDEBUG=1 to test", 2 unless $DEBUG;

  my $TESTURL = 'http://www.google.com/index.html';
  
  ok $rio = Bio::Root::IO->new(-url=>$TESTURL), 'default -url method';
  
  if ($Bio::Root::IO::HAS_LWP) {
    $Bio::Root::IO::HAS_LWP = 0;
    ok $rio = Bio::Root::IO->new(-url=>$TESTURL), 'non-LWP -url method';
  } 
  else {
    ok 1, 'non-LWP -url method not needed as non-LWP was default';
  }
}  
  
1;
