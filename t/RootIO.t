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
	use lib 't', '.';
    }
    use Test;    
    plan tests => 14;
}

use Bio::Root::IO;

my $obj = new Bio::Root::IO();
ok defined($obj) && $obj->isa('Bio::Root::IO');

eval { $obj->throw('Testing throw') };
ok $@ =~ /Testing throw/;# 'throw failed';

$obj->verbose(-1);
eval { $obj->throw('Testing throw') };
ok $@=~ /Testing throw/;# 'verbose(-1) throw did not work properly' . $@;

eval { $obj->warn('Testing warn') };
ok !$@;

$obj->verbose(1);
eval { $obj->throw('Testing throw') };
ok $@ =~ /Testing throw/;# 'verbose(1) throw did not work properly' . $@;

my @stack = $obj->stack_trace();
ok scalar @stack, 2;

my $verbobj = new Bio::Root::IO(-verbose=>1,-strict=>1);
ok $verbobj->verbose(), 1;

ok $obj->verbose(-1);

#<tests for handle read and write abilities>
my($handle1,$file1) = $obj->tempfile;
my($handle2,$file2) = $obj->tempfile;

ok open(I,"$file1");
ok open(O,">$file2");

my $iotest;

$iotest = Bio::Root::IO->new(-file=>$file1);
ok $iotest->mode eq 'r';

$iotest = Bio::Root::IO->new(-fh=>\*I);
ok $iotest->mode eq 'r';

$iotest = Bio::Root::IO->new(-file=>">$file2");
ok $iotest->mode eq 'w';

$iotest = Bio::Root::IO->new(-fh=>\*O);
ok $iotest->mode eq 'w';
#</tests for handle read and write abilities>

1;


