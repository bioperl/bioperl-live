# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;


BEGIN {
    my $NUMTESTS = 26;
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't', '.';
    }
    use Test;    
    plan tests => $NUMTESTS;

}

$| = 1;

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

#############################################
# <tests for handle read and write abilities>
#############################################
my($handle,$file) = $obj->tempfile;

ok open(I,"t/data/test.waba");
ok open(O,">$file");

my $rio;
my $wio;

#test with files
ok $rio = Bio::Root::IO->new(-file=>"t/data/test.waba");
ok $wio = Bio::Root::IO->new(-file=>">$file");

ok $rio->mode eq 'r';

ok $wio->mode eq 'w';

#test with handles
ok $rio = Bio::Root::IO->new(-fh=>\*I);
ok $wio = Bio::Root::IO->new(-fh=>\*O);
ok $rio->mode eq 'r';
ok $wio->mode eq 'w';

##############################################
# </tests for handle read and write abilities>
##############################################

##############################################
# <tests _pushback for multi-line buffering>
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

##############################################
# </tests _pushback for multi-line buffering>
##############################################

ok close(I);
ok close(O);

##############################################
# <tests http retrieval>
##############################################

eval {
  $rio = Bio::Root::IO->new(-file=>'http://www.google.com/index.html');
};

if($@){
  skip("couldn't get google.com, network down?",1);  
} else {
  ok(1);
}

##############################################
# </tests http retrieval>
##############################################


1;
