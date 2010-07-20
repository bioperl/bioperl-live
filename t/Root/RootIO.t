# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;
use warnings;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 67);
	
    use_ok('Bio::Root::IO');
}

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
is scalar @stack, 2, 'stack_trace()';

my $verbobj = Bio::Root::IO->new(-verbose=>1,-strict=>1);
is $verbobj->verbose(), 1, 'set verbosity to 1';

ok $obj->verbose(-1);

#############################################
# tests for finding executables
#############################################
ok(my $io = Bio::Root::IO->new());
# An executable file
my $test_file = 'test_file.txt';
open FILE, '>', $test_file or die "Could not write file '$test_file': $!\n";
print FILE 'test';
close FILE;
chmod 0777, $test_file or die "Could not change permission of file '$test_file': $!\n";
ok ($obj->exists_exe($test_file), 'executable file');
# A not executable file
chmod 0444, $test_file or die "Could not change permission of file '$test_file': $!\n";
ok (! $obj->exists_exe($test_file), 'non-executable file');
unlink $test_file or die "Could not delete file '$test_file': $!\n";
# An executable dir
my $test_dir = 'test_dir';
mkdir $test_dir or die "Could not write dir '$test_dir': $!\n";
chmod 0777, $test_dir or die "Could not change permission of dir '$test_dir': $!\n";
ok (! $obj->exists_exe($test_dir), 'executable dir');
rmdir $test_dir or die "Could not delete dir '$test_dir': $!\n";


#############################################
# tests for handle read and write abilities
#############################################

ok my $TESTINFILE = Bio::Root::IO->catfile(qw(t data test.waba));

my($handle,$file) = $obj->tempfile;
ok $handle;
ok $file;

#test with files

ok my $rio = Bio::Root::IO->new(-file=>$TESTINFILE);
is $rio->file, $TESTINFILE;
is $rio->mode, 'r', 'filename, read';

ok my $wio = Bio::Root::IO->new(-file=>">$file");
is $wio->file, ">$file";
is $wio->mode, 'w', 'filename, write';

# test with handles

ok open(my $I, $TESTINFILE);
ok open(my $O, '>', $file);

ok $rio = Bio::Root::IO->new(-fh=>$I);
is $rio->_fh, $I;
is $rio->mode, 'r', 'handle, read';

ok $wio = Bio::Root::IO->new(-fh=>$O);
is $wio->_fh, $O;
is $wio->mode, 'w', 'handle, write';

SKIP: {
    my $warn;
    local $SIG{__WARN__} = sub { $warn = shift };
    my $tempfile = eval { require File::Temp; File::Temp->new }
	or skip 'could not create File::Temp object, maybe your File::Temp is 10 years old', 3;
    my $temp_io;
    ok $temp_io = Bio::Root::IO->new( -fh => $tempfile );
    is $temp_io->mode, 'w', 'is a write handle';
    # wish i could just use Test::Warn.  but then there's ... THE DEPENDENCY HOBGOBLIN! (TM)
    $temp_io->close;
    ok !$warn, 'no warnings';
}

##############################################
# tests _pushback for multi-line buffering
##############################################

my $rio1 = 

my $line1 = $rio->_readline;
my $line2 = $rio->_readline;

ok $rio->_pushback($line1);
ok $rio->_pushback($line2);

my $line3 = $rio->_readline;
my $line4 = $rio->_readline;
my $line5 = $rio->_readline;

is $line1, $line3;
is $line2, $line4;
isnt $line5, $line4;

ok close($I);
ok close($O);

##############################################
# test _print and _insert
##############################################

ok my $fio = Bio::Root::IO->new(-file=>">$file");
ok $fio->_print("line 1\n"), '_print';
ok $fio->_print("line 2\n");
ok $fio->_insert("insertion at line 2\n",2), '_insert at middle of file';
ok $fio->_print("line 3\n");
ok $fio->_print("line 4\n");
$fio->close;

open my $checkio, '<', $file;
my @content = <$checkio>;
close $checkio;
is_deeply \@content, ["line 1\n","insertion at line 2\n","line 2\n","line 3\n","line 4\n"];

ok $fio = Bio::Root::IO->new(-file=>">$file");
ok $fio->_insert("insertion at line 1\n",1), '_insert in empty file';
$fio->close;

open $checkio, '<', $file;
@content = <$checkio>;
close $checkio;
is_deeply \@content, ["insertion at line 1\n"];


##############################################
# test Win vs UNIX line ending
##############################################
{

my $unix_rio = Bio::Root::IO->new(-file => test_input_file('U71225.gb.unix'));
my $win_rio = Bio::Root::IO->new(-file => test_input_file('U71225.gb.win'));

for (1..5) {
    is($unix_rio->_readline, $win_rio->_readline);
}

}

##############################################
# test Win vs UNIX line ending using PerlIO::eol
##############################################

SKIP: {
    test_skip(-tests => 5, -requires_module => 'PerlIO::eol');
    local $Bio::Root::IO::HAS_EOL = 1;
    my $unix_rio = Bio::Root::IO->new(-file => test_input_file('U71225.gb.unix'));
    my $win_rio = Bio::Root::IO->new(-file => test_input_file('U71225.gb.win'));
    for (1..5) {
        is($unix_rio->_readline, $win_rio->_readline);
    }
}


##############################################
# tests http retrieval
##############################################

SKIP: {
  test_skip(-tests => 2, -requires_networking => 1);

  my $TESTURL = 'http://www.google.com/index.html';
  
  ok $rio = Bio::Root::IO->new(-url=>$TESTURL), 'default -url method';
  
  $Bio::Root::IO::HAS_LWP = 0;
  lives_ok {$rio = Bio::Root::IO->new(-url=>$TESTURL)};
}

##############################################
# test -string
##############################################

my $teststring = "Foo\nBar\nBaz";
ok $rio = Bio::Root::IO->new(-string =>$teststring), 'default -string method';

$line1 = $rio->_readline;
is($line1, "Foo\n");

$line2 = $rio->_readline;
is($line2, "Bar\n");
$rio->_pushback($line2);

$line3 = $rio->_readline;
is($line3, "Bar\n");
$line3 = $rio->_readline;
is($line3, "Baz");
