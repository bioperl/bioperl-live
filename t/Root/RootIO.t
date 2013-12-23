# -*-Perl-*- Test Harness script for Bioperl

use strict;
use warnings;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    test_begin(-tests => 84);
    use_ok('Bio::Root::IO');
    use_ok('Bio::SeqIO');
    use_ok('Bio::Assembly::IO');
}


ok my $obj = Bio::Root::IO->new();
isa_ok $obj, 'Bio::Root::IO';



#############################################
# tests for exceptions/debugging/verbosity
#############################################

throws_ok { $obj->throw('Testing throw') } qr/Testing throw/, 'Throw';

$obj->verbose(-1);
throws_ok { $obj->throw('Testing throw') } qr/Testing throw/;

eval { $obj->warn('Testing warn') };
ok !$@, 'Warn';

$obj->verbose(1);
throws_ok { $obj->throw('Testing throw') } qr/Testing throw/;

ok my @stack = $obj->stack_trace(), 'Stack trace';
is scalar @stack, 2;

ok my $verbobj = Bio::Root::IO->new( -verbose => 1, -strict => 1 ), 'Verbosity';
is $verbobj->verbose(), 1;

ok $obj->verbose(-1);


#############################################
# tests for finding executables
#############################################

ok my $io = Bio::Root::IO->new();

# An executable file
my $out_file = 'test_file.txt';
my $out_fh;
open  $out_fh, '>', $out_file or die "Could not write file '$out_file': $!\n";
print $out_fh 'test';
close $out_fh;
chmod 0777, $out_file or die "Could not change permission of file '$out_file': $!\n";
ok $obj->exists_exe($out_file), 'executable file';

# A not executable file
chmod 0444, $out_file or die "Could not change permission of file '$out_file': $!\n";
ok (! $obj->exists_exe($out_file), 'non-executable file');
unlink $out_file or die "Could not delete file '$out_file': $!\n";

# An executable dir
my $out_dir = 'test_dir';
mkdir $out_dir or die "Could not write dir '$out_dir': $!\n";
chmod 0777, $out_dir or die "Could not change permission of dir '$out_dir': $!\n";
ok (! $obj->exists_exe($out_dir), 'executable dir');
rmdir $out_dir or die "Could not delete dir '$out_dir': $!\n";


#############################################
# tests for handle read and write abilities
#############################################

ok my $in_file = Bio::Root::IO->catfile(qw(t data test.waba));
my $in_fh;

($out_fh, $out_file) = $obj->tempfile;
ok $out_fh;
ok $out_file;

# test with files

ok my $rio = Bio::Root::IO->new( -file => $in_file ), 'Read from file';
is $rio->file, $in_file;
is $rio->mode, 'r';

ok my $wio = Bio::Root::IO->new( -file => ">$out_file" ), 'Write to file';
is $wio->file, ">$out_file";
is $wio->mode, 'w';

# test with handles

ok open($in_fh, $in_file);
ok open($out_fh, '>', $out_file);

ok $rio = Bio::Root::IO->new( -fh => $in_fh ), 'Read from handle';
is $rio->_fh, $in_fh;
is $rio->mode, 'r';

ok $wio = Bio::Root::IO->new( -fh => $out_fh ), 'Write to handle';
is $wio->_fh, $out_fh;
is $wio->mode, 'w';

ok close $in_fh;
ok close $out_fh;

SKIP: {
    my $tempfile = eval { require File::Temp; File::Temp->new }
       or skip 'could not create File::Temp object, maybe your File::Temp is 10 years old', 3;

    my $temp_io = Bio::Root::IO->new( -fh => $tempfile );
    isa_ok $temp_io, 'Bio::Root::IO';
    is $temp_io->mode, 'w', 'is a write handle';
    warnings_like sub { $temp_io->close }, '', 'no warnings in ->close call';
}


##############################################
# tests _pushback for multi-line buffering
##############################################

ok $rio = Bio::Root::IO->new( -file => $in_file ), 'Pushback';

ok my $line1 = $rio->_readline;
ok my $line2 = $rio->_readline;

ok $rio->_pushback($line2);
ok $rio->_pushback($line1);

ok my $line3 = $rio->_readline;
ok my $line4 = $rio->_readline;
ok my $line5 = $rio->_readline;

is $line1, $line3;
is $line2, $line4;
isnt $line5, $line4;

ok $rio->close;



##############################################
# test _print and _insert
##############################################

ok my $fio = Bio::Root::IO->new( -file => ">$out_file" );
ok $fio->_print("line 1\n"), '_print';
ok $fio->_print("line 2\n");
ok $fio->_insert("insertion at line 2\n",2), '_insert at middle of file';
ok $fio->_print("line 3\n");
ok $fio->_print("line 4\n");
$fio->close;

open my $checkio, '<', $out_file;
my @content = <$checkio>;
close $checkio;
is_deeply \@content, ["line 1\n","insertion at line 2\n","line 2\n","line 3\n","line 4\n"];

ok $fio = Bio::Root::IO->new(-file=>">$out_file");
ok $fio->_insert("insertion at line 1\n",1), '_insert in empty file';
$fio->close;

open $checkio, '<', $out_file;
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
    test_skip(-tests => 7, -requires_module => 'PerlIO::eol');
    local $Bio::Root::IO::HAS_EOL = 1;
    ok my $unix_rio = Bio::Root::IO->new(-file => test_input_file('U71225.gb.unix'));
    ok my $win_rio  = Bio::Root::IO->new(-file => test_input_file('U71225.gb.win'));
    for (1..5) {
        is $unix_rio->_readline, $win_rio->_readline;
    }
}


##############################################
# test Path::Class support
##############################################

SKIP: {
    test_skip(-tests => 2, -requires_module => 'Path::Class');
    my $f = sub { Bio::Root::IO->new( -file => Path::Class::file(test_input_file('U71225.gb.unix') ) ) };
    lives_ok(sub { $f->() } , 'Bio::Root::IO->new can handle a Path::Class object');
    isa_ok($f->(), 'Bio::Root::IO');
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
ok $rio = Bio::Root::IO->new(-string => $teststring), 'default -string method';

ok $line1 = $rio->_readline;
is $line1, "Foo\n";

ok $line2 = $rio->_readline;
is $line2, "Bar\n";
ok $rio->_pushback($line2);

ok $line3 = $rio->_readline;
is $line3, "Bar\n";
ok $line3 = $rio->_readline;
is $line3, 'Baz';
