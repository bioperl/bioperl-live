# -*-Perl-*- Test Harness script for Bioperl

use strict;
use warnings;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    test_begin(-tests => 154);
    use_ok 'Bio::Root::IO';
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
# -X test file will fail in Windows regardless of chmod,
# because it looks for the executable suffix (like ".exe")
if ($^O =~ m/mswin/i) {
    # An executable file
    my $exec_file = 'test_exec.exe';
    open my $exe_fh, '>', $exec_file or die "Could not write file '$exec_file': $!\n";
    close $exe_fh;
    ok $obj->exists_exe($exec_file), 'executable file';
    unlink $exec_file or die "Could not delete file '$exec_file': $!\n";

    # A not executable file
    ok (! $obj->exists_exe($out_file), 'non-executable file');
    unlink $out_file  or die "Could not delete file '$out_file': $!\n";
}
else {
    # An executable file
    chmod 0777, $out_file or die "Could not change permission of file '$out_file': $!\n";
    ok $obj->exists_exe($out_file), 'executable file';

    # A not executable file
    chmod 0444, $out_file or die "Could not change permission of file '$out_file': $!\n";
    ok (! $obj->exists_exe($out_file), 'non-executable file');
    unlink $out_file or die "Could not delete file '$out_file': $!\n";
}

# An executable dir
my $out_dir = 'test_dir';
mkdir $out_dir or die "Could not write dir '$out_dir': $!\n";
chmod 0777, $out_dir or die "Could not change permission of dir '$out_dir': $!\n";
ok (! $obj->exists_exe($out_dir), 'executable dir');
rmdir $out_dir or die "Could not delete dir '$out_dir': $!\n";


#############################################
# tests for handle read and write abilities
#############################################

# Test catfile

ok my $in_file = Bio::Root::IO->catfile(qw(t data test.waba));
is $in_file, test_input_file('test.waba');

ok my $in_file_2 = Bio::Root::IO->catfile(qw(t data test.txt));

$out_file = test_output_file();


# Test with files

ok my $rio = Bio::Root::IO->new( -input => $in_file ), 'Read from file';
is $rio->file, $in_file;
is_deeply [$rio->cleanfile], [undef, $in_file];
is $rio->mode, 'r';
ok $rio->close;

ok $rio = Bio::Root::IO->new( -file => '<'.$in_file );
is $rio->file, '<'.$in_file;
is_deeply [$rio->cleanfile], ['<', $in_file];
1 while $rio->_readline; # read entire file content
is $rio->mode, 'r';
ok $rio->close;

ok my $wio = Bio::Root::IO->new( -file => ">$out_file" ), 'Write to file';
is $wio->file, ">$out_file";
is_deeply [$wio->cleanfile], ['>', $out_file];
is $wio->mode, 'w';
ok $wio->close;

ok $rio = Bio::Root::IO->new( -file => "+>$out_file" ), 'Read+write to file';
is $rio->file, "+>$out_file";
is_deeply [$rio->cleanfile], ['+>', $out_file];
is $rio->mode, 'rw';
ok $rio->close;


# Test with handles

my $in_fh;
open $in_fh , '<', $in_file  or die "Could not read file '$in_file': $!\n", 'Read from GLOB handle';
ok $rio = Bio::Root::IO->new( -fh => $in_fh );
is $rio->_fh, $in_fh;
is $rio->mode, 'r';
close $in_fh;

open $out_fh, '>', $out_file or die "Could not write file '$out_file': $!\n", 'Write to GLOB handle';
ok $wio = Bio::Root::IO->new( -fh => $out_fh );
is $wio->_fh, $out_fh;
is $wio->mode, 'w';
close $out_fh;

SKIP: {
    eval { require File::Temp; }
       or skip 'could not create File::Temp object, maybe your File::Temp is 10 years old', 4;

    $out_fh = File::Temp->new;
    ok $wio = Bio::Root::IO->new( -fh => $out_fh ), 'Read from File::Temp handle';
    isa_ok $wio, 'Bio::Root::IO';
    is $wio->mode, 'rw', 'is a write handle';
    warnings_like sub { $wio->close }, '', 'no warnings in ->close()';
    ok $wio->close;
}


# Exclusive arguments
open $in_fh , '<', $in_file  or die "Could not read file '$in_file': $!\n", 'Read from GLOB handle';
throws_ok {$rio = Bio::Root::IO->new( -input => $in_file, -fh     => $in_fh     )} qr/Providing both a file and a filehandle for reading/, 'Exclusive arguments';
throws_ok {$rio = Bio::Root::IO->new( -input => $in_file, -file   => $in_file_2 )} qr/Input file given twice/;
throws_ok {$rio = Bio::Root::IO->new( -input => $in_file, -string => 'abcedf'   )} qr/File or filehandle provided with -string/;
throws_ok {$rio = Bio::Root::IO->new( -fh    => $in_fh  , -file   => $in_file   )} qr/Providing both a file and a filehandle for reading/;
throws_ok {$rio = Bio::Root::IO->new( -fh    => $in_fh  , -string => 'abcedf'   )} qr/File or filehandle provided with -string/;
throws_ok {$rio = Bio::Root::IO->new( -file  => $in_file, -string => 'abcedf'   )} qr/File or filehandle provided with -string/;
close $in_fh;

lives_ok  {$rio = Bio::Root::IO->new( -input => $in_file, -file   => $in_file   )} 'Same file';


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
ok $fio->close;

open my $checkio, '<', $out_file or die "Could not read file '$out_file': $!\n";
my @content = <$checkio>;
close $checkio;
is_deeply \@content, ["line 1\n","insertion at line 2\n","line 2\n","line 3\n","line 4\n"];

ok $fio = Bio::Root::IO->new(-file=>">$out_file");
ok $fio->_insert("insertion at line 1\n",1), '_insert in empty file';
ok $fio->close;

open $checkio, '<', $out_file or die "Could not read file '$out_file': $!\n";
@content = <$checkio>;
close $checkio;
is_deeply \@content, ["insertion at line 1\n"];


##############################################
# test Win vs UNIX line ending
##############################################

{
    ok my $unix_rio = Bio::Root::IO->new(-file => test_input_file('U71225.gb.unix'));
    ok my $win_rio  = Bio::Root::IO->new(-file => test_input_file('U71225.gb.win' ));
    ok my $mac_rio  = Bio::Root::IO->new(-file => test_input_file('U71225.gb.mac' ));

    my $expected = "LOCUS       U71225                  1164 bp    DNA     linear   VRT 27-NOV-2001\n";
    is $unix_rio->_readline, $expected;
    is $win_rio->_readline , $expected;
    like $mac_rio->_readline, qr#^LOCUS.*//\n$#ms;
    # line spans entire file because lines end with "\r" but $/ is "\n"

    $expected = "DEFINITION  Desmognathus quadramaculatus 12S ribosomal RNA gene, partial\n";
    is $unix_rio->_readline, $expected;
    is $win_rio->_readline , $expected;
    is $mac_rio->_readline , undef;

    $expected = "            sequence; tRNA-Val gene, complete sequence; and 16S ribosomal RNA\n";
    is $unix_rio->_readline, $expected;
    is $win_rio->_readline , $expected;
    is $mac_rio->_readline , undef;

    $expected = "            gene, partial sequence, mitochondrial genes for mitochondrial RNAs.\n";
    is $unix_rio->_readline, $expected;
    is $win_rio->_readline , $expected;
    is $mac_rio->_readline , undef;

    $expected = "ACCESSION   U71225\n";
    is $unix_rio->_readline, $expected;
    is $win_rio->_readline , $expected;
    is $mac_rio->_readline , undef;

    # In Windows the "-raw" parameter has no effect, because Perl already discards
    # the '\r' from the line when reading in text mode from the filehandle
    # ($line = <$fh>), and put it back automatically when printing
    if ($^O =~ m/mswin/i) {
        is $win_rio->_readline( -raw => 1) , "VERSION     U71225.1  GI:2804359\n";
    }
    else {
        is $win_rio->_readline( -raw => 1) , "VERSION     U71225.1  GI:2804359\r\n";
    }
    is $win_rio->_readline( -raw => 0) , "KEYWORDS    .\n";
}


##############################################
# test Win vs UNIX line ending using PerlIO::eol
##############################################

SKIP: {
    test_skip(-tests => 20, -requires_module => 'PerlIO::eol');

    local $Bio::Root::IO::HAS_EOL = 1;
    ok my $unix_rio = Bio::Root::IO->new(-file => test_input_file('U71225.gb.unix'));
    ok my $win_rio  = Bio::Root::IO->new(-file => test_input_file('U71225.gb.win' ));
    ok my $mac_rio  = Bio::Root::IO->new(-file => test_input_file('U71225.gb.mac' ));

    my $expected = "LOCUS       U71225                  1164 bp    DNA     linear   VRT 27-NOV-2001\n";
    is $unix_rio->_readline, $expected;
    is $win_rio->_readline , $expected;
    is $mac_rio->_readline , $expected;

    $expected = "DEFINITION  Desmognathus quadramaculatus 12S ribosomal RNA gene, partial\n";
    is $unix_rio->_readline, $expected;
    is $win_rio->_readline , $expected;
    is $mac_rio->_readline , $expected;

    $expected = "            sequence; tRNA-Val gene, complete sequence; and 16S ribosomal RNA\n";
    is $unix_rio->_readline, $expected;
    is $win_rio->_readline , $expected;
    is $mac_rio->_readline , $expected;

    $expected = "            gene, partial sequence, mitochondrial genes for mitochondrial RNAs.\n";
    is $unix_rio->_readline, $expected;
    is $win_rio->_readline , $expected;
    is $mac_rio->_readline , $expected;

    $expected = "ACCESSION   U71225\n";
    is $unix_rio->_readline, $expected;
    is $win_rio->_readline , $expected;
    is $mac_rio->_readline , $expected;

    # $HAS_EOL ignores -raw
    is $win_rio->_readline( -raw => 1) , "VERSION     U71225.1  GI:2804359\n";
    is $win_rio->_readline( -raw => 0) , "KEYWORDS    .\n";
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
# test -string
##############################################

my $teststring = "Foo\nBar\nBaz";
ok $rio = Bio::Root::IO->new(-string => $teststring), 'Read string';

is $rio->mode, 'r';

ok $line1 = $rio->_readline;
is $line1, "Foo\n";

ok $line2 = $rio->_readline;
is $line2, "Bar\n";
ok $rio->_pushback($line2);

ok $line3 = $rio->_readline;
is $line3, "Bar\n";
ok $line3 = $rio->_readline;
is $line3, 'Baz';


##############################################
# test tempfile()
##############################################
{
ok my $obj = Bio::Root::IO->new(-verbose => 0);

isa_ok $obj, 'Bio::Root::IO';

my $TEST_STRING = "Bioperl rocks!\n";

my ($tfh,$tfile);

eval {
    ($tfh, $tfile) = $obj->tempfile();
    isa_ok $tfh, 'GLOB';
    print $tfh $TEST_STRING;
    close $tfh;
    open my $IN, '<', $tfile or die "Could not read file '$tfile': $!\n";
    my $val = join '', <$IN>;
    is $val, $TEST_STRING;
    close $IN;
    ok -e $tfile;
    undef $obj;
};
undef $obj;
if ( $@ ) {
    ok 0;
} else {
    ok ! -e $tfile, 'auto UNLINK => 1';
}

$obj = Bio::Root::IO->new();

eval {
    my $tdir = $obj->tempdir(CLEANUP=>1);
    ok -d $tdir;
    ($tfh, $tfile) = $obj->tempfile(dir => $tdir);
    close $tfh;
    ok -e $tfile;
    undef $obj; # see Bio::Root::IO::_io_cleanup
};

if ( $@ ) {
    ok 0;
} else {
    ok ! -e $tfile, 'tempfile deleted';
}

eval {
    $obj = Bio::Root::IO->new(-verbose => 0);
    ($tfh, $tfile) = $obj->tempfile(UNLINK => 0);
    isa_ok $tfh, 'GLOB';
    close $tfh;
    ok -e $tfile;
    undef $obj; # see Bio::Root::IO::_io_cleanup
};

if ( $@ ) {
   ok 0;
} else {
   ok -e $tfile, 'UNLINK => 0';
}

ok unlink( $tfile) == 1 ;


ok $obj = Bio::Root::IO->new;

# check suffix is applied
my ($fh1, $fn1) = $obj->tempfile(SUFFIX => '.bioperl');
isa_ok $fh1, 'GLOB';
like $fn1, qr/\.bioperl$/, 'tempfile suffix';
ok close $fh1;

# check single return value mode of File::Temp
my $fh2 = $obj->tempfile;
isa_ok $fh2, 'GLOB';
ok $fh2, 'tempfile() in scalar context';
ok close $fh2;
}
