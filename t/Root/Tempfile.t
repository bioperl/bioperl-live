# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 18);
	
	use_ok('Bio::Root::IO');
}

ok my $obj = Bio::Root::IO->new(-verbose => 0);

isa_ok($obj, 'Bio::Root::IO');

my $TEST_STRING = "Bioperl rocks!\n";

my ($tfh,$tfile);

eval {
    ($tfh,$tfile) = $obj->tempfile();
    print $tfh $TEST_STRING; 
    close($tfh);
    open(my $IN, $tfile) or die("cannot open $tfile");    
    my $val = join("", <$IN>) ;
    ok( $val eq $TEST_STRING );
    close $IN;
    ok( -e $tfile );
    undef $obj; 
};
undef $obj;
if( $@ ) {
    ok(0);
} else { 
   ok( ! -e $tfile, 'auto UNLINK => 1' );
}

$obj = Bio::Root::IO->new();

eval {
    my $tdir = $obj->tempdir(CLEANUP=>1);
    ok( -d $tdir );
    ($tfh, $tfile) = $obj->tempfile(dir => $tdir);
    close $tfh;
    ok( -e $tfile );
    undef $obj; # see Bio::Root::IO::_io_cleanup
};

if( $@ ) { ok(0); } 
else { ok( ! -e $tfile, 'tempfile deleted' ); }

eval {
    $obj = Bio::Root::IO->new(-verbose => 0);
    ($tfh, $tfile) = $obj->tempfile(UNLINK => 0);
    close $tfh;
    ok( -e $tfile );   
    undef $obj; # see Bio::Root::IO::_io_cleanup
};

if( $@ ) { ok(0) }
else { ok( -e $tfile, 'UNLINK => 0') }

ok unlink( $tfile) == 1 ;


ok $obj = Bio::Root::IO->new;

# check suffix is applied
my($fh1, $fn1) = $obj->tempfile(SUFFIX => '.bioperl');
ok $fh1;
like $fn1, qr/\.bioperl$/, 'tempfile suffix';
ok close $fh1;

# check single return value mode of File::Temp
my $fh2 = $obj->tempfile;
ok $fh2, 'tempfile() in scalar context';
ok close $fh2;
