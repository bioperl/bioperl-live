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
	use lib 't';
    }
    use Test;    
    plan tests => 6;
}

use Bio::Root::IO;

my $obj = new Bio::Root::IO();

ok defined($obj) && $obj->isa('Bio::Root::IO');

# doesn't work in perl 5.00405
my ($tfh,$tfile,$tdir,$val);
eval {
    ($tfh,$tfile) = $obj->tempfile();
    print $tfh ("test1"); 
    close($tfh);
    open(IN, $tfile) or die("cannot open $tfile");    
    $val = join("", <IN>) ;
    close IN;
    ok( -e $tfile );
    undef $obj;
};
if( $@ ) {
    ok(0);
} else { 
    ok( ! -e $tfile );
}

$obj = new Bio::Root::IO();

eval {
    ($tdir) = $obj->tempdir(CLEANUP=>1);
    ($tfh, $tfile) = $obj->tempfile(dir => $tdir);
    close $tfh;
    ok( -e $tfile );
    ok( -d $tdir );
    undef $obj;
};

if( $@ ) { ok(0); } 
else { ok( ! -e $tfile ); }

1;


