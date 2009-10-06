# -*-perl-*-
# $Id$

BEGIN {
    use lib '../..';
    use Bio::Root::Test;
    @AnyDBM_File::ISA = qw( Bio::DB::SQLite_File );
    test_begin( -tests => 79,
		-requires_module => 'DBD::SQLite' );
    use Fcntl qw(O_CREAT O_RDWR O_RDONLY);
    use_ok( AnyDBM_File );
}
use vars qw( $DB_HASH $DB_BTREE $DB_RECNO
             &R_DUP &R_CURSOR &R_FIRST &R_LAST
             &R_NEXT &R_PREV &R_IAFTER &R_IBEFORE
             &R_NOOVERWRITE &R_SETCURSOR );	       

use Bio::DB::AnyDBMImporter qw(:bdb);

my %db;
my $flags = O_CREAT | O_RDWR;
ok tie( %db, 'AnyDBM_File', 'my.db', $flags, 0666, $DB_BTREE, 0), "tie";
my $db = tied %db;

ok @db{qw( a b c d e f g h )} = qw( 1 2 3 4 5 6 7 8 ), "set";

# test: put, get, del, seq 
# all flags, boundaries
my ($key, $value, $ret);
$key = 'd';
$value = 40;

ok !$db->put($key, $value), "put replace";
ok !$db->get($key, $ret), "get";
is($ret, 40, "correct put/get");

ok !$db->seq($key, $value, R_FIRST), "seq R_FIRST";
is ($key, 'a', "key correct");
is ($value, 1, "value correct");

ok !$db->seq($key, $value, R_NEXT), "seq R_NEXT";
is ($key, 'b', "key correct");
is ($value, 2, "value correct");

ok !$db->del($key, R_CURSOR), "del R_CURSOR";
ok !$db->seq($key, $value, R_CURSOR), "seq R_CURSOR";
is ($key, 'c', "key correct on cursor update from del");
is ($value, 3, "value correct on cursor update from del");

ok delete $db{d}, "create a stray undef with the other API";
ok !$db->seq($key, $value, R_NEXT), "run into the undef with seq/R_NEXT";
is ($key, 'e', "key correct seq R_NEXT");
is ($value, 5, "value correct seq R_NEXT");
ok !$db->seq($key, $value, R_PREV), "seq R_PREV";
is ($key, 'c', "key correct seq R_PREV");
is ($value, 3, "value correct seq R_PREV");

ok !$db->seq($key, $value, R_LAST), "seq R_LAST";
is ($key, 'h', "key correct seq R_LAST");
is ($value, 8, "value correct seq R_LAST");
# note following, NO "!"
ok $db->seq($key, $value, R_NEXT), "check fail for R_NEXT";
!$db->seq($key, $value, R_FIRST);
ok $db->seq($key, $value, R_PREV), "check fail for R_PREV";
# put with flags
$key = 'Y';
$value = 101;
ok !$db->put($key, $value, R_SETCURSOR), "put R_SETCURSOR";
ok !$db->seq($key, $value, R_CURSOR), "what's at the cursor";
is ($value, 101, "cursor floated");
$key = 'Z';
$value = 102;
ok !$db->put($key, $value), "put no flags";
ok !$db->seq($key, $value, R_NEXT), "what's after the cursor?";
is($value, 102, "what we just put");
ok !$db->seq($key, $value, R_LAST), "what's at the bottom?";
is($value, 102, "what we just put");
# no '!'
ok $db->put($key, $value, R_NOOVERWRITE), "put R_NOOVERWRITE";


undef $db;
ok untie %db;

# array tests 

my @db;
ok tie( @db, 'AnyDBM_File', undef, $flags, 0666, $DB_RECNO), "tied array";
$db = tied @db;

@db = qw( a b c d e f g h );

$key = 3;
$value = 'D';
ok !$db->put($key, $value), "put replace";
ok !$db->get($key, $ret), "get";
is($ret, 'D', "correct put/get");

ok !$db->seq($key, $value, R_FIRST), "seq R_FIRST";
is ($key, 0, "key correct");
is ($value, 'a', "value correct");

ok !$db->seq($key, $value, R_NEXT), "seq R_NEXT";
is ($key, 1, "key correct");
is ($value, 'b', "value correct");

ok !$db->del($key, R_CURSOR), "del R_CURSOR";
ok !$db->seq($key, $value, R_CURSOR), "seq R_CURSOR";
is ($key, 1, "key correct on cursor update from del");
is ($value, 'c', "value correct on cursor update from del");

ok delete $db[2], "create a stray undef with the other API";
ok !$db->seq($key, $value, R_NEXT), "run into the undef with seq/R_NEXT";
is ($key, 2, "key correct seq R_NEXT");
is ($value, 'e', "value correct seq R_NEXT");
ok !$db->seq($key, $value, R_PREV), "seq R_PREV";
is ($key, 1, "key correct seq R_PREV");
is ($value, 'c', "value correct seq R_PREV");

ok !$db->seq($key, $value, R_LAST), "seq R_LAST";
is ($key, 5, "key correct seq R_LAST");
is ($value, 'h', "value correct seq R_LAST");
# note following, NO "!"
ok $db->seq($key, $value, R_NEXT), "check fail for R_NEXT";
!$db->seq($key, $value, R_FIRST);
ok $db->seq($key, $value, R_PREV), "check fail for R_PREV";
# put with flags
$key = 100; # key will be ignored for arrays
$value = 'Z';
ok !$db->put($key, $value, R_IBEFORE), "put R_IBEFORE";
is($key, 0, "key set to recno");
ok !$db->seq($key, $value, R_CURSOR), "what's at the cursor";
is($value, 'a', "cursor didn't float");
ok !$db->seq($key, $value, R_FIRST), "what's at the top";
is($value, 'Z', "correct put R_IBEFORE");

$key = 101;
$value = 'X';
ok !$db->put($key, $value, R_SETCURSOR), "put R_SETCURSOR";
ok !$db->seq($key, $value, R_CURSOR), "what's at the cursor";
is ($value, 'X', "cursor floated");
$key = 102;
$value = 'Y';
ok !$db->put($key, $value, R_IAFTER), "put R_IAFTER";
ok !$db->seq($key, $value, R_NEXT), "what's after the cursor?";
is($value, 'Y', "what we just put");
ok !$db->seq($key, $value, R_LAST), "what's at the bottom?";
is($value, 'Y', "what we just put");

undef $db;
ok untie(@db);

1;
