# -*-perl-*-
# $Id$

BEGIN {
    use lib '../..';
    use Bio::Root::Test;
    @AnyDBM_File::ISA = qw( Bio::DB::SQLite_File );
    test_begin( -tests => 20,
		-requires_module => 'DBD::SQLite' );
    use Fcntl qw(O_CREAT O_RDWR O_RDONLY);
    use_ok( AnyDBM_File );
}

my %db;
my $flags = O_CREAT | O_RDWR;
ok tie( %db, 'AnyDBM_File', 'my.db', $flags, 0666, 0), "tie";
ok $db{'butcher'} = 1, "set";
ok $db{'baker'} = 2, "set";
ok $db{'candlestick maker'} = 3, "set";
ok $db{'ooh, isnt this a very very very very very very very very very very long key, my goodness gracious me'} = 4, "set";
ok $a = $db{'baker'}, "access";
is($a, 2), "value correct";
ok exists $db{'baker'}, "exists";
ok delete $db{'baker'}, "delete";
ok !exists $db{'baker'}, "delete deletes";
ok my @k = keys %db, "iterate (keys)";
is_deeply( [ sort @k ] , ['butcher','candlestick maker', 'ooh, isnt this a very very very very very very very very very very long key, my goodness gracious me'] , "all keys retrieved");
my $f = 1;
while (my ($id, $val) = each %db) {
    1;
    $f *= $val;
}
is($f, 12, "iterate (each)");

my $file = (tied %db)->file;
ok ( -e $file, "now you see it" );
untie %db;
ok ( ! -e $file, "now you don't");

tie( %db, 'AnyDBM_File', 'my.db', $flags, 0666, 1);
ok ( -e 'my.db', "now you see it" );
untie %db;
ok ( -e 'my.db', "now you still see it");
ok ( !(tied %db), "but tied obj is gone" );

ok ( unlink('my.db'), "now you don't");

1;
