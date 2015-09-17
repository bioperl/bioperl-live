# -*-Perl-*- Test Harness script for Bioperl
# $Id$


use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 56);

    use_ok('Bio::Root::Utilities');
}

# Object creation
my $u = Bio::Root::Utilities->new();
isa_ok($u, 'Bio::Root::Utilities') ;

# month2num() and num2month()

my @month = qw(XXX Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
for my $i (1 .. 12) {
  is $u->month2num($month[$i]), $i;
  is $u->num2month($i), $month[$i];
}

# untaint()

is $u->untaint(''), '';
is $u->untaint('nice string'), 'nice string';
is $u->untaint('bad *?&^$! string'), 'bad ';
is $u->untaint( q{100% relaxed&;''\"|*?!~<>^()[]{}$}, 1 ), '100% relaxed';

# mean_stdev()

my($mu,$sd);

($mu,$sd) = $u->mean_stdev();
is $mu, undef;
is $sd, undef;

($mu,$sd) = $u->mean_stdev(42);
is $mu, 42;
is $sd, undef;

($mu,$sd) = $u->mean_stdev(-1,0,1);
is $mu, 0;
is $sd, 1;

# file_date(), file_flavor(), date_format()

my $file = test_input_file('test.txt');
my $file2 = test_input_file('test 2.txt');
my $fdate = $u->file_date($file);
like $fdate ,  qr/\d{4}-\d{2}-\d{2}/, 'file_date()';
ok $u->file_flavor($file), 'unix (\n or 012 or ^J)';

my $date = $u->date_format();
like $date, qr/\d{4}-\d{2}-\d{2}/, 'date format';
my $date2 = $u->date_format('yyyy-mmm-dd', $date);
like $date2 , qr/\d{4}-[a-z]{3}-\d{2}/i, 'date format';
my $date3 = $u->date_format('mdhms');
like $date3 , qr/[a-z]{3}\d{1,2} \d{1,2}:\d{1,2}:\d{1,2}/, 'date format';
my $date4 = $u->date_format('d-m-y', '11/22/60');
like $date4 , qr/\d{1,2}-[a-z]{3}-\d{4}/i, 'date format';
my $date5 = $u->date_format('mdy', '1/5/01');
like $date5 , qr/[a-z]{3} \d{1,2}, \d{4}/i, 'date format';

# External executable-related functions.

my $exe = $u->find_exe('some-weird-thing-no-one-will-have');
ok ! defined $exe ;

# compress() and uncompress() using gzip.
SKIP: {
    my $gzip = $u->find_exe('gzip');
    skip "gzip not found, skipping gzip tests", 12 unless $gzip;
    ok -x $gzip;

    # test compression/decompression of a simple file
    my $zfile = $u->compress($file);

    # In Windows, the folder separator '\' may brake
    # the following qr{}, so change it to '/'
    $zfile =~ s'\\'/'g;
    $file  =~ s'\\'/'g;

    like $zfile, qr/$file.gz|tmp.bioperl.gz/;
    ok -s $zfile;
    if ($zfile =~ /tmp.bioperl.gz/) {
        ok -e $file;
    }
    else {
        ok ! -e $file;
    }
    my $unzfile = $u->uncompress($zfile);
    ok ! -e $zfile;
    ok -e $file;

    # test compression/decompression of a filename with spaces keeping the original intact
    my $zfile2 = $file2.'.gz';
    my $return = $u->compress(-file => $file2, -outfile => $zfile2, -tmp => 1);
    is $return, $zfile2;
    ok -e $zfile2;
    ok -e $file2;
    unlink $file2 or die "Problem deleting $file2: $!\n";
    $return = $u->uncompress(-file => $zfile2, -outfile => $file2, -tmp => 1);
    is $return, $file2;
    ok -e $file2;
    ok -e $zfile2;
    unlink $zfile2 or die "Problem deleting $zfile2: $!\n";
}

# send_mail()

# $u->send_mail(-to=>'sac@bioperl.org',  # <--- your address here!
#               -subj=>'Root-Utilities.t',
#               -msg=>'Hey, your send_mail() method works!');
