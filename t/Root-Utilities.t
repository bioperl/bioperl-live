#-*-Perl-*-

# Test file for Bio::Root::Utilities
# Author: Torsten Seemann

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use lib './';
BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;

    plan tests => 46;
}
use Bio::Root::Utilities;
use FindBin qw($Bin);

# Test

ok(1);

# Object creation

my $u = Bio::Root::Utilities->new();
ok($u);

# month2num() and num2month()

my @month = qw(XXX Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
for my $i (1 .. 12) {
  ok $u->month2num($month[$i]), $i;
  ok $u->num2month($i), $month[$i];
}

# untaint()

ok $u->untaint(''), '';
ok $u->untaint('nice string'), 'nice string';
ok $u->untaint('bad *?&^$! string'), 'bad ';
ok $u->untaint( q{100% relaxed&;''\"|*?!~<>^()[]{}$}, 1 ), '100% relaxed';

# mean_stdev()

my($mu,$sd);

($mu,$sd) = $u->mean_stdev();
ok $mu, undef;
ok $sd, undef;

($mu,$sd) = $u->mean_stdev(42);
ok $mu, 42;
ok $sd, undef;

($mu,$sd) = $u->mean_stdev(-1,0,1);
ok $mu, 0;
ok $sd, 1;

# file_date(), file_flavor(), date_format()

my $data = "$Bin/data";
my $file = "$data/test.txt";
ok -d $data, 1;
ok -f $file, 1;

my $fdate = $u->file_date($file);
ok $fdate =~ /\d{4}-\d{2}-\d{2}/, 1;
ok $u->file_flavor($file), 'unix (\n or 012 or ^J)';

my $date = $u->date_format();
ok $date =~ /\d{4}-\d{2}-\d{2}/, 1;
my $date2 = $u->date_format('yyyy-mmm-dd', $date);
ok $date2 =~ /\d{4}-[a-z]{3}-\d{2}/i, 1;
my $date3 = $u->date_format('mdhms');
ok $date3 =~ /[a-z]{3}\d{1,2} \d{1,2}:\d{1,2}:\d{1,2}/, 1;
my $date4 = $u->date_format('d-m-y', '11/22/60');
ok $date4 =~ /\d{1,2}-[a-z]{3}-\d{4}/i, 1;
my $date5 = $u->date_format('mdy', '1/5/01');
ok $date5 =~ /[a-z]{3} \d{1,2}, \d{4}/i, 1;

#print "date: $date\ndate2: $date2\ndate3: $date3\ndate4: $date4\ndate5: $date5\n";

# External executable-related functions.
# Can't reliably test methods that depend on access to binaries:
# e.g. find_exe(), compress(), uncompress(), send_mail()

my $exe = $u->find_exe('some-weird-thing-no-one-will-have');
ok $exe, undef;

# compress() and uncompress() using gzip.
# Uncomment if you have gzip.

#my $gzip = $u->find_exe('gzip');
#ok $gzip, '/usr/bin/gzip';        # <--- modify for your gzip location

#my $zfile = $u->compress($file);
#ok $zfile, "$file.gz";
#my $unzfile = $u->uncompress($zfile);

#print "zfile: $zfile\n";
#print "unzfile: $unzfile\n";

# send_mail()
# Uncomment and edit to test sending mail

# $u->send_mail(-to=>'sac@bioperl.org',  # <--- your address here!
#               -subj=>'Root-Utilities.t',
#               -msg=>'Hey, your send_mail() method works!');
