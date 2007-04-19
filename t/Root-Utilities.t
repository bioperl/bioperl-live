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
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;

    plan tests => 47;
}
use_ok('Bio::Root::Utilities');
use_ok('FindBin') ;
use FindBin qw($Bin);

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

my $data = "$Bin/data";
my $file = "$data/test.txt";
is -d $data, 1;
is -f $file, 1;

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
# Can't reliably test methods that depend on access to binaries:
# e.g. find_exe(), compress(), uncompress(), send_mail()

my $exe = $u->find_exe('some-weird-thing-no-one-will-have');
ok ! defined $exe ;

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
