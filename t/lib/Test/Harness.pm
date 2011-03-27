# -*- Mode: cperl; cperl-indent-level: 4 -*-

package Test::Harness;

require 5.00405;
use Test::Harness::Straps;
use Test::Harness::Assert;
use Exporter;
use Benchmark;
use Config;
use strict;


use vars qw(
    $VERSION 
    @ISA @EXPORT @EXPORT_OK 
    $Verbose $Switches $Debug
    $verbose $switches $debug
    $Columns
    $Timer
    $ML $Last_ML_Print
    $Strap
    $has_time_hires
);

BEGIN {
    eval q{use Time::HiRes 'time'};
    $has_time_hires = !$@;
}

=head1 NAME

Test::Harness - Run Perl standard test scripts with statistics

=head1 VERSION

Version 2.64

=cut

$VERSION = '2.64';

# Backwards compatibility for exportable variable names.
*verbose  = *Verbose;
*switches = *Switches;
*debug    = *Debug;

$ENV{HARNESS_ACTIVE} = 1;
$ENV{HARNESS_VERSION} = $VERSION;

END {
    # For VMS.
    delete $ENV{HARNESS_ACTIVE};
    delete $ENV{HARNESS_VERSION};
}

my $Files_In_Dir = $ENV{HARNESS_FILELEAK_IN_DIR};

# Stolen from Params::Util
sub _CLASS {
    (defined $_[0] and ! ref $_[0] and $_[0] =~ m/^[^\W\d]\w*(?:::\w+)*$/s) ? $_[0] : undef;
}

# Strap Overloading
if ( $ENV{HARNESS_STRAPS_CLASS} ) {
    die 'Set HARNESS_STRAP_CLASS, singular, not HARNESS_STRAPS_CLASS';
}
my $HARNESS_STRAP_CLASS  = $ENV{HARNESS_STRAP_CLASS} || 'Test::Harness::Straps';
if ( $HARNESS_STRAP_CLASS =~ /\.pm$/ ) {
    # "Class" is actually a filename, that should return the
    # class name as its true return value.
    $HARNESS_STRAP_CLASS = require $HARNESS_STRAP_CLASS;
    if ( !_CLASS($HARNESS_STRAP_CLASS) ) {
        die "HARNESS_STRAP_CLASS '$HARNESS_STRAP_CLASS' is not a valid class name";
    }
}
else {
    # It is a class name within the current @INC
    if ( !_CLASS($HARNESS_STRAP_CLASS) ) {
        die "HARNESS_STRAP_CLASS '$HARNESS_STRAP_CLASS' is not a valid class name";
    }
    eval "require $HARNESS_STRAP_CLASS";
    die $@ if $@;
}
if ( !$HARNESS_STRAP_CLASS->isa('Test::Harness::Straps') ) {
    die "HARNESS_STRAP_CLASS '$HARNESS_STRAP_CLASS' must be a Test::Harness::Straps subclass";
}

$Strap = $HARNESS_STRAP_CLASS->new;

sub strap { return $Strap };

@ISA = ('Exporter');
@EXPORT    = qw(&runtests);
@EXPORT_OK = qw(&execute_tests $verbose $switches);

$Verbose  = $ENV{HARNESS_VERBOSE} || 0;
$Debug    = $ENV{HARNESS_DEBUG} || 0;
$Switches = '-w';
$Columns  = $ENV{HARNESS_COLUMNS} || $ENV{COLUMNS} || 80;
$Columns--;             # Some shells have trouble with a full line of text.
$Timer    = $ENV{HARNESS_TIMER} || 0;

=head1 SYNOPSIS

  use Test::Harness;

  runtests(@test_files);

=head1 DESCRIPTION

B<STOP!> If all you want to do is write a test script, consider
using Test::Simple.  Test::Harness is the module that reads the
output from Test::Simple, Test::More and other modules based on
Test::Builder.  You don't need to know about Test::Harness to use
those modules.

Test::Harness runs tests and expects output from the test in a
certain format.  That format is called TAP, the Test Anything
Protocol.  It is defined in L<Test::Harness::TAP>.

C<Test::Harness::runtests(@tests)> runs all the testscripts named
as arguments and checks standard output for the expected strings
in TAP format.

The F<prove> utility is a thin wrapper around Test::Harness.

=head2 Taint mode

Test::Harness will honor the C<-T> or C<-t> in the #! line on your
test files.  So if you begin a test with:

    #!perl -T

the test will be run with taint mode on.

=head2 Configuration variables.

These variables can be used to configure the behavior of
Test::Harness.  They are exported on request.

=over 4

=item C<$Test::Harness::Verbose>

The package variable C<$Test::Harness::Verbose> is exportable and can be
used to let C<runtests()> display the standard output of the script
without altering the behavior otherwise.  The F<prove> utility's C<-v>
flag will set this.

=item C<$Test::Harness::switches>

The package variable C<$Test::Harness::switches> is exportable and can be
used to set perl command line options used for running the test
script(s). The default value is C<-w>. It overrides C<HARNESS_PERL_SWITCHES>.

=item C<$Test::Harness::Timer>

If set to true, and C<Time::HiRes> is available, print elapsed seconds
after each test file.

=back


=head2 Failure

When tests fail, analyze the summary report:

  t/base..............ok
  t/nonumbers.........ok
  t/ok................ok
  t/test-harness......ok
  t/waterloo..........dubious
          Test returned status 3 (wstat 768, 0x300)
  DIED. FAILED tests 1, 3, 5, 7, 9, 11, 13, 15, 17, 19
          Failed 10/20 tests, 50.00% okay
  Failed Test  Stat Wstat Total Fail  List of Failed
  ---------------------------------------------------------------
  t/waterloo.t    3   768    20   10  1 3 5 7 9 11 13 15 17 19
  Failed 1/5 test scripts, 80.00% okay. 10/44 subtests failed, 77.27% okay.

Everything passed but F<t/waterloo.t>.  It failed 10 of 20 tests and
exited with non-zero status indicating something dubious happened.

The columns in the summary report mean:

=over 4

=item B<Failed Test>

The test file which failed.

=item B<Stat>

If the test exited with non-zero, this is its exit status.

=item B<Wstat>

The wait status of the test.

=item B<Total>

Total number of tests expected to run.

=item B<Fail>

Number which failed, either from "not ok" or because they never ran.

=item B<List of Failed>

A list of the tests which failed.  Successive failures may be
abbreviated (ie. 15-20 to indicate that tests 15, 16, 17, 18, 19 and
20 failed).

=back


=head1 FUNCTIONS

The following functions are available.

=head2 runtests( @test_files )

This runs all the given I<@test_files> and divines whether they passed
or failed based on their output to STDOUT (details above).  It prints
out each individual test which failed along with a summary report and
a how long it all took.

It returns true if everything was ok.  Otherwise it will C<die()> with
one of the messages in the DIAGNOSTICS section.

=cut

sub runtests {
    my(@tests) = @_;

    local ($\, $,);

    my ($tot, $failedtests,$todo_passed) = execute_tests(tests => \@tests);
    print get_results($tot, $failedtests,$todo_passed);

    my $ok = _all_ok($tot);

    assert(($ok xor keys %$failedtests), 
           q{ok status jives with $failedtests});

    if (! $ok) {
        die("Failed $tot->{bad}/$tot->{tests} test programs. " .
            "@{[$tot->{max} - $tot->{ok}]}/$tot->{max} subtests failed.\n");
    }

    return $ok;
}

# my $ok = _all_ok(\%tot);
# Tells you if this test run is overall successful or not.

sub _all_ok {
    my($tot) = shift;

    return $tot->{bad} == 0 && ($tot->{max} || $tot->{skipped}) ? 1 : 0;
}

# Returns all the files in a directory.  This is shorthand for backwards
# compatibility on systems where C<glob()> doesn't work right.

sub _globdir {
    local *DIRH;

    opendir DIRH, shift;
    my @f = readdir DIRH;
    closedir DIRH;

    return @f;
}

=head2 execute_tests( tests => \@test_files, out => \*FH )

Runs all the given C<@test_files> (just like C<runtests()>) but
doesn't generate the final report.  During testing, progress
information will be written to the currently selected output
filehandle (usually C<STDOUT>), or to the filehandle given by the
C<out> parameter.  The I<out> is optional.

Returns a list of two values, C<$total> and C<$failed>, describing the
results.  C<$total> is a hash ref summary of all the tests run.  Its
keys and values are this:

    bonus           Number of individual todo tests unexpectedly passed
    max             Number of individual tests ran
    ok              Number of individual tests passed
    sub_skipped     Number of individual tests skipped
    todo            Number of individual todo tests

    files           Number of test files ran
    good            Number of test files passed
    bad             Number of test files failed
    tests           Number of test files originally given
    skipped         Number of test files skipped

If C<< $total->{bad} == 0 >> and C<< $total->{max} > 0 >>, you've
got a successful test.

C<$failed> is a hash ref of all the test scripts that failed.  Each key
is the name of a test script, each value is another hash representing
how that script failed.  Its 