use warnings; use strict;
use Test::More;
use File::Basename qw(dirname basename); use File::Spec;

# Funky terminals like xterm on cygwin can mess up output comparison.
$ENV{'TERM'}='dumb';

package Helper;
use English qw( -no_match_vars ) ;
use Config;
use File::Basename qw(dirname basename); use File::Spec;
require Exporter;
our (@ISA, @EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(run_bio_program run_bio_program_nocheck test_file_name
             test_no_arg_opts test_one_arg_opts);

my $debug = $^W;

my $dirname = dirname(__FILE__);

sub test_file_name($)
{
    File::Spec->catfile($dirname, '..', 'test-files', shift)
}

# Runs bio program in a subshell. 0 is returned if everything went okay.
# nonzero if something went wrong. We check the results against
# $check_filename. Filenames are expanded to the proper locations
# in directories.
sub run_bio_program($$$$;$)
{
    my ($bio_program, $data_filename, $run_opts, $check_filename,
	$other_opts) = @_;
    $other_opts = {} unless defined $other_opts;
    $other_opts->{do_test} = 1 unless exists $other_opts->{do_test};

    my $full_data_filename = '';
    if ($data_filename ne '/dev/null') {
	$full_data_filename = test_file_name($data_filename);
    }

    my $full_check_filename = File::Spec->catfile($dirname, 'check-data',
						  "${bio_program}-${check_filename}");
    my $full_bio_progname = File::Spec->catfile($dirname, '..', 'bin', $bio_program);

    my $ext_file = sub {
        my ($ext) = @_;
        my $new_fn = $full_check_filename;
        $new_fn =~ s/\.right\z/.$ext/;
        return $new_fn;
    };

    my $err_filename = $ext_file->('err');

    my $cmd = "$EXECUTABLE_NAME $full_bio_progname $run_opts $full_data_filename" .
	" 2>$err_filename";
    print $cmd, "\n" if $debug;
    my $output = `$cmd`;
    print "$output\n" if $debug;
    my $rc = $CHILD_ERROR >> 8;
    my $test_rc = $other_opts->{exitcode} || 0;
    if ($other_opts->{do_test}) {
	Test::More::note("testing " . $other_opts->{note}) if $other_opts->{note};
	Test::More::note( "running $bio_program $run_opts $full_data_filename" );
	Test::More::is($rc, $test_rc,
		       "command ${bio_program} executed giving exit code $test_rc\n$cmd\n");
    }
    return $rc if $rc;

    open(RIGHT_FH, "<$full_check_filename") ||
	die "Cannot open $full_check_filename for reading - $OS_ERROR";
    undef $INPUT_RECORD_SEPARATOR;
    my $right_string = <RIGHT_FH>;
    ($output, $right_string) = $other_opts->{filter}->($output, $right_string)
	if $other_opts->{filter};
    my $got_filename;
    $got_filename = $ext_file->('got');
    # TODO : Perhaps make sure we optionally use eq_or_diff from
    # Test::Differences here.
    my $equal_output = $right_string eq $output;
    Test::More::ok($right_string eq $output, 'Output comparison')
	if $other_opts->{do_test};
    unlink $err_filename if -z $err_filename;
    if ($equal_output) {
        unlink $got_filename;
	return 0;
    } else {
        open (GOT_FH, '>', $got_filename)
            or die "Cannot open '$got_filename' for writing - $OS_ERROR";
        print GOT_FH $output;
        close GOT_FH;
        Test::More::diag("Compare $got_filename with $check_filename:");
	# FIXME use a better diff test.
	if ($OSNAME eq 'MSWin32') {
	    # Windows doesn't do diff.
	    Test::More::diag("Got:\n", $output, "Need:\n", $right_string);
	} else {
	    my $output = `diff -au $full_check_filename $got_filename 2>&1`;
	    my $rc = $? >> 8;
	    # GNU diff returns 0 if files are equal, 1 if different and 2
	    # if something went wrong. We also should take care of the
	    # case where diff isn't installed. So although we expect a 1
	    # for GNU diff, we'll also take accept 0, but any other return
	    # code means some sort of failure.
	    $output = `diff $check_filename $got_filename 2>&1`
		if ($rc > 1) || ($rc < 0) ;
	    Test::More::diag($output);
	    return 1;
	}
    }
}

# Runs a bioprogram but skips output checking
sub run_bio_program_nocheck($$$;$)
{
    my ($bio_program, $data_filename, $run_opts, $other_opts) = @_;
    $other_opts = {} unless defined $other_opts;
    $other_opts->{do_test} = 1 unless exists $other_opts->{do_test};

    my $full_data_filename = '';
    if ($data_filename ne '/dev/null') {
	$full_data_filename = test_file_name($data_filename);
    }

    my $full_bio_progname = File::Spec->catfile($dirname, '..', 'bin', $bio_program);

    my $err_filename = "$$.err";

    my $cmd = "$EXECUTABLE_NAME $full_bio_progname $run_opts $full_data_filename" .
	" 2>$err_filename";
    print $cmd, "\n" if $debug;
    my $output = `$cmd`;
    print "$output\n" if $debug;
    my $rc = $CHILD_ERROR >> 8;
    my $test_rc = $other_opts->{exitcode} || 0;
    if ($other_opts->{do_test}) {
	Test::More::note("testing " . $other_opts->{note}) if $other_opts->{note};
	Test::More::note( "running $bio_program $run_opts $full_data_filename" );
	Test::More::is($rc, $test_rc, "command ${bio_program} executed giving exit code $test_rc");
    }
    return $rc;
}

sub test_no_arg_opts($$$) {
    my ($bio_program, $data_filename, $notes) = @_;
    Test::More::note( "Testing ${bio_program} single-letter options on ${data_filename}" );
    for my $opt (keys %$notes) {
	run_bio_program($bio_program, $data_filename, "--${opt}", "opt-${opt}.right",
			{note=>$notes->{$opt}});
    }
}


sub test_one_arg_opts($$$) {
    my ($bio_program, $data_filename, $opts) = @_;

    for my $tup (@$opts) {
	my ($opt, $arg, $note) = @$tup;
	Test::More::note( "Testing ${bio_program} option-value options on ${data_filename}" );

	run_bio_program($bio_program, $data_filename, "--$opt $arg",
			"opt-$opt.right", {note=>$note});
    }
}



# Demo program
unless(caller) {
    run_bio_program('bioaln', 'test-bioaln.cds', '-a', 'opt-a.right');
    Test::More::done_testing();
}
1;
