package Bio::Root::Test;
use strict;
use warnings;
# According to Ovid, 'use base' can override signal handling, so use
# old-fashioned way. This should be a Test::Builder::Module subclass
# for consistency (as are any Test modules)
use Test::Most;
use Test::Builder;
use Test::Builder::Module;
use File::Temp qw(tempdir);
use File::Spec;

our @ISA = qw(Test::Builder::Module);

# ABSTRACT: a common base for all Bioperl test scripts
# AUTHOR:   Sendu Bala <bix@sendu.me.uk>
# OWNER:    Sendu Bala
# LICENSE:  Perl_5

# CONTRIBUTOR: Chris Fields <cjfields@bioperl.org>

=head1 SYNOPSIS

  use lib '.'; # (for core package tests only)
  use Bio::Root::Test;

  test_begin(-tests => 20,
             -requires_modules => [qw(IO::String XML::Parser)],
             -requires_networking => 1);

  my $do_network_tests = test_network();
  my $output_debugging = test_debug();

  # Bio::Root::Test rewraps Test::Most, so one can carry out tests with
  # Test::More, Test::Exception, Test::Warn, Test::Deep, Test::Diff syntax

  SKIP: {
    # these tests need version 2.6 of Optional::Module to work
    test_skip(-tests => 10, -requires_module => 'Optional::Module 2.6');
    use_ok('Optional::Module');

    # 9 other optional tests that need Optional::Module
  }

  SKIP: {
    test_skip(-tests => 10, -requires_networking => 1);

    # 10 optional tests that require internet access (only makes sense in the
    # context of a script that doesn't use -requires_networking in the call to
    # &test_begin)
  }

  # in unix terms, we want to test with a file t/data/input_file.txt
  my $input_file = test_input_file('input_file.txt');

  # we want the name of a file we can write to, that will be automatically
  # deleted when the test script finishes
  my $output_file = test_output_file();

  # we want the name of a directory we can store files in, that will be
  # automatically deleted when the test script finishes
  my $output_dir = test_output_dir();

=head1 DESCRIPTION

This provides a common base for all BioPerl test scripts. It safely handles the
loading of Test::Most, itself a simple wrapper around several highly used test
modules: Test::More, Test::Exception, Test::Warn, Test::Deep, and Test::Diff. It
also presents an interface to common needs such as skipping all tests if
required modules aren't present or if network tests haven't been enabled. See
test_begin().

In the same way, it allows you to skip just a subset of tests for those same
reasons, in addition to requiring certain executables and environment variables.
See test_skip().

It also has two further methods that let you decide if network tests should be
run, and if debugging information should be printed. See test_network() and
test_debug().

Finally, it presents a consistent way of getting the path to input and output
files. See test_input_file(), test_output_file() and test_output_dir().

=cut

# TODO: Evil magic ahead; can we clean this up?

{
    my $Tester = Test::Builder->new;

    no warnings 'redefine';
    sub Test::Warn::_canonical_got_warning {
        my ($called_from, $msg) = @_;
        my $warn_kind = $called_from eq 'Carp' ? 'carped' : ($called_from =~ /Bio::/ ? 'Bioperl' : 'warn');

        my $warning;
        if ($warn_kind eq 'Bioperl') {
            ($warning) = $msg =~ /\n--------------------- WARNING ---------------------\nMSG: (.+)\n---------------------------------------------------\n$/m;
            $warning ||= $msg; # shouldn't ever happen
        }
        else {
            my @warning_stack = split /\n/, $msg;   # some stuff of uplevel is included
            $warning = $warning_stack[0];
        }

        return {$warn_kind => $warning}; # return only the real message
    }

    sub Test::Warn::_diag_found_warning {
        my @warns = @_;
        foreach my $warn (@warns) {
            if (ref($warn) eq 'HASH') {
                   ${$warn}{carped}  ? $Tester->diag("found carped warning: ${$warn}{carped}")
                : (${$warn}{Bioperl} ? $Tester->diag("found Bioperl warning: ${$warn}{Bioperl}")
                : $Tester->diag("found warning: ${$warn}{warn}"));
            } else {
                $Tester->diag( "found warning: $warn" );
            }
        }
        $Tester->diag( "didn't find a warning" ) unless @warns;
    }

    sub Test::Warn::_cmp_got_to_exp_warning {
        my ($got_kind, $got_msg) = %{ shift() };
        my ($exp_kind, $exp_msg) = %{ shift() };
        return 0 if ($got_kind eq 'warn') && ($exp_kind eq 'carped');

        my $cmp;
        if ($got_kind eq 'Bioperl') {
            $cmp = $got_msg =~ /^\Q$exp_msg\E$/;
        }
        else {
            $cmp = $got_msg =~ /^\Q$exp_msg\E at \S+ line \d+\.?$/;
        }

        return $cmp;
    }
}

our @EXPORT = (@Test::Most::EXPORT,
               #@Bio::Root::Test::Warn::EXPORT,
               # Test::Warn method wrappers

               # BioPerl-specific
               qw(
                test_begin
                test_skip
                test_output_file
                test_output_dir
                test_input_file
                test_network
                test_email
                test_debug
                float_is
             ));

our $GLOBAL_FRAMEWORK = 'Test::Most';
our @TEMP_FILES;

=head2 test_begin

 Title   : test_begin
 Usage   : test_begin(-tests => 20);
 Function: Begin your test script, setting up the plan (skip all tests, or run
           them all)
 Returns : True if tests should be run.
 Args    : -tests               => int (REQUIRED, the number of tests that will
                                        be run)
           -requires_modules    => []  (array ref of module names that are
                                        required; if any don't load, all tests
                                        will be skipped. To specify a required
                                        version of a module, include the version
                                        number after the module name, separated
                                        by a space)
           -requires_module     => str (as above, but for just one module)
           -requires_networking => 1|0 (default 0, if true all tests will be
                                        skipped if network tests haven't been
                                        enabled in Build.PL)
           -requires_email      => 1   (if true the desired number of tests will
                                        be skipped if either network tests
                                        haven't been enabled in Build.PL or an
                                        email hasn't been entered)
           -excludes_os         => str (default none, if OS suppied, all tests
                                        will skip if running on that OS (eg.
                                        'mswin'))
           -framework           => str (default 'Test::Most', the Test module
                                        to load. NB: experimental, avoid using)

           Note, supplying -tests => 0 is possible, allowing you to skip all
           tests in the case that a test script is testing deprecated modules
           that have yet to be removed from the distribution

=cut

sub test_begin {
    my ($skip_all, $tests, $framework) = _skip(@_);
    $GLOBAL_FRAMEWORK = $framework;

    if ($framework eq 'Test::Most') {
        # ideally we'd delay loading Test::Most until this point, but see BEGIN
        # block

        if ($skip_all) {
            eval "plan skip_all => '$skip_all';";
        }
        elsif (defined $tests && $tests == 0) {
            eval "plan skip_all => 'These modules are now probably deprecated';";
        }
        elsif ($tests) {
            eval "plan tests => $tests;";
        }

        return 1;
    }
    # go ahead and add support for other frameworks here
    else {
        die "Only Test::Most is supported at the current time\n";
    }

    return 0;
}

=head2 test_skip

 Title   : test_skip
 Usage   : SKIP: {
                   test_skip(-tests => 10,
                             -requires_module => 'Optional::Module 2.01');
                   # 10 tests that need v2.01 of Optional::Module
           }
 Function: Skip a subset of tests for one of several common reasons: missing one
           or more optional modules, network tests haven't been enabled, a
           required binary isn't present, or an environmental variable isn't set
 Returns : n/a
 Args    : -tests               => int (REQUIRED, the number of tests that are
                                        to be skipped in the event one of the
                                        following options isn't satisfied)
           -requires_modules    => []  (array ref of module names that are
                                        required; if any don't load, the desired
                                        number of tests will be skipped. To
                                        specify a required version of a module,
                                        include the version number after the
                                        module name, separated by a space)
           -requires_module     => str (as above, but for just one module)
           -requires_executable => Bio::Tools::Run::WrapperBase instance
                                       (checks WrapperBase::executable for the
                                        presence of a binary, skips if absent)
           -requires_env        => str (checks %ENV for a specific env. variable,
                                        skips if absent)
           -excludes_os         => str (default none, if OS suppied, desired num
                                        of tests will skip if running on that OS
                                        (eg. 'mswin'))
           -requires_networking => 1   (if true the desired number of tests will
                                        be skipped if network tests haven't been
                                        enabled in Build.PL)
           -requires_email      => 1   (if true the desired number of tests will
                                        be skipped if either network tests
                                        haven't been enabled in Build.PL or an
                                        email hasn't been entered)

=cut

sub test_skip {
    my ($skip, $tests, $framework) = _skip(@_);
    $tests || die "-tests must be a number greater than 0";

    if ($framework eq 'Test::Most') {
        if ($skip) {
            eval "skip('$skip', $tests);";
        }
    }
    # go ahead and add support for other frameworks here
    else {
        die "Only Test::Most is supported at the current time\n";
    }
}

=head2 test_output_file

 Title   : test_output_file
 Usage   : my $output_file = test_output_file();
 Function: Get the full path of a file suitable for writing to.
           When your test script ends, the file will be automatically deleted.
 Returns : string (file path)
 Args    : none

=cut

sub test_output_file {
    die "test_output_file takes no args\n" if @_;

    # RT 48813
    my $tmp = File::Temp->new();
    push(@TEMP_FILES, $tmp);
    close($tmp); # Windows needs this
    return $tmp->filename;
}

=head2 test_output_dir

 Title   : test_output_dir
 Usage   : my $output_dir = test_output_dir();
 Function: Get the full path of a directory suitable for storing temporary files
           in.
           When your test script ends, the directory and its contents will be
           automatically deleted.
 Returns : string (path)
 Args    : none

=cut

sub test_output_dir {
    die "test_output_dir takes no args\n" if @_;

    return tempdir(CLEANUP => 1);
}

=head2 test_input_file

 Title   : test_input_file
 Usage   : my $input_file = test_input_file();
 Function: Get the path of a desired input file stored in the standard location
           (currently t/data), but correct for all platforms.
 Returns : string (file path)
 Args    : list of strings (ie. at least the input filename, preceded by the
           names of any subdirectories within t/data)
           eg. for the file t/data/in.file pass 'in.file', for the file
           t/data/subdir/in.file, pass ('subdir', 'in.file')

=cut

sub test_input_file {
    return File::Spec->catfile('t', 'data', @_);
}

=head2 test_network

 Title   : test_network
 Usage   : my $do_network_tests = test_network();
 Function: Ask if network tests should be run.
 Returns : boolean
 Args    : none

=cut

sub test_network {
    return $ENV{AUTHOR_TESTING} || $ENV{RELEASE_TESTING};
}

=head2 test_email

 Title   : test_email
 Usage   : my $do_network_tests = test_email();
 Function: Ask if email address provided
 Returns : boolean
 Args    : none

=cut

sub test_email {
    return $ENV{AUTHOR_TESTING} || $ENV{RELEASE_TESTING};
}

=head2 test_debug

 Title   : test_debug
 Usage   : my $output_debugging = test_debug();
 Function: Ask if debugging information should be output.
 Returns : boolean
 Args    : none

=cut

sub test_debug {
    return $ENV{'BIOPERLDEBUG'} || 0;
}

=head2 float_is

 Title   : float_is
 Usage   : float_is($val1, $val2);
 Function: test two floating point values for equality
 Returns : Boolean based on test (can use in combination with diag)
 Args    : two scalar values (floating point numbers) (required via prototype)
           test message (optional)

=cut

sub float_is ($$;$) {
    my ($val1, $val2, $message) = @_;
    # catch any potential undefined values and directly compare
    if (!defined $val1 || !defined $val2) {
        is($val1, $val2 ,$message);
    } else {
        is(sprintf("%g",$val1), sprintf("%g",$val2),$message);
    }
}

=head2 _skip

Decide if should skip and generate skip message
=cut

sub _skip {
    my %args = @_;

    # handle input strictly
    my $tests = $args{'-tests'};
    #(defined $tests && $tests =~ /^\d+$/) || die "-tests must be supplied and be an int\n";
    delete $args{'-tests'};

    my $req_mods = $args{'-requires_modules'};
    delete $args{'-requires_modules'};
    my @req_mods;
    if ($req_mods) {
        ref($req_mods) eq 'ARRAY' || die "-requires_modules takes an array ref\n";
        @req_mods = @{$req_mods};
    }
    my $req_mod = $args{'-requires_module'};
    delete $args{'-requires_module'};
    if ($req_mod) {
        ref($req_mod) && die "-requires_module takes a string\n";
        push(@req_mods, $req_mod);
    }

    my $req_net = $args{'-requires_networking'};
    delete $args{'-requires_networking'};

    my $req_email = $args{'-requires_email'};
    delete $args{'-requires_email'};

    my $req_env = $args{'-requires_env'};
    delete $args{'-requires_env'};

    # strip any leading $ in case someone passes $FOO instead of 'FOO'
    $req_env =~ s{^\$}{} if $req_env;

    my $req_exe = $args{'-requires_executable'};
    delete $args{'-requires_executable'};

    if ($req_exe && (!ref($req_exe) || !$req_exe->isa('Bio::Tools::Run::WrapperBase'))) {
        die "-requires_exe takes an argument of type Bio::Tools::Run::WrapperBase";
    }

    my $os = $args{'-excludes_os'};
    delete $args{'-excludes_os'};

    my $framework = $args{'-framework'} || $GLOBAL_FRAMEWORK;
    delete $args{'-framework'};

    # catch user mistakes
    while (my ($key, $val) = each %args) {
        die "unknown argument '$key' supplied, did you mistake 'required...' for 'requires...'?\n";
    }

    # test user requirments and return
    if ($os) {
        if ($^O =~ /$os/i) {
            return ('Not compatible with your Operating System', $tests, $framework);
        }
    }

    foreach my $mod (@req_mods) {
        my $skip = _check_module($mod);
        if ($skip) {
            return ($skip, $tests, $framework);
        }
    }

    if ($req_net && ! test_network()) {
        return ('Network tests have not been requested', $tests, $framework);
    }

    if ($req_email && ! test_email()) {
        return ('Valid email not provided; required for tests', $tests, $framework);
    }

    if ($req_exe) {
        my $eval = eval {$req_exe->executable};
        if ($@ or not defined $eval) {
            my $msg = 'Required executable for '.ref($req_exe).' is not present';
            diag($msg);
            return ($msg, $tests, $framework);
        }
    }

    if ($req_env && !exists $ENV{$req_env}) {
        my $msg = 'Required environment variable $'.$req_env. ' is not set';
        diag($msg);
        return ($msg, $tests, $framework);
    }

    return ('', $tests, $framework);
}

=head2 _check_module

=cut

sub _check_module {
    my $mod = shift;

    my $desired_version;
    if ($mod =~ /(\S+)\s+(\S+)/) {
        $mod = $1;
        $desired_version = $2;
    }

    eval "require $mod;";

    if ($@) {
        if ($@ =~ /Can't locate/) {
            return "The optional module $mod (or dependencies thereof) was not installed";
        }
        else {
            return "The optional module $mod generated the following error: \n$@";
        }
    }
    elsif ($desired_version) {
        no strict 'refs';
        unless (defined ${"${mod}::VERSION"}) {
            return "The optional module $mod didn't have a version, but we want v$desired_version";
        }
        elsif (${"${mod}::VERSION"} < $desired_version) {
            return "The optional module $mod was out of date (wanted v$desired_version)";
        }
    }

    return;
}

1;
