#
# BioPerl module for Bio::Root::Test
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Sendu Bala
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Root::Test - A common base for all Bioperl test scripts.

=head1 SYNOPSIS

  use lib '.'; # (for core package tests only)
  use Bio::Root::Test;

  test_begin(-tests => 20,
             -requires_modules => [qw(IO::String XML::Parser)],
             -requires_networking => 1);

  my $do_network_tests = test_network();
  my $output_debugging = test_debug();

  # carry out tests with Test::More, Test::Exception and Test::Warn syntax

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
loading of Test::More, Test::Exception and Test::Warn (actually, a subclass
compatible with Bioperl warnings) prior to tests being run. It also presents an
interface to common needs such as skipping all tests if required modules aren't
present or if network tests haven't been enabled. See test_begin().

In the same way, it allows you to skip just a subset of tests for those same
reasons, in addition to requiring certain executables and environment variables.
See test_skip().

It also has two further methods that let you decide if network tests should be
run, and if debugging information should be printed. See test_network() and
test_debug().

Finally, it presents a consistent way of getting the path to input and output
files. See test_input_file(), test_output_file() and test_output_dir().

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Root::Test;

use strict;
use warnings;

use File::Temp qw(tempdir);
use File::Spec;
use Exporter qw(import);

BEGIN {
    # For prototyping reasons, we have to load Test::More's methods now, even
    # though theoretically in future the user may use a different Test framework
    
    # We want to load Test::More, Test::Exception and Test::Warn. Preferably the
    # users own versions, but if they don't have them, the ones in t/lib.
    # However, this module is in t/lib so t/lib is already in @INC so Test::* in
    # t/lib will be used first, which we don't want: get rid of t/lib in @INC
    no lib 't/lib';
    eval { require Test::More;
           require Test::Exception;
           require Test::Warn; };
    if ($@) {
        eval "use lib 't/lib';";
    }
    eval "use Test::More;
          use Test::Exception;";
    die "$@\n" if $@;
    
    # now that the users' Test::Warn has been loaded if they had it, we can
    # use Bio::Root::TestWarn
    eval "use Bio::Root::Test::Warn;";
    die "$@\n" if $@;
}

# re-export Test::More, Test::Exception and Test::Warn methods and export our own
our @EXPORT = qw(ok use_ok require_ok
                 is isnt like unlike is_deeply
                 cmp_ok
                 skip todo todo_skip
                 pass fail
                 eq_array eq_hash eq_set
                 $TODO
                 plan
                 can_ok isa_ok
                 diag
                 BAIL_OUT
                 
                 dies_ok
                 lives_ok
                 throws_ok
                 lives_and
                 
                 warning_is
                 warnings_are
                 warning_like
                 warnings_like
                 
                 test_begin
                 test_skip
                 test_output_file
                 test_output_dir
                 test_input_file
                 test_network
                 test_email
                 test_debug
                 float_is
                 );

if (Test::More->can('done_testing')) {
    push @EXPORT, 'done_testing';
}

our $GLOBAL_FRAMEWORK = 'Test::More';
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
           -framework           => str (default 'Test::More', the Test module
                                        to load. NB: experimental, avoid using)
           
           Note, supplying -tests => 0 is possible, allowing you to skip all
           tests in the case that a test script is testing deprecated modules
           that have yet to be removed from the distribution

=cut

sub test_begin {
    my ($skip_all, $tests, $framework) = _skip(@_);
    $GLOBAL_FRAMEWORK = $framework;
    
    if ($framework eq 'Test::More') {
        # ideally we'd delay loading Test::More until this point, but see BEGIN
        # block
        
        if ($skip_all) {
            eval "plan skip_all => '$skip_all';";
        }
        elsif (defined $tests && $tests == 0) {
            eval "plan skip_all => 'All tests are being skipped, probably because the module(s) being tested here are now deprecated';";
        }
        elsif ($tests) {
            eval "plan tests => $tests;";
        }
        
        return 1;
    }
    # go ahead and add support for other frameworks here
    else {
        die "Only Test::More is supported at the current time\n";
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
    
    if ($framework eq 'Test::More') {
        if ($skip) {
            eval "skip('$skip', $tests);";
        }
    }
    # go ahead and add support for other frameworks here
    else {
        die "Only Test::More is supported at the current time\n";
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
    require Module::Build;
    my $build = Module::Build->current();
    return $build->notes('Network Tests');
}

=head2 test_email

 Title   : test_email
 Usage   : my $do_network_tests = test_email();
 Function: Ask if email address provided
 Returns : boolean
 Args    : none

=cut

sub test_email {
    require Module::Build;
    my $build = Module::Build->current();
    # this should not be settable unless the network tests work
    return $build->notes('email');
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

# decide if should skip and generate skip message
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

    if ($req_exe && !$req_exe->executable) {
        my $msg = 'Required executable for '.ref($req_exe).' is not present';
        diag($msg);
        return ($msg, $tests, $framework);
    }
    
    if ($req_env && !exists $ENV{$req_env}) {
        my $msg = 'Required environment variable $'.$req_env. ' is not set';
        diag($msg);
        return ($msg, $tests, $framework);
    }
    
    return ('', $tests, $framework);
}

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
