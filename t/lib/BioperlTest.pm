# $Id$
#
# BioPerl module for BioperlTest
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Sendu Bala
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

BioperlTest - A common base for all Bioperl test scripts.

=head1 SYNOPSIS

  use lib 't/lib';
  use BioperlTest;
  
  test_begin(-tests => 20,
             -requires_modules => [qw(IO::String XML::Parser)],
             -requires_networking => 1);

  my $do_network_tests = test_network();
  my $output_debugging = test_debug();

  # carry out tests in Test::More syntax
  
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

=head1 DESCRIPTION

This provides a common base for all Bioperl test scripts. It presents an
interface to common needs such as skipping all tests if required modules aren't
present or if network tests haven't been enabled. See test_begin().

In the same way, it allows you to skip just a subset of tests for those same
reasons. See test_skip().

It also has two further methods that let you decide if network tests should be
run, and if debugging information should be printed. See test_network() and
test_debug().

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package BioperlTest;

use strict;
use warnings;

use File::Temp ();
use File::Spec;
use Exporter qw(import);

BEGIN {
    # For prototyping reasons, we have to load Test::More's methods now, even
    # though theoretically in future the user may use a different Test framework
    
    # We want to load Test::More. Preferably the users own version, but if they
    # don't have it, the one in t/lib. However, this module is in t/lib so
    # t/lib is already in @INC so Test::More in t/lib will be used first, which
    # we don't want: get rid of t/lib in @INC
    no lib 't/lib';
    eval { require Test::More; };
    if ($@) {
        eval "use lib 't/lib';";
    }
    eval "use Test::More;";
}

# re-export Test::More methods and export our own
our @EXPORT = qw(ok use_ok require_ok
                 is isnt like unlike is_deeply
                 cmp_ok
                 skip todo todo_skip
                 pass fail
                 eq_array eq_hash eq_set
                 $TODO
                 plan
                 can_ok  isa_ok
                 diag
                 BAIL_OUT
                 
                 test_begin
                 test_skip
                 test_output_file
                 test_input_file
                 test_network
                 test_debug);

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
        elsif ($tests == 0) {
            eval "plan skip_all => 'All tests are being skipped, probably because the module(s) being tested here are now deprecated';";
        }
        else {
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
                             -requires_modules => ['Optional::Module']);

                   # 10 tests that need Optional::Module
           }
 Function: Skip a subset of tests for one of two common reasons: missing one or
           more optional modules, or network tests haven't been enabled.
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
           -excludes_os         => str (default none, if OS suppied, desired num
                                        of tests will skip if running on that OS
                                        (eg. 'mswin'))
           -requires_networking => 1   (if true the desired number of tests will
                                        be skipped if network tests haven't been
                                        enabled in Build.PL)

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
    
    my $tmp = File::Temp->new();
    push(@TEMP_FILES, $tmp);
    
    return $tmp->filename;
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
    return $build->notes('network');
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


# decide if should skip and generate skip message
sub _skip {
    my %args = @_;
    
    # handle input strictly
    my $tests = $args{'-tests'};
    (defined $tests && $tests =~ /^\d+$/) || die "-tests must be supplied and be an int\n";
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
        return "The optional module $mod (or dependencies thereof) was not installed";
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
