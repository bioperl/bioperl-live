# BioPerl Installation

There are two principle ways to get BioPerl onto your system. One (and
the still most common) is to install it locally on your machine. To do
this, refer to and follow the instructions below for installing BioPerl on
Unix, Linux, and Mac OS X. For installing BioPerl on Windows, please
read INSTALL.WIN. The other way is to use a Docker image with a
BioPerl build, whether one we support or one you build yourself.

# Installing BioPerl locally

Note that this documentation is for Unix, Linux, and MacOSX. For installing BioPerl on Windows, please read INSTALL.WIN.

## System Requirments

 * `Perl 5.6.1 or higher` Version 5.10 or higher is highly
   recommended. Modules are tested against version 5.14 and
   above.
 * `make` For Mac OS X, this requires installing the Xcode Developer
   Tools.

## Preliminary Preparation

These are optional, but regardless of your subsequent choice of
installation method, it will help to carry out the following steps.
They will increase the likelihood of installation success
(especially of the optional dependencies).

* Upgrade CPAN:

```
perl -MCPAN -e shell
```

* Install or upgrade `Module::Build`, and make it your preferred installer:

```
cpan>install Module::Build
cpan>o conf prefer_installer MB
cpan>o conf commit
```

* Install the `expat` library by whatever method is appropriate for your system (e.g. `apt`, `yum`, `homebrew`).

## Installing BioPerl the Easy Way

We highly recommend using
[cpanminus](https://metacpan.org/pod/distribution/App-cpanminus/bin/cpanm) for
installing BioPerl and its dependencies. We also highly recommend (if possible)
using a tool like [perlbrew](https://perlbrew.pl) to locally install a modern
version of perl (a version that is higher than perl 5.16).  The linked
pages describe how to install each tool; make sure if you install perlbrew that
you follow up with installing `cpanm`:

```
perlbrew install-cpanm
```

Then, you can install BioPerl:

```
cpanm Bio::Perl
```

You can also use the older CPAN shell to install BioPerl. For example:

```
perl -MCPAN -e shell
```

Or you might have the `cpan` alias installed:

```
cpan
```

Then find the name of the latest BioPerl package:

```
cpan>d /bioperl/
 ....

 Distribution    C/CJ/CJFIELDS/BioPerl-1.007001.tar.gz
 Distribution    C/CJ/CJFIELDS/BioPerl-1.007001.tar.gz
 Distribution    C/CJ/CJFIELDS/BioPerl-1.007001.tar.gz
```

And install the most recent:

```
cpan>install C/CJ/CJFIELDS/BioPerl-1.007001.tar.gz
```

If you've installed everything perfectly and all the network
connections are working then you will pass all the tests run in the
`./Build test` phase. Sometimes you may see a failed test. Remember that
there are over 900 modules in BioPerl and the test suite is running more
than 12000 individual tests, a failed test may not affect your usage
of BioPerl.

If there's a failed test and you think that the failed test will not
affect how you intend to use BioPerl then do:

```
cpan>force install C/CJ/CJFIELDS/BioPerl-1.007001.tar.gz
```

If you're concerned about a failed test and need assistance or advice
then contact bioperl-l@bioperl.org, and provide us the detailed
results of the failed install.

## Installing BioPerl from Github

**NOTE:** We generally do not recommend installing the latest code from Github
unless you absolutely need the latest bug fixes.

The very latest version of Bioperl is at github.com. If you want this
version then download it from https://github.com/bioperl/bioperl-live as a
`tar.gz` or `zip` file, or retrieve it using the command line:

```
git clone https://github.com/bioperl/bioperl-live.git
cd bioperl-live
```

### Using cpanm

If you have `cpanm`, you can install within the checkout directory by simply using:

```
cpanm --interactive .
```

to run interative installation, or you can leave out the `--interactive` flag to accept the defaults.

### Using the Installation Script

Issue the build commands:

```
perl Build.PL
```

You will be asked a few questions about installing BioPerl scripts
and running various test suites, hit `return` to accept the defaults.

Test:

```
./Build test
```

Install:

```
./Build install
```

You may need root permissions in order to run `./Build install`, so you
will want to talk to your systems manager if you don't have the necessary
privileges. Or you can install the package in your own home
directory, see INSTALLING BIOPERL USING `local::lib`.

## Installing Bioperl using `local::lib`

If you lack permission to install Perl modules into the standard
system directories you can install them in your home directory
using [local::lib](https://metacpan.org/pod/local::lib). The instructions for first installing
`local::lib` are found in the link.

Once `local::lib` is installed you can install BioPerl using a
command like this:

```
perl -MCPAN -Mlocal::lib -e 'CPAN::install(C/CJ/CJFIELDS/BioPerl-1.007001.tar.gz)'
```

## Installing BioPerl Scripts

BioPerl comes with a set of production-quality scripts that are
kept in the scripts/ directory. You can install these scripts if you'd
like, simply answer the questions during `perl Build.PL`.
The installation directory can be specified by:

```
perl Build.PL ./Build install --install_path script=/foo/scripts
```

By default they install to `/usr/bin` or similar, depending on platform.

## The Test System

The BioPerl test system is located in the `t/` directory and is
automatically run whenever you execute the `./Build test` command.

The tests have been organized into groups
based upon the specific task or class the module being tested belongs
to. If you want to investigate the behavior of a specific test such as
the Seq test you would type:

```
./Build test --test_files t/Seq/Seq.t --verbose
```

The `--test_files` argument can be used multiple times to try a set of test
scripts in one go. The `--verbose` arguement outputs the detailed test results, instead of just the summary you see during `./Build test`.

The `--test-files` argument can also work as a glob. For instance, to run tests on all SearchIO modules, use the following:

```
./Build test --test_files t/SearchIO* --verbose
```

You can also use the command-line tool `prove` to run tests as well, which
is quite useful if you are developing code:

```
prove -lrv t/SearchIO*
```

If you are trying to learn how to use a module, often the test suite
is a good place to look. All good extreme programmers try and write a
test BEFORE they write the module to insure that their module behaves
the way they expect. You'll notice some `ok` and `skip` commands in a
test, this is part of the Perl test suite that signifies a passed test
with an 'ok N', where N is the test number. Alternatively you can tell
Perl to skip tests. This is useful when, for example, your test
detects that the network is not present and thus should skip, not
fail, any tests that require a network connection.

The core developers have indicated that future releases of BioPerl
will require that new modules come with a test suite with some minimal
tests.  Modules that lack adequate tests or could otherwise be
considered 'unstable' will be moved into a separate developer
distribution until adequate tests are added and the API stablizes.

[how to install Docker]: https://docs.docker.com/engine/installation/
[bioperl/bioperl]: https://hub.docker.com/r/bioperl/bioperl/
[bioperl/bioperl-deps]: https://hub.docker.com/r/bioperl/bioperl-deps/

# Using BioPerl via Docker

If you don't have Docker installed already, instructions for [how to install Docker] on Linux, MacOSX, and Windows are available online.

We officially support several builds (latest, stable, and releases)
hosted in the [bioperl/bioperl] repo on Docker Hub. These images do not
have a pre-defined entrypoint. If you have a BioPerl script in the
current directory, you can run it as simple as this:

```
docker run -t --rm -v `pwd`:/work -w /work bioperl/bioperl perl my-script.pl
```

Or run an interactive shell:

```
docker run -ti --rm -v `pwd`:/work -w /work bioperl/bioperl bash
```

You can also build your own Docker image of BioPerl, using the same
base image and pre-built dependencies that we use. Simply build off of
the [bioperl/bioperl-deps] image.

