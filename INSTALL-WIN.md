# Installing BioPerl on Windows

# Introduction

This installation guide was written by Barry Moore, Nathan Haigh
and other BioPerl authors based on the original work of Paul Boutros. The
guide has been updated by Paul Cantalupo and Francisco Ossandon.
*Note* For Windows it is recommended to install BioPerl version 1.6.924 or
newer, since many Windows-specific bugs from previous versions were fixed.

Please report problems and/or fixes to the BioPerl mailing list.

# Perl Options

There are a couple of ways of installing Perl on a Windows machine. One is
to get the most recent build from [Strawberry Perl](http://strawberryperl.com/)
and the other is to get it from [ActiveState](http://www.activestate.com/).
Both are software companies that provide free builds of Perl for Windows
users.

Strawberry Perl is recommended since is more CPAN-friendly and
because it includes a compiler (gcc), related tools and other external
libraries. These installation steps were verified on December 2015 using
ActivePerl 5.20.2.2002 (5.22 has been released but the MinGW package is not
available for it yet) and Strawberry Perl 5.22.0.1 from a clean install.

*Note* Only ActivePerl **5.18 or greater** is supported by the BioPerl team.
This is because the necessary MinGW package needed for CPAN installations
is only available for both 32-bit and 64-bit Windows since this version
(32-bit was available on previous versions but only in the Business edition,
see [ActivePerl MinGW PPM webpage](http://www.activestate.com/activeperl/downloads)).

*Note* Support for installation through ActivePerl Perl Package Manager has been dropped in favor of CPAN.

# Installing Perl on Windows

1. Download the [Strawberry Perl MSI](http://strawberryperl.com/releases.html) or [ActivePerl MSI](http://www.activestate.com/activeperl/downloads).
2. Run the Installer (accepting all defaults is fine).

You can also build Perl yourself (which requires a C compiler). The Perl
source for building it yourself is available from [CPAN](http://www.cpan.org/).
This approach is not recommended unless you have specific reasons for doing so
and know what you're doing.

[Cygwin](http://en.wikipedia.org/wiki/Cygwin) is a [UNIX](http://en.wikipedia.org/wiki/UNIX)
emulation environment for Windows and comes with its own copy of Perl. Information on Cygwin and BioPerl is found below.

### Installation using the ActiveState Perl Package Manager

This installation is not supported nor recommended anymore, because the last
BioPerl package produced for it was the old version 1.6.1 (2009) and many
Windows-specific bugs were fixed in more recent versions 1.6.923 and 1.6.924.
Please install using the CPAN instructions below.

# CPAN for Strawberry Perl and ActivePerl

Both CPAN and manual methods ultimately need the accessory compiling program
MinGW, which incorporates the necessary tools like `dmake` and `gcc`. MinGW comes
by default with Strawberry Perl, but must it be installed through PPM for
ActivePerl (CPAN will display a warning if not present). Also, CPAN need to be
upgraded to >= v1.81, [Module::Build](https://metacpan.org/pod/Module::Build)
to be upgraded (>= v0.2805) and [Test::Harness](https://metacpan.org/pod/Test::Harness) to be upgraded to >= v2.62.

### Dmake for ActivePerl

1) Install MinGW package through PPM: Using a cmd window type

`ppm install MinGW`.

It is `important` to check if ActiveState provides the
[MinGW package](http://code.activestate.com/ppm/MinGW/) for your ActivePerl version,
since each version have to wait its own release. For example,
although MinGW has been available since ActivePerl version 5.18,
the release for newest version 5.22 it's still not available (December 2015).

### Start the install with CPAN

1. Open a cmd window by going to `Start >> Run` and typing `cmd`
into the box and pressing return.

2. Type `cpan` to enter the CPAN shell.

3. At the `cpan>` prompt, type `install CPAN` to upgrade to the
latest version.

4. Quit (by typing 'q') and reload CPAN. You may be asked some
configuration questions; accepting defaults is fine.

5. At the `cpan>` prompt, type `o conf prefer_installer MB` to tell
CPAN to prefer to use `Build.PL` scripts for installation, and the type `o conf commit` to save that choice.

6. At the `cpan>` prompt, type `install Module::Build`.

7. At the `cpan>` prompt, type `install Test::Harness`.

8. At the `cpan>` prompt, type `install Test::Most`.

### Finish the install with CPAN

You can now follow the instructions INSTALLING BioPerl THE EASY WAY
USING CPAN in the INSTALL file.

### Finish the install with BioPerl from GitHub

For the bleeding edge version install manually using a ZIP file from the
GitHub repository:

1. Go to [GitHub](https://github.com/BioPerl/BioPerl-live) and press the
`Download ZIP` button.

2. Extract the archive in the normal way.

3. In a cmd window `cd` to the directory you extracted to. Eg. if
you extracted to directory 'BioPerl-live', `cd BioPerl-live`

4. Type `perl Build.PL` and answer the questions appropriately.

5. Type `perl Build test`. All the tests should pass, but if they don't,
[let us know](https://github.com/BioPerl/BioPerl-live/issues). Your usage of
BioPerl may not be affected by the failure, so you can choose to continue
anyway.

6. Type `perl Build install` to install BioPerl.

# BioPerl in Cygwin

Cygwin is a Unix emulator and shell environment available free at
http://www.cygwin.com. Some users claim that installation of BioPerl is
easier within Cygwin than within Windows, but these may be users with UNIX
backgrounds. A note on Cygwin: it doesn't write to your Registry, it doesn't
alter your system or your existing files in any way, it doesn't create
partitions, it simply creates a `cygwin/` directory and writes all of its
files to that directory. To uninstall Cygwin just delete that directory.

To get BioPerl running first install the basic Cygwin package as well as
the Cygwin `perl`, `make`, `binutils`, and `gcc` packages. Clicking the View
button in the upper right of the installer window enables you to see
details on the various packages. Then start up Cygwin and follow the
BioPerl installation instructions for UNIX in BioPerl's `INSTALL` file.

## Cygwin paths

If you're trying to use some application or resource outside of the Cygwin
directory and you're having a problem remember that Cygwin's path syntax
may not be the correct one. Cygwin understands `/home/jacky` or
`/cygdrive/e/cygwin/home/jacky` (when referring to the E: drive) but the
external resource may want `E:/cygwin/home/jacky`. So your `*rc` files may end
up with paths written in these different syntaxes, depending.

For example, here's how to set the environmental variable TMPDIR, programs
like BLAST and clustalw need a place to create temporary files:

```
setenv TMPDIR e:/cygwin/tmp    # csh, tcsh
export TMPDIR=e:/cygwin/tmp    # sh, bash
```
