#!/usr/bin/perl
#
=head1 NAME

pod.pl - check the POD documentation syntax in modules and scripts

=head1 SYNOPSIS

B<pod.pl> [B<-d|--dir> path ] [B<-v|--verbose>] B<-b|--blankline>
    [B<-?|-h|--help>]

=head1 DESCRIPTION

Checks Plain Old Documentation (POD) with highest possible stringency
in every bioperl module and script in CVS modules 'core' and 'run'.

Amounts to same as running

  podchecker -warnings -warnings

on every file.

=head2 Results

The results are written into file '/tmp/bioperl_pod_check' and
displayed after the run. The output is filtered not to show
confirmations of correct syntax. The result file is not removed.

The aim is to have as few warnings, and no errors, as possible.  Links
to web URLs give a warning but that seems to be spurious, so they are
filtered out.  Currently there are a few cases of "multiple occurrence
of link target" in several modules which are harmless.

=head1 SEE ALSO

L<podchecker>, L<Pod::Checker>

=cut

use File::Find;
use Pod::Checker;
use Getopt::Long;
use strict;

sub podcheck;
sub blankline;

#
## Directories to check
#
my @dirs = qw( ../Bio/ ../scripts . );

# command line options
my ($verbose, $blankline, $dir, $help) = (0, undef, undef, undef);
GetOptions(
           'v|verbose' => \$verbose,
           'dir:s' => \$dir,
           'blankline' => \$blankline,
	   'h|help|?' => sub{ exec('perldoc',$0); exit(0) }
	   );

# setup
my $tmpfile = '/tmp/bioperl_pod_check';
our %POD_CHECKER_OPTIONS = ( '-warnings' => 2 );
our %FIND_OPTIONS = ( wanted => \&podcheck, no_chdir => 1 );

# run
open (F, ">$tmpfile") || die "can't open file $tmpfile: $!";
$FIND_OPTIONS{wanted} = \&blankline if $blankline;

if ($dir) {
    find \%FIND_OPTIONS, $dir;
} else {
    find \%FIND_OPTIONS, @dirs;
}
close F;
open (F, "grep -v OK $tmpfile|") || die "can't open file $tmpfile: $!";
while (<F>) { print unless /http/ and /non-escaped/ }


# this is where the action is
sub podcheck {
    return unless /\.PLS$/ or /\.p[ml]$/ ;
    return unless -e $_;
    print "$_\n" if $verbose;
    my $checker = Pod::Checker->new( %POD_CHECKER_OPTIONS );
    $checker->parse_from_file($_, \*F);
    print "$_\tno POD\n" if $checker->num_errors() < 0;
}

=head1 OPTIONS

=over 3

=item B<-d | --dir> path

Overides the default directories to check by one directory 'path' and
all its subdirectories.

=item B<-b | --blankline>

Checks POD command paragraphs (lines starting with '=' character) for
preceding nonblank lines. These lines are printed out with '++'.

Also, if verbose is turned on, it will report on lines whitespace
characters which prevent paragrafs to be recognised by older POD
parsers (marked with '+'). Modern perlpod parsers (5.6.0 and later, I
suppose) allow for whitespace lines surrounding command lines, but
since bioperl still supports older versions, these lines should be
cleaned to contain only '\n' and no space or tab characters.


See: L<perlpodspec>


=cut

sub blankline {
    return unless /\.PLS$/ or /\.p[ml]$/ ;
    return unless -e $_;
    my $file = $_;
    open (F, $_) or warn "can't open file $_: $!" && return;
    local $/="";
    while (<F>) {
        print "$file: +|$1|\n" if /[ \t]\n(=[a-z][^\n]+$)/m and $verbose;
        print "$file: ++|$1|\n" if /\w\n(=[a-z][^\n]+$)/m and $verbose;
        print "$file:|$1|+\n" if /(^=[a-z][^\n]+)\n[\t ]/m;
        #print "$file:|$1|++\n" if /(^=[^\n]+)\n\w/m;
    }
    close F;
}

__END__

=item B<-v | --verbose>

Show the progress through files during the POD checking.

=item B<-? | -h  | --help>

This help text.

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Heikki Lehvaslaiho

Email heikki-at-bioperl-dot-org

=cut


# find . -name '*.pm' -print | xargs  perl -e '$/=""; while (<>) {$n = $1 if /^package\s+([\w:]+)/; print "$n:|$1|"  if  /(\s\s^=[^\n]+$)/m ; }'  ;

# find . -name '*.pm' -print | xargs  perl -e '$/=""; while (<>) {$n = $1 if /^package\s+([\w:]+)/; print "$n:|$1|\n"  if /(^=[^\n]+\n[\t ])/m ; }'  ;
