#!/usr/bin/perl -w

=head1 NAME

pod.pl - check the POD documentation syntax in modules and scripts

=head1 SYNOPSIS

B<pod.pl> [B<-d|--dir> path ] [B<-v|--verbose>] [B<-?|-h|--help>]

=head1 DESCRIPTION

Checks Plain Old Documentation (POD) with highest possible stringency
in every bioperl module and script in CVS modues core and run.

Amounts to same as running

  podchecker -warnings -warnings

on every file.

You are expected to have checked out CVS module 'bioperl_all'.
Otherwise, bioperl-run module is not found.


=head2 Results

The results are written into file '/tmp/bioperl_pod_check' and
displayed after the run. The outut is filtered not to show
conformations of correct syntax. The result file is not removed.

The aim is to have as few warnings, and no errors, as possible.  Links
to web URLs give a waring but that seems to be spurious, so they are
filtered out.  Currently there are a few cases of "multiple occurrence
of link target" in several modlues which do no harm.

=head1 SEE ALSO

L<podchecker>, L<Pod::Checker>

=cut

use File::Find;
use Pod::Checker;
use Getopt::Long;
use strict;

#
## Directories to check
#
my @dirs = qw( ../Bio/ ../../run/Bio  ../scripts ../../run/scripts . );

# command line options
my ($verbose, $dir, $help) = (0, undef, undef);
GetOptions(
           'v|verbose' => \$verbose,
           'dir:s' => \$dir,
	   'h|help|?' => sub{ exec('perldoc',$0); exit(0) }
	   );

# setup
my $tmpfile = '/tmp/bioperl_pod_check';
our %POD_CHECKER_OPTIONS = ( '-warnings' => 2 );
our %FIND_OPTIONS = ( wanted => \&podcheck, no_chdir => 1 );

# run
open (F, ">$tmpfile") || die "can't open file $tmpfile: $!";
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
    my $checker = new Pod::Checker %POD_CHECKER_OPTIONS;
    $checker->parse_from_file($_, \*F);
}

__END__

=head1 OPTIONS

=over 2

=item B<-d | --dir> path

Overides the default directories to check by one directory 'path' and
all its subdirectories.

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

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Heikki Lehvaslaiho

Email heikki@ebi.ac.uk

=cut
