#!/usr/bin/perl

=head1 NAME

check_NAMEs.pl - check NAME in module POD has fully qualified object name

=head1 SYNOPSIS

B<check_NAMEs.pl> [B<-d|--dir> path] [B<-v|--verbose>] [B<-?|-h|--help>]

=head1 DESCRIPTION

This script is designed to find all Bioperl modules which don't
have the fully qualified object name with correct capitalization
in the "NAME" section of the POD. 

The full name is required for the PDOC POD to HTML script to
correctly render the module documentation.

=cut

use strict;
use File::Find;
use Getopt::Long;

#
# command line options
#

my ($verbose, $dir, $help) = (0, '../Bio/', undef);
GetOptions(
    'v|verbose' => \$verbose,
    'd|dir:s' => \$dir,
    'h|help|?' => sub{ exec('perldoc',$0); exit(0) }
);

#
# globals
#

my $num_found = 0;

#
# find all modules
#

print STDERR "Searching for incorrect NAME POD section of all modules in: $dir\n";
find( \&find_modules, $dir );
print STDERR "$num_found found.\n";

# this is where the action is

sub find_modules {
    # only want files with .pm
    return unless m/\.pm$/;
    return unless -f $_;
    
    my $fname = $_;
    my $pm = $File::Find::name;
    $pm =~ s{.*?/(?=Bio/)}{};  # remove up to first slash before Bio/
    $pm =~ s{\.pm$}{};         # remove .pm suffix
    $pm =~ s{/}{::}g;          # convert / to ::
    
    print STDERR "# $File::Find::name\n" if $verbose;
    
    # slurp in the file
    my $text = do { local( @ARGV, $/ ) = $fname ; <> } ;

    # check if the NAME section has the _full_ module name in it
    if ($text !~ m/^=head1\s+NAME.*?^$pm/xms) {
      print "$pm\n";
      $num_found++;
    }
}


=head1 OPTIONS

=over 3

=item B<-d | --dir> path

Overides the default directory to recursively look for .pm file 
(Default is '../Bio')

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

=head1 AUTHOR - Torsten Seemann

Email: torsten-dot-seemann-at-infotech-dot-monash-dot-edu-dot-au

=cut
