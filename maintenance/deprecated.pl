#!/usr/bin/perl

use strict;
use warnings;

use version;
use Bio::Root::Version;
use File::Find;
use Getopt::Long;
use Perl6::Form;
use Carp;

#
# command line options
#

my ($verbose, $dir, $depfile, $help, $new, $outfile, $write, $version) =
(0, undef, "../DEPRECATED", undef, [], '../DEPRECATED.NEW', 0, $Bio::Root::Version::VERSION);
GetOptions(
        'v|verbose' => \$verbose,
		'b|bp_version:s' => \$version,
        'dir:s' => \$dir,
        'depfile:s' => \$depfile,
        'n|new=s@' => \$new,
        'o|outfile:s' => \$outfile,
        'w|write' => \$write,
        'h|help|?' => sub{ exec('perldoc',$0); exit(0) }
	   );

# Default directories to check
my @dirs = qw(../Bio/ );

# use version to consolidate old vs new versioning schemes
my $base_version = version->new( $version );

print "Version: $base_version\n";

my %deprecated;
my %removed;
my @dep_data;

# parse DEPRECATED file

open my $DFILE, '<', $depfile or die "Could not read file '$depfile': $!\n";
my $seen_top;
while (my $data = <$DFILE>) {
    if ($data =~ /^-+$/) {
        $seen_top = 1;
        next;
    }
    next unless $seen_top;
    chomp $data;
    my ($module, $dep, $rem, $note) = split(/\s+/,$data,4);
    next unless $module;
    my $d = version->new($dep);
    my $r = version->new($rem);
	print "$module Dep: $d Rem: $r\n" if $verbose;
    if ($rem <= $base_version) {
        $removed{$module}++;
    } elsif ($dep <= $base_version) {
        $deprecated{$module}++;
    }
    push @dep_data, {module => $module,
                     dep => $dep,
                     remove => $rem,
                     note => $note}
}
close $DFILE;

for my $new (@$new) {
    my ($module, $dep, $rem, $note) = split(',',$new,4);
    last if !$module || !$dep || !$rem;
    if ($module !~ /Bio/) {
        croak "Can only deprecate BioPerl modules, not $module"
    }
    push @dep_data, {module => $module,
                 dep => $dep,
                 remove => $rem,
                 note => $note}
}

# run through all files in core (checks to see if anything is still present)

if ($dir) {
    find {wanted => \&parse_core, no_chdir => 1}, $dir;
} else {
    find {wanted => \&parse_core, no_chdir => 1}, @dirs;
}

#
# results
#

# uses Perl6::Form

if ($write || @$new) {

open my $NEWDEP, '>', $outfile or croak "Could not write file '$outfile': $!\n";

print $NEWDEP <<HEAD;
# These are modules which are deprecated and later removed from the toolkit
# See http://www.bioperl.org/wiki/Deprecated_modules for the latest details

HEAD

# may replace with better formatting, but it needs to be round-tripped

print $NEWDEP form
"Deprecated                     Version    Version                                             ",
"Module                        Deprecated  Removed  Notes                                      ",
"----------------------------------------------------------------------------------------------";

for my $datum (@dep_data) {
    my ($mod, $dep, $rem, $note) = map {$datum->{$_}} qw (module dep remove note);
    
print $NEWDEP form
"{[[[[[[[[[[[[[[[[[[[[[[[[[[[[[} {|||||}  {|||||}   {[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[}",
    $mod,     $dep,        $rem,       $note;
}

}

#
##
### end main
##
#

#
# this is where the action is
#

sub parse_core {
    my $file = $_;
    return unless $file =~ /\.PLS$/ || $file =~ /\.p[ml]$/ ;
    return unless -e $file;
    open my $F, '<', $file or die "Could not read file '$file': $!\n";
    while (my $line = <$F>) {
        if ($line =~ /(?:['"])?\b(use|require)\s+([A-Za-z0-9:_\.\(\)]+)\s*([^;'"]+)?(?:['"])?\s*;/) {
            my ($use, $mod) = ($1, $2);
			if (exists $removed{$mod}) {
                print "$File::Find::name: Line $.: $mod is removed\n";
            } elsif (exists $deprecated{$mod}) {
                print "$File::Find::name: Line $.: $mod is deprecated\n";
            } 
        }
    }
    close $F;
}

# $Id: deprecated.pl 10084 2006-07-04 22:23:29Z mauricio $
#
=head1 NAME

deprecated.pl - Check modules and scripts for use of deprecated modules and
methods, indicates presence in a file to STDERR. Optionally accepts new modules
and adds them to a newly formatted deprecation file.

=head1 SYNOPSIS

B<deprecated.pl> [B<-d|--dir> path ] [B<-v|--verbose>] [B<-a|--depfile>]
    [B<-n|--new>] [B<-w|--write>] [B<-o|--outfile>]
    [B<-?|-h|--help>]

=head1 OPTIONS

=over 3

=item B<-d | --dir> path

Overides the default directories to check by one directory 'path' and
all its subdirectories.

=item B<-a | --depfile>

path from working directory that contains the DEPRECATED file.

=item B<-n | --new>

New addition to the deprecation list; this should be in the form of
'Module,dep_release,remove_release,notes'. Notes should only be 40 chars long.

=item B<-b | --bp_version>

BioPerl version.  This only appears to work correctly when using numerical
versions (1.5.2 instead of 1.005002)

=item B<-w | --write>

Write out new deprecation file to $outfile.  If --new is used this is assumed.

=item B<-o | --outfile>

Name of output file to write deprecation table to. DEPRECATED.NEW is the default
name

=item B<-v | --verbose>

Show the progress through files during the checking.

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

=head1 AUTHOR - Chris Fields

Email cjfields-at-bioperl-dot-org

=cut
