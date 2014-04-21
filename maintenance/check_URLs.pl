#!/usr/bin/perl

=head1 NAME

check_URLs.pl - validate URLs located in module code and POD

=head1 SYNOPSIS

B<check_URLs.pl> [B<-d|--dir> path] [B<-v|--verbose>] [B<-?|-h|--help>]
    [B<-o|--outfile> filename]

=head1 DESCRIPTION

Checks code and POD of all bioperl-live modules for URLs, and validates them.

Output is a series of lines containing two fields, tab separated.
The first field is the file with the bad URL, the second is the URL itself.

The whole URL is not retrieved, only a HTTP "HEAD" request is done
to see if the URL exists on the server. The request is done using 
B<LWP::Simple> so the B<http_proxy> environmental variable will be
honoured.

The URL parsing may not be perfect - although I use the B<Regexp::Common::URI>
module, I have to manually clean up some URLs which are embedded in Perl
strings to convert the matched URL to a more probable real world URL,
e.g. most URLs don\'t end in "'," or ")" :-)

=cut

use strict;
use Data::Dumper;
use File::Find;
use Getopt::Long;
use Regexp::Common qw(URI);
use LWP::Simple qw($ua head);

$ua->timeout(15);

#
# command line options
#

my ($verbose, $dir, $help) = (0, '../Bio/', undef);
my $file;
GetOptions(
    'v|verbose' => \$verbose,
    'd|dir:s' => \$dir,
    'o|outfile:s' => \$file,
    'h|help|?' => sub{ exec('perldoc',$0); exit(0) }
);

my $fh;

if (defined $file) {
    open $fh, '>', $file or die "Could not write file '$file': $!\n";
} else {
    $fh = \*STDOUT;
}

#
# find all modules
#

find( \&find_modules, $dir );

#
# validate unique URLs and print fail cases to stdout
#

my %cached_urls;

sub check_url {
    my ($url, $file) = @_;
    if (exists $cached_urls{$url}) {
        print STDERR "$url checked in ".$cached_urls{$url}[0].":".$cached_urls{$url}[1]."\n" if $verbose;
        print $fh "$file\t$url\n" if $cached_urls{$url}[1] ne 'ok';
        return;
    }
    print STDERR "Checking $url in $file... " if $verbose;
    my $ok = head($url);
    my $status = $ok ? 'ok' : 'FAIL!';
    print STDERR "$status!\n" if $verbose;
    print $fh "$file\t$url\n" if !$ok;
    $cached_urls{$url} = [$file, $status];
}

close $fh if $file; # don't close STDOUT

#
# this is where the action is
#

sub find_modules {
    # only want files with .pm
    return unless m/\.pm$/;
    return unless -f $_;
    
    my $fname = $_;
    
    # slurp in the file
    my $text = do { local( @ARGV, $/ ) = $fname ; <> } ;
    
    # keep track of URLs
    while ($text =~ m/$RE{URI}{HTTP}{-keep}/g) {
        my $url = $1 or next;
        # remove Perl code if URL was embedded in string and other stuff
        $url =~ s/\s*[.,;'")]*\s*$//;
        check_url($url, $File::Find::name);
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
