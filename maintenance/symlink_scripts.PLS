#!perl
use Config;
use File::Basename qw(&basename &dirname);
use Cwd;

$origdir = cwd;
chdir dirname($0);
$file = basename( $0, '.PL', '.PLS' );
$file .= $^O eq 'VMS' ? '.com' : '.pl';

open OUT, ">$file" or die "Can't create $file: $!";

print "Extracting $file (with variable substitutions)\n";

print OUT "$Config{startperl}\n";

print OUT <<'!NO!SUBS!';
use strict;

my %symlink_scripts = ('bp_bulk_load_gff.pl' => 'bp_pg_bulk_load_gff.pl');

my $dir = $ARGV[0];

foreach my $target ( keys ( %symlink_scripts ) ) {
    unlink "$dir/".$symlink_scripts{$target} if -e "$dir/".$symlink_scripts{$target};
    # place symlink in eval to catch error on systems that don't allow symlinks
    eval { symlink( "$dir/$target", "$dir/".$symlink_scripts{$target} ); 1}
        or print STDERR "Cannot create symbolic link named $dir/"
            . $symlink_scripts{$target}
            . " on your system for $dir/$target\n";
}

!NO!SUBS!
close OUT or die "Can't close $file: $!";
chmod 0755, $file or die "Can't reset permissions for $file: $!\n";
exec("$Config{'eunicefix'} $file") if $Config{'eunicefix'} ne ':';
chdir $origdir;
__END__

=head1 NAME

symlink_scripts.pl - install script to create symbolic links

=head1 SYNOPSIS

  perl symlink_scripts.pl

=head1 DESCRIPTION

Used during "make install".

This script will create a symlink in the 'installscript' directory (as defined
during install) to a script in that same directory.  It was written to create a
symlink with the name 'bp_pg_bulk_load_gff.pl' that targeted
'bp_bulk_load_gff.pl' but can be extended by adding files to the
%symlink_scripts hash.

Perl function 'symlink' is used to keep the script from crashing on systems
that don't allow symbolic linking.

=head1 SEE ALSO

=cut

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

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Ben Faga

Email faga@cshl.edu

=cut
