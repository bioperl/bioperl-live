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
!NO!SUBS!

print OUT 'my $dir = "' . $Config{'installscript'} . '";';

print OUT <<'!NO!SUBS!';

foreach my $target ( keys ( %symlink_scripts ) ) {
    unlink "$dir/".$symlink_scripts{$target} if -e "$dir/".$symlink_scripts{$target};
    eval { symlink( "$dir/$target", "$dir/".$symlink_scripts{$target} ); 1} 
        or print STDERR "Cannot create symbolic link named $dir/"
            . $symlink_scripts{$target}
            . " on you system for $dir/$target\n";
}

!NO!SUBS!
close OUT or die "Can't close $file: $!";
chmod 0755, $file or die "Can't reset permissions for $file: $!\n";
exec("$Config{'eunicefix'} $file") if $Config{'eunicefix'} ne ':';
chdir $origdir;
