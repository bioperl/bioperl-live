#!/usr/bin/perl
use Module::Build;
use strict;
use warnings;

my $build = Module::Build->current;

my %symlink_scripts = ('bp_bulk_load_gff.pl' => 'bp_pg_bulk_load_gff.pl');

#my $blib_dir = File::Spec->catdir($build->blib, 'script');
# using blib prior to installation, post build, always 'works', but the
# installation process installs the symlink as the actual file, so we may as
# well have just done a copy

my $install_dir = $build->install_destination('script');
$build->log_info("Will try to install symlinks to $install_dir\n");
my $orig_dir = $build->cwd;
chdir($install_dir);

while (my ($source, $destination) = each %symlink_scripts) {
	if ($^O !~ /Win32/) {
		eval { symlink($source, $destination) };
		$build->log_warn("Cannot create symbolic link named $destination on your system for $source in $install_dir\n") if $@;
	} else {
		# Win32 perl does not implement symlink(), as it would not work on all filesystems.
		require File::Copy;
		eval { File::Copy::copy($source, $destination) };
		$build->log_warn("Cannot create copy of script named $destination on your system for $source in $install_dir\n") if $@;
	}
}

chdir($orig_dir);

exit;

__END__

=head1 NAME

symlink_script.pl - install script to create symbolic links

=head1 SYNOPSIS

  perl Build.pl
  ./Build install

=head1 DESCRIPTION

Used during "./Build install". Only works if the script installation directory
used during "perl Build.pl" matches that used for the actual installation during
"./Build install". So if you install to a special place, do

  perl Build.pl --install_base /home/me
  ./Build install

not

  perl Build.pl
  ./Build install --install_base /home/me

This script will create a symlink to a script in that same directory. It was
written to create a symlink with the name 'bp_pg_bulk_load_gff.pl' that targeted
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

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=cut
