#!/usr/bin/perl -w

use strict;
use Config;
use File::Copy;
use IO::File;
use ExtUtils::MakeMaker;

# this script should be run from the bioperl top level
use constant SCRIPTS => './scripts';
use constant MODE    => 0555;         # -r-xr-xr-x

my $install_dir = shift || $Config{installscript};
my $prompt_mode = shift || 'n';

my $usage = "

usage:

  perl install_bioperl_scripts.pl  install_dir  prompt_mode

prompt_mode = 'a' (all) | 'i' (interactive) | 'n' (none)
Note:  this script must be run from the parent directory of /scripts

";

die $usage unless $install_dir;


$prompt_mode = 'all'  if $prompt_mode =~ /^a/i;
$prompt_mode = 'none' if $prompt_mode =~ /^n/i;
$prompt_mode = 'some' if $prompt_mode =~ /^i/i;
die $usage if $prompt_mode eq 'none';

print "\n** BioPerl Scripts Install** \n\n";

chdir SCRIPTS  or die "Can't chdir to ",SCRIPTS,": $!";
opendir(F,".") or die "Can't opendir ",SCRIPTS,": $!";
while (my $file_or_dir = readdir(F)) {
  next if $file_or_dir =~ /^\./;
  next if $file_or_dir eq 'CVS';
  next unless -d $file_or_dir;
  next unless prompt_to_install($file_or_dir);
  print "Installing scripts in $file_or_dir...\n";
  install_contents($file_or_dir,$install_dir);
}
closedir F;

sub prompt_to_install {
  my $f = shift;
  print "\n* Script Directory $f *\n";
  if (-e "$f/TAG" && (my $g = IO::File->new("$f/TAG"))) {
    print while <$g>;
  }
  return 1 if $prompt_mode eq 'all';
  my $result = prompt("Install scripts in $f? y/n",'n');
  return $result =~ /^[yY]/;
}

sub install_contents {
  my $dir  = shift;
  my $dest = shift;
  my $bangline = $Config{startperl};
  opendir (D,$dir) or die "Can't open $dir: $!\n";
  while (my $script = readdir(D)) {
    next unless $script =~ /\.PLS$/;
    my $in  = IO::File->new("$dir/$script")    or die "Can't open $dir/$script: $!";
    $script =~ s/\.PLS$/\.pl/;                   # change from .PLS to .pl
    $script =~ s/^/bp_/ unless $script =~ /^bp/; # add the "bp" prefix
    print "\tInstalling $script....\n";
    unlink "$dest/$script" if -e "$dest/$script";
    my $out = IO::File->new(">$dest/$script")  or die "Can't open $dest/$script: $!";
    my $doneit;
    while (<$in>) {
      next if $doneit;
      if (s/^\#\!\S+/$bangline/) {
	$doneit++;
      }
    } continue {
      print $out $_;
    }
    close $in;
    close $out;
    chmod MODE,"$dest/$script" or die "Can't change mode of $script to ",MODE,": $!";
  }
  closedir D;
}

1;

__END__

