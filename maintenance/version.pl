#!/usr/bin/perl
#
=head1 version

This script is to add or modify version declaration for each bioperl pm.

[Currently, it just add version. Later I will update it to modify version.]

=head1 USAGE

  perl version.pl <module directory> <version>

=cut

use strict;

if(@ARGV < 2) {
    die "USAGE: perl version.pl <module directory> <version>\n";
}
my $dir=shift || "$ENV{HOME}/src/bioperl-live/";
my $version=shift || '1.4';

sub traveral_dir {
    my ($dir, )=@_;
    opendir DIR, $dir;
    my @allfiles= grep{$_ ne '.' and $_ ne '..'}readdir DIR;
    closedir DIR;
    my @full_path = map{"$dir/$_"} @allfiles;
    my @out = grep -f, @full_path;
    foreach(grep -d, @full_path){
        push @out, traveral_dir($_);
    }
    return @out;
}

my @pm=sort grep /\.pm$/, traveral_dir($dir);

use ExtUtils::MakeMaker;

map {
    my $f=$_; 
    my $v = MM->parse_version($f);
    print "$v\t$f\n";
    my $ep ='s/^(package\s+[\w:]+;\r?)$/$1\nour \$VERSION="'. $version.'";/';

    if((not defined $v) or $v eq 'undef'){ # This is strange on parse_version. 
    # It return scalar 'undef', not the undef can be detected by defined.
        `perl -p -i -e '$ep' $f`;
    }

} @pm;

