#!/usr/bin/perl -w
use strict;

my $base       = '/home/websites/bioperl.org';
my $srcdir     = "$base/src/git";
my $deob_index = "$base/src/Deobfuscator/bin/deob_index.pl";

my @modules = qw(
    bioperl-corba-client
    bioperl-corba-server
    bioperl-db
    bioperl-dev
    bioperl-ext
    bioperl-gui
    bioperl-live
    bioperl-microarray
    bioperl-network
    bioperl-pedigree
    bioperl-pipeline
    bioperl-pise
    bioperl-run
);

chdir $srcdir;
for my $module (@modules) {
    system("/usr/bin/perl $deob_index -s $module $srcdir/$module/Bio $srcdir/$module");
}
exit();
