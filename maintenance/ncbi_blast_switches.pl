#!/usr/bin/perl

# This script determines all the valid command line switches
# from the four main NCBI BLAST tools, and produces Perl code
# to put into Bio/Tools/Run/StandAloneBlast.pm
#
# Torsten Seemann
# 27 June 2006


my @exe = qw(blastall blastpgp rpsblast bl2seq);

for my $exe (@exe) {
  open(HELP, "$exe - |") or die $!;
  my @switch;
  while (<HELP>) {
    next unless m/^\s*-(\w)\s/;
    push @switch, $1;
  }
  close(HELP);
  print "\t\@",uc($exe),"_PARAMS = qw(", join(q{ }, sort @switch), ");\n";
}

