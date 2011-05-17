#!/usr/bin/perl
# prosite2perl -- convert Prosite patterns to Perl regular expressions
#
# Jordan Dimov (jdimov@cis.clarion.edu)
#
# Submitted to bioperl scripts project 2001/08/03 
#
# Description: 
# Prosite patterns to Perl regular expressions.
# The prositeRegEx($) sub accepts a string
# containing a Prosite pattern and returns a
# string containing a valid Perl regex.  The code
# is self-explanatory.

sub prositeRegEx($);

while (<>) {
  chomp ($_);
  print prositeRegEx ($_), "\n";
}

sub prositeRegEx ($) {
  my $regex = shift;
  $regex =~ s/[\-\.]//g;    
  $regex =~ s/\{/[^/g; 
  $regex =~ tr/x()<>}/.{}^$]/;
  return ($regex);
}
