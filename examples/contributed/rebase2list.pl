#!/usr/bin/perl
# Generate an enzyme list for RestrictionEnzyme.pm from rebase
# From Ryan Brinkman

my $strider = $ARGV[0]; #commercial_version_rebase_strider_format

open my $FILEIN, '<', $strider or die "Could not read file '$strider': $!\n";

while (<$FILEIN>){
   chomp;
   if ( /^[A-Z]\S+,\S+/ ){
      ($enzyme,$cutsite)=split(',');
      if ($cutsite =~ m#/#){
	 $match=$-[0];
      }
      ($seqfixed=$cutsite) =~ s/\///g;
      $seqfixed=uc $seqfixed;
      print " \'$enzyme\'\t=> \'".$seqfixed." ".$match."\'\,\n";
   }
}
close $FILEIN;
