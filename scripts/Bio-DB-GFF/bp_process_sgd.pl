#!/usr/bin/perl

# This script will convert from SGD format to GFF format
# See http://db.yeastgenome.org/schema/Schema.html

use strict;
use warnings;

# hard-coded length data that I couldn't get directly
my %CHROMOSOMES = (I => 230_203,
		   II => 813_139,
		   III => 316_613,
		   IV  => 1_531_929,
		   V   => 576_869,
		   VI => 270_148,
		   VII => 1_090_937,
		   VIII => 562_639,
		   IX => 439_885,
		   X => 745_444,
		   XI => 666_445,
		   XII => 1_078_173,
		   XIII => 924_430,
		   XIV => 784_328,
		   XV  => 1_091_284,
		   XVI => 948_061,
		   Mit => 85_779);
my @ROMAN = qw(I II III IV V VI VII VIII IX X
	       XI XII XIII XIV XV XVI Mit);

if ($ARGV[0] =~ /^--?h/) {
  die <<USAGE;
 Usage: $0 <SGD features file>

This script massages the SGD sequence annotation flat files located at
ftp://genome-ftp.stanford.edu/pub/yeast/data_dump/feature/chromosomal_features.tab
into a version of the GFF format suitable for display by the generic
genome browser.

To use this script, get the SGD chromosomal_features.tab file from the
FTP site listed above, and run the following command:

  % process_sgd.pl chromosomal_features.tab > yeast.gff

The yeast.gff file can then be loaded into a Bio::DB::GFF database
using the following command:

  % bulk_load_gff.pl -d <databasename> yeast.gff

USAGE
;
}

# first print out chromosomes
# We hard coded the lengths because they are not available in the features table.
for my $chrom (sort keys %CHROMOSOMES) {
  print join("\t",$chrom,'chromosome','Component',1,$CHROMOSOMES{$chrom},'.','.','.',qq(Sequence "$chrom")),"\n";
}

# this is hard because the SGD idea of a feature doesn't really map onto the GFF idea.
while (<>) {
  chomp;
  my($id,$gene,$aliases,$type,$chromosome,$start,$stop,$strand,$sgdid,$sgdid2,$description,$date) = split "\t";
  my $ref = $ROMAN[$chromosome-1];
  $description =~ s/"/\\"/g;
  $description =~ s/;/\\;/g;

  $strand = $strand eq 'W' ? '+' : '-';
  ($start,$stop) = ($stop,$start) if $strand eq '-';
  die "Strand logic is messed up" if $stop < $start;

  if ($gene) {
     my @aliases = split(/\|/,$aliases);
     my $aliases = join " ; ",map {qq(Alias "$_")} @aliases;
     my $group = qq(Gene "$gene" ; Note "$description");
     $group .= " ; $aliases" if $aliases;
     print join("\t",$ref,'sgd','gene',$start,$stop,'.',$strand,'.',$group),"\n";
     $description .= "\\; AKA @aliases" if @aliases;
  }

  print join("\t",$ref,'sgd',$type,$start,$stop,'.',$strand,'.',qq($type "$id" ; Note "$description")),"\n";
}

__END__

=head1 NAME

bp_process_sgd.pl - Massage SGD annotation flat files into a version suitable for the Generic Genome Browser

=head1 SYNOPSIS

  % bp_process_sgd.pl chromosomal_features.tab > yeast.gff

=head1 DESCRIPTION

This script massages the SGD sequence annotation flat files located at
ftp://genome-ftp.stanford.edu/pub/yeast/data_dump/feature/chromosomal_features.tab
into a version of the GFF format suitable for display by the generic
genome browser.

To use this script, get the SGD chromosomal_features.tab file from the
FTP site listed above, and run the following command:

  % bp_process_sgd.pl chromosomal_features.tab > yeast.gff

The yeast.gff file can then be loaded into a Bio::DB::GFF database
using the following command:

  % bulk_load_gff.pl -d <databasename> yeast.gff

=head1 SEE ALSO

L<Bio::DB::GFF>, L<bulk_load_gff.pl>, L<load_gff.pl>

=head1 AUTHOR

Lincoln Stein, lstein@cshl.org

Copyright (c) 2002 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut


