#!/usr/bin/perl

# $Id$
# This script will convert from SGD format to GFF format
# See http://genome-www4.stanford.edu/Saccharomyces/SGD/doc/db_specifications.html

=head1 NAME

sgd_to_gff.pl - Massage SGD's feature dump format into a form suitable for Bio::DB::GFF

=head2 SYNOPSIS

   perl sgd_to_gff.pl chromosomal_features.tab > sgd.gff

=head2 DESCRIPTION

This script massages the SGD yeast sequence feature file located at
ftp://genome-ftp.stanford.edu/pub/yeast/data_dump/feature/chromosomal_feature.tab
into a GFF format suitable for use with Bio::DB::GFF.  This lets you
view the yeast annotations with the generic genome browser
(http://www.gmod.org).

To use this script, get the SGD features file at the above URL.  Then
run this command:

  sgd_to_gff.pl chromosomal_feature.tab > sgd.gff

The resulting database will have the following feature types
(represented as "method:source"):

  Component:chromosome       A chromosome
  gene:sgd                   A named gene
  rRNA:sgd		     A ribosomal RNA
  ARS:sgd		     An origin of replication
  CEN:sgd		     Centromere
  snRNA:sgd		     Small nuclear RNA
  RNA:sgd		     An RNA gene
  ORF:sgd		     An open reading frame
  ORF|Pseudogene:sgd	     A probably pseudogene
  LTR:sgd		     A long terminal repeat
  Ty ORF:sgd		     ??
  Transposon:sgd	     A transposon
  Pseudogene|Ty ORF:sgd	     ??
  snoRNA:sgd		     Small nucleolar RNA
  tRNA:sgd		     Transfer RNA

=head2 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>

=cut

use strict;

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

Converts "chromosomal_features.tab" file from SGD into a GFF file
suitable for loading into Bio::DB::GFF.  You can get this file at
ftp://genome-ftp.stanford.edu/pub/yeast/data_dump/feature/chromosomal_feature.tab

The output of running this script is suitable for loading using
load_gff.pl or bulk_load_gff.pl.
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
    my @genes = ($gene,split/\|/,$aliases);
    foreach (@genes) {
      print join("\t",$ref,'sgd','gene',$start,$stop,'.',$strand,'.',qq(Gene "$_" ; Note "$description")),"\n";
    }
    $description = "$gene\\; $description";
  }

  print join("\t",$ref,'sgd',$type,$start,$stop,'.',$strand,'.',qq($type "$id" ; Note "$description")),"\n";
}
