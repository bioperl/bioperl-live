#!/usr/bin/perl
if ($ARGV[0]=~/^-?-h/ || @ARGV < 1) {
die <<'USAGE';

This script massages the RELEASE 3 Flybase/Gadfly GFF files located at
http://www.fruitfly.org/sequence/release3download.shtml into the
"correct" version of the GFF format.

To use this script, download the whole genome FASTA file and save it
to disk.  (The downloaded file will be called something like
"na_whole-genome_genomic_dmel_RELEASE3.FASTA", but the link on the
HTML page doesn't give the filename.)  Do the same for the whole
genome GFF annotation file (the saved file will be called something
like "whole-genome_annotation-feature-region_dmel_RELEASE3.GFF".)  If
you wish you can download the ZIP compressed versions of these files.

Next run this script on the two files, indicating the name of the
downloaded FASTA file first, followed by the gff file:

 % process_gadfly.pl na_whole-genome_genomic_dmel_RELEASE3.FASTA whole-genome_annotation-feature-region_dmel_RELEASE3.GFF > fly.gff

The gadfly.gff file and the fasta file can now be loaded into a Bio::DB::GFF database
using the following command:

  % bulk_load_gff.pl -d fly -fasta na_whole-genome_genomic_dmel_RELEASE3.FASTA fly.gff 

(Where "fly" is the name of the database.  Change it as appropriate.
The database must already exist and be writable by you!)

The resulting database will have the following feature types
(represented as "method:source"):

  Component:arm              A chromosome arm
  Component:scaffold	     A chromosome scaffold (accession #)
  Component:gap	             A gap in the assembly
  clone:clonelocator         A BAC clone
  gene:gadfly                A gene accession number
  transcript:gadfly          A transcript accession number
  translation:gadfly         A translation
  codon:gadfly               Significance unknown
  exon:gadfly                An exon
  symbol:gadfly              A classical gene symbol
  similarity:blastn          A BLASTN hit
  similarity:blastx          A BLASTX hit
  similarity:sim4            EST->genome using SIM4
  similarity:groupest        EST->genome using GROUPEST
  similarity:repeatmasker    A repeat

IMPORTANT NOTE: This script will *only* work with the RELEASE3 gadfly
files and will not work with earlier releases.

USAGE
;
}

use strict;
use warnings;

foreach (@ARGV) {
  $_ = "gunzip -c $_ |" if /\.gz$/;
}

if ($ARGV[0] =~ /fasta/i) {
  process_fasta();
} else {
  die "call as process_gadfly.pl \"release3_dna.FASTA\" \"release3_features.GFF\"";
}

while (<>) {
  next if /^\#/;
  chomp;
  my ($ref,$csource,$cmethod,$start,$stop,$cscore,$strand,$cphase,$cgroup) = split "\t";
  next if $start > $stop;  # something wrong. Don't bother fixing it.

  my $fixed_group = fix_group($csource,$cmethod,$cgroup);
  print join("\t",$ref,$csource,$cmethod,$start,$stop,$cscore,$strand,$cphase,$fixed_group),"\n";
  dump_symbol($ref,$csource,$cmethod,$start,$stop,$cscore,$strand,$cphase,$cgroup) if $cgroup =~ /symbol/i;
}

sub fix_group {
  my ($source,$method,$group) = @_;
  my (@group,$gene);
  push @group,"Transcript $1" if $group =~ /transgrp=([^; ]+)/;
  push @group,"Gene $1"       if $method eq 'gene' && $group =~ /genegrp=([^; ]+)/;

  $gene ||= qq(Note "FlyBase $1")  if $group =~ /dbxref=FlyBase:(\w+)/;
  $gene ||= qq(Note "GadFly $1")   if $group =~ /genegrp=([^; ]+)/;
  push @group,qq(Note "Symbol $1") if $group =~ /symbol=([^; ]+)/ && "Gene $1" ne $group[0];
  push @group,$gene;
  return join ' ; ',@group;
}

# called when we encounter a gene symbol
sub dump_symbol {
  my ($ref,$csource,$cmethod,$start,$stop,$cscore,$strand,$cphase,$cgroup) = @_;
  my ($symbol) = $cgroup=~/symbol=([^;]+)/;
  my ($gene)   = $cgroup=~/genegrp=([^;]+)/;
  return if $symbol eq $gene;
  $cmethod = 'symbol';
  print join("\t",$ref,$csource,$cmethod,$start,$stop,$cscore,$strand,$cphase,qq(Symbol "$symbol")),"\n";
}

sub process_fasta {
  my $file = shift @ARGV;
  open my $F, '<', $file or die "Could not read file '$file': $!\n";
  print STDERR "Reading big FASTA file, please be patient...\n";
  my ($current_id,%lengths);
  while (my $line = <$F>) {
    if ($line =~ /^>(\S+)/) {
      $current_id = $1;
      next;
    }
    die "this doesn't look like a fasta file to me" unless $current_id;
    chomp $line;
    $lengths{$current_id} += length $line;
  }
  close $F;

  foreach my $id (sort keys %lengths) {
    print join("\t", $id, 'arm', 'Component', 1, $lengths{$id},
                     '.', '+', '.', qq(Sequence "$id")
               ), "\n";
  }
}

__END__

=head1 NAME

bp_process_gadfly.pl - Massage Gadfly/FlyBase GFF files into a version suitable for the Generic Genome Browser

=head1 SYNOPSIS

  % bp_process_gadfly.pl ./RELEASE2 > gadfly.gff

=head1 DESCRIPTION

This script massages the RELEASE 3 Flybase/Gadfly GFF files located at
http://www.fruitfly.org/sequence/release3download.shtml into the "correct"
version of the GFF format.

To use this script, download the whole genome FASTA file and save it
to disk.  (The downloaded file will be called something like
"na_whole-genome_genomic_dmel_RELEASE3.FASTA", but the link on the
HTML page doesn't give the filename.)  Do the same for the whole
genome GFF annotation file (the saved file will be called something
like "whole-genome_annotation-feature-region_dmel_RELEASE3.GFF".)  If
you wish you can download the ZIP compressed versions of these files.

Next run this script on the two files, indicating the name of the
downloaded FASTA file first, followed by the gff file:

 % bp_process_gadfly.pl na_whole-genome_genomic_dmel_RELEASE3.FASTA whole-genome_annotation-feature-region_dmel_RELEASE3.GFF > fly.gff

The gadfly.gff file and the fasta file can now be loaded into a Bio::DB::GFF database
using the following command:

  % bulk_load_gff.pl -d fly -fasta na_whole-genome_genomic_dmel_RELEASE3.FASTA fly.gff 

(Where "fly" is the name of the database.  Change it as appropriate.
The database must already exist and be writable by you!)

The resulting database will have the following feature types
(represented as "method:source"):

  Component:arm              A chromosome arm
  Component:scaffold	     A chromosome scaffold (accession #)
  Component:gap	             A gap in the assembly
  clone:clonelocator         A BAC clone
  gene:gadfly                A gene accession number
  transcript:gadfly          A transcript accession number
  translation:gadfly         A translation
  codon:gadfly               Significance unknown
  exon:gadfly                An exon
  symbol:gadfly              A classical gene symbol
  similarity:blastn          A BLASTN hit
  similarity:blastx          A BLASTX hit
  similarity:sim4            EST->genome using SIM4
  similarity:groupest        EST->genome using GROUPEST
  similarity:repeatmasker    A repeat

IMPORTANT NOTE: This script will *only* work with the RELEASE3 gadfly
files and will not work with earlier releases.

=head1 SEE ALSO

L<Bio::DB::GFF>, L<bulk_load_gff.pl>, L<load_gff.pl>

=head1 AUTHOR

Lincoln Stein, lstein@cshl.org

Copyright (c) 2002 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
