#!/usr/bin/perl


if ($ARGV[0]=~/^-?-h/i) {
die <<USAGE;

This script massages the Flybase/Gadfly GFF files located at
ftp://ftp.fruitfly.org/pub/genomic/gadfly/ into the "correct" version
of the GFF format.

To use this script, download the Gadfly GFF distribution archive which
are organized by chromosome arm (e.g. "RELEASE2GFF.2L.tar.gz").
Unpack them will yield a directory named after the release,
e.g. RELEASE2, containing a directory named after the chromosome arm.
Do this repeatedly in order to create a directory that contains each
of the chromosome arms, i.e.:

   RELEASE2/gff/X
   RELEASE2/gff/2L
   RELEASE2/gff/2R
   ...

Give the release directory as the argument to this script, and capture
the script's output to a file:

  % process_gadfly.pl ./RELEASE2 > fly.gff

The gadfly.gff file can then be loaded into a Bio::DB::GFF database
using the following command:

  % bulk_load_gff.pl -d fly fly.gff

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

NOTE: Loading the DNA:

To load the fly DNA, download the FASTA format file for the
corresponding release in chromosome arm format
(e.g. ftp://ftp.fruitfly.org/pub/genomic/fasta/na_arms.dros.RELEASE2.Z),
and uncompress the file with the "uncompress" command.  This file will
also need to be converted to fix the names of the chromosome arms.
Use the following Perl command to do this:

 % perl -pi -e 's/^>Chromosome_arm_(\S+)/>$1/' na_arms.dros.RELEASE2

Now use bulk_load_gff.pl with the -fasta option in order to load the
DNA as well as the GFF data:

  % bulk_load_gff.pl -d fly -fasta na_arms.dros.RELEASE2 fly.gff

USAGE
;
}

use strict;
my $dir = shift || './RELEASE2';
$dir   .= "/gff";
my $map = read_map($dir) || die;
for my $scaffold (sort {
                       $map->{$a}[0] cmp $map->{$b}[0]
			 ||
		       $map->{$a}[1] <=> $map->{$b}[1]
		     } keys %$map) {
  convert_scaffold($dir,$scaffold,$map->{$scaffold});
}


sub read_map {
  my $dir = shift;
  my ($segments) = <$dir/*/SEGMENTS*.gff>;
  my %arms;
  $segments or die "Can't find SEGMENTS file";
  open (F,$segments) or die "Can't open $segments: $!";
  my %position;
  while (<F>) {
    chomp;
    my ($ref,$source,$method,$start,$stop,undef,$strand,undef,$group) = split "\t";
    $group =~ /name=(\w+)/ or next;
    $position{$1}=[$ref,$start,$stop,$strand];
    $arms{$ref}{min} = $start if !defined($arms{$ref}{min}) || $arms{$ref}{min} > $start;
    $arms{$ref}{max} = $stop  if !defined($arms{$ref}{max}) || $arms{$ref}{max} < $stop;
  }
  for my $ref (keys %arms) {
    print join("\t",$ref,'arm','Component',$arms{$ref}{min},$arms{$ref}{max},'.','+','.',qq(Sequence "$ref")),"\n";
  }
  return \%position;
  close F;
}

sub convert_scaffold {
  my ($dir,$scaffold,$pos) = @_;
  my ($ref,$rstart,$rstop,$rstrand) = @$pos;
  my ($file) = <$dir/*/$scaffold*.gff>;
  unless ($file && -r $file) {
    warn "$scaffold: Can't find corresponding GFF file.  Skipping.\n";
    return;
  }

  print join("\t",$ref,'scaffold','Component',$rstart,$rstop,'.',$rstrand,'.',qq(Sequence "$scaffold")),"\n";
  open (F,$file) or die "Can't open $file: $!";

  my @resultset  = ();
  my $oldresultset = '';
  while (<F>) {
    next if /^#/;
    chomp;
    my ($cref,$csource,$cmethod,$cstart,$cstop,$cscore,$cstrand,$cphase,$cgroup) = split "\t";
    next if $cstart > $cstop;  # something wrong. Don't bother fixing it.

    if ($cgroup =~ /resultset_id=([^ ;]+)/) { # put aside to deal with later
      my $resultset = $1;
      if ($resultset ne $oldresultset && @resultset) {
	dump_resultset($scaffold,$pos,\@resultset);
	@resultset = ();
      }
      push @resultset,[$cref,$csource,$cmethod,$cstart,$cstop,$cscore,$cstrand,$cphase,$cgroup];
      $oldresultset = $resultset;
      next;
    }

    else {
      (my $scref = $cref) =~ s/\.\d+$//;  # get rid of version number, if there is one
      unless ($scref eq $scaffold) {
	warn "$scaffold: feature uses $cref as reference.  Skipping";
	next;
      }
    }

    my ($start,$stop,$strand) = fix_coordinates($rstart,$rstop,$rstrand,$cstart,$cstop,$cstrand);
    my $fixed_group = fix_group($csource,$cmethod,$cgroup);
    print join("\t",$ref,$csource,$cmethod,$start,$stop,$cscore,$strand,$cphase,$fixed_group),"\n";
    dump_symbol($ref,$csource,$cmethod,$start,$stop,$cscore,$strand,$cphase,$cgroup),"\n" if $cgroup =~ /symbol/i;
  }
  close F;
  dump_resultset($scaffold,$pos,\@resultset);
}

# called when a full resultset is found
sub dump_resultset {
  my ($scaffold,$pos,$results) = @_;
  return unless @$results;
  local $^W = 0;
  my ($rref,$rstart,$rstop,$rstrand) = @$pos;

  my $genomic;
  for my $d (@$results,[]) {

    if (index($d->[0],$scaffold) >= 0) {  # bit of the genomic sequence
      $genomic = $d;
      next;
    }

    # if we get here, we're in the target
    next unless $genomic; # don't know what to do with the thing
    my ($tref,$tsource,$tmethod,$tstart,$tstop,$tscore,$tstrand,$tphase,$tgroup) = @$d;
    my ($qref,$qsource,$qmethod,$qstart,$qstop,$qscore,$qstrand,$qphase,$qgroup) = @$genomic;

    # Comments in the reference field are a no-no
    $tref =~ s/\s+\(\d+ total\)$//;

    # While we're at it, might as well change the EST nomenclature so that the
    # 5', 3' pairs will automatically match in the viewer...
    $tstrand = '-' if $tref =~ /revcomp/;

    $tref =~ s/prime(_revcomp)?$//;

    my ($method,$group);
    ($tstart,$tstop) = ($tstop,$tstart) if $tstrand eq '-';

    if ($qmethod eq 'alignment') {
      $method = 'similarity';
      $group = qq(Target "Sequence:$tref" $tstart $tstop);
    }

    elsif ($qmethod eq 'HSP' && $qsource eq 'blastx') {
      $method = 'similarity';
      $group = qq(Target "Protein:$tref" $tstart $tstop);
    }

    elsif ($qmethod eq 'HSP' && $qsource eq 'blastn') {
      $method = 'similarity';
      $group  = qq(Target "Sequence:$tref" $tstart $tstop);
    }

    elsif ($qsource eq 'clonelocator') {
      $method  = 'clone';
      $group   = qq(Target "Clone:$tref" $tstart $tstop);
    }

    elsif ($qsource eq 'gap') {
      $method = 'Component';
      if ($qgroup =~ /span_id=(:\w+)/) {
	$group = "Gap $1";
      }
    }

    elsif (defined $tstart && defined $tstop) {
      $method = 'similarity';
      $group = qq(Target "Sequence:$tref" $tstart $tstop);
    }

    $method ||= $qmethod;
    next unless $group;

    my ($start,$stop,$strand) = fix_coordinates($rstart,$rstop,$rstrand,$qstart,$qstop,$qstrand);

    print join("\t",$rref,$qsource,$method,$start,$stop,$qscore,$qstrand,$qphase,$group),"\n";
    undef $genomic;
  }
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

sub fix_coordinates {
  my ($rstart,$rstop,$rstrand,$cstart,$cstop,$cstrand) = @_;  
  # Fix the coordinates.  We are going to accomodate (-) strand scaffolds, even though they
  # don't seem to occur
  if ($rstrand eq '+') {
    $cstart += ($rstart - 1);
    $cstop  += ($rstart - 1);
  } else {
    $cstart += ($rstop - $cstart + 1);
    $cstop  += ($rstop - $cstart + 1);
    $cstrand = $cstrand eq '+' ? '-' : '+';
  }
  return ($cstart,$cstop,$cstrand);
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

__END__

=head1 NAME

process_gadfly.pl - Massage Gadfly/FlyBase GFF files into a version suitable for the Generic Genome Browser

=head1 SYNOPSIS

  % process_gadfly.pl ./RELEASE2 > gadfly.gff

=head1 DESCRIPTION

This script massages the Flybase/Gadfly GFF files located at
ftp://ftp.fruitfly.org/pub/genomic/gadfly/ into the "correct" version
of the GFF format.

To use this script, get the Gadfly GFF distribution archive which is
organized by GenBank accession unit (e.g. "RELEASE2GFF.tar.gz").
Unpacking it will yield a directory named after the release,
e.g. RELEASE2.

Give that directory as the argument to this script, and capture the
script's output to a file:

  % process_gadfly.pl ./RELEASE2 > gadfly.gff

The gadfly.gff file can then be loaded into a Bio::DB::GFF database
using the following command:

  % bulk_load_gff.pl -d <databasename> gadfly.gff

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

=head1 SEE ALSO

L<Bio::DB::GFF>, L<bulk_load_gff.pl>, L<load_gff.pl>

=head1 AUTHOR

Lincoln Stein <lstein@cshl.org>.

Copyright (c) 2002 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut


