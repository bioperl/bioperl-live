#!/usr/bin/perl

use constant ACEDB => 'sace://aceserver.cshl.org:2005'; 
use strict;
use warnings;
use Ace;

my @framework = qw(mex-3 spe-15 lin-17 unc-11 dhc-1 unc-40 smg-5
		   unc-13 unc-29 eat-16 lin-11 spe-9 par-6 unc-59 unc-54 mab-9 lin-42
		   sri-71 smu-2 vab-1 bli-2 dpy-10 him-14 mig-5 unc-4 bli-1 sqt-1 rol-1
		   his-14 unc-52 unc-45 par-2 let-805 sel-8 mab-21 daf-4 sma-3 lin-39
		   unc-32 tax-4 ced-9 tra-1 nob-1 daf-1 ced-2 lin-1 unc-17 dpy-13 unc-5
		   smg-7 dif-1 lin-49 elt-1 daf-14 dpy-20 dpy-26 unc-30 tra-3 sup-24
		   rho-1 egl-8 unc-60 srh-36 apx-1 unc-62 let-418 dpy-11 let-413 sel-9
		   unc-42 egl-9 sma-1 sqt-3 odr-3 hda-1 unc-76 gcy-20 skr-5 par-4 unc-51
		   egl-17 lim-6 fox-1 fax-1 lon-2 unc-97 unc-6 unc-18 mec-10 sop-1 mab-18
		   sdc-2 odr-7 unc-9 unc-3 gas-1 ace-1);
my %framework = map {$_=>1} @framework;
my %framework_seen = ();

my $USAGE = <<USAGE;
This script massages the Wormbase GFF files located at
ftp://www.wormbase.org/pub/wormbase/GENE_DUMPS into a version of the
GFF format suitable for display by the generic genome browser.  It
mainly adds comments to the annotations and designates certain
well-spaced genetic loci as framework landmarks.

This script requires the AcePerl distribution, which is available on
CPAN (look for the "Ace" module).

To use this script, get the WormBase GFF files from the FTP site
listed above and place them in a directory.  It might be a good idea
to name the directory after the current release, such as WS61.  You do
not need to uncompress the files.

Then give that directory as the argument to this script and capture
the script's output to a file:

  % process_wormbase.pl ./WS61 > wormbase.gff

It may take a while before you see output from this script, since it
must first fetch gene and protein database from the remote AceDB
running at www.wormbase.org.
The wormbase.gff file can then be loaded into a Bio::DB::GFF database
using the following command:

  % bulk_load_gff.pl -d <databasename> wormbase.gff
USAGE
;
#'

die $USAGE if $ARGV[0]=~/^-?-h/i;

my $db = Ace->connect(-url=>ACEDB,
		      -query_timeout=>500) or die "Can't open ace database:",Ace->error;

if (-d $ARGV[0]) {
  @ARGV = <$ARGV[0]/*.gff.gz>;
}

@ARGV || die $USAGE;

foreach (@ARGV) { # GFF FILES
  $_ = "gunzip -c $_ |" if /\.gz$/;
}

my (%NOTES,%LOCUS,%GENBANK,%CONFIRMED,%ORFEOME);
get_confirmed($db,\%CONFIRMED);
get_genbank($db,\%GENBANK);
get_loci($db,\%LOCUS);
get_notes($db,\%NOTES);
get_orfeome($db,\%ORFEOME);

while (<>) {
  chomp;
  next if /^\#/;
  my ($ref,$source,$method,$start,$stop,$score,$strand,$phase,$group) = split /\t/;
  next if $source eq 'assembly_tag';  # don't want 'em, don't need 'em
  $ref    =~ s/^CHROMOSOME_//;
  $group  =~ s/CHROMOSOME_//;

  $source ='' if $source eq '*UNKNOWN*';

  if ($method eq 'Sequence' && ($source eq 'curated' || $source eq 'RNA') && $group =~ /Sequence "(\w+\.\d+[a-z]?)"/) {
    my @notes;
    push @notes,map { qq(Note "$_")        } @{$NOTES{$1}}     if $NOTES{$1};
    push @notes,map { qq(Note "$_")        } @{$LOCUS{$1}}     if $LOCUS{$1};
    push @notes,qq(Confirmed_by "$CONFIRMED{$1}")              if $CONFIRMED{$1};
    $group = join ' ; ',$group,@notes;
    if (my $loci = $LOCUS{$1}) {
      foreach (@$loci) {
        print join("\t",$ref,$source,'gene',$start,$stop,$score,$strand,$phase,"Locus $_"),"\n";
        print join("\t",$ref,'framework','gene',$start,$stop,$score,$strand,$phase,"Locus $_"),"\n" 
          if $framework{$_} && !$framework_seen{$_}++;
      }
    }
  }

  if ($method eq 'Sequence' && $source eq 'Genomic_canonical' && $group =~ /Sequence "(\w+)"/) {
    if (my $accession = $GENBANK{$1}) {
      $group .= qq( ; Note "Genbank $accession");
      print join("\t",$ref,'Genbank',$method,$start,$stop,$score,$strand,$phase,"Genbank \"$accession\""),"\n";
    }
  }

  if ($method eq 'reagent' && $source eq 'Orfeome_project' && $group =~ /PCR_product "([^\"]+)"/) {
    my $amp = $ORFEOME{$1};
    $group .= qq( ; Amplified $amp) if defined $amp;
  }

  # fix variant fields: Variant "T" => Note "T"
  $group =~ s/(?:Variant|Insert) "(\w+)"/Note "$1"/;

  # fix UTR fields
  if ($group =~ /UTR "([35])_UTR:(\S+)"/) {
    $method = 'UTR';
    $source = "$1_UTR";
    $group = qq(Sequence "$2");
  }

  print join("\t",$ref,$source,$method,$start,$stop,$score,$strand,$phase,$group),"\n";
}

sub get_loci {
  my ($db,$hash) = @_;  # hash keys are predicted gene names, values are one or more loci names
  my @genes = $db->fetch(-query=>'find Locus Genomic_sequence',-filltag=>'Genomic_sequence');
  foreach my $obj (@genes) {
    my @genomic = $obj->Genomic_sequence or next;
    foreach (@genomic) {
      push @{$hash->{$_}},$obj;
    }
  }
}

sub get_notes {
  my ($db,$hash) = @_;  # hash keys are predicted gene names, values are one or more brief identifications
  my @genes = $db->fetch(-query=>'find Sequence Brief_identification',-filltag=>'Brief_identification');
  foreach my $obj (@genes) {
    my @notes = $obj->Brief_identification or next;
    $hash->{$obj} = \@notes;
  }
}

sub get_genbank {
  my ($db,$hash) = @_;   # hash keys are cosmid names, values are genbank accessions (1 to 1)
  my @cosmids = $db->fetch(-query=>'find Genome_Sequence Database',-filltag=>'Database');
  for my $cosmid (@cosmids) {
    my ($database,undef,$accession) = $cosmid->Database(1)->row;
    next unless $accession;
    $hash->{$cosmid} = $accession;
  }
}

sub get_confirmed {
  my ($db,$hash) = @_;  # hash keys are predicted gene names, values are confirmation type
  my @confirmed = $db->fetch(-query=>'find Sequence Confirmed_by',-filltag=>'Confirmed_by');
  foreach my $obj (@confirmed) {
    my $confirmed_by = $obj->Confirmed_by || 'Unknown';
    $hash->{$obj} = $confirmed_by;
  }
}

sub get_orfeome {
  my ($db,$hash) = @_;
  my @mv_primers = $db->fetch(-query=>'find PCR_Product mv*',-filltag=>'Amplified');
  for my $obj (@mv_primers) {
    my $amplified = $obj->Amplified;
    $hash->{$obj} = $amplified;
  }
}

__END__

=head1 NAME

bp_process_wormbase.pl - Massage WormBase GFF files into a version suitable for the Generic Genome Browser

=head1 SYNOPSIS

  % bp_process_wormbase.pl ./WS61 > wormbase.gff

=head1 DESCRIPTION

This script massages the Wormbase GFF files located at
ftp://www.wormbase.org/pub/wormbase/GENE_DUMPS into a version of the
GFF format suitable for display by the generic genome browser.  It
mainly adds comments to the annotations and designates certain
well-spaced genetic loci as framework landmarks.

This script requires the AcePerl distribution, which is available on
CPAN (look for the "Ace" module).

To use this script, get the WormBase GFF files from the FTP site
listed above and place them in a directory.  It might be a good idea
to name the directory after the current release, such as WS61.  You do
not need to uncompress the files.

Then give that directory as the argument to this script and capture
the script's output to a file:

  % bp_process_wormbase.pl ./WS61 > wormbase.gff

It may take a while before you see output from this script, since it
must first fetch gene and protein database from the remote AceDB
running at www.wormbase.org.
The wormbase.gff file can then be loaded into a Bio::DB::GFF database
using the following command:

  % bulk_load_gff.pl -d <databasename> wormbase.gff

=head1 SEE ALSO

L<Bio::DB::GFF>, L<bulk_load_gff.pl>, L<load_gff.pl>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>

Copyright (c) 2002 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

