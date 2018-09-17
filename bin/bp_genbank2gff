#!/usr/bin/perl

use lib '.';

use strict;
use warnings;
use Bio::DB::GFF;
use Getopt::Long;

=head1 NAME

bp_genbank2gff.pl - Load a Bio::DB::GFF database from GENBANK files.

=head1 SYNOPSIS

  % bp_genbank2gff.pl -d genbank -f localfile.gb
  % bp_genbank2gff.pl -d genbank --accession AP003256
  % bp_genbank2gff.pl --accession AP003256 --stdout

=head1 DESCRIPTION

This script loads a Bio::DB::GFF database with the features contained
in a either a local genbank file or an accession that is fetched from
genbank.  Various command-line options allow you to control which
database to load and whether to allow an existing database to be
overwritten.

The database must already have been created and the current user must
have appropriate INSERT and UPDATE privileges.  The --create option
will initialize a new database with the appropriate schema, deleting
any tables that were already there.

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options.
e.g. -d instead of --database.

   --create                 Force creation and initialization of database
   --dsn       <dsn>        Data source (default dbi:mysql:test)
   --user      <user>       Username for mysql authentication
   --pass      <password>   Password for mysql authentication
   --proxy     <proxy>      Proxy server to use for remote access
   --stdout                 direct output to STDOUT
   --adaptor   <adaptor>    adaptor to use (eg dbi::mysql, dbi::pg, dbi::oracle)   --viral                  the genome you are loading is viral (changes tag
                                 choices)
   --source    <source>     source field for features ['genbank']
    EITHER --file           Arguments that follow are Genbank/EMBL file names
    OR --gb_folder          What follows is a folder full of gb files to process    OR --accession          Arguments that follow are genbank accession numbers
                                 (not gi!)
    OR --acc_file           Accession numbers (not gi!) in a file (one per line,                                 no punc.)
    OR --acc_pipe           Accession numbers (not gi!) from a STDIN pipe (one
                                 per line)


=head1 SEE ALSO

L<Bio::DB::GFF>, L<bulk_load_gff.pl>, L<load_gff.pl>

=head1 AUTHOR

Scott Cain, cain@cshl.org

Copyright (c) 2003 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

package Bio::DB::GFF::Adaptor::biofetch_to_stdout;
use CGI 'escape';
use Bio::DB::GFF::Util::Rearrange;
use Bio::DB::GFF::Adaptor::biofetch;
use vars '@ISA';
@ISA = 'Bio::DB::GFF::Adaptor::biofetch';

sub load_gff_line {
  my ($self,$options) = @_;
  # synthesize GFF3-compatible line
  my @attributes;
  if (my $id = $options->{gname}) {
    my $parent = $id;
    $parent    =~ s/\..\d+$//                  if     $options->{method} =~ /^(mRNA|transcript|exon|gene)$/;
    push @attributes,"Parent=".escape($parent) if     $options->{method} =~ /^(variation|exon|CDS|transcript|mRNA|coding)$/;
    push @attributes,"ID=".escape($id)         unless $options->{method} =~ /^(exon|CDS)$/;
  }
  if (my $tstart = $options->{tstart}) {
    my $tstop    = $options->{tstop};
    my $target   = escape($options->{gname});
    push @attributes,"Target=$target+$tstart+$tstop";
  }
  my %a;
  if (my $attributes = $options->{attributes}) {
    for my $a (@$attributes) {
      my ($tag,$value) = @$a;
      push @{$a{escape($tag)}},escape($value);
    }
    for my $a (keys %a) {
       push @attributes,"$a=".join(',',@{$a{$a}});
    }
  }
  ${$options}{'score'} = "." unless ${$options}{'score'};
  ${$options}{'strand'} = "." unless ${$options}{'strand'};
  ${$options}{'phase'} = "." unless ${$options}{'phase'};
  my $last_column = join ';',@attributes;
  if ($options->{method} eq 'origin') {
     print "##sequence-region $options->{gname} $options->{start} $options->{stop}\n";
  }
  print join("\t",@{$options}{qw(ref source method start stop score strand phase)},$last_column),"\n";
}

sub load_sequence_string {
  my $self = shift;
  my ($acc,$seq)  = @_;
  return unless $seq;
  $seq =~ s/(.{1,60})/$1\n/g;
  print ">$acc\n\L$seq\U\n";
}

sub setup_load {
   my $self = shift;
   print "##gff-version 3\n";
}

sub finish_load { }

1;

package main;

my $USAGE = <<USAGE;

Usage: $0 [options] [<gff file 1> <gff file 2>] ...
Load a Bio::DB::GFF database from GFF files.

 Options:
   --create                 Force creation and initialization of database
   --dsn       <dsn>        Data source (default dbi:mysql:test)
   --user      <user>       Username for mysql authentication
   --pass      <password>   Password for mysql authentication
   --proxy     <proxy>      Proxy server to use for remote access
   --stdout                 direct output to STDOUT
   --adaptor   <adaptor>    adaptor to use (eg dbi::mysql, dbi::pg, dbi::oracle)
   --viral                  the genome you are loading is viral (changes tag
                                 choices)
   --source    <source>     source field for features ['genbank']
    EITHER --file           Arguments that follow are Genbank/EMBL file names
    OR --gb_folder          What follows is a folder full of gb files to process
    OR --accession          Arguments that follow are genbank accession numbers
                                 (not gi!)
    OR --acc_file           Accession numbers (not gi!) in a file (one per line,
                                 no punc.) 
    OR --acc_pipe           Accession numbers (not gi!) from a STDIN pipe (one
                                 per line)   


This script loads a Bio::DB::GFF database with the features contained
in a either a local genbank file or an accession that is fetched from
genbank.  Various command-line options allow you to control which
database to load and whether to allow an existing database to be
overwritten.

USAGE
;

my ($DSN,$ADAPTOR,$CREATE,$USER,$VIRAL,$PASSWORD,$gbFOLDER,
    $FASTA,$ACC,$accFILE, $accPIPE, $FILE,$PROXY,$STDOUT,$SOURCE);


GetOptions (
            'dsn:s'       => \$DSN,
            'user:s'      => \$USER,
            'password:s'  => \$PASSWORD,
            'adaptor:s'   => \$ADAPTOR,
            'accession'   => \$ACC,
            'file'        => \$FILE,
            'viral'       => \$VIRAL,
            'acc_file'    => \$accFILE,
            'acc_pipe'    => \$accPIPE,
	    'source:s'    => \$SOURCE,
            'gb_folder=s' => \$gbFOLDER,
            'proxy:s'     => \$PROXY,
            'stdout'      => \$STDOUT,
            'create'      => \$CREATE) or die $USAGE;


die $USAGE unless ($DSN || $STDOUT);  # at a minimum we need to have a place to write to!

# some local defaults
$DSN     ||= 'dbi:mysql:test';
$ADAPTOR ||= $STDOUT ? 'memory' : 'dbi::mysql';

# Ensure that biofetch inherits from the "right" adaptor.
# This is a horrible hack and should be fixed.
eval "use Bio::DB::GFF::Adaptor::${ADAPTOR}";
local @Bio::DB::GFF::Adaptor::biofetch::ISA = "Bio::DB::GFF::Adaptor::${ADAPTOR}";

my $biofetch = $STDOUT ? 'biofetch_to_stdout' : 'biofetch';
my @dsn      = $STDOUT ? () : (-dsn => $DSN);

my @auth;
push @auth,(-user=>$USER)     if defined $USER;
push @auth,(-pass=>$PASSWORD) if defined $PASSWORD;
push @auth,(-proxy=>$PROXY)   if defined $PROXY;

my %preferred_tags = (
		      strain        => 10,
		      organism      => 20,
		      protein_id    => 40,
		      locus_tag     => 50,
		      locus         => 60,
		      gene          => 70,
		      standard_name => 80,
                     );
$preferred_tags{'product'} = 90 if $VIRAL; # added this to the default list for viral genomes
       # since most functions come from post-translational processing, so the default labels are c**p!

my $db = Bio::DB::GFF->new(-adaptor=>$biofetch,
			   @dsn,
			   @auth,
			   -preferred_tags => \%preferred_tags,
			   -source=> $SOURCE || 'Genbank')
  or die "Can't open database: ",Bio::DB::GFF->error,"\n";

if ($CREATE) {
  $db->initialize(1);
}

die "you must specify either an accession to retrieve from\nembl or a local file containing data in embl format\n" if (($FILE || $ACC) && !scalar(@ARGV));

if ($ACC) {
  while ($_ = shift) {
    status(loading => $_);
    my $result = $db->load_from_embl(/^NC_/?'refseq':'embl' => $_);
    status(done    => $result);
  }
  exit 1;
}

elsif ($FILE) {
  while ($_ = shift) {
    status('loading' => $_);
    my $result = $db->load_from_file($_);
    status (done => $result);
  }
  exit 1;
}

elsif ($accFILE){
    my $filename = shift;
    die "you must supply a filename after the --accFILE command line flag\n" unless $filename;
    die "file $filename does not exist\n" unless (-e $filename && !(-d $filename));
    open my $IN, '<', $filename or die "Could not read file '$filename' for reading accession numbers: $!\n";
    while (my $line = <$IN>){
        chomp $line;
	status(loading => $line);
        my $result = $db->load_from_embl(/^NC_/?'refseq':'embl' => $line);
	status(done => $result);
    }
    close $IN;
    exit 1;
}

elsif ($gbFOLDER){
    my $dir = $gbFOLDER;
    die "folder $dir does not exist\n" unless (-e $dir && -d $dir);
    opendir DIR, "$dir" || die "can't open directory $dir for reading: $!\n";
    my @files = readdir DIR;
    foreach my $file(@files){
        if (!(-e "$gbFOLDER/$file") || (-d "$gbFOLDER/$file")){
            print STDERR " $gbFOLDER/$file is not a filename!  Skipping...\n";
            next
        }
        my $result = $db->load_from_file("$gbFOLDER/$file");
        print STDERR $result ? "ok\n" : "failed\n";        
    }
} elsif ($accPIPE){
    my @accessions = <STDIN>;
    chomp @accessions;
    foreach (@accessions){
      status(loading => $_);
      my $result = $db->load_from_embl(/^NC_/?'refseq':'embl' => $_);
      status(done => $result);
    }
    exit 1;
}

else {
  my $done;
  while ($_ = shift) {
    $done = 1;
    status(loading => $_);
    my $result = $db->load_from_file($_);
    status(done => $result);
  }

  $done || die "\n\nno source of data provided\n\n";
  exit 1;
}

sub status {
  my ($state,$msg) = @_;
  return if $STDOUT;
  if ($state eq 'loading') {
    print STDERR "Loading $msg...";
  } elsif ($state eq 'done') {
    print STDERR $msg ? "ok\n" : "failed\n";
  }
}
