#!/usr/bin/perl -w

use lib '.';

# $Id: bp_genbank2gff.PLS,v 1.5 2003-07-30 15:17:06 mwilkinson Exp $
use strict;
use Bio::DB::GFF;
use Getopt::Long;

=head1 NAME

bp_genbank2gff.pl - Load a Bio::DB::GFF database from GENBANK files.

=head1 SYNOPSIS

  % bp_genbank2gff.pl -d genbank -f localfile.gb
  % bp_genbank2gff.pl -d genbank -a AP003256

=head1 DESCRIPTION

This script loads a Bio::DB::GFF database with the features contained
in a either a local genbank file or an accession that is fetched from
genbank.  Various command-line options allow you to control which
database to load and whether to allow an existing database to be
overwritten.

This script currently only uses MySQL, though it is a proof-of-
principle and could easily be extended to work with other RDMS
that are supported by GFF through adaptors.

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options.
e.g. -d instead of --database.

   --create                  Force creation and initialization of database
   --dsn       <dsn>         Data source (default dbi:mysql:test)
   --user      <user>        Username for mysql authentication
   --pass      <password>    Password for mysql authentication
   --stdout                  Write converted GFF file to stdout rather than loading
   --proxy     <proxy>       Proxy server to use for remote access
    EITHER --file            Arguments that follow are Genbank/EMBL file names
    OR --accession           Arguments that follow are genbank accession numbers
    OR --acc_file            accession numbers in a file (one per line, no punc.) 
    OR --acc_pipe            accession numbers from a STDIN pipe (one per line)   

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

sub new {
  my $class = shift;
  my $self  = bless {},$class;
  my ($proxy) = rearrange(['PROXY'],@_);
  if ($proxy) {
    my @args = ref($proxy) ? @$proxy : eval $proxy;
    $self->{_proxy} = \@args if @args;
  }
  $self;
}

sub load_gff_line {
  my ($self,$options) = @_;
  # synthesize GFF3-compatible line
  my @attributes;
  if (my $parent = $options->{gname}) {
     push @attributes,"Parent=".escape($parent) unless $options->{method} =~ /^(mRNA|gene)$/;
     push @attributes,"ID=".escape($parent);
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
   --create                  Force creation and initialization of database
   --dsn       <dsn>         Data source (default dbi:mysql:test)
   --user      <user>        Username for mysql authentication
   --pass      <password>    Password for mysql authentication
   --proxy     <proxy>       Proxy server to use for remote access
   --stdout                  direct output to STDOUT
    EITHER --file            Arguments that follow are Genbank/EMBL file names
    OR --gb_folder           What follows is a folder full of gb files to process
    OR --accession           Arguments that follow are genbank accession numbers (not gi!)
    OR --acc_file            accession numbers (not gi!) in a file (one per line, no punc.) 
    OR --acc_pipe            accession numbers (not gi!) from a STDIN pipe (one per line)   


This script loads a Bio::DB::GFF database with the features contained
in a either a local genbank file or an accession that is fetched from
genbank.  Various command-line options allow you to control which
database to load and whether to allow an existing database to be
overwritten.

This script currently only uses MySQL, though it is a proof-of-
principle and could easily be extended to work with other RDMS
that are supported by GFF through adaptors.

USAGE
;

my ($DSN,$ADAPTOR,$CREATE,$USER,$PASSWORD,$gbFOLDER,$FASTA,$ACC,$accFILE, $accPIPE, $FILE,$PROXY,$STDOUT);


GetOptions (
            'dsn:s'       => \$DSN,
            'user:s'      => \$USER,
            'password:s'  => \$PASSWORD,
            'adaptor:s'   => \$ADAPTOR,
            'accession'   => \$ACC,
            'file'        => \$FILE,
            'acc_file'    => \$accFILE,
            'acc_pipe'    => \$accPIPE,
            'gb_folder=s' => \$gbFOLDER,
            'proxy:s'     => \$PROXY,
            'stdout'      => \$STDOUT,
            'create'      => \$CREATE) or die $USAGE;


die $USAGE unless ($DSN || $STDOUT);  # at a minimum we need to have a place to write to!

# some local defaults
$DSN     ||= 'dbi:mysql:test';
#$ADAPTOR = $STDOUT ? 'biofetch_to_stdout' : 'biofetch';
$ADAPTOR = $STDOUT ? 'biofetch_to_stdout' : $ADAPTOR;

my @auth;
push @auth,(-user=>$USER)     if defined $USER;
push @auth,(-pass=>$PASSWORD) if defined $PASSWORD;
push @auth,(-proxy=>$PROXY)   if defined $PROXY;

my $db = Bio::DB::GFF->new(-adaptor=>$ADAPTOR,-dsn => $DSN,@auth)
  or die "Can't open database: ",Bio::DB::GFF->error,"\n";

if ($CREATE) {
  $db->initialize(1);
}

die "you must specify either an accession to retrieve from\nembl or a local file containing data in embl format\n" if (($FILE || $ACC) && !scalar(@ARGV));

if ($ACC) {
  while ($_ = shift) {
    print STDERR "Loading $_...";
    my $result = $db->load_from_embl(/^NC_/?'refseq':'embl' => $_);
    print STDERR $result ? "ok\n" : "failed\n";
  }
  exit 1;
} elsif ($FILE) {
  while ($_ = shift) {
    print STDERR "Loading $_...\n";
    my $result = $db->load_from_file($_);
    print STDERR $result ? "ok\n" : "failed\n";
  }
  exit 1;
} elsif ($accFILE){
    my $filename = shift;
    die "you must supply a filename after the --accFILE command line flag\n" unless $filename;
    die "file $filename does not exist\n" unless (-e $filename && !(-d $filename));
    open IN, "$filename" || die "Can't open file $filename for reading accession numbers: $!\n";
    while (<IN>){
        chomp;
        print STDERR "Loading $_...";   
        my $result = $db->load_from_embl(/^NC_/?'refseq':'embl' => $_);
        print STDERR $result ? "ok\n" : "failed\n";
    }
    exit 1;
} elsif ($gbFOLDER){
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
        print STDERR "Loading $_...";   
        my $result = $db->load_from_embl(/^NC_/?'refseq':'embl' => $_);
        print STDERR $result ? "ok\n" : "failed\n";        
    }
    exit 1;
} else {
    my $done;
    while ($_ = shift) {
        $done = 1;
        print STDERR "Loading $_...\n";
        my $result = $db->load_from_file($_);
        print STDERR $result ? "ok\n" : "failed\n";
    }
    $done || die "\n\nyou didn't provide any source of data you silly sausage!\ndying from boredom...Aaaaaaaaaahhhhhhh\n\n";
    exit 1;
}
