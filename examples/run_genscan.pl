#!/usr/bin/perl -w
# Brian Osborne
# script to run genscan on all nucleotide sequences in a fasta file
# and save results as fasta, creates <file>.gs.pept and <file>.gs.cds

use Bio::SeqIO;
use Getopt::Long;
use Bio::Tools::Genscan;
use Carp;
use strict;

# directory with GENSCAN matrices
my $genscanDir = "/dbs/genscan";
# GENSCAN location
my $binDir = "/usr/local/bin";
# GENSCAN matrix
my $matrix = "HumanIso.smat";

my ($file,$in);

GetOptions( "f|file=s" => \$file );
usage() if ( !$file );

# create output files for predicted protein and CDS
open PEPT, ">>$file.gs.pept" or croak "Error opening $file.gs.pept: $!\n";
open DNA, ">>$file.gs.cds" or croak "Error opening $file.gs.cds: $!\n";

$in = Bio::SeqIO->new(-file => $file , -format => 'Fasta');

while ( my $seq = $in->next_seq() ) {
    croak "Input sequence is protein\n" if ( $seq->moltype eq 'protein' );
    my $str = $seq->seq;
    my $id = $seq->id;
    $id = $1 if ( $id =~ /gi\|(\d+)/ );
    # create temp file, input to GENSCAN
    open OUT,">$id.temp.fa" or croak "Error opening $id.temp.fa: $!\n";
    print OUT ">$id.fa\n$str\n\n";

    # the contents of the *raw file will be parsed
    system "genscan $genscanDir/$matrix $id.temp.fa -cds > $id.gs.raw";
    unlink "$id.temp.fa";
    my $genscan = Bio::Tools::Genscan->new( -file => "$id.gs.raw");
    while ( my $gene = $genscan->next_prediction() ) {
	my $prt = $gene->predicted_protein;
	my $cds = $gene->predicted_cds;

	if ( defined $cds  ) {
	    $cds->display_id =~ /predicted_(CDS_\d+)/;
	    print DNA ">" . $id . "_" . $1 . " " . $cds->display_id . "\n"
	      . $cds->seq . "\n";
	}
	if ( defined $prt ) {
	    $prt->display_id =~ /predicted_(peptide_\d+)/;
	    print PEPT ">" . $id . "_" . $1 . " " . $prt->display_id . "\n"
	      . $prt->seq . "\n";
	}
    }
    $genscan->close();
    unlink "$id.gs.raw";
}


sub usage {
    print "
Usage    : $0 -f <file>
Function : run genscan on all nucleotide sequences in a multiple fasta file
Output   : <file>.gs.pept and <file>.gs.cds

";
    exit;
}
