#!/usr/bin/perl
# Author: cckim@stanford.edu

# Description: This program designs primers for constructing knockouts
# of genes by transformation of PCR products (ref: Datsenko & Wanner,
# PNAS 2000).  A tab-delimited file containing ORF START STOP is read,
# and primers flanking the start & stop coordinates are designed based
# on the user-designated sequence file.  In addition, primers flanking
# the knockout regions are chosen for PCR screening purposes once the
# knockout is generated.  The script uses Bioperl in order to
# determine the primer sequences, which requires getting subsequences
# and reverse complementing some of the objects.

# make_primers.pl
# Purpose: Design primers for the Wanner method of PCR product-based knockouts
# Input: FASTA sequence file, tab-delimited coordinates file
# Output: Primer output file
# July 4, 2001
# Charles C. Kim

###########
# MODULES #
###########
use Bio::Seq;
use Getopt::Std;

#############
# VARIABLES #
#############
$upgap = 0; # the number of nt upstream of the 5' end to include in the deletion
$downgap = 0; # the number of nucleotides downstream of the 3' end to include
              # in the deletion
$oligolength = 40; # the length of the homologous region on each primer
$seqfile = '';   # don't specify these filenames unless you want to run
$coordfile = ''; # the program on these filenames exclusively
$outfile = '';   #
%fiveprime_primers = (
		      "P1" => "GTGTAGGCTGGAGCTGCTTC",
		      );
%threeprime_primers = (
		       "P2" => "CATATGAATATCCTCCTTAG",
		       "P4" => "ATTCCGGGGATCCGTCGACC",
		       );

#########
# FILES #
#########
getopts('s:c:o:');  # sequence file, coordinates file, output file

$seqfile = $opt_s if $opt_s;
$coordfile = $opt_c if $opt_c;
$outfile = $opt_o if $opt_o;

&open_readfile(*SEQFILE, 'sequence', $seqfile);
&open_readfile(*COORDFILE, 'coordinate', $coordfile);
&open_writefile(*PRIMERFILE, 'output', $outfile);

########
# MAIN #
########

$seq = '';
$count = 0;
while (<SEQFILE>) {
    if (/>/) {
	$count++;
	if ($count > 1) {
	    die "More than one sequence present in the input file\n";
	}
	next;
    }
    chomp($_);
    $_ =~ tr/gatc/GATC/;
    $seq .= $_;
}
close SEQFILE;

$seq = Bio::Seq-> new('-seq'=>$seq );

while (<COORDFILE>) {
    chomp($_);
    next if !$_;
    (my $name, my $start, my $stop) = split(/\t/, $_);
    if ($start < $stop) {
	$upprimer = $seq->subseq($start-$oligolength-$upgap, $start-1-$upgap);
	$downprimer = $seq->subseq($stop+1+$downgap,$stop+$oligolength+$downgap);
	$downprimer = Bio::Seq->new('-seq'=>$downprimer);
	$downprimer = $downprimer->revcom();
	$downprimer = $downprimer->seq();
	$uppcr = $seq->subseq($start-$oligolength-$upgap-20,$start-1-$upgap-$oligolength);
	$downpcr = $seq->subseq($stop+1+$downgap+$oligolength,$stop+$oligolength+$downgap+20);
	$downpcr = Bio::Seq->new('-seq'=>$downpcr);
	$downpcr = $downpcr->revcom();
	$downpcr = $downpcr->seq();
    }
    elsif ($start > $stop) {
	$upprimer = $seq->subseq($start+$upgap+1,$start+$oligolength+$upgap);
	$downprimer = $seq->subseq($stop-$oligolength-$downgap, $stop-1-$downgap);
	$upprimer = Bio::Seq->new('-seq'=>$upprimer);
	$upprimer = $upprimer->revcom();
	$upprimer = $upprimer->seq();
	$uppcr = $seq->subseq($start+$oligolength+$upgap+1,$start+$oligolength+$upgap+20);
	$downpcr = $seq->subseq($stop-$oligolength-$downgap-20,$stop-1-$downgap-$oligolength);
	$uppcr = Bio::Seq->new('-seq'=>$uppcr);
	$uppcr = $uppcr->revcom();
	$uppcr = $uppcr->seq();
    }
    else { die "Problem with start and stop coordinates\n"; }
    print PRIMERFILE "$name\n";
    print PRIMERFILE "5'pcr\t$uppcr\n";
    print PRIMERFILE "3'pcr\t$downpcr\n";
    print PRIMERFILE "\tExpected wildtype product size: ",abs($start-$stop)+121," bp\n";
    foreach $entry (sort keys %fiveprime_primers) {
	print PRIMERFILE "5'+$entry\t$upprimer$fiveprime_primers{$entry}\n";
    }
    foreach $entry (sort keys %threeprime_primers) {
	print PRIMERFILE "3'+$entry\t$downprimer$threeprime_primers{$entry}\n";
    }
    print PRIMERFILE "\n";
    $upprimer = '';
    $downprimer = '';
    $uppcr = '';
    $downpcr = '';
}


###############
# SUBROUTINES #
###############

sub open_readfile {
    my $filehandle = $_[0];
    my $filetype = $_[1] if $_[1];
    my $filename = $_[2] if $_[2];
    unless ($filename) {
	print "Enter $filetype filename: ";
	chomp ($filename=<STDIN>);
    }
    unless (-e $filename) { die "$filename not found\n"; }
    open $filehandle,'<', $filename or die "Could not read file '$filename': $!\n";
    $filehandle = '';
    $filetype = '';
    $filename = '';
}

sub open_writefile {
    my $filehandle = $_[0];
    my $filetype = $_[1] if $_[1];
    my $filename = $_[2] if $_[2];
    unless ($filename) {
	print "Enter $filetype filename: ";
	chomp ($filename=<STDIN>);
    }
    if (-e $filename) {
	print "$filename already exists!  Overwrite (Y/N)? ";
	chomp ($_ = <STDIN>);
	while (/[^yn]/i) {
	    print 'Y or N, please: ';
	    chomp ($_ = <STDIN>);
	}
	if (/n/i) { die "$filename not overwritten.\n"; }
	else { open $filehandle, '>', $filename or die "Could nott write file '$filename':·$!\n"; }
    }
    else { open $filehandle, '>', $filename or die "Could not write file '$filename': $!\n"; }
    $filehandle = '';
    $filetype = '';
    $filename = '';
}
