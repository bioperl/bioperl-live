#!/usr/bin/perl -w

# $Id$
#
# Author: Charles Kim <cckim@stanford.edu>
# Submitted to bioperl scripts project 2001/08/06

# Note that this script could be run much faster by
# utilizing the compseq program as part of EMBOSS (jason stajich)
#

# oligomer_freq.pl
# We use this to determine what primers are useful for frequent priming of 
# nucleic acid for random labeling
# Input: Sequence file, oligomer length
# Output: Tab-delimited text file of oligomer frequencies
# Written July 2, 2001
# Charles C. Kim

###########
# MODULES #
###########
use Bio::Seq;
use Bio::SeqIO;
use Getopt::Std;

#########################
# VARIABLES & FILENAMES #
#########################

getopts('s:o:l:'); # sequence file, output file, length

if ($opt_s) {
    $infile = $opt_s;
}
else {
    print 'Enter your concatenated FASTA sequence filename: ';
    chomp ($infile=<STDIN>);
}
unless (-e $infile) { die "$infile not found\n"; }

if ($opt_o) {
    $outfile = $opt_o;
    if (-e $outfile) {
	print "$outfile already exists!  Overwrite (Y/N)? ";
	chomp ($_ = <STDIN>);
	while (/[^yn]/i) {
	    print 'Y or N, please: ';
	    chomp ($_ = <STDIN>);
	}
	if (/n/i) { die "$outfile not overwritten.\n"; }
    }
}

else {
    print 'Enter an output filename: ';
    chomp ($outfile=<STDIN>);
    if (-e $outfile) {
	print "$outfile already exists!  Overwrite (Y/N)? ";
	chomp ($_ = <STDIN>);
	while (/[^yn]/i) {
	    print 'Y or N, please: ';
	    chomp ($_ = <STDIN>);
	}
	if (/n/i) { die "$outfile not overwritten.\n"; }
    }
}

if ($opt_l) {
    if ($opt_l !~ /\d/) { die "Specified length is non-numeric\n"; }
    $oligomerlength = $opt_l;
}
else {
    while () {
	print 'Enter an oligomer length to count: ';
	chomp($oligomerlength=<STDIN>);
	if ($oligomerlength !~ /\d/) {
	    print "Value is non-numeric!\n";
	}
	else {last;}
    }
}


########
# MAIN #
########

if ($oligomerlength >= 9) {
    print "An oligomer length of $oligomerlength will generate ";
    print 4 ** $oligomerlength, " combinations,\nwhich could cause ";
    print "an out of memory error.  Proceed? (y/n) ";
    chomp($_=<STDIN>);
    if (/y/i) { ; }
    else { die "Program terminated\n"; }
}
@oligoseqs = &generate_all_oligos($oligomerlength);
%oligos = ();
foreach $entry (@oligoseqs) {
    $oligos{$entry} = 0;
}

$in = Bio::SeqIO->new( -file => "$infile",
		       -format => 'Fasta');
$seqnumber = 0;
$oligocounts = 0;
while (my $seq = $in->next_seq() ) {
    my $len = $seq->length();
    my $position = 1;
    if ($position+$oligomerlength > $len) {
	$exception = 2;
	next;
    }
    while ($position + $oligomerlength-1 <= $len) {
	my $chunk = $seq->subseq($position, $position+$oligomerlength-1);
	my $chunkseq = Bio::Seq->new( -seq => $chunk );
	$oligoseq = $chunkseq->seq();
	$oligoseq =~ tr/gatc/GATC/;
	if ($oligoseq =~ /[^GATC]/) { $exception = 1;}
	$oligos{$oligoseq}++;
	$position++;
	if ($position%250000 == 0) {print "$position\n";}
    }
    $oligocounts += $position-1;
    $seqnumber++;
}

open(OUTFILE, ">$outfile") or die "Can't open $outfile\n";
print OUTFILE "$seqnumber sequences analyzed\n";
print OUTFILE "$oligocounts total $oligomerlength-mers counted\n";
print OUTFILE "$oligomerlength-mer\tNumber\tFrequency\n";
foreach $key (sort keys %oligos) {
    print OUTFILE "$key\t$oligos{$key}\t", $oligos{$key}/$oligocounts, "\n";
}

if ($exception) {
    if ($exception == 1) {
	print "Non-standard (non-GATC) bases were found in sequence\n";
    }
    if ($exception == 2) {
	print "Oligomer length greater than sequence length\n";
    }
}

#&notify();

###############
# SUBROUTINES #
###############

sub generate_all_oligos {
    my $oligolength = $_[0];
    my $iter = 1;
    my @newarray = (A,C,G,T);
    my @bases = (A,C,G,T);

    while ($iter < $oligolength) {
	my @oldarray = @newarray;
	@newarray = ();
	foreach my $oligoseq (@oldarray) {
	    foreach my $newbase (@bases) {
		my $newoligo = $oligoseq . $newbase;
		push @newarray, $newoligo;
	    }
	}
	$iter++;
    }   
    return @newarray;
}

# if you wanted to be notified about status of running
#my $EMAILADDRESS = undef;
#die("Must change script to a valid email addres for notification") unless( defined $EMAILADDRESS );

#sub notify {
#    $address = $EMAILADDRESS;
#    $address = $_[0] if $_[0];
#    open(SENDMAIL, "|/usr/lib/sendmail -oi -t") or die "Can't fork for sendmail: $!\n";
#    print SENDMAIL <<"EOF";
#From: Computer
#To: $address
#Subject: Program Finished
#	  
#EOF
#    close(SENDMAIL) or warn "sendmail didn't close nicely";
#}
