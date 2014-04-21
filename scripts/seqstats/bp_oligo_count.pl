#!perl
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
use Getopt::Long;

#########################
# VARIABLES & FILENAMES #
#########################

use strict;
use warnings;

my ($format, $infile, $help, $outfile, $oligomerlength) = ('fasta');
GetOptions(
           'f|format:s'            => \$format,
           'i|in|s|sequence:s'     => \$infile,
           'h|help|?'              => \$help,
           'o|out:s'               => \$outfile,
           'length:i'              => \$oligomerlength
          );

my $USAGE = "Usage:\toligo_count [-h/--help] [-l/--length OLIGOLENGTH]\n".
    "\t[-f/--format SEQFORMAT] [-i/--in/-s/--sequence SEQFILE]\n".
    "\t[-o/--out OUTFILE]\n".
    "\tDefault SEQFORMAT is fasta\n";

print $USAGE and exit if $help;

unless ($infile ) {
    print 'Enter your concatenated FASTA sequence filename: ';
    chomp ($infile=<STDIN>);
}
unless (-e $infile) { die "$infile not found\n"; }

if ($outfile) {
    if (-e $outfile) {
	print "$outfile already exists!  Overwrite (Y/N)? ";
	chomp ($_ = <STDIN>);
	while (/[^yn]/i) {
	    print 'Y or N, please: ';
	    chomp ($_ = <STDIN>);
	}
	if (/n/i) { die "$outfile not overwritten.\n"; }
    }
#} else {
#    print 'Enter an output filename: ';
#    chomp ($outfile=<STDIN>);
#    if (-e $outfile) {
#	print "$outfile already exists!  Overwrite (Y/N)? ";
#	chomp ($_ = <STDIN>);
#	while (/[^yn]/i) {
#	    print 'Y or N, please: ';
#	    chomp ($_ = <STDIN>);
#	}
#	if (/n/i) { die "$outfile not overwritten.\n"; }
#    }
}

unless ($oligomerlength) {
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
my @oligoseqs = &generate_all_oligos($oligomerlength);
my %oligos = ();
foreach  (@oligoseqs) {
    $oligos{$_} = 0;
}

my $in = Bio::SeqIO->new( -file => $infile,
		       -format => $format);
my $seqnumber = 0;
my $oligocounts = 0;
my $exception;
while (my $seq = $in->next_seq() ) {
    my $len = $seq->length();
    my $position = 1;
    if ($position+$oligomerlength > $len) {
	$exception = 2;
	next;
    }
    $seq = uc $seq->seq; #string
    $exception = 1 if $seq =~ /[^GATC]/;

    while ($position + $oligomerlength-1 <= $len) {
	$oligos{substr $seq, $position-1, $oligomerlength}++;
	$position++;
	if ($position%250000 == 0) {print "$position\n";}
    }
    $oligocounts += $position-1;
    $seqnumber++;
}

my $OUTFILE;
if ($outfile) {
    open $OUTFILE, '>', $outfile or die "Could not open file '$outfile': $!\n";
} else {
    open $OUTFILE, '>-'; # STDOUT
}
print $OUTFILE "$seqnumber sequences analyzed\n";
print $OUTFILE "$oligocounts total $oligomerlength-mers counted\n";
print $OUTFILE "$oligomerlength-mer\tNumber\tFrequency\n";
foreach my $key (sort keys %oligos) {
    print $OUTFILE "$key\t$oligos{$key}\t", $oligos{$key}/$oligocounts, "\n";
}
close $OUTFILE;

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
    my @newarray = qw{A C G T};
    my @bases = qw{A C G T};

    while ($iter < $oligolength) {
	my @oldarray = @newarray;
	@newarray = ();
	foreach my $oligoseq (@oldarray) {
	    foreach my $newbase (@bases) {
		push @newarray, $oligoseq . $newbase;
	    }
	}
	$iter++;
    }
    return @newarray;
}

# if you wanted to be notified about status of running
#my $EMAILADDRESS = undef;
#die("Must change script to a valid email addres for notification") 
#    unless( defined $EMAILADDRESS );

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

__END__

=head1 NAME

bp_oligo_count - oligo count and frequency

=head1 SYNOPSIS

  Usage:  bp_oligo_count [-h/--help] [-l/--length OLIGOLENGTH]
          [-f/--format SEQFORMAT] [-i/--in/-s/--sequence SEQFILE]
          [-o/--out OUTFILE]

=head1 DESCRIPTION

This scripts counts occurrence and frequency for all oligonucleotides
of given length.

It can be used to determine what primers are useful for
frequent priming of nucleic acid for random labeling.

Note that this script could be run by utilizing the compseq
program which is part of EMBOSS.

=head1 OPTIONS

The default sequence format is fasta. If no outfile is given, the
results will be printed to standard out. All other options can entered
interactively.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Charles C. Kim

Email cckim@stanford.edu

=head1 HISTORY

Written July 2, 2001

Submitted to bioperl scripts project 2001/08/06

E<gt>E<gt> 100 x speed optimization by Heikki Lehvaslaiho

=cut
