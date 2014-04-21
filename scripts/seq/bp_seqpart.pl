#!/usr/bin/env perl
use strict;
use warnings;

=head1 NAME

I<bp_seqpart.pl> - Takes one or more sequence files and splits them into a number of load balanced files.

=head1 USAGE

 bp_seqpart.pl -n <NUM_PARTS> [-h, -p <PREFIX>, -f <FORMAT>, -o <OUT_DIR>] <FILES...>

   -n number of files to create through partitioning
   -h this help message
   -p prefix for all FASTA file names output, files are of the form <outdir>/<prefix>#.<format>
   -f format of the files, defaults to FASTA but you can specify anything supported by SeqIO from BioPerl
   -o output directory where to dump the split sequence files

=head1 DESCRIPTION

Script wrapping SeqIO that allows partitioning of multiple sequence files into near equal sized parts for later parallel processing. Even if you have 10 input files outputting to 10 files will balance the files to contain similar total length of sequence. ID's are ignored when deciding on how to balance each sequence.

=head1 AUTHOR

B<Matt Oates> - I<Matt.Oates@bristol.ac.uk>

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 EDIT HISTORY

2012-04-03 - Matt Oates
	First features added.
=cut

=head1 DEPENDANCY
B<Getopt::Long> Used to parse command line options.
B<Pod::Usage> Used for usage and help output.
B<Bio::SeqIO> Used to cut up sequences and parse FASTA.
=cut
use Getopt::Long;                     #Deal with command line options
use Pod::Usage;                       #Print a usage man page from the POD comments after __END__
use Bio::SeqIO;                       #Deal with sequence parsing, format and file IO

# Command Line Options
my $help;               #Same again but this time should we output the POD man page defined after __END__
my $prefix = 'part';    #Name each part
my $format = 'fasta';   #Sequence format we are using, default to fasta
my $outdir = '.';       #Use the current directory as default
my $num_splits;         #Number of files to split into
my @partitions;         #Details of each partition for the split

#Set command line flags and parameters.
GetOptions("help|h!" => \$help,
           "prefix|p=s" => \$prefix,
           "format|f=s" => \$format,
           "num-splits|n=i" => \$num_splits,
           "outdir|o=s" => \$outdir,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;

pod2usage(-exitstatus => 0, -verbose => 1, -msg => 'Please specify the number of split parts with -n <N>') 
    unless defined $num_splits;

#Setup a bunch of empty partitions including some SeqIO file handles to write to
@partitions = map { 
    $_ = { length => 0, 
           size => 0,
           file => Bio::SeqIO->new(
                    -file   => ">$outdir/$prefix$_.$format",
                    -format => $format,
               )
         } 
    } 1..$num_splits;

#Get sequences from all the files specified.
foreach my $file (@ARGV) {
    #Open each input file in turn for reading
    my $in  = Bio::SeqIO->new(
        -file => "<$file",
        -format => $format
    );
    #While there are still sequences to consume
    while ( my $seq = $in->next_seq() ) {
        #Sort the partitions on how full they are
        @partitions = sort {$a->{size} <=> $b->{size}} @partitions;
        #Add the length of the current seq to the smallest partition size
        my $length = $seq->length;
        $partitions[0]{size} += $length;
        #Increase the length of the partition
        $partitions[0]{length}++;
        #Write this sequence to the partitions file
        $partitions[0]{file}->write_seq($seq);
    }
}

#Report some basic statistics after the job
my $part = 1;
foreach my $partition (@partitions) {
    print STDERR "$outdir/$prefix$part.$format\n";
    print STDERR "\tSequence count = $partition->{length}\n";
    print STDERR "\tSequence characters = $partition->{size}\n";
    $part++;
}

1;
__END__

