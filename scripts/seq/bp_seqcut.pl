#!/usr/bin/env perl
use strict;
use warnings;

=head1 NAME

I<bp_seqcut.pl>

=head1 USAGE

 bp_seqcut.pl [options -h,-s,-e,-f,-w] <FILES>...
 bp_seqcut.pl [options -h,-f,-w] s-e <FILES>...

   -h this help message
   -s which residue to start cutting on
   -e which residue to finish cutting on
   -f format of the files, defaults to FASTA but you can specify anything supported by SeqIO from BioPerl
   -w hard wrap width, this might not be supported depending on which format you are using

=head1 Description

A script to cut FASTA sequences with a given range `fastacut -s 1 -e 10 *.fasta` or `fastacut 1-10 *.fasta`.
This is just a convenience wrapper around the Bio::SeqIO module. Useful if you wish to trim out a section of an
alignment to build a profile of a specific region of sequence.

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

2010-11-22 - Matt Oates
	First features added.
=cut



# Includes
=head1 DEPENDANCY
B<Getopt::Long> Used to parse command line options.
B<Pod::Usage> Used for usage and help output.
B<Bio::SeqIO> Used to cut up sequences and parse FASTA.
=cut
use Getopt::Long;                     #Deal with command line options
use Pod::Usage;                       #Print a usage man page from the POD comments after __END__
use Bio::SeqIO;

# Command Line Options
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $format = "Fasta";
my $start;
my $end;
my $width = 72; #Default for Jalview output
my $outfile = '/dev/stdout';

#Set command line flags and parameters.
GetOptions("help|h!" => \$help,
           "start|s=s" => \$start,
           "format|f=s" => \$format,
           "end|e=s" => \$end,
           "width|w=s" => \$width,
           "outfile|o=s" => \$outfile,
        ) or die "Fatal Error: Problem parsing command-line ".$!;

#Get other command line arguments that weren't optional flags.
($start,$end) = split (/-/, shift) unless ($start and $end);
my @files = @ARGV;

#Print out some help if it was asked for or if no arguments were given.
pod2usage(-exitstatus => 0, -verbose => 2) if $help;

pod2usage(-exitstatus => 0, -verbose => 1, -msg => 'Please specify the sequence files you wish to cut.') 
    unless scalar @files;

pod2usage(-exitstatus => 0, -verbose => 1, -msg => 'Please specify the region you wish to cut -s 1 -e 10 or 1-10.') 
    unless defined $end;

my $out = Bio::SeqIO->newFh(-file => ">$outfile", -format => $format) or die "Couldn't open selected output sequence file.";

#Open and iterate over all sequence in all files
foreach my $file (@files) {
	my $in  = Bio::SeqIO->new(-file => $file, -format => $format);
	while ( my $seq = $in->next_seq() ) {
            #Alter the ID to be postfixed with '/s-e'
            $seq->display_id($seq->display_id."/$start-$end");
            #Edit the sequence we have cut out
	        my $sequence = $seq->subseq($start,$end);
	        $sequence =~ s/([^\n]{0,$width})/$1\n/gi;
	        chomp $sequence;
            $seq->seq($sequence);
            #Print the sequence back out
        	print $out $seq;
	}
}

1;
__END__

