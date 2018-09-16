#!/usr/bin/perl
#
################################################################################
#11-17-2001
#Jianwen Fang (jwfang1999@yahoo.com)
#
#THis program returns reverse complement sequences of all sequences in the current directory
#and save them in the same directory, using the same name with extension ".rev"
###############################################################################


use strict;
use Bio::Seq;
use Bio::SeqIO;

my @files = ();
my $folder = '.';
my $inputFormat;
my $outputFormat;
my $numSeq;

   #Fasta       FASTA format
   #EMBL        EMBL format
   #GenBank     GenBank format
   #GCG         GCG format
   #raw         Raw format (one sequence per line, no ID)

my @format = ('Fasta', 'EMBL', 'GenBank', 'GCG', 'Raw');

print("\nWhat is the format of the original sequence files?\n");
print("type 0 for Fasta; 1 for EMBL; 2 for GenBank; 3 for GCG; 4 for Raw\n");
$inputFormat = <STDIN>;
chomp ($inputFormat);

print("\nWhat is the format of the reverse complement sequence files you want?\n");
print("type 0 for Fasta; 1 for EMBL; 2 for GenBank; 3 for GCG; 4 for Raw\n");
$outputFormat = <STDIN>;
chomp ($outputFormat);

unless(opendir(FOLDER, $folder))
{
	print "cannot open folder $folder!\n";
	exit;
}
	
@files = grep(!/^\.\.?$/, readdir(FOLDER));

foreach my $file (@files)
	{
	   if($file =~ /seq/i)
	    {
		    getRevcom($file);
	        $numSeq++;
	    }
	}
	
print "$numSeq reverse complement sequences have been saved in current directory\n";
exit;

############################################################################
#subroutine getRevcom take an backward sequence file name(should with .seq extension) as parameter
#return its revcom sequence using the same name with the extension replaced with rev
############################################################################
sub getRevcom
{
	my $seqFile = $_[0];
	my $in = Bio::SeqIO->new('-file'=>$seqFile, '-format'=>$format[$inputFormat]);
	my $seq = $in->next_seq();
	my $revcomSeq = $seq->revcom();
	my @outSeqFile = split (/\./, $seqFile);
	pop @outSeqFile;
	push(@outSeqFile, 'rev');
	my $outSeqFile = join('.', @outSeqFile);
	print "$outSeqFile\n";
	my $out = Bio::SeqIO->new('-file'=>">$outSeqFile", '-format'=>$format[$outputFormat]);
	$out->write_seq($revcomSeq);
}


