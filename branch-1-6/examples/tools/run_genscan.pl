#!/usr/bin/perl
# Brian Osborne
# script to run genscan on all nucleotide sequences in a fasta file
# and save results as the fasta files <file>.gs.pept and <file>.gs.cds,
# and <file>.gs.exons

use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long;
use Bio::Tools::Genscan;
use strict;

# GENSCAN matrix
my $matrix = "/home/bosborne/src/genscan/HumanIso.smat";

my ($file,$i);

GetOptions( "f|file=s" => \$file );
usage() if ( !$file );

my $pept_out = Bio::SeqIO->new(-file   => ">$file.gs.pept",
			       -format => "fasta");
my $cds_out = Bio::SeqIO->new(-file   => ">$file.gs.cds",
			      -format => "fasta");
my $exons_out = Bio::SeqIO->new(-file   => ">$file.gs.exons",
				-format => "fasta");

my $in = Bio::SeqIO->new(-file => $file , -format => 'Fasta');

while ( my $seq = $in->next_seq() ) {
   die "Input sequence is protein\n" if ( $seq->alphabet eq 'protein' );

   # create temp file, input to GENSCAN
   my $temp_out = Bio::SeqIO->new(-file   => ">temp.fa",
				   -format => "fasta");
   $temp_out->write_seq($seq);

   my $file_id = $seq->display_id;
   $file_id =~ s/\|/-/g;

   system "genscan $matrix temp.fa -cds > $file_id.gs.raw";
   unlink "temp.fa";

   my $genscan = Bio::Tools::Genscan->new( -file => "$file_id.gs.raw");
   while ( my $gene = $genscan->next_prediction() ) {
      $i++;
      my $pept = $gene->predicted_protein;
      my $cds = $gene->predicted_cds;
      my @exon_arr = $gene->exons;

      if ( defined $cds  ) {
	 my $cds_seq = Bio::Seq->new(-seq => $cds->seq,
				     -display_id => $cds->display_id);
	 $cds_out->write_seq($cds_seq);
      }

      if ( defined $pept ) {
	 my $pept_seq = Bio::Seq->new(-seq => $pept->seq,
				      -display_id => $pept->display_id);
	 $pept_out->write_seq($pept_seq);
      }

      for my $exon (@exon_arr) {
	 my $desc = $exon->strand . " " . $exon->start . "-" . $exon->end .
	   " " . $exon->primary_tag . " " . "GENSCAN_predicted_$i";
	 my $exon_seq = Bio::Seq->new(-seq => $seq->subseq($exon->start,
							   $exon->end),
				      -display_id => $seq->display_id,
				      -desc => $desc );
	 $exons_out->write_seq($exon_seq);
      }
   }
   $genscan->close();
   unlink "$file_id.gs.raw";
}

sub usage {
    print "
Usage    : $0 -f <file>
Function : run genscan on all nucleotide sequences in a multiple fasta file
Output   : <file>.gs.pept, <file>.gs.cds, <file>.gs.exons

";
    exit;
}
