#!/usr/bin/perl -w
# Example that shows values returned by Bio::SearchIO::Blast.
# Note that some methods here will return objects or arrays, not text.
# For example, $hsp->get_aln will return a Bio::SimpleAlign object,
# not the alignment in a printable form.
# Brian Osborne

use strict;
use Bio::SearchIO;
use Bio::SimpleAlign;
use Bio::AlignIO;

my $file = shift or die "No BLAST input\n";
my $in = new Bio::SearchIO(-format => 'blast',
			   -file   => $file );
while ( my $result = $in->next_result ) {
   print "algorithm is " . $result->algorithm . "\n";
   print "algorithm_version is " . $result->algorithm_version . "\n";
   print "query_name is " . $result->query_name  . "\n";
   print "query_accession is " . $result->query_accession  . "\n";
   print "query_length is " . $result->query_length  . "\n";
   print "query_description is " . $result->query_description  . "\n";
   print "database_name is " . $result->database_name  . "\n";
   print "database_letters is " . $result->database_letters  . "\n";
   print "database_entries is " . $result->database_entries  . "\n";
   my @stats = $result->available_statistics;
   print "available_statistics is @stats" . "\n";
   my @params = $result->available_parameters;
   print "available_parameters is @params" . "\n";
   print "num_hits is " . $result->num_hits . "\n";
   print "hits is " . $result->hits . "\n\n";
   while ( my $hit = $result->next_hit ) {
      print "name is " . $hit->name  . "\n";
      print "accession is " . $hit->accession  . "\n";
      print "description is " . $hit->description  . "\n";
      print "length is " . $hit->length  . "\n";
      print "algorithm is " . $hit->algorithm  . "\n";
      print "raw_score is " . $hit->raw_score  . "\n";
      print "significance is " . $hit->significance  . "\n";
      print "bits is " . $hit->bits  . "\n"; 
      print "hsps is " . $hit->hsps  . "\n";
      print "num_hsps is " . $hit->num_hsps  . "\n";
      print "ambiguous_aln is " . $hit->ambiguous_aln  . "\n"; 
      print "overlap is " . $hit->overlap  . "\n"; 
      print "n is " . $hit->n  . "\n"; 
      print "logical_length is " . $hit->logical_length  . "\n";
      print "length_aln is " . $hit->length_aln  . "\n";
      print "gaps is " . $hit->gaps  . "\n";
      print "matches is " . $hit->matches  . "\n";
      print "frac_identical is " . $hit->frac_identical  . "\n";
      print "frac_conserved is " . $hit->frac_conserved  . "\n";
      print "frac_aligned_query is " . $hit->frac_aligned_query  . "\n";
      print "frac_aligned_hit is " . $hit->frac_aligned_hit  . "\n";
      print "frac_aligned_sbjct is " . $hit->frac_aligned_sbjct  . "\n";
      print "num_unaligned_sbjct is " . $hit->num_unaligned_sbjct  . "\n";
      print "num_unaligned_hit is " . $hit->num_unaligned_hit  . "\n";
      print "num_unaligned_query is " . $hit->num_unaligned_query  . "\n";
      print "seq_inds is " . $hit->seq_inds  . "\n";
      print "strand is " . $hit->strand  . "\n";
      print "frame is " . $hit->frame  . "\n"; 
      print "rank is " . $hit->rank . "\n";
      print "locus is " . $hit->locus  . "\n";
      print "each_accession_number is " . $hit->each_accession_number  . "\n";
      print "tiled_hsps is " . $hit->tiled_hsps  . "\n\n";
      while ( my $hsp = $hit->next_hsp ) {
	 print "algorithm is " . $hsp->algorithm . "\n";
	 print "evalue is " . $hsp->evalue  . "\n";
	 print "frac_identical is " . $hsp->frac_identical  . "\n";
	 print "frac_conserved is " . $hsp->frac_conserved  . "\n";
	 print "gaps is " . $hsp->gaps         . "\n";
	 print "query_string is " . $hsp->query_string . "\n";
	 print "hit_string is " . $hsp->hit_string . "\n";
	 print "homology_string is " . $hsp->homology_string . "\n";
	 print "length is " . $hsp->length  . "\n";
	 print "hsp_length is " . $hsp->hsp_length  . "\n";
	 print "frame is " . $hsp->frame  . "\n";
	 print "num_conserved is " . $hsp->num_conserved . "\n";
	 print "num_identical is " . $hsp->num_identical . "\n";
	 print "rank is " . $hsp->rank  . "\n";
	 print "seq_inds is " . $hsp->seq_inds . "\n";
	 print "score is " . $hsp->score  . "\n";
	 print "bits is " . $hsp->bits  . "\n";
	 print "range is " . $hsp->range  . "\n";
	 print "percent_identity is " . $hsp->percent_identity . "\n";
	 print "strand is " . $hsp->strand . "\n";
	 print "start is " . $hsp->start . "\n";
	 print "end is " . $hsp->end . "\n";
	 print "matches is " . $hsp->matches . "\n";
	 # use AlignIO to get alignment as text...
	 my $aln = $hsp->get_aln;
	 my $alnIO = Bio::AlignIO->new(-format=>"selex");
	 print "alignment is\n";
	 print $alnIO->write_aln($aln),"\n\n";
      }
   }
}
