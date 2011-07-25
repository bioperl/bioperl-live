#!/usr/bin/perl
# Example that shows values returned by Bio::SearchIO::Blast.
# Note that some methods will return objects or arrays, not text.
# For example, $hsp->get_aln will return a Bio::SimpleAlign object,
# not the alignment in a printable form.
# This script was used to create the table in the SearchIO HOWTO,
# found at http://bioperl.open-bio.org/wiki/HOWTO:SearchIO
# Brian Osborne

use strict;
use Bio::SearchIO;
use Bio::SimpleAlign;
use Bio::AlignIO;

my $file = shift or die "Usage: $0 <BLAST-report-file>\n";
my $in = new Bio::SearchIO(-format => 'blast',
			   # comment out the next line to read STDIN
			   -file   => $file );

while ( my $result = $in->next_result ) {
   my @stats =  $result->available_statistics;
   my @params = $result->available_parameters;
   print "Result\tavailable_statistics\t@stats\n";
   print "Result\tavailable_parameters\t@params\n";
   print "Result\talgorithm\t" .         $result->algorithm . "\n";
   print "Result\talgorithm_version\t" . $result->algorithm_version . "\n";
   print "Result\tquery_name\t" .        $result->query_name . "\n";
   print "Result\tquery_accession\t" .   $result->query_accession . "\n";
   print "Result\tquery_length\t" .      $result->query_length . "\n";
   print "Result\tquery_description\t" . $result->query_description . "\n";
   print "Result\tdatabase_name\t" .     $result->database_name . "\n";
   print "Result\tdatabase_letters\t" .  $result->database_letters . "\n";
   print "Result\tdatabase_entries\t" .  $result->database_entries . "\n";
   print "Result\tnum_hits\t" .          $result->num_hits . "\n";
   print "Result\thits\t" .              $result->hits . "\n";
   while ( my $hit = $result->next_hit ) {
      my $id =         $hit->matches('id');
      my $cons =       $hit->matches('cons');
      my @accs =       $hit->each_accession_number;
      my @qidentical = $hit->seq_inds('query','identical');
      my @qconserved = $hit->seq_inds('query','conserved');
      my @hidentical = $hit->seq_inds('hit','identical');
      my @hconserved = $hit->seq_inds('hit','conserved');
      print "Hit\tseq_inds('query','identical')\t@qidentical\n";
      print "Hit\tseq_inds('query','conserved')\t@qconserved\n";
      print "Hit\tseq_inds('hit','identical')\t@hidentical\n";
      print "Hit\tseq_inds('hit','conserved')\t@hconserved\n";
      print "Hit\teach_accession_number\t@accs\n";
      print "Hit\tmatches('id')\t" . $id . "\n";
      print "Hit\tmatches('cons')\t" . $cons . "\n";
      print "Hit\tname\t" .                $hit->name . "\n";
      print "Hit\taccession\t" .           $hit->accession . "\n";
      print "Hit\tdescription\t" .         $hit->description . "\n";
      print "Hit\tlength\t" .              $hit->length . "\n";
      print "Hit\talgorithm\t" .           $hit->algorithm . "\n";
      print "Hit\traw_score\t" .           $hit->raw_score . "\n";
      print "Hit\tsignificance\t" .        $hit->significance . "\n";
      print "Hit\tbits\t" .                $hit->bits . "\n";
      print "Hit\thsps\t" .                $hit->hsps . "\n";
      print "Hit\tnum_hsps\t" .            $hit->num_hsps . "\n";
      print "Hit\tambiguous_aln\t" .       $hit->ambiguous_aln . "\n";
      print "Hit\toverlap\t" .             $hit->overlap . "\n";
      print "Hit\tn\t" .                   $hit->n . "\n"; 
      print "Hit\tlogical_length\t" .      $hit->logical_length . "\n";
      print "Hit\tlength_aln\t" .          $hit->length_aln . "\n";
      print "Hit\tgaps\t" .                $hit->gaps . "\n";
      print "Hit\tfrac_identical\t" .      $hit->frac_identical . "\n";
      print "Hit\tfrac_conserved\t" .      $hit->frac_conserved . "\n";
      print "Hit\tfrac_aligned_query\t" .  $hit->frac_aligned_query . "\n";
      print "Hit\tfrac_aligned_hit\t" .    $hit->frac_aligned_hit . "\n";
      print "Hit\tnum_unaligned_sbjct\t" . $hit->num_unaligned_sbjct . "\n";
      print "Hit\tnum_unaligned_hit\t" .   $hit->num_unaligned_hit . "\n";
      print "Hit\tnum_unaligned_query\t" . $hit->num_unaligned_query . "\n";
      print "Hit\tstrand\t" .              $hit->strand . "\n";
      print "Hit\tframe\t" .               $hit->frame . "\n";
      print "Hit\trank\t" .                $hit->rank . "\n";
      print "Hit\tlocus\t" .               $hit->locus . "\n";
      while ( my $hsp = $hit->next_hsp ) {
	 my ($qid,$qcons) = $hsp->matches('hit');
	 my ($id,$cons) = $hsp->matches('query');
	 @qidentical = $hsp->seq_inds('query','identical');
	 @qconserved = $hsp->seq_inds('query','conserved');
	 @hidentical = $hsp->seq_inds('hit','identical');
	 @hconserved = $hsp->seq_inds('hit','conserved');
	 my @hrange =  $hsp->range('hit');
	 my @qrange =  $hsp->range('query');
	 my $aln =     $hsp->get_aln;
	 my $alnIO = Bio::AlignIO->new(-format=>"clustalw");
	 print "HSP\trange('hit')\t@hrange\n";
	 print "HSP\trange('query')\t@qrange\n";
	 print "HSP\tseq_inds('hit','identical')\t@hidentical\n";
	 print "HSP\tseq_inds('query','conserved')\t@qconserved\n";
	 print "HSP\tseq_inds('query','identical')\t@qidentical\n";
	 print "HSP\tseq_inds('hit','conserved')\t@hconserved\n";
	 print "HSP\tmatches('hit')\t" .   $qid . " " . $qcons . "\n";
	 print "HSP\tmatches('query')\t" . $id . " " . $cons . "\n";
	 print "HSP\talgorithm\t" .        $hsp->algorithm . "\n";
	 print "HSP\tevalue\t" .           $hsp->evalue . "\n";
	 print "HSP\tfrac_identical\t" .   $hsp->frac_identical . "\n";
	 print "HSP\tfrac_conserved\t" .   $hsp->frac_conserved . "\n";
	 print "HSP\tgaps\t" .             $hsp->gaps     . "\n";
	 print "HSP\tquery_string\t" .     $hsp->query_string . "\n";
	 print "HSP\thit_string\t" .       $hsp->hit_string . "\n";
	 print "HSP\thomology_string\t" .  $hsp->homology_string . "\n";
	 print "HSP\tlength('total')\t" .  $hsp->length('total') . "\n";
	 print "HSP\tlength('hit')\t" .    $hsp->length('hit') . "\n";
	 print "HSP\tlength('query')\t" .  $hsp->length('query') . "\n";
	 print "HSP\thsp_length\t" .       $hsp->hsp_length . "\n";
	 print "HSP\tframe\t" .            $hsp->frame . "\n";
	 print "HSP\tnum_conserved\t" .    $hsp->num_conserved . "\n";
	 print "HSP\tnum_identical\t" .    $hsp->num_identical . "\n";
	 print "HSP\trank\t" .             $hsp->rank . "\n";
	 print "HSP\tscore\t" .            $hsp->score . "\n";
	 print "HSP\tbits\t" .             $hsp->bits . "\n";
	 print "HSP\tpercent_identity\t" . $hsp->percent_identity . "\n";
	 print "HSP\tstrand()\t" .         $hsp->strand() . "\n";
	 print "HSP\tstart('hit')\t" .     $hsp->start('hit') . "\n";
	 print "HSP\tstart('query')\t" .   $hsp->start('query') . "\n";
	 print "HSP\tend('hit')\t" .       $hsp->end('hit') . "\n";
	 print "HSP\tend('query')\t" .     $hsp->end('query') . "\n";
	 print "HSP\talignment\n";
	 print $alnIO->write_aln($aln),"\n\n";
      }
   }
}

__END__
