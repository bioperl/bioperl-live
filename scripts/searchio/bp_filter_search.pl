#!perl
use strict;
use warnings;


=head1 NAME

bp_filter_search - filters searchio results, outputting a tab delimited summary

=head1 SYNOPSIS

  #bp_filter_search -format blast -score 200 < search.bl > search.tab

=head1 DESCRIPTION 

This script filters searchio results allowing a number of different
filters to be applied before outputting to stdout in a tab delimited
format.

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

=head1 AUTHOR

  Ewan Birney <birney@ebi.ac.uk>

=cut

use Bio::SearchIO;
use Getopt::Long;

my ($format,$score);

$format  = 'blast';
$score   = 150;

GetOptions(
	   'format:s'  => \$format,
	   'score:s'  => \$score,
	);


my $searchin = Bio::SearchIO->new( -format => $format);


while( (my $result = $searchin->next_result()) ) { 
  while( (my $hit = $result->next_hit())) {

      if( $score ) {
        if( $hit->raw_score < $score ) {
           next;
        }
       }


      foreach my $hsp ( $hit->hsps() ) {
         print $result->query_name,"\t",$hit->score,"\t",$hsp->start,"\t",$hsp->end,"\t",$hsp->strand,"\t",$hsp->hseq_id,"\t",$hsp->hstart,"\t",$hsp->hend,"\t",$hsp->strand,"\n";
      }
    }
}


