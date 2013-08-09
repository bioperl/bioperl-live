#! /usr/bin/perl
#####################################################
#                 Fasta
#                     |
#                     Align
#
#                     By
#                 Antony Vincent 
#          (a.vincent.0312@gmail.com)
#           	
# FastAlign is a perl script which uses the heuristic method
# of tfasty36 for find similarity between a query sequence
# in amino acids and a sequence in nucleotides. It provides
# a more intuitive output to find exon-intron junctions. 
# The query string is in amino acids and the hit string is
# in nucleotides. There are extra nucleotides at the end of
# the hit string (option -diff and by default = 10), that 
# allow to verify if the intron start with common rules 
# (5'-GTGCGA-... for group II intron and after an exonic T 
# for group I intron).
# 
# The FASTA version can be changed by the user by changing
# the line with tfasty36 for tfastyXX.
#
# If you have Emboss, you can genarate a graphic with option
# -graph 1.
#
# For complete help: type perl fastalign.pl -help
#              Last Update: 01/06/13
#######################################################


=head1

Copyright (C) 2013  Antony Vincent

Licence:

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

use strict; 
use Bio::SearchIO;
use Bio::SeqIO;
use Getopt::Long;
use Bio::SeqUtils;


## Set the default variables
	my $identity     = 75;
	my $length      = 50;
	my $diff      = 10;
	my $out      = 'output';
	my $graphic      = 10;
	my $query;
	my $library;
	my $help;

GetOptions(
    'seq=s'    => \$query,
    'db=s'     => \$library, 
    'graph=s'     => \$graphic,
    'i=i'     => \$identity,
    'l=i'     => \$length,
    'diff=s'    => \$diff,
    'out=s'     => \$out,   
    'help!'     => \$help,
) 
or die "Incorrect usage! Try perl fastalign.pl -help for an exhaustif help.\n";
###
	if( $help ) 
{ # if start
    print "\n";
    print "Two options are required:\n";
    print "	-seq: Your sequence in amino acids\n";
    print "	-db: The sequence in nucleotides (Could be a whole genome or a partial sequence...)\n";
    print "\n";
    print "There are few miscellaneous options:\n";
    print "	-i: Minimum identity percentage (default = 75)\n";
    print "	-l: Minimum match length (default = 50)\n";
    print "	-diff: Difference between the hit string and the alignement (default = 10)\n";
    print "	-out: Name of the output file (default = output.txt)\n";
    print "	-graph: If this option = 1, a graph will be generated (default = no)\n";
} # if help

	else 
{ # else start
my $date = `date`;
open (WRITE, ">>$out.txt"); ## Open the output file
print WRITE "		Fasta\n";
print WRITE "		    |\n";
print WRITE "		    Align\n\n";
print WRITE "Date:", $date, "\n";
print WRITE "PARAMETERS\n";
print WRITE "Minimum identity =", $identity, "\n";
print WRITE "Minimum length =", $length, "\n";
print WRITE "Diff =", $diff, "\n\n";

	if ( $graphic == 1 )
{
open (WRITE, ">>$out.txt"); ## Open the output file
open (WRITE2, ">>lindna.lnp"); ## Open the lindna config file

## start the lindna header
print WRITE2 "start";
print WRITE2 "\t";
print WRITE2 "1";
print WRITE2 "\n";
print WRITE2 "End";
print WRITE2 "\t";
my $seqio_obj = Bio::SeqIO->new(-file => "$library", -format => "fasta" );
my $seq_obj = $seqio_obj->next_seq;
my $count_obj = $seq_obj->length;
print WRITE2 "$count_obj";
print WRITE2 "\n\n";
print WRITE2 "group";
print WRITE2 "\n";
}
	else
{
	print "No graphic generated \n";
}
## run tfasty36
my $fh;                   
my $fasta   = "tfasty36"; # <------ You can change this line for newest fasta algorithm

my $command = "$fasta $query $library";
 
open $fh,"$command |" || die("cannot run fasta cmd of $command: $!\n");
 
my $searchio  = Bio::SearchIO->new(-format => 'fasta', -fh => $fh);
print WRITE "Fasta algorithm:", $fasta, "\n\n";
			## start the parsing part of the script

while( my $result = $searchio->next_result ) {
  ## $result is a Bio::Search::Result::ResultI compliant object
  while( my $hit = $result->next_hit ) {
    ## $hit is a Bio::Search::Hit::HitI compliant object
    while( my $hsp = $hit->next_hsp ) {
      ## $hsp is a Bio::Search::HSP::HSPI compliant object
      if( $hsp->length('total') > $length ) {
        if ( $hsp->percent_identity >= $identity ) {

          		## Generals informations
print WRITE "Rank=", $hsp->rank, "\n",
            "Query=",   $result->query_name, "\n",
            "Hit=",        $hit->name, "\n" ,
            "Length=",     $hsp->length('total'), "\n",
            "Percent_id=", $hsp->percent_identity, "\n",
	    "Strand=", $hsp->strand('hit'), "\n";

	     
			print WRITE "\n";

		
  	 
		
			## Search for nucleotide sequences
			print WRITE "\n";
	my $start_hit = $hsp->start('hit'), "\n";
	my $end_hit = $hsp->end('hit') , "\n";	
   	my $in  = Bio::SeqIO->new(-file => "$library" , '-format' => 'fasta');

    while ( my $seq = $in->next_seq() ) {#1
    	  
		## looking for query position
			my $start_query = $hsp->start('query'), "\n";
			my $end_query = $hsp->end('query') , "\n";
		## aa_to_3aa
        my $seqobj = Bio::PrimarySeq->new ( -seq => $hsp->query_string);
	my $polypeptide_3char = Bio::SeqUtils->seq3($seqobj);

                ## modify the homology string
      my $homo = $hsp->homology_string;
         $homo =~ s/:/|||/g;
         $homo =~ s/\./***/g;
         $homo =~ s/ /XXX/g;
	
	
			## HSP
           
	print WRITE "Query($start_query,$end_query)\n";
	print WRITE "Hit($start_hit,$end_hit)\n\n";
	print WRITE $polypeptide_3char, "\n";
	print WRITE $homo, "\n";
	print WRITE $seq->subseq($start_hit,$end_hit+$diff), "\n";


        if ( $graphic == 1)
{ ## if start
	 ## write in lindna file
	print WRITE2 "label", "\n", "Block", "\t", "$start_hit", "\t",
                      "$end_hit", "\t", "3", "\t", "H", "\n";
        print WRITE2 "Exon", $hsp->rank, "\n";
        print WRITE2 "endlabel";
        print WRITE2 "\n\n";
} ## if end
	else
{print "No lindna file generated\n";}
} #1
print WRITE "\n";

        }
      }
    }  
  }
}
        if ( $graphic == 1)
{ ## if start
print WRITE2 "endgroup";
system ("lindna -infile lindna.lnp -ruler y -blocktype filled -graphout svg");
system ("rm *.lnp");
} ## if end
} # else end
