#!/usr/bin/perl -w

# Simon Chan, SFU Co-op Student
# Bioinformatics, Xenon Genetics Inc.
# http://www.xenongenetics.com
# schan@xenongenetics.com
# skchan@sfu.ca
# BioPerlBlast.pl

# Well, here's an oversimplified description of the program:
# This program takes a query file containing sequences in fasta format and 
# blasts them against a fasta database created by formatdb.

# Here's a much more detailed description:
# Some Vocab - Query: sequences in the file.  Subject: the sequences in the database.

# Each sequence in the file is blasted against the database.  The best hit for this sequence is determined.  Then the program
# asks,

# Question: My current query hits subject X.  Have I hit this subject before?
#
#  Yes: Was that previous P value better or this current P value better 
#       ("better" is defined as lower)?
#
#      Yes: Previous P value better than current P value
#           Get the next best subject hit for this query and repeat until you 
#           have no more hits or a current P value is better
#           than a previous P value
#
#      No:  Current P value is better than previous P value
#           So now, the best hit for this subject is the current query.
#           Take that query and re-blast it and ask the Question again.
#
#  No : My current query has not hit subject X before.  This query best matches 
#       this subject (so far, anyways).  Do next query.
#       Repeat as necessary.

# *************************************************************************************************************************************

use strict;

use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use Bio::Tools::BPlite;
use Getopt::Long;

my ($query_file, $database) = '';
&GetOptions("-q=s"=>\$query_file , "-d=s"=>\$database);

unless($query_file && $database){

    print "You must enter the path to the query file and database!\n";
    print "prompt> perl BioPerlBlast.pl -q path_to_query_file -d path_to_database\n";
    exit;

}


my $blast_factory = Bio::Tools::Run::StandAloneBlast->new( 'database' => "$database", 'program' => "blastn");



# SeqIO
my $seq_in = Bio::SeqIO->new( -file => "$query_file", -format => "fasta");

# Set some max P value.  that is, anything above a this certain P value will be considered to be a "not-good" hit.
my $signif = '1e-10'; # 1e-10 or whatever value you want.

# Temp array for storing new exons.  See *** below
my @new_exons = ();



my (%best_pig_hit, %best_hit, %seen_before, %exon, %hit) = ();

# %best_pig_hit  key: $subject => value: best P value;
# %best_hit key: $subject => value: $query
# %seen_before key: $subject => value: 1.  Used to see if the subject has been hit before.
# %exon key: $query => value: best P value for that $query
# %hit key: $query => value: 1.  Used to see if the query hits ANYTHING at all.
# NOTE: we call the ids for the sequences in the database (the "subjects") 'PIGs'


### the subroutine blast_the_sucker is what "drives" the program.
### Everything below the while loop is where the output is modified and displayed.
### That is, once the flow exits the while loop, all the blasting is done.


while (my $seq = $seq_in->next_seq()){

                 # blast the sequence
                 blast_the_sucker($seq);

} # close next_seq


### ***
### Below is a check to see if any exon that produced a blast hit "accidentally" ended up in the @new_exons array.
### So, I first get all the exons from the best_hit hash and mark them as seen with the seen_exon hash.

my %seen_exon = ();
foreach my $subject (keys %best_hit){


        my $exon = $best_hit{$subject};
        $seen_exon{$exon} = 1;


        # Yes, I know.  I could have just had one line: $seen_before{$best_hit{$subject}} = 1;
        # but that would have been hard to read and confusing if you forgot what the values of those hashes
        # represent.  Sometimes, you have to sacrifice the "one-liners" for readability. ;-)
        # Same goes for using "foreach my $subject" instead of the quicker "for my $subject"

}


### Then, iterate through @new_exons and push the exon into another array if the exon has not been marked.


my @NEW_EXONS = ();
foreach my $exon (@new_exons){push @NEW_EXONS, $exon unless $seen_exon{$exon};}


### All done!
# @NEW_EXONS now contains all exons that did not produce a blast hit.
# process %best_pig_hit, %best_hit, %exon in anyway way you wish.  Add the necessary code/subroutines to accomplish your
# task.  Hope this script has made your life a little less PAINFUL! :)



# print a friendly, 50's nostalgic-like, good-bye. ;-)
print "\n\nFinished.  Have a nice day!\n\n";

exit;

########################################################################################################################################
################################################ Subroutines Below: ####################################################################
########################################################################################################################################

sub blast_the_sucker {

           my ($seq) = @_;

           # BPlite object
           my $blast_report = $blast_factory->blastall($seq);

           # if the query does not produce a blast hit, it will not
           # enter the SUBJECT while loop.
           # therefore, $hit{$blast_report->query} will not be defined
           # before the end of this sub, if $hit is not defined for the query,
           # it will be pushed into @new_exons


           SUBJECT: while( my $sbjct = $blast_report->nextSbjct ){

                 # determine if the query even produces a blast hit:
                 $hit{$blast_report->query} = 1;


                 my $query;

                 # used to determine if any hit occurs...
                 my $good_hit_counter = 0;

                 HSP: while ( my $hsp = $sbjct->nextHSP ) {

                       $query = '';
                       $query = $blast_report->query;
                       $query =~ s/\s.*//; # get rid of anything after the query id.

                       my $p_value = $hsp->P;

                       my $subject = $hsp->hit->seqname;
                       $subject =~ s/\s.*//; # get rid of anything after the subject id.

                       print "Query $query hits Subject $subject with a P value of $p_value\n";


                       # %exon is populated in the sub no_I_have_not_seen_this_subject_before
                       # the first time through this sub, %exon will be undefined.  Thus,
                       # the else block is run.

                       # if the current P value is greater than the significant value, get the next HSP, if any.
                       # if no more HSPs, it's a new exon!

	               if ($exon{$query}){

                                  # if the current P value is greater than $signif, or the current P value is greater than the
                                  # exon's best P value, get the next HSP.
		                  if ($p_value > $signif || $p_value > $exon{$query}){

                                                 next HSP;

		                  } else {

                                                 ++$good_hit_counter;
                                                 have_I_seen_this_subject_before($query, $subject, $p_value);
                                  }



                        ### This else block is only run the first time around for the each exon
	                } else {

                                  # if the current P value is greater than $signif, or the current P value is greater than the
                                  # exon's best P value, get the next HSP.
		                  if ($p_value > $signif){

                                                 next HSP;

		                  } else {

                                                 ++$good_hit_counter;
                                                 have_I_seen_this_subject_before($query, $subject, $p_value);
                                  }
                         }

		 } # close nextHSP


                      if ($good_hit_counter == 0){

                        my $tag = 0;
		        foreach my $new_exon (@new_exons){

		            if ($new_exon eq $query){
                               ++$tag;
                               last;

		            }

		         }

                         # hmmm....doubles in @new_exons...
                         # Well, that can be easily solved.
                         # before you analyze the results, 
                         # mark each of the keys in %best_hit
                         # then iterate through @new_exons and push
                         # the exon into another array unless you've seen it before.
                         # (that is, find the difference between 2 arrays)


                         push @new_exons, $query if $tag== 0;

		      } # close if($good_hit_counter == 0){

	   } # close nextSbjct

           print "\n\n----------\n\n";

	   unless($hit{$blast_report->query}){

                 print "Query ", $blast_report->query, " did not produce hits!  New Exon!\n";
                 my $blast_query = $blast_report->query;
                 $blast_query =~ s/\s.*//;
                 push @new_exons, $blast_query;
	   }

           %hit = ();

} # close sub blast_the_sucker


########################################################################################################################################

# if this subject was seen before, run the sub yes_I_have_seen_this_subject_before
# else, run the sub no_I_have_not_seen_this_subject_before

# yes, loooong subroutine names.

sub have_I_seen_this_subject_before {

my ($best_query, $best_subject, $best_p_value) = @_;

if ($seen_before{$best_subject}){


       yes_I_have_seen_this_subject_before($best_query, $best_subject, $best_p_value);

} else {

       no_I_have_not_seen_this_subject_before($best_query, $best_subject, $best_p_value);

}

return;
} # close sub have_I_seen_this_subject_before 

########################################################################################################################################

# so, no, this subject has not been hit before.

# get the info, populate the hashes and mark this subject as seen before.

sub no_I_have_not_seen_this_subject_before {

my ($current_query, $current_subject, $current_p_value) = @_;


  print "The best hit for $current_query is $current_subject with a P value of $current_p_value\n";

  # mark this subject as seen before.
  $seen_before{$current_subject}  = 1;

  # key: the subject.  the value: the best p value so far.
  $best_pig_hit{$current_subject} = $current_p_value;

  # key: the subject.  the value: the query which produces the best blast hit so far.
  $best_hit{$current_subject}     = $current_query;

  # key: the query.    the value: the best p value it produces.
  $exon{$current_query}           = $current_p_value;


 return;
} # close sub no_I_have_not_seen_this_subject_before

########################################################################################################################################

# so the subject hit has been seen before.

# is the current the p value better or the previous better (by better, we mean lower) ?

sub yes_I_have_seen_this_subject_before {

my ($current_query, $subject, $current_p_value) = @_;

my $past_p_value    = $best_pig_hit{$subject};
my $past_query      = $best_hit{$subject};

# determine which is better: the current P value or the past one:

if ($past_p_value > $current_p_value){

   # we have a new best p value for this subject!

   $best_pig_hit{$subject} = $current_p_value;

   $best_hit{$subject}     = $current_query;

   $exon{$current_query}   = $current_p_value;

   print "*** new best query for subject $subject is $best_hit{$subject} with P value $best_pig_hit{$subject}\n";

   # stick the sequence of the $past_query back into blast queue....
   # fetch the sequnce for $past_query,
   # use Bio::Seq to make seq object.

   my $seq_in = Bio::SeqIO->newFh(-file => "$query_file", -format => "fasta");

   my ($Seq,$Desc)='';

   while (my $line = <$seq_in>){

     if ($line->display_id eq $past_query){

        $Seq = $line->seq;
        $Desc= $line->desc;
        last;
     }

   }

   my $seq_obj = Bio::Seq->new(-id => "$past_query", -seq => "$Seq");

   print "Re-blast $past_query ************************************************\n";
   sleep(3);

   blast_the_sucker($seq_obj);



} elsif ($current_p_value > $past_p_value){

# no, the current query P value is not better than the past P value.
print "$current_query has a p value of $current_p_value for $subject. Not better than $past_p_value of query $past_query! Getting next hit!\n"; sleep (3);


return; # get the next HSP

} else {

# do something if they are equal.


}


return;
} # close sub yes_I_have_seen_this_subject_before


#######################################################################################################################################

