#!/usr/bin/perl

# PROGRAM  : standaloneblast.pl
# PURPOSE  : Demonstrate possible uses of Bio::Tools::StandAloneBlast.pm
# AUTHOR   : Peter Schattner schattner@alum.mit.edu
# CREATED  : Nov 01 2000
# REVISION : $Id$
#
# INSTALLATION
#

# You will need to enable Blast to find the Blast program. This can be done
# in (at least) three ways:
#  1. Modify your $PATH variable to include your Blast directory as in (for Linux):
#	export PATH=$PATH:/home/peter/blast   or
#  2. define an environmental variable blastDIR:
#	export BLASTDIR=/home/peter/blast   or
#  3. include a definition of an environmental variable BLASTDIR in every script that will
#     use StandAloneBlast.pm.
#	BEGIN {$ENV{BLASTDIR} = '/home/peter/blast/'; }
#
#  We also need to select the database to be used
my $amino_database = 'swissprot';


# 
#  We are going to demonstrate 3 possible applications of StandAloneBlast.pm:
#	1. Test effect of varying choice of substitution matrices	
#	2. Test effect of varying choice of gap penalty 
#	3. Comparison of results of psiblast depending on whether psiblast itself is used
#	to identify an alignment to use for blasting or whether an external alignment is given to 
#	psiblast
#

use Getopt::Long;
use Bio::SimpleAlign;
use Bio::Tools::StandAloneBlast;
use Bio::Tools::BPlite::Sbjct;
use Bio::AlignIO;
use Bio::SeqIO;
use strict;

# set some default values
my $queryseq = 't/cysprot1.fa';
my $executable = 'blastpgp';
my $queryaln = 't/cysprot.msf';
my @params = ('database' => $amino_database);
my $do_only = ''; 	# string listing examples to be executed. Default is to execute
			# all tests (ie 1,2 and 3)
my $example1param = 'MATRIX';  # parameter to be varied in example 1
my $example2param = 'GAP';  # parameter to be varied in example 1
my $example1values = [ 'BLOSUM62', 'BLOSUM80', 'PAM70']; # MATRIX values to be tried
my $example2values = [ 7, 9, 25]; # GAP values to be tried
my $queryalnformat = 'msf';
my  $jiter = 2;
my  $masktype = 'all';
my $helpflag = 0;   # Flag to show usage info.

# get user options
my @argv       = @ARGV;  # copy ARGV before GetOptions() massacres it.
my $paramvalstring;

&GetOptions("h!" => \$helpflag, "help!" => \$helpflag,
	"in=s" => \$queryseq,
	"inaln=s" => \$queryaln,
	"alnfmt=s" => \$queryalnformat,
	"param=s" => \$example1param,
	"exec=s" => \$executable,
	"paramvals=s" => \$paramvalstring,
	"do=i" =>  \$do_only,
	"mask=s" => \$masktype,
	"iter=i" =>  \$jiter,
	) ;

if ($paramvalstring) { @$example1values = split (":", $paramvalstring); }

 if ($helpflag) { &example_usage(); exit 0;}

# create factory & set user-specified global blast parameters
foreach my $argv (@argv) {
	unless ($argv =~ /^(.*)=>(.*)$/) { next;}
	push (@params, $1 => $2);
}
my  $factory = Bio::Tools::StandAloneBlast->new(@params);
	
# If "do" variable not set, do all three examples
if ( ! $do_only)  {
	&vary_params($queryseq, $example1param, $example1values);
# To compare gap penalties of 7, 9 and 25 we need to set the scoring matrix to BLOSUM62
#  and extension penalty to 2 (these are limitations of BLAST)

	$factory->MATRIX('BLOSUM62');  

	$factory->EXTENSION(2);  
	&vary_params($queryseq, $example2param, $example2values);
# For the psiblast tests we want to restore gap opening and extension values to their defaults

	$factory->GAP(11);  

	$factory->EXTENSION(1);  
	&aligned_blast($queryseq, $queryaln, $queryalnformat, $jiter, $masktype);
} elsif ($do_only  == 1) {
	&vary_params($queryseq,$example1param, $example1values);
} elsif ($do_only  == 3 ) {
	&aligned_blast ($queryseq, $queryaln, $queryalnformat, $jiter, $masktype);
}  else {
	&example_usage();
}

exit 0;

##########
## End of "main"




#################################################
#   compare_runs(): Prints out display of which hits were found by different methods
#	Various methods are labeled by "tags" found in array @runtags
#
#  args: 
#	$typetag  -  label describing type of "tags"
#	$runtags  -  reference to array @runtags
#	$hashhits  - reference to hash of all the hits found by all runs (%hashhits) 
#		value for each hit is string which is the concatenation of all the "tags" of
#		runs that found that hit
#  returns: nothing  

sub compare_runs {
my $typetag = shift;
my $runtags = shift;

my $hashhits = shift;   

my ($tag, @taghits);

print "Comparing BLAST results... \n";

# Get total number of hits found by any method
my $numhits = keys %$hashhits ;   # scalar context to get total number of hits by all methods
print  "Total number of hits found: $numhits \n";

# Get total number of hits found by every method
my $alltags =  join ( "" ,  @$runtags );
my  @alltaghits = grep  $$hashhits{$_} =~ /$alltags/  ,  keys %$hashhits;
print  " Number of hits found by every method / parameter-value: " ,   scalar(@alltaghits), "\n";

# If one desires to see the hits found by all methods, uncomment next 2 lines
#print  " Hits were found all methods / parameters: \n";
#print   join ( "\n", @alltaghits ) ,  "\n";

# For each method/parameter-value (labeled by type) display  hits found 
# exclusively by that method
  foreach $tag (@$runtags)  {
	@taghits = grep  $$hashhits{$_} =~ /^$tag$/  ,  keys %$hashhits;
	print  " Hits found only when $typetag was $tag: \n";
	print   join ( "\n", @taghits ) ,  "\n";
  }
return 1;
}


#################################################

#   vary_params(): Example demonstrating varying of parameter

#

#  args: 
#	$queryseq  - query sequence (can be filename (fasta),  or Bio:Seq object) 
#	$param  - name of parameter to be varied 
#	$values  - reference to array of values to be used for the parameter 
#  returns: nothing  



sub vary_params {

my $queryseq = shift;   
my $param = shift;   
my $values = shift;  


print "Beginning $param parameter-varying example... \n";

# Now we'll perform several blasts, 1 for each value of the selected parameter.
# We also compute a simple consensus string for each alignment.
# In the first default case, we vary the MATRIX substitution parameter,  
# creating 3 BLAST reports, using MATRIX values of  BLOSUM62, BLOSUM80 or PAM70.
# In the second default case, we vary the GAP penalty parameter,  
# creating 3 BLAST reports, using GAP penalties of 7, 9 and 25.
# In either case we then automatically parse the resulting report to identify which hits
# are found with any of the parameter values and which with only some of them.
#
# To test the BLAST results to some other parameter it is only necessary to change the 
# parameters passed to the script on the commandline.  The only tricky part is that the BLAST
# program itself only supports a limited range of parameters.  See the BLAST documentation.

my ($report, $sbjct, $paramvalue);

my $hashhits = { };     # key is hit id, value is string of param values for which hit was found

foreach $paramvalue (@$values)  {

	$factory->$param($paramvalue);  # set parameter value

            print "Performing BLAST with $param = $paramvalue \n";

	$report = $factory->$executable($queryseq);
	while ($sbjct = $report->nextSbjct) {
		$hashhits->{$sbjct->name} .= "$paramvalue";			
	}
}

&compare_runs( $param , $values , $hashhits);  
 
return 1;

}


#################################################
#  aligned_blast ():
#
#
#  args: 
#	$queryseq  - query sequence (can be filename (fasta),  or Bio:Seq object) 
#	$queryaln  - file containing alignment to be used to "jumpstart" psiblast in "-B mode"
#			$queryaln *must contain $queryseq with the same name and length
#				(psiblast is very picky)
#	$queryalnformat  - format of alignment (can = "fasta", "msf", etc)
#	$jiter  - number of iterations in psiblast run
#	$masktype  - key to indicate what sort of mask to use in creating jumpstart alignment file
#		Currently only options are: "all" => use position specific matrix at all residues,  or
#						   "none" => use default (eg BLOSUM) at all residues
#  returns: nothing  





# For this example, we'll compare the results of psiblast depending on whether psiblast itself is 

#  used to identify an alignment to use for blasting or whether an external alignment is given to 
#  psiblast

sub aligned_blast {


my     $queryseq  =  shift; 
my	$queryaln  =  shift; 
my	$queryalnformat  =  shift;
my	$jiter  =  shift;
my	$masktype  =  shift;

my $hashhits = { };
my ($sbjct, $id);

print "\nBeginning aligned blast example... \n";

# First we do a "plain" psiblast multiple-iteration search

print "\nBeginning multiple-iteration psiblast ... \n";

my $tag1 = 'iterated';
$factory->j($jiter);    # 'j' is blast parameter for # of iterations


my $report1 = $factory->blastpgp($queryseq);
my $total_iterations = $report1->number_of_iterations;
my $last_iteration = $report1->round($total_iterations);


 while($sbjct = $last_iteration->nextSbjct) {
		$hashhits->{$sbjct->name} .= "$tag1";			
	}


# Then we do a  single-iteration psiblast search but with a specified alignment to 
#  "jump start" psiblast


print "\nBeginning jump-start psiblast ... \n";


my $tag2 = 'jumpstart'; 

$factory->j('1');    # perform single iteration

# Get the alignment file
my $str = Bio::AlignIO->new(-file=> "$queryaln", '-format' => "$queryalnformat", );
my $aln = $str->next_aln();


# Create the proper mask for 'jumpstarting'
my $mask = &create_mask($aln, $masktype);


my $report2 = $factory->blastpgp($queryseq, $aln, $mask);
while($sbjct = $report2->nextSbjct) {
		$hashhits->{$sbjct->name} .= "$tag2";			
}

my $tagtype = 'iterated_or_jumpstart'; 
my $values = [ $tag1, $tag2];

&compare_runs( $tagtype , $values , $hashhits);  
 
return 1;

}


#################################################
# create_mask(): creates a mask for the psiblast jumpstart alignment
#	that determines at what residues position-specific scoring matrices (PSSMs)
#	are used and at what residues default scoring matrices (eg BLOSUM)
#	are used. See psiblast documentation for more details,
#	Currently only two very simple masks are implemented here:
#	Either *all* positions use the position specific matrix or *none* do
#
#  args: 
#	$aln  -  SimpleAlign object with alignment
#	$masktype  -  label describing type of "tags"
#  returns: actual mask, ie a string of 0's and 1's which is the same length as each
#		sequence in the alignment and has a "1" where (PSSMs) are to be used
#		and a "0" at  

sub create_mask {
my $aln = shift;
my $masktype = shift;
my $mask = "";

unless ( $aln->is_flush() )  { die "psiblast jumpstart requires all sequences to be same length \n"; }
my $len = $aln->length_aln();
if ($masktype =~ /all/i  ) {  $mask =   '1' x $len; }
if ($masktype =~ /none/i  ) {  $mask =   '0' x $len; }
return $mask ;
}


#----------------
sub example_usage {
#----------------

#-----------------------
# Prints usage information for general parameters.

    print STDERR <<"QQ_PARAMS_QQ";

 Command-line accessible script variables and commands:
 -------------------------------
 -h 		:  Display this usage info and exit.
 -in <str>	:  File containing input sequences in fasta format (default = $queryseq) .
 -inaln <str>	:  File containing input alignment for example 3 (default = $queryaln) .
 -alnfmt <str>	:  Format of input alignment for example 3, eg "msf", "fasta", "pfam".
		   (default = $queryalnformat) .
 -do <int>	:  Number of test to be executed ("1" => vary parameters,
		   "3" => compare iterated & jumpstart psiblast.) If omitted,
		   three default tests performed.
 -exec <str>  	:  Blast executable to be used in example 1.  Can be "blastall" or
		   "blastpgp" (default is "blastpgp").
 -param <str>  	:  Parameter to be varied in example 1. Any blast parameter
		   which takes inteer values can be varied (default = 'MATRIX')
 -paramvals <str>:  String containing parameter values in example 1, separated
		   by ":"'s. (default = 'BLOSUM62:BLOSUM80:PAM70')
 -iter <int>    :  Maximum number of iterations in psiblast in example 3 (default = 2)
 -mask <str>	:  Type of alignment mask to be used in example 3. can equal "all"
		   or "none". (default = 'all')

In addition, any valid Blast parameter can be set using the syntax "parameter=>value" as in "database=>swissprot"

So some typical command lines might be:
 >standaloneblast.pl -do 1 -param MATRIX -paramvals 'BLOSUM62:BLOSUM80'
or
 >standaloneblast.pl -do 1 -exec blastall -param q -paramvals '-1:-7' -in='t/dna1.fa' "pr=>blastn" "d=>ecoli.nt"
or
 >standaloneblast.pl -do 3 -mask none -iter 3

QQ_PARAMS_QQ
}
