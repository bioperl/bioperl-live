#!/usr/bin/perl
# PROGRAM  : clustalw.pl
# PURPOSE  : Demonstrate possible uses of Bio::Tools::Run::Alignment::Clustalw.pm
# AUTHOR   : Peter Schattner schattner@alum.mit.edu
# CREATED  : Oct 06 2000
#
# INSTALLATION
#
# You will need to have installed clustalw and to ensure that Clustalw.pm can find it.
# This can be done in different ways (bash syntax):
#	   export PATH=$PATH:/home/peter/clustalw1.8
#  or
#     define an environmental variable CLUSTALDIR:
#	   export CLUSTALDIR=/home/peter/clustalw1.8   
#  or
#     include a definition of an environmental variable CLUSTALDIR in every
#     script that will use Clustal.pm.
#	   BEGIN {$ENV{CLUSTALDIR} = '/home/peter/clustalw1.8/'; }
#
#  We are going to demonstrate 3 possible applications of Clustalw.pm:
#	1. Test effect of varying clustalw alignment parameter(s) on resulting alignment
#	2. Test effect of changing the order that sequences are added to the alignment
#		on the resulting alignment
#	3. Test effect of incorporating an "anchor point" in the alignment process
#
#  Before we can do any tests, we need to set up the environment, create the factory
#  and read in the unaligned sequences.
#

#BEGIN {
#	$ENV{CLUSTALDIR} = '/home/peter/clustalw1.8/';
#}

use Getopt::Long;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::SeqIO;
use strict;

# set some default values
my $infile = 't/data/cysprot1a.fa';
my @params = ('quiet' => 1 );
my $do_only = '123';   # string listing examples to be executed. Default is to
			              # execute all tests (ie 1,2 and 3)
my $param = 'ktuple';  # parameter to be varied in example 1
my $startvalue = 1;    # initial value for parameter $param
my $stopvalue = 3;     # final value for parameter $param
my $regex = 'W[AT]F';  # regular expression for 'anchoring' alignment in example 3
my $extension = 30; 	  # distance regexp anchor should be extended in each direction
			              # for local alignment in example 3
my $helpflag = 0;      # Flag to show usage info.

# get user options
my @argv = @ARGV;  # copy ARGV before GetOptions() massacres it.

&GetOptions("h!" => \$helpflag, "help!" => \$helpflag,
				"in=s" => \$infile,
				"param=s" => \$param,
				"do=s" =>  \$do_only,
				"start=i" =>  \$startvalue,
				"stop=i" =>  \$stopvalue,
				"ext=i" =>  \$extension,
				"regex=s" =>  \$regex,) ;

if ($helpflag) { &clustalw_usage(); exit 0;}

# create factory & set user-specified global clustalw parameters
foreach my $argv (@argv) {
	unless ($argv =~ /^(.*)=>(.*)$/) { next;}
	push (@params, $1 => $2);
}
my  $factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);
	

# put unaligned sequences in a Bio::Seq array
my $str = Bio::SeqIO->new(-file=> $infile, '-format' => 'Fasta');
my ($paramvalue, $aln, $subaln, @consensus, $seq_num, $string, $strout, $id);
my @seq_array =();
while ( my $seq = $str->next_seq() ) { push (@seq_array, $seq) ;}

# Do each example that has digit present in variable $do_only
$_= $do_only;
/1/ && &vary_params();
/2/ && &vary_align_order();
/3/ && &anchored_align();

## End of "main"

#################################################
#   vary_params(): Example demonstrating varying of clustalw parameter
#

sub vary_params {

	print "Beginning parameter-varying example... \n";

	# Now we'll create several alignments, 1 for each value of the selected
	# parameter. We also compute a simple consensus string for each alignment.
	# (In the default case, we vary the "ktuple" parameter,  creating 3
	# alignments using ktuple values from 1 to 3.)

	my $index =0;
	for ($paramvalue = $startvalue; $paramvalue < ($stopvalue + 1); $paramvalue++) {
		$factory->$param($paramvalue);  # set parameter	value
		print "Performing alignment with $param = $paramvalue \n";
		$aln = $factory->align(\@seq_array);
		$string = $aln->consensus_string(); # Get consensus of alignment
		# convert '?' to 'X' at non-consensus positions
		$string =~ s/\?/X/g;
		$consensus[$index] = Bio::Seq->new(-id=>"$param=$paramvalue",-seq=>$string);
		$index++;
	}
	# Compare consensus strings for alignments with different $param values by
	# making an alignment of the different consensus strings
	# $factory->ktuple(1);  # set ktuple parameter	
	print "Performing alignment of $param consensus sequences \n";
	$aln = $factory->align(\@consensus);
	$strout = Bio::AlignIO->newFh('-format' => 'msf');
	print $strout $aln;

	return 1;
}


#################################################
#   vary_align_order():
#
# For our second example, we'll test the effect of changing the order
# that sequences are added to the alignment

sub vary_align_order {

	print "\nBeginning alignment-order-changing example... \n";

	@consensus = ();  # clear array
	for ($seq_num = 0; $seq_num < scalar(@seq_array); $seq_num++) {
		my $obj_out = shift @seq_array;  # remove one Seq object from array and save
		$id = $obj_out->display_id;
		# align remaining sequences
		print "Performing alignment with sequence $id left out \n";
		$subaln = $factory->align(\@seq_array);
		# add left-out sequence to subalignment
		$aln = $factory->profile_align($subaln,$obj_out);
		$string = $aln->consensus_string(); # Get consensus of alignment
		# convert '?' to 'X' for non-consensus positions
		$string =~ s/\?/X/g;
		$consensus[$seq_num] = Bio::Seq->new(-id=>"$id left out",-seq=>$string);
		push @seq_array, $obj_out;  # return Seq object for next (sub) alignment
	}

	# Compare consensus strings for alignments created in different orders
	# $factory->ktuple(1);  # set ktuple parameter	
	print "\nPerforming alignment of consensus sequences for different reorderings \n";
	print "Each consensus is labeled by the sequence which was omitted in the initial alignment\n";
	$aln = $factory->align(\@consensus);
	$strout = Bio::AlignIO->newFh('-format' => 'msf');
	print $strout $aln;

	return 1;
}

#################################################
#   anchored_align()
#
# For our last example, we'll test a way to perform a local alignment by
# "anchoring" the alignment to a regular expression.  This is similar
# to the approach taken in the recent dbclustal program.
# In principle, we could write a script to search for a good regular expression
# to use. Instead, here we'll simply choose one manually after looking at the
# previous alignments.

sub anchored_align {

	my @local_array = ();
	my @seqs_not_matched = ();

	print "\n Beginning anchored-alignment example... \n";

	for ($seq_num = 0; $seq_num < scalar(@seq_array); $seq_num++) {
		my $seqobj = $seq_array[$seq_num];
		my $seq =  $seqobj->seq();
		my $id =  $seqobj->id();
		# if $regex is not found in the sequence, save sequence id name and set
		# array value =0 for later
		unless ($seq =~/$regex/) {
			$local_array[$seq_num] = 0;
			push (@seqs_not_matched, $id) ;
			next;
		}
		# find positions of start and of subsequence to be aligned
		my $match_start_pos = length($`);
		my $match_stop_pos = length($`) + length($&);
		my	$start =  ($match_start_pos - $extension) > 1 ? 
		  ($match_start_pos - $extension) +1 : 1;
		my	$stop =  ($match_stop_pos + $extension) < length($seq) ?
		  ($match_stop_pos + $extension) : length($seq);
		my $string = $seqobj->subseq($start, $stop);

		$local_array[$seq_num] = Bio::Seq->new(-id=>$id, -seq=>$string);
	}
	@local_array = grep $_ , @local_array; # remove array entries with no match

	# Perform alignment on the local segments of the sequences which match "anchor"
	$aln = $factory->align(\@local_array);
	my $consensus  = $aln->consensus_string(); # Get consensus of local alignment

	if (scalar(@seqs_not_matched) ) {
		print " Sequences not matching $regex : @seqs_not_matched \n"
	} else {
		print " All sequences match $regex : @seqs_not_matched \n"
	}
	print "Consensus sequence of local alignment: $consensus \n";

	return 1;
}

#----------------
sub clustalw_usage {
#----------------

#-----------------------
# Prints usage information for general parameters.

    print STDERR <<"QQ_PARAMS_QQ";

 Command-line accessible script variables and commands:
 -------------------------------
 -h 		       :  Display this usage info and exit.
 -in <str>	    :  File containing input sequences in fasta format (default = $infile) .
 -do <str>	    :  String listing examples to be executed. Default is to execute
		             all tests (ie default = '123')
 -param <str>   :  Parameter to be varied in example 1. Any clustalw parameter
		             which takes inteer values can be varied (default = 'ktuple')
 -start <int>   :  Initial value for varying parameter in example 1 (default = 1)
 -stop <int>    :  Final value for varying parameter (default = 3)
 -regex   <str> :  Regular expression for 'anchoring' alignment in example 3
                   (default = $regex)
 -ext <int>     :  Distance regexp anchor should be extended in each direction
		             for local alignment in example 3   (default = 30)

In addition, any valid Clustalw parameter can be set using the syntax 
"parameter=>value" as in "ktuple=>3"

So a typical command lines might be:
 > clustalw.pl -param=pairgap -start=2 -stop=3 -do=1 "ktuple=>3"
or
 > clustalw.pl -ext=10 -regex='W[AST]F' -do=23 -in='t/cysprot1a.fa'

QQ_PARAMS_QQ

}
