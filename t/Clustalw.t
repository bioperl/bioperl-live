## Bioperl Test Harness Script for Modules
##


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..9\n";
	use vars qw($loaded); }
# Modify following line as required to point to clustalw program directory on your system
BEGIN { $ENV{CLUSTALDIR} = '/home/peter/clustalw1.8/' if( !defined $ENV{CLUSTALDIR} ); }

END {print "not ok 1\n" unless $loaded;}

use Bio::Tools::Alignment::Clustalw;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::SeqIO;
use strict;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 


## Create clustalw alignment factory
my @params = ('ktuple' => 2, 'matrix' => 'BLOSUM');
my  $factory = Bio::Tools::Alignment::Clustalw->new(@params);

if( $factory ) {
    print "ok 2\n";
} else {
    print "not ok 2 , couldn't create alignment factory \n";	
}


## Set parameter test

my $ktuple = 3;
$factory->ktuple($ktuple);

my $new_ktuple = $factory->ktuple();

if( $new_ktuple == 3) {
    print "ok 3\n";
} else {
    print "not ok 3 , couldn't set factory parameter\n";	
}


## Get parameter test

my $what_matrix = $factory->matrix();

if( $what_matrix eq 'BLOSUM') {
    print "ok 4\n";
} else {
    print "not ok 4 , couldn't get factory parameter\n";	
}


## Alignment test (from fasta file)
my $bequiet = 1;
$factory->quiet($bequiet);  # Suppress clustal messages to terminal

my $inputfilename = 't/cysprot.fa';
my $aln;

# If the clustalw program isn't found and executable at the expected location,
# there is no point in executing the remaining tests...

my $clustal_present = Bio::Tools::Alignment::Clustalw->exists_clustal();
unless ($clustal_present) {
	warn "Clustalw program not found. Skipping tests 5 to 9.\n";
    	print "ok 5\n", "ok 6\n", "ok 7\n", "ok 8\n", "ok 9\n";
	exit 0;
}
$aln = $factory->align($inputfilename);
if( $aln->{order}->{'0'} eq 'CATH_HUMAN-1-335' ) {
    print "ok 5\n";
} else {
    print "not ok 5 , failed clustalw alignment using input file\n";	
}

## Alignment test (from BioSeq array)

my $str = Bio::SeqIO->new(-file=> 't/cysprot.fa', '-format' => 'Fasta');
my @seq_array =();

while ( my $seq = $str->next_seq() ) {
	push (@seq_array, $seq) ;
    }

my $seq_array_ref = \@seq_array;

$aln = $factory->align($seq_array_ref);

if( $aln->{order}->{'0'} eq 'CATH_HUMAN-1-335' ) {
    print "ok 6\n";
} else {
    print "not ok 6 , failed clustalw alignment using BioSeq array\n";	
}



## Profile alignment test (from alignment files)

	
my $profile1 = 't/cysprot1a.msf';
my $profile2 = 't/cysprot1b.msf';
$aln = $factory->profile_align($profile1,$profile2);

if( $aln->{order}->{'1'} eq 'CATH_HUMAN-1-335' ) {
    print "ok 7\n";
} else {
    print "not ok 7 , failed clustalw profile alignment using input file\n";	
}

## Profile alignment test (from SimpleAlign objects)

my $str1 = Bio::AlignIO->new(-file=> 't/cysprot1a.msf');
my $aln1 = $str1->next_aln();
my $str2 = Bio::AlignIO->new(-file=> 't/cysprot1b.msf');
my $aln2 = $str2->next_aln();

$aln = $factory->profile_align($aln1,$aln2);

if( $aln->{order}->{'1'} eq 'CATH_HUMAN-1-335' ) {
    print "ok 8\n";
} else {
    print "not ok 8 , failed clustalw profile alignment using SimpleAlign input\n";	
}


## Test aligning (single) new sequence to an alignment

$str1 = Bio::AlignIO->new(-file=> 't/cysprot1a.msf');
$aln1 = $str1->next_aln();
$str2 = Bio::SeqIO->new(-file=> 't/cysprot1b.fa');
my $seq = $str2->next_seq();
$aln = $factory->profile_align($aln1,$seq);


if( $aln->{order}->{'1'} eq 'CATH_HUMAN-1-335' ) {
    print "ok 9\n";
} else {
    print "not ok 9 , failed adding new sequence to alignment\n";	
}





