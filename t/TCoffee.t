# -*-Perl-*-
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
# Modify following line as required to point to tcoffee program directory on your system

END {print "not ok 1\n" unless $loaded;
    unlink( qw(cysprot.dnd cysprot1a.dnd) );
}

use Bio::Tools::Run::Alignment::TCoffee;
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


sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

## Create tcoffee alignment factory
my @params = ('ktuple' => 2, 'matrix' => 'BLOSUM');
my  $factory = Bio::Tools::Run::Alignment::TCoffee->new(@params);

test 2, $factory, " couldn't create alignment factory";

## Set parameter test

my $ktuple = 3;
$factory->ktuple($ktuple);

my $new_ktuple = $factory->ktuple();
test 3, $new_ktuple == 3, " couldn't set factory parameter";

## Get parameter test

my $what_matrix = $factory->matrix();
test 4, $what_matrix eq 'BLOSUM', "couldn't get factory parameter";

## Alignment test (from fasta file)
my $bequiet = 1;
$factory->quiet($bequiet);  # Suppress tcoffee messages to terminal

my $inputfilename = 't/cysprot.fa';
my $aln;

# If the tcoffee program isn't found and executable at the expected location,
# there is no point in executing the remaining tests...

my $coffee_present = Bio::Tools::Run::Alignment::TCoffee->exists_tcoffee();
unless ($coffee_present) {
	warn "tcoffee program not found. Skipping tests 5 to 9.\n";
    	print "ok 5\n", "ok 6\n", "ok 7\n", "ok 8\n", "ok 9\n";
	exit 0;
}
$aln = $factory->align($inputfilename);

test 5, $aln->{order}->{'0'} eq 'CYS1_DICDI-1-343', "failed tcoffee alignment using input file";
## Alignment test (from BioSeq array)

my $str = Bio::SeqIO->new(-file=> 't/cysprot.fa', '-format' => 'Fasta');
my @seq_array =();

while ( my $seq = $str->next_seq() ) {
    push (@seq_array, $seq) ;
}

my $seq_array_ref = \@seq_array;

$aln = $factory->align($seq_array_ref);
	
test 6,$aln->{order}->{'0'} eq 'CYS1_DICDI-1-343', "failed tcoffee alignment using BioSeq array ";

## Profile alignment test (from alignment files)

	
my $profile1 = 't/cysprot1a.msf';
my $profile2 = 't/cysprot1b.msf';
$aln = $factory->profile_align($profile1,$profile2);

test 7, $aln->{order}->{'1'} eq 'CATL_HUMAN-1-333', " failed tcoffee profile alignment using input file ". $aln->{order}->{'1'} ;

## Profile alignment test (from SimpleAlign objects)

my $str1 = Bio::AlignIO->new(-file=> 't/cysprot1a.msf');
my $aln1 = $str1->next_aln();
my $str2 = Bio::AlignIO->new(-file=> 't/cysprot1b.msf');
my $aln2 = $str2->next_aln();

$aln = $factory->profile_align($aln1,$aln2);
test 8, $aln->{order}->{'1'} eq 'CATL_HUMAN-1-333', "failed tcoffee profile alignment using SimpleAlign input ". $aln->{order}->{'1'};

## Test aligning (single) new sequence to an alignment

$str1 = Bio::AlignIO->new(-file=> 't/cysprot1a.msf');
$aln1 = $str1->next_aln();
$str2 = Bio::SeqIO->new(-file=> 't/cysprot1b.fa');
my $seq = $str2->next_seq();
$aln = $factory->profile_align($aln1,$seq);

test 9, $aln->{order}->{'1'} eq 'CATH_RAT-1-333', "failed adding new sequence to alignment ". $aln->{order}->{'1'};

END {
opendir(CDIR, ".");
foreach my $file  ( readdir(CDIR)) {
  if( $file =~ /\.dnd$/ ) {
	unlink($file);
  }
}
}
