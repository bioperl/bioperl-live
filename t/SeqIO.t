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
BEGIN { $| = 1; print "1..18\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::SeqIO;
use Bio::SeqIO::MultiFile;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

$str = Bio::SeqIO->new(-file=> 't/test.fasta', '-format' => 'Fasta');

if( $str ) {
    print "ok 2\n";
} else {
    print "not ok 2\n";	
}

$seq = $str->next_seq();

if( $seq->id eq 'roa1_drome' ) {
    print "ok 3\n";
} else {
    print "not ok 3\n";
}

if ($seq->length == 358) {
    print "ok 4\n";
} else {
    print "not ok 4\n";
}

#####
## ChrisDag -- testing out Bio::SeqIO::Raw & SeqIO::GCG
##
## We open a file, test.raw which has 2 raw lines of
## sequence. No formatting at all. Raw sequences are delimited
## simply by a newline. This code tests to make sure we can
## create 2 sequential bioseq objects out of the raw file without
## breaking or getting confused.
##
## Not tested yet: ability to write a raw formatted stream
## Not tested yet: ability to write a GCG formatted stream

$str = Bio::SeqIO->new(-file=> 't/test.raw', '-format' => 'Raw');

if( $str ) {
    print "ok 5\n";
} else {
    print "not ok 5 , unable to open stream from raw sequence DB\n";	
}

if($seq = $str->next_seq()) { print "ok 6\n";
 } else { print "not ok 6 , failed to read 1st raw sequence from stream,\n"; }
print "Sequence 1 of 2 from Raw stream:\n", $seq->seq;

print "\n\n";

if($seq = $str->next_seq()) { print "ok 7\n";
 } else { print "not ok 7 , failed to read 2nd raw sequence from stream.\n"; }
print "Sequence 2 of 2 from Raw stream:\n", $seq->seq;
print $seq->seq;
print "\n";


## Now we test Bio::SeqIO::GCG

$str = Bio::SeqIO->new(-file=> 't/test.gcg', '-format' => 'GCG');

if( $str ) {
    print "ok 8\n";
} else {
    print "not ok 8 , unable to open stream from GCG sequence file\n";	
}

if($seq = $str->next_seq()) { print "ok 9\n";
 } else { print "not ok 9 , failed to read GCG sequence from stream,\n"; }
print "Sequence 1 of 1 from GCG stream:\n", $seq->seq, "\n";

## Now we test Bio::SeqIO::GCG output writing

$str = Bio::SeqIO->new(-file=> '>t/gcg.out', '-format' => 'GCG');

$str->write_seq($seq);

print "ok 10\n";

#####
## End of ChrisDag's SeqIO tests.
#####

## Now we test Bio::SeqIO::GenBank
$str = Bio::SeqIO->new(-file=> 't/test.genbank', '-format' => 'GenBank');


if( $str ) {
    print "ok 11\n";
} else {
    print "not ok 11 , unable to open stream from GenBank sequence file\n";	
}

if($seq = $str->next_seq()) { print "ok 12\n";
 } else { print "not ok 12 , failed to read GenBank sequence from stream,\n"; }
print "Sequence 1 of 1 from GenBank stream:\n", $seq->seq, "\n";


$str = Bio::SeqIO->new(-file=> '>t/genbank.out', '-format' => 'GenBank');

$str->write_seq($seq);

print "ok 13\n";

# please leave this as the last line:
$str = undef;


# EMBL format

$ast = Bio::SeqIO->new( '-format' => 'embl' , -file => 't/roa1.dat');

while( my $as = $ast->next_seq() ) {
       if( ! defined $as->seq ) {
	   print "not ok 14\n";
	   }
      }


print "ok 14\n";

$ast = Bio::SeqIO->new( '-format' => 'GenBank' , -file => 't/roa1.genbank');

while( my $as = $ast->next_seq() ) {
       if( ! defined $as->seq ) {
	   print "not ok 15\n";
	   }
      }


print "ok 15\n";

$mf = Bio::SeqIO::MultiFile->new( '-format' => 'Fasta' , -files => ['t/multi_1.fa','t/multi_2.fa']);

print "ok 16\n";

# read completely to the end

while( $seq = $mf->next_seq() ) {
	$temp = $seq->display_id;
}
$temp = undef;
print "ok 17\n";


$ast = Bio::SeqIO->new( '-format' => 'Swiss' , -file => 't/roa1.swiss');

while( my $as = $ast->next_seq() ) {
   if( ! defined $as->seq || $as->id ne 'ROA1_HUMAN' ) {
	print "not ok 18\n";
	print STDERR "id is ".$as->id."\n";
   } else {
     print "ok 18\n";
   }
}



