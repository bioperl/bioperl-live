## Bioperl Test Harness Script for Modules
##
# CVS Version
# $Id$


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
BEGIN { $| = 1; print "1..16\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::SeqFeature::Generic;
use Bio::SeqFeature::FeaturePair;

$loaded = 1;
print "ok 1\n";    # 1st test passes.

$feat = new Bio::SeqFeature::Generic ( -start => 40,
				       -end => 80,
				       -strand => 1,
				       -primary => 'exon',
				       -source  => 'internal',
				       -tag => {
					   silly => 20,
					   new => 1
					   }
				       );

if( $feat->start == 40 ) { 
    print "ok 2\n";
} else {
    print "not ok 2\n";
}


if( $feat->end == 80 ) { 
    print "ok 3\n";
} else {
    print "not ok 3\n";
}


if( $feat->primary_tag eq 'exon' ) { 
    print "ok 4\n";
} else {
    print "not ok 4\n";
}


if( $feat->source_tag eq 'internal' ) { 
    print "ok 5\n";
} else {
    print "not ok 5\n";
}




$str = $feat->gff_string();
$str = ""; # shut up -w

# we need to figure out the correct mapping of this stuff
# soon

#if( $str ne "SEQ\tinternal\texon\t40\t80\t1\t.\t." ) {
#    print "not ok 3\n";
#} else {
#    print "ok 3\n";
#}

print "ok 6\n";

$pair = new Bio::SeqFeature::FeaturePair();

print "ok 7\n";

$feat2 = new Bio::SeqFeature::Generic ( -start => 400,
				       -end => 440,
				       -strand => 1,
				       -primary => 'other',
				       -source  => 'program_a',
				       -tag => {
					   silly => 20,
					   new => 1
					   }
				       );

$pair->feature1($feat);
$pair->feature2($feat2);

print "ok 8\n";

if( $pair->start == 40 ) { 
    print "ok 9\n";
} else {
    print "not ok 9\n";
}


if( $pair->end == 80 ) { 
    print "ok 10\n";
} else {
    print "not ok 10\n";
}


if( $pair->primary_tag eq 'exon' ) { 
    print "ok 11\n";
} else {
    print "not ok 11\n";
}


if( $pair->source_tag eq 'internal' ) { 
    print "ok 12\n";
} else {
    print "not ok 12\n";
}



if( $pair->hstart == 400 ) { 
    print "ok 13\n";
} else {
    print "not ok 13\n";
}


if( $pair->hend == 440 ) { 
    print "ok 14\n";
} else {
    print "not ok 14\n";
}


if( $pair->hprimary_tag eq 'other' ) { 
    print "ok 15\n";
} else {
    print "not ok 15\n";
}


if( $pair->hsource_tag eq 'program_a' ) { 
    print "ok 16\n";
} else {
    print "not ok 16\n";
}

$pair->invert;

if( $pair->start == 440 ) {
    print "ok 17\n";
} else {
    print "not ok 17\n";
}




