# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS);

my $error;

BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    $error = 0;
    if( $@ ) {
	use lib 't';
    }
    use Test;

    $NUMTESTS = 12;
    plan tests => $NUMTESTS;

    eval { require DB_File; };
    if( $@ ) {
	print STDERR "DB_File not installed. This means the SeqFeatCollection wont work\n";
	for( 1..$NUMTESTS ) {
	    skip("DB_File",1);
	}
       $error = 1; 
    }

}

if( $error ==  1 ) {
    exit(0);
}
my $testnum;
my $verbose = $ENV{'BIOPERLDEBUG'} || 0;

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 


#First of all we need to create an flat db
use Bio::SeqFeature::Collection;
use Bio::Root::IO;
use Bio::Location::Simple;
use Bio::Tools::GFF;
use Bio::SeqIO;

my $simple = new Bio::SeqIO(-format => 'genbank',
			    -file   =>  Bio::Root::IO->catfile
			    ("t","data","AB077698.gb"));

my @features;
my $seq = $simple->next_seq();
@features = $seq->top_SeqFeatures();
ok(scalar @features, 11);

my $col = new Bio::SeqFeature::Collection(-verbose => $verbose);

ok($col);
ok($col->add_features( \@features), 11);
my @feat = $col->features_in_range(-range => ( new Bio::Location::Simple
						   (-start => 100,
						    -end   => 300,
						    -strand => 1) ),
				       -contain => 0);
ok(scalar @feat, 5);
if( $verbose ) {    
    foreach my $f ( @feat ) {
	print $f->location->to_FTstring(), "\n";    	
    }
}

ok(scalar $col->features_in_range(-range => ( new Bio::Location::Simple
						   (-start => 100,
						    -end   => 300,
						    -strand => -1) ),
				      -strandmatch => 'ignore',
				      -contain => 1), 2);

@feat = $col->features_in_range(-start => 79,
				-end   => 1145,
				-strand => 1,
				-strandmatch => 'strong',
				-contain => 1);
ok(scalar @feat, 5);
if( $verbose ) {    
    foreach my $f ( sort { $a->start <=> $b->start} @feat ) {
	print $f->primary_tag, " ", $f->location->to_FTstring(), "\n";
    }
}

ok($feat[0]->primary_tag, 'CDS');
ok($feat[0]->has_tag('gene'));

$verbose = 0;
# specify input via -fh or -file
my $gffio = Bio::Tools::GFF->new(-file => Bio::Root::IO->catfile
				 ("t","data","myco_sites.gff"), 
				 -gff_version => 2);
@features = ();
# loop over the input stream
while(my $feature = $gffio->next_feature()) {
    # do something with feature
    push @features, $feature;
}
$gffio->close();

ok(scalar @features, 412);
$col = new Bio::SeqFeature::Collection(-verbose => $verbose,
				       -usefile => 1);

ok($col);

ok($col->add_features( \@features), 412);

my $r = new Bio::Location::Simple(-start => 67700,
				  -end   => 150000,
				  -strand => 1);

@feat = $col->features_in_range(-range => $r,
				-strandmatch => 'ignore',
				-contain => 0);

ok(scalar @feat, 56);

if( $verbose ) {
    my @fts =  sort { $a->start <=> $b->start}  grep { $r->overlaps($_,'ignore') } @features;
    
    if( $verbose ) {
	foreach my $f ( @fts ) {
	    print $f->primary_tag, "    ", $f->location->to_FTstring(), "\n";
	}
	print "\n";
    }

    my %G = map { ($_,1) } @feat; 
    my $c = 0;
    foreach my $A ( @fts ) {
	if( ! $G{$A} ) {
	    print "missing ", $A->primary_tag, " ", $A->location->to_FTstring(), "\n";
	} else { 
	    $c++;
	}
    }
    print "Number of features correctly retrieved $c\n";
    foreach my $f ( sort { $a->start <=> $b->start} @feat ) {
	print $f->primary_tag, "    ", $f->location->to_FTstring(), "\n";
    }
}

