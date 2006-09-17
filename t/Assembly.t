# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

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

    $NUMTESTS = 19;
    plan tests => $NUMTESTS;
    eval { require DB_File };
    if( $@ ) {
		 print STDERR "DB_File not installed. This means the Assembly modules are not available.  Skipping tests.\n";
		 for( 1..$NUMTESTS ) {
			 skip("DB_File not installed",1);
		 }
		 $error = 1; 
    }
}


END {
	foreach ( $Test::ntest..$NUMTESTS) {
		skip('unable to run all of the DB tests',1);
	}
}

if( $error ==  1 ) {
    exit(0);
}

#syntax test

require Bio::Assembly::IO;
require Bio::Assembly::Scaffold;
require Bio::Assembly::Contig;
require Bio::Assembly::ContigAnalysis;

use Data::Dumper;

ok 1;

#
# Testing IO
#

# -file => ">".Bio::Root::IO->catfile("t","data","primaryseq.embl")

ok my $in = Bio::Assembly::IO->new
    (-file=>Bio::Root::IO->catfile
     ("t","data","consed_project","edit_dir","test_project.phrap.out"));

ok my $sc = $in->next_assembly;
#print Dumper $sc;

#
# Testing Scaffold
#


ok $sc->id, "NoName";
ok $sc->id('test'), "test";

ok $sc->annotation;
skip "no annotations in Annotation collection?", $sc->annotation->get_all_annotation_keys, 0;
ok $sc->get_nof_contigs, 1;
ok $sc->get_nof_sequences_in_contigs, 2;
skip "should return a number", $sc->get_nof_singlets, 0;
skip $sc->get_seq_ids, 2;
skip $sc->get_contig_ids, 1;
skip "nothing to test", $sc->get_singlet_ids;


#
# Testing Contig
#

#
# Testing ContigAnalysis
#

#
# Testing Ace 
#

my $aio = Bio::Assembly::IO->new(
    -file=>Bio::Root::IO->catfile
     ("t","data","consed_project","edit_dir","test_project.fasta.screen.ace.2"),
    -format=>'ace',
);

my $assembly = $aio->next_assembly();
my @contigs = $assembly->all_contigs();

my $direction = $contigs[0]->strand;
ok $direction, 1;

my $features =  $contigs[0]->get_features_collection;
my @contig_features = $features->get_all_features;
ok @contig_features, 8;

my @annotations = grep {$_->primary_tag eq 'Annotation'} @contig_features;
ok @annotations, 2;
my $had_tag = 0;
foreach my $an (@annotations) {
	if ($an->has_tag('extra_info')) {
		$had_tag++;
		ok (($an->get_tag_values('extra_info'))[0], "contig extra\ninfo\n");
	}
	elsif ($an->has_tag('comment')){
		$had_tag++;
		ok (($an->get_tag_values('comment'))[0], "contig tag\ncomment\n");
	}
}
ok $had_tag, 2;