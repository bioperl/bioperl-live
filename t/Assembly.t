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

    $NUMTESTS = 13;
    plan tests => $NUMTESTS;

}

#syntax test

use Bio::Assembly::IO;
use Bio::Assembly::Scaffold;
use Bio::Assembly::Contig;
use Bio::Assembly::ContigAnalysis;

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
skip "should return a number", $sc->get_nof_contigs, undef ; 
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
