# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error;
use strict;
use vars qw($DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'};

BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    use vars qw($TESTCOUNT);
    $TESTCOUNT = 23;
    plan tests => $TESTCOUNT;
    
    $error  = 0;
    eval { require XML::Parser::PerlSAX; };
    if( $@ ) {
	print STDERR "XML::Parser::PerlSAX not loaded. This means game test cannot be executed. Skipping\n";
	foreach ( $Test::ntest..$TESTCOUNT ) {
	    skip('XML::Parser::PerlSAX installed',1);
	}
	$error = 1;
    } 
    # make sure we can load it, assuming that the prerequisites are really met

    if( $error == 0 ) {
	eval { require Bio::SeqIO::game; };
	if( $@ ) {
	    print STDERR "game.pm not loaded. This means game test cannot be executed. Skipping\n";
	    foreach ( $Test::ntest..$TESTCOUNT ) {
		skip('game.pm not loaded because XML::Writer not loaded',1);
	    }
	    $error = 1;
	} 
    }
}

if( $error == 1 ) {
    exit(0);
}

END{ 
    unlink('testgameout.game')
}
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqIO::MultiFile;
use Bio::Root::IO;
my $verbose = $DEBUG ? 1 : -1;
my $str = Bio::SeqIO->new('-file'=> Bio::Root::IO->catfile("t","data","test.game"), 
			  '-format' => 'game',
			  '-verbose' => $verbose);
ok ($str);
my $seq = $str->next_primary_seq();
ok($seq);

ok ($seq->display_id(), 'AE003417' );
ok ($seq->id(), 'AE003417' );
ok ($seq->alphabet, 'dna');
ok ($seq->desc, 'E003417|Drosophila melanogaster genomic scaffold 142000013386054 section 1 of 35, complete sequence.|AE003417.1 GI:7290018');
my $str2 = Bio::SeqIO->new(-file=> Bio::Root::IO->catfile("t","data","test.game"), 
			   '-format' => 'game',
			   '-verbose' => $verbose,
			   );
ok ($str2);

$seq = $str2->next_seq();
ok $seq;
my @features = ( [qw(exon FBan0003038 fCG3038:1 2565 2154 -1)],
		 [qw(exon FBan0003038 fCG3038:2 2077 1181 -1)] );
foreach my $f ( $seq->all_SeqFeatures() ) {
    last unless @features;
    my $index = 0;
    ok($f->primary_tag, $features[0]->[$index++]);
    ok(($f->each_tag_value('annotation_id'))[0], $features[0]->[$index++]);
    ok(($f->each_tag_value('id'))[0],$features[0]->[$index++]);
    ok($f->end,    $features[0]->[$index++]);
    ok($f->start,  $features[0]->[$index++]);
    ok($f->strand, $features[0]->[$index++]);
    shift @features;
}

$str2 = Bio::SeqIO->new('-format'  => 'game', 
			'-file'    => '>testgameout.game', 
			'-verbose' => $verbose);
$str2->write_seq($seq);

ok ( $seq->id, 'AE003417');

my @feats = $seq->all_SeqFeatures();

ok @feats, 5;
ok( $feats[0]->primary_tag, 'exon');

