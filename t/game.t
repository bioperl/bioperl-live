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
use Bio::SeqIO;
use Bio::Root::IO;
my $verbose = $DEBUG ? 1 : -1;
my $str = Bio::SeqIO->new('-file'=> Bio::Root::IO->catfile("t","data","test.game"), 
			  '-format' => 'game',
			  '-verbose' => $verbose);
ok ($str);
my $seq = $str->next_seq();
ok($seq);

# exercise game parsing
$str = new Bio::SeqIO(
    -format =>'game',
    -file => Bio::Root::IO->catfile ( qw(t data test.game))
		      );
$seq = $str->next_seq;
ok(defined $seq);
ok(defined $seq->seq);
ok($seq->alphabet, 'dna');
ok($seq->display_id, 'L16622');
ok($seq->length, 28735);
ok($seq->species->binomial, 'Caenorhabditis elegans');
my @feats = $seq->get_SeqFeatures;
ok(scalar(@feats), 7);
my $source = grep { $_->primary_tag eq 'source' } @feats;
ok($source);
my @genes = grep { $_->primary_tag eq 'gene' } @feats;
ok(scalar(@genes), 3);
ok($genes[0]->has_tag('gene'));
my $gname;
if ( $genes[0]->has_tag('gene') ) {
    ($gname) = $genes[0]->get_tag_values('gene');
}
ok($gname, 'C02D5.3');
ok($genes[0]->strand, 1);
my $cds   = grep { $_->primary_tag eq 'CDS' } @feats;
ok($cds, 3);

# make sure we can read what we write
# test XML-writing
my $testfile = "testgameout.game";
# map argument is require to write a <map_position> element
my $out = new Bio::SeqIO(-format => 'game', -file => ">$testfile", -map => 1);
$out->write_seq($seq);
$out->close();

$str = new Bio::SeqIO(-format =>'game', -file => $testfile);
$seq = $str->next_seq;
ok(defined $seq);
ok(defined $seq->seq);
ok($seq->alphabet, 'dna');
ok($seq->display_id, 'L16622');
ok($seq->length, 28735);
ok($seq->species->binomial, 'Caenorhabditis elegans');

my $genes = grep { $_->primary_tag eq 'gene' } @feats;
$cds   = grep { $_->primary_tag eq 'CDS' } @feats;
ok($genes, 3);
ok($cds, 3);
unlink $testfile;

