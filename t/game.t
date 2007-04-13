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
    eval { require Test::More; };
    if( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
    use vars qw($TESTCOUNT);
    $TESTCOUNT = 26;
    
    eval { require XML::Parser::PerlSAX; };
    if( $@ ) {
        plan skip_all => "XML::Parser::PerlSAX not loaded. This means game test cannot be executed. Skipping";
    } else {
        plan tests => $TESTCOUNT;
    }
    # make sure we can load it, assuming that the prerequisites are really met
	use_ok('Bio::SeqIO::game');
    use_ok('Bio::SeqIO');
    use_ok('Bio::Root::IO');
}

END{ 
    unlink('testgameout.game')
}

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
is($seq->alphabet, 'dna');
is($seq->display_id, 'L16622');
is($seq->length, 28735);
is($seq->species->binomial, 'Caenorhabditis elegans');
my @feats = $seq->get_SeqFeatures;
is(scalar(@feats), 7);
my $source = grep { $_->primary_tag eq 'source' } @feats;
ok($source);
my @genes = grep { $_->primary_tag eq 'gene' } @feats;
is(scalar(@genes), 3);
ok($genes[0]->has_tag('gene'));
my $gname;
if ( $genes[0]->has_tag('gene') ) {
    ($gname) = $genes[0]->get_tag_values('gene');
}
is($gname, 'C02D5.3');
is($genes[0]->strand, 1);
my $cds   = grep { $_->primary_tag eq 'CDS' } @feats;
is($cds, 3);

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
is($seq->alphabet, 'dna');
is($seq->display_id, 'L16622');
is($seq->length, 28735);
is($seq->species->binomial, 'Caenorhabditis elegans');

my $genes = grep { $_->primary_tag eq 'gene' } @feats;
$cds   = grep { $_->primary_tag eq 'CDS' } @feats;
is($genes, 3);
is($cds, 3);
unlink $testfile;

