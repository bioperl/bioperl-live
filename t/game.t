# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    plan tests => 9;
    
    eval { require XML::Parser::PerlSAX; };
    if( $@ ) {
	print STDERR "XML::Parser::PerlSAX not loaded. This means game test cannot be executed. Skipping\n";
	foreach ( 1..9 ) {
	    skip(1,1);
	}
	exit(0);
    } 
}

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqIO::MultiFile;
use XML::Parser::PerlSAX;
use vars qw($DEBUG);

my $str = Bio::SeqIO->new('-file'=> 't/test.game', 
		       '-format' => 'game');
ok ($str);
my $seq = $str->next_primary_seq();
ok($seq);

ok ($seq->display_id(), 'AE003417' );
ok ($seq->id(), 'AE003417' );

my $str2 = Bio::SeqIO->new(-file=> 't/test.game', '-format' => 'game');
ok ($str2);

$seq = $str2->next_seq();
ok $seq;

$str2->write_seq($seq) if( $DEBUG);

ok ( $seq->id, 'AE003417', 
     'id was not AE003417 it was ' .$seq->id);

my @feats = $seq->all_SeqFeatures();

ok @feats, 5;
ok( $feats[0]->primary_tag, 'exon', 
    'primary tag was not exon it was ' . $feats[0]->primary_tag);

