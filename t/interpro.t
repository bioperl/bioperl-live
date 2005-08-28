# -*-Perl-*-
# Bioperl Test Harness Script for Modules
# $Id$
#
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error = 0;

use strict;
BEGIN {
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test; };
	if ( $@ ) {
		use lib 't';
	}
   # interpro uses XML::DOM
   eval {require XML::DOM::XPath};
   if ( $@ ) {
      $error = 1;
      warn "XML::DOM::XPath not found - skipping interpro tests\n";
   }
	use Test;
	plan tests => 9;
}

if ( $error == 1 ) {
	exit(0);
}

use Bio::SeqIO;

my $verbose = $ENV{'BIOPERLDEBUG'};
ok(1);

my $t_file = Bio::Root::IO->catfile("t","data","test.interpro");
my $a_in = Bio::SeqIO->new( -file => $t_file,
									 -verbose => $verbose,
									 -format => 'interpro');

my $seq = $a_in->next_seq();
ok($seq);
ok($seq->isa('Bio::Seq::RichSeq'));
ok(scalar( $seq->get_SeqFeatures() ) == 6);

my($feat) = $seq->get_SeqFeatures();
ok($feat->isa('Bio::SeqFeature::Generic'));

ok($feat->display_name eq 'Retinoblastoma-associated protein, B-box');

ok($seq = $a_in->next_seq());
ok(scalar( $seq->get_SeqFeatures() ) == 40);

ok(!($seq = $a_in->next_seq()));
