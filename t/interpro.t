# -*-Perl-*-
# Bioperl Test Harness Script for Modules
# $Id$
#
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN {
  use vars qw($error $NUMTESTS);
  $error = 0;
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
  $NUMTESTS = 17;
  use Test;
  plan tests => $NUMTESTS;
}

END { 
  foreach ( $Test::ntest..$NUMTESTS) {
    skip('Unable to run interpro tests because XML::DOM::XPath is not installed',1);
  }
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

# Bug 1908 (enhancement)
$t_file = Bio::Root::IO->catfile("t","data","interpro_ebi.xml");
my $b_in = Bio::SeqIO->new( -file => $t_file,
									 -verbose => $verbose,
									 -format => 'interpro');
$seq = $b_in->next_seq();
ok($seq);

my @features = $seq->get_SeqFeatures;
ok scalar @features,2;
ok $features[0]->primary_tag, 'region';
ok $features[0]->display_name,'Protein of unknown function DUF1021';
ok $features[0]->location->end,78;

my @dblinks = $features[0]->annotation->get_Annotations('dblink');
ok (scalar @dblinks,3);
ok $dblinks[1]->primary_id,'IPR009366';
ok $dblinks[2]->primary_id,'PF06257.1';

__END__

