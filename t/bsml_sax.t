# -*-Perl-*-
# Bioperl Test Harness Script for Modules
# $Id$
#
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($error $NUMTESTS);
BEGIN {
  $NUMTESTS = 15;
  $error = 0;
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if ( $@ ) {
	use lib 't';
    }
    # bsml_sax uses XML::SAX
    eval {require XML::SAX;
	  require XML::SAX::Writer;
	  require XML::SAX::Base;
	  1;
      };
    if ( $@ ) {
	$error = 1;
	warn "XML::SAX::Base or XML::SAX or XML::SAX::Writer not found - skipping bsml_sax tests\n";
   } 
    use Test;
    plan tests => $NUMTESTS;
}

END { 
   foreach ( $Test::ntest..$NUMTESTS) {
      skip('Unable to run BSML_sax tests because XML::SAX is not installed',1);
   }
}


if ( $error == 1 ) {
  exit(0);
}


use Bio::SeqIO;
use Bio::Root::IO;

my $verbose = $ENV{'BIOPERLDEBUG'};
ok(1);

my $str = Bio::SeqIO->new(-format => 'bsml_sax',
			  -verbose => $verbose,
			  -file => Bio::Root::IO->catfile
			  (qw(t data U83300.bsml) ));
ok(my $seq = $str->next_seq);
my @refs = $seq->annotation->get_Annotations('reference');
ok(@refs, 2);
ok($seq->display_id,'MIVN83300');
ok($seq->molecule ,'dna');
ok(! $seq->is_circular);
ok($seq->get_dates,2);
ok($seq->accession_number, 'U83300');
ok($seq->seq_version,1);
my @feats = $seq->get_SeqFeatures;
ok(@feats, 2);
ok($feats[1]->start, 1);
ok($feats[1]->end, 946);
ok($feats[1]->get_tag_values('db_xref'), 3);
ok($seq->annotation->get_Annotations('reference'),2);
ok($seq->annotation->get_Annotations('dblink'),2);

