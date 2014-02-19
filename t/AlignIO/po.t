# -*-Perl-*- Test Harness script for Bioperl
# $Id: po.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 11);
	
	use_ok('Bio::AlignIO::po');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# POA
# just skip on perl 5.6.0 and earlier as it causes a crash on 
# default perl with OS X 10.2
# fink perl 5.6.0 does not seem to have the problem
# can't figure out what it is so just skip for now
SKIP: {
	skip("skipping due to bug in perl 5.6.0 that comes with OS X 10.2", 10) unless ($^O ne 'darwin' || $] > 5.006);
	
	$str = Bio::AlignIO->new(
			  -file   => test_input_file('testaln.po'),
			  -format => 'po',
			  );
	isa_ok($str, 'Bio::AlignIO');
	$aln = $str->next_aln();
	isa_ok($aln,'Bio::Align::AlignI');
	is $aln->num_sequences, 6;
	
	# output is? i.e. does conversion from clustalw to po give the same alignment?
	$str = Bio::AlignIO->new(
		  '-file'   => test_input_file('testaln.clustalw'),
		  '-format' => 'clustalw');
	isa_ok($str,'Bio::AlignIO');
	$aln = $str->next_aln();
	isa_ok($aln,'Bio::Align::AlignI');
	$strout = Bio::AlignIO->new(
		 '-file'   => ">" . test_output_file(),
		 '-format' => 'po');
	$status = $strout->write_aln($aln);
	is $status, 1, "po output test";
	
	$str = Bio::AlignIO->new(
		 '-file'   => test_input_file('testaln.po'),
		 '-format' => 'po');
	isa_ok($str,'Bio::AlignIO');
	my $aln2 = $str->next_aln();
	isa_ok($aln2,'Bio::Align::AlignI');
	is $aln2->num_sequences, $aln->num_sequences;
	is $aln2->length, $aln->length;
}
