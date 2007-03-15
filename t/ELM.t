# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'
use strict;
use vars qw($NUMTESTS $DEBUG $ERROR);

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test::More; };

    $DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;

    $NUMTESTS = 17;

    eval {
	require IO::String; 
	require LWP::UserAgent;
	require HTML::HeadParser
    }; 
	if ($@) {
		plan skip_all => "IO::String or LWP::UserAgent not installed. This means that the module is not usable. Skipping tests";
	}
    elsif (!$DEBUG) {
		plan skip_all => 'Skipping all tests since they require network access, set BIOPERLDEBUG=1 to test';
	}
	else {
		plan tests => $NUMTESTS;
	}	
	require_ok('Bio::Tools::Analysis::Protein::ELM');
	use_ok('Bio::SeqIO');
	use_ok('Bio::PrimarySeq');
	require_ok('Bio::WebAgent');
}

use Data::Dumper;

ok 1;

my $verbose = 0;
$verbose = 1 if $DEBUG;

ok my $tool = Bio::WebAgent->new(-verbose =>$verbose);

my $seqio=new Bio::SeqIO( -verbose => $verbose,
                  -format => 'swiss',
                  -file   => Bio::Root::IO->catfile('t','data', 'swiss.dat'));

my $seq = $seqio->next_seq();
ok $tool = Bio::Tools::Analysis::Protein::ELM->new( 
					-seq=>$seq->primary_seq);
ok $tool->compartment(['golgi', 'er']);
ok my $cmp = $tool->compartment();
is $cmp->[1], 'GO:0005783';
ok $tool->species(9606);
is $tool->species, 9606;

ok $tool->run ();
exit if $tool->status eq 'TERMINATED_BY_ERROR';
ok my $raw = $tool->result('');
print $raw if $verbose;
ok my $parsed = $tool->result('parsed');
is $parsed->{'CLV_NDR_NDR_1'}{'locus'}[0], '54-56';
ok my @res = $tool->result('Bio::SeqFeatureI');
