# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id: Signalp2.t,v 1.1 cjfields 2007/03/14 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use Data::Dumper;
my $DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
BEGIN {
	eval { require Test::More; };
	if( $@ ) {
		use lib 't/lib';
	}
	use Test::More;
	plan tests => 12;
    use_ok('Bio::Tools::Signalp');
    use_ok('Bio::Root::IO');
}

# global setting

my $verbose = $ENV{BIOPERLDEBUG} || 0;

# shared variables

my $infile;
my $parser;
my @feat;

# negative example without "YES" features

ok $infile = File::Spec->catfile(qw(t data signalp.negative.out));
ok $parser = Bio::Tools::Signalp->new(-file=>$infile, -verbose=>$verbose);

while ( my $feat = $parser->next_result ) {
  push @feat, $feat;
}
is @feat, 0;
is $parser->_seqname, 'my_fasta_id';
is $parser->_fact1,   'NO';

# positive example with "YES" features

ok $infile = File::Spec->catfile(qw(t data signalp.positive.out));
ok $parser = Bio::Tools::Signalp->new(-file=>$infile, -verbose=>$verbose);

#
#  The current module does NOT parse stuff properly
#  It is probably from version 2 but version 3 is used today
#  This has to be investigated!!!! --Torsten
#  FIXME / TODO / BUG / *** 
# 

while ( my $feat = $parser->next_result ) {
  push @feat, $feat;
}
is @feat , 1;
is $parser->_seqname, 'my_fasta_id';
is $parser->_fact1,   'YES';

