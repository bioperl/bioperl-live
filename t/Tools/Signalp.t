# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 11);
	
    use_ok('Bio::Tools::Signalp');
}

# global setting

my $verbose = test_debug();

# shared variables

my $infile;
my $parser;
my @feat;

# negative example without "YES" features

ok $infile = test_input_file('signalp.negative.out');
ok $parser = Bio::Tools::Signalp->new(-file=>$infile, -verbose=>$verbose);

while ( my $feat = $parser->next_result ) {
  push @feat, $feat;
}
is @feat, 0;
is $parser->_seqname, 'my_fasta_id';
is $parser->_fact1,   'NO';

# positive example with "YES" features

ok $infile = test_input_file('signalp.positive.out');
ok $parser = Bio::Tools::Signalp->new(-file=>$infile, -verbose=>$verbose);

#
#  The current module does NOT parse stuff properly
#  It is probably from version 2 but version 3 is used today
#  This has to be investigated!!!! --Torsten
#  FIXME / TODO? / BUG / *** 
# 

while ( my $feat = $parser->next_result ) {
  push @feat, $feat;
}
is @feat , 1;
is $parser->_seqname, 'my_fasta_id';
is $parser->_fact1,   'YES';
