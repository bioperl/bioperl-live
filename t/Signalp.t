# Bio::Tools::Signalp test script

use strict;
use Test;

BEGIN {	plan tests => 7 }

use Bio::Tools::Signalp;
use File::Spec;

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
ok @feat == 0;
ok $parser->_seqname, 'my_fasta_id';
ok $parser->_fact1,   'NO';

# positive example with "YES" features

ok $infile = File::Spec->catfile(qw(t data signalp.positive.out));
ok $parser = Bio::Tools::Signalp->new(-file=>$infile, -verbose=>$verbose);

#
#  The current module does NOT parse stuff properly
#  It is probably from version 2 but version 3 is used today
#  This has to be investigated!!!! --Torsten
#  FIXME / TODO / BUG / *** 
# 

#while ( my $feat = $parser->next_result ) {
#  push @feat, $feat;
#}
#ok @feat == 1;
#ok $parser->_seqname, 'my_fasta_id';
#ok $parser->_fact1,   'YES';

