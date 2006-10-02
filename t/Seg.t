# Bio::Tools::Seg test script

use strict;
use Test;

BEGIN {	plan tests => 15 }

use Bio::Tools::Seg;
use File::Spec;

ok my $infile = File::Spec->catfile(qw(t data seg.out));
ok my $parser = Bio::Tools::Seg->new(-file=>$infile);

my @feat;
while ( my $feat = $parser->next_result ) {
  push @feat, $feat;
}

ok @feat == 3;

#>LBL_0012(32-46) complexity=2.47 (12/2.20/2.50)
#gdggwtfegwggppe

ok $feat[0]->seq_id, 'LBL_0012';
ok $feat[0]->start,  32;
ok $feat[0]->end,    46;
ok $feat[0]->score,  2.47;

#>LBL_0012(66-80) complexity=2.31 (12/2.20/2.50)
#kfssrasakavakks

ok $feat[1]->seq_id, 'LBL_0012';
ok $feat[1]->start,  66;
ok $feat[1]->end,    80;
ok $feat[1]->score,  2.31;

#>LBL_0012(123-138) complexity=2.31 (12/2.20/2.50)
#svivsqsqgvvkgvgv

ok $feat[2]->seq_id, 'LBL_0012';
ok $feat[2]->start,  123;
ok $feat[2]->end,    138;
ok $feat[2]->score,  2.31;
