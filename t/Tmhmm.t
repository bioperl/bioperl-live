use strict;
use Test;

BEGIN {	plan tests => 12 }

use Bio::Root::IO;
use Bio::Tools::Tmhmm;

ok my $infile = Bio::Root::IO->catfile(qw(t data tmhmm.out));
ok my $parser = Bio::Tools::Tmhmm->new(-file=>$infile);

my @feat;
while ( my $feat = $parser->next_result ) {
  push @feat, $feat;
}

ok @feat==3;

ok $feat[0]->seq_id,      'my_sequence_id';
ok $feat[0]->source_tag,  'TMHMM2.0';
ok $feat[0]->primary_tag, 'transmembrane';

ok $feat[0]->start,  54;
ok $feat[0]->end,    76;

ok $feat[1]->start,  116;
ok $feat[1]->end,    138;

ok $feat[2]->start,  151;
ok $feat[2]->end,    173;

  
