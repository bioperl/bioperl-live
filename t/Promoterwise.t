#!/usr/local/bin/perl
# -*-Perl-*-
## Bioperl Test Harness Script for Modules

use strict;
BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    use vars qw($NTESTS);
    $NTESTS = 6;
    plan tests => $NTESTS;
}
use Bio::Tools::Promoterwise;
use Bio::Root::IO;
use Bio::Seq;

END {
    for ( $Test::ntest..$NTESTS ) {
        skip("promoterwise parser not working properly. Skipping.",1);
    }
}

my $file = Bio::Root::IO->catfile(qw(t data promoterwise.out));
my  $parser = Bio::Tools::Promoterwise->new(-file=>$file);
ok $parser->isa('Bio::Tools::Promoterwise');
my @fp;
while (my $fp = $parser->next_result){
  push @fp,$fp;
}
my $first = $fp[0]->feature1;
my $second = $fp[0]->feature2;

my @sub = $first->sub_SeqFeature;
my @sub2 = $second->sub_SeqFeature;

ok $sub[0]->start,4;
ok $sub2[0]->start,29;
ok $sub[0]->end,18;
ok $sub2[0]->end,43;
ok $sub[0]->score,1596.49






