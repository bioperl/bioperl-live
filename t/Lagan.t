use strict;
use constant NUMTESTS => 6;
BEGIN {
  eval { require Test; };
  if ($@) {
     use lib 't';
  }
  use Test;
  plan tests => NUMTESTS;
}
END {
  foreach ( $Test::ntest..NUMTESTS) {
	skip('Lagan or env variables not installed correctly',1);
  }
}

use Bio::Tools::Run::Alignment::Lagan;
use Bio::SeqIO;

my $verbose = -1;

ok(1);

my @params = (  'order' => "\"-gs -7 -gc -2 -mt 2 -ms -1\"",
		'match' => 12,
		'mismatch' => -8,
                'gapstart' => -50,
		'gapend' => -50,
		'gapcont' => -2  );

my $lagan = new Bio::Tools::Run::Alignment::Lagan(@params);

ok $lagan;

my $lagan_installed = $lagan->executable('lagan.pl');
if ( ! $lagan_installed ) {
   skip ("Lagan not installed", 1);
   exit;
} else {
   ok($lagan_installed);
}

my $mlagan_installed = $lagan->executable('mlagan');
if ( ! $mlagan_installed ) {
   skip ("mlagan not installed", 1);
   exit;
} else {
   ok($mlagan_installed);
}

my $in = Bio::SeqIO->new(-file=> Bio::Root::IO->catfile("t","data","dna2.fa"), '-format' => 'Fasta');

my $seq1 = $in->next_seq();
my $seq2 = $in->next_seq();

my $report_out = $lagan->lagan($seq1, $seq2);

ok($report_out);

my $tree = "((Test1 Test2) Test1)";

my @seq;
push @seq, $seq1;
push @seq, $seq2;
push @seq, $seq1;
my $seq_ref = \@seq;
bless $seq_ref, "ARRAY";

my $report_out = $lagan->mlagan($seq_ref, $tree);

ok($report_out);
