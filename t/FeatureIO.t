#$Id$
use strict;
use vars qw($DEBUG $TESTCOUNT);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

BEGIN {
  eval { require Test; };
  if( $@ ) {
    use lib 't';
  }
  use Test;

  $TESTCOUNT = 22;
  plan tests => $TESTCOUNT;
}

use Bio::FeatureIO;
use Data::Dumper;
ok(1);

my $io;
my $f;
my $s;
my $fcount;
my $scount;

################################################################################
#
# use FeatureIO::gff to read a FASTA file.
#
$fcount = 0;
$scount = 0;

ok( $io = Bio::FeatureIO->new( -file => Bio::Root::IO->catfile('t','data','dna1.fa') ) );

#read features
while($f = $io->next_feature()){
warn $f;
  $fcount++;
}
ok($fcount == 0);

#then try to read sequences again.  should get seqs now
while($s = $io->next_seq()){
  $scount++;
}
ok($scount == 1);

################################################################################
#
# use FeatureIO::gff to read a GFF3 file.
#
$fcount = 0;
$scount = 0;

ok( $io = Bio::FeatureIO->new( -file => Bio::Root::IO->catfile('t','data','knownGene.gff3') ) );

#try to read sequences first.  should be undef
while($s = $io->next_seq()){
  $scount++;
}
ok($scount == 0);

#then read features
while($f = $io->next_feature()){
  $fcount++;
}
ok($fcount == 14);

#then try to read sequences again.  should still be undef
while($s = $io->next_seq()){
  $scount++;
}
ok($scount == 0);

################################################################################
#
# use FeatureIO::gff to read a GFF3 file w/ directivized FASTA tail
#
$fcount = 0;
$scount = 0;

ok( $io = Bio::FeatureIO->new( -file => Bio::Root::IO->catfile('t','data','hybrid1.gff3') ) );

#try to read sequences first.  should be undef
while($s = $io->next_seq()){
  $scount++;
}
ok($scount == 0);

#then read features
while($f = $io->next_feature()){
  $fcount++;
}
ok($fcount == 6);

#then try to read sequences again.
while($s = $io->next_seq()){
  $scount++;
}
ok($scount == 1);

################################################################################
#
# use FeatureIO::gff to read a GFF3 file w/ non-directivized FASTA tail
#
$fcount = 0;
$scount = 0;

ok( $io = Bio::FeatureIO->new( -file => Bio::Root::IO->catfile('t','data','hybrid2.gff3') ) );

#try to read sequences first.  should be undef
while($s = $io->next_seq()){
  $scount++;
}
ok($scount == 0);

#then read features
while($f = $io->next_feature()){
  $fcount++;
}
ok($fcount == 6);

################################################################################
#
# use FeatureIO::gff to read a GFF3 file of directives
#
$fcount = 0;
$scount = 0;

ok( $io = Bio::FeatureIO->new( -file => Bio::Root::IO->catfile('t','data','directives.gff3') ) );

#read features
while($f = $io->next_feature()){
  $fcount++;
}
ok($fcount == 1); #sequence-region

################################################################################
#
# use FeatureIO::gff to read a GFF3 file as aggregated feature groups
#
$fcount = 0;
$scount = 0;

ok( $io = Bio::FeatureIO->new( -file => Bio::Root::IO->catfile('t','data','hybrid1.gff3') ) );

#try to read sequences first.  should be undef
while($s = $io->next_seq()){
  $scount++;
}
ok($scount == 0);

#read feature groups
$f = $io->next_feature_group();
ok($f == 1);
$f = $io->next_feature_group();
ok($f == 0);

#then try to read sequences again.
while($s = $io->next_seq()){
  $scount++;
}
ok($scount == 1);
