# -*-Perl-*-
# $Id$
use strict;
use vars qw($DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
use constant NUMTESTS => 34;
my $error;
BEGIN {
  eval { require Test::More; };
  if( $@ ) {
	  use lib 't/lib';
  }
  use Test::More;
  eval {require Graph;};
  if ($@) {
	plan skip_all => "Graph is not installed. Bio::FeatureIO::gff cannot be run";
  } else {
	plan tests => NUMTESTS;
  }
  require_ok('Bio::FeatureIO');
  use_ok('Bio::Root::IO');
}

use Data::Dumper;

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
is($fcount, 0);

#then try to read sequences again.  should get seqs now
while($s = $io->next_seq()){
  $scount++;
}
is($scount,  1);

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
is($scount,0);

#then read features
while($f = $io->next_feature()){
  $fcount++;
}
is($fcount, 15);

#then try to read sequences again.  should still be undef
while($s = $io->next_seq()){
  $scount++;
}
is($scount,0);

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
is($scount , 0);

#then read features
while($f = $io->next_feature()){
  $fcount++;
}
is($fcount , 6);

#then try to read sequences again.
while($s = $io->next_seq()){
  $scount++;
}
is($scount , 1);

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
is($scount , 0);

#then read features
while($f = $io->next_feature()){
  $fcount++;
}
is($fcount , 6);

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
is($fcount , 1); #sequence-region

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
is($scount , 0);

#read feature groups
$f = $io->next_feature_group();
is($f , 1);
$f = $io->next_feature_group();
is($f , 0);

#then try to read sequences again.
while($s = $io->next_seq()){
  $scount++;
}
is($scount , 1);

################################################################################
#
# use FeatureIO::gff to read a PTT file.
#
$fcount = 0;

my $ptt_in = Bio::FeatureIO->new(
  -file => Bio::Root::IO->catfile('t','data','test.ptt'), 
  -format => 'ptt',
);
ok($ptt_in);

while (my $f = $ptt_in->next_feature) {
  $fcount++;
  if ($fcount==2) {
    # 2491..3423  + 310 24217063  metF  LB002 - COG0685E  5,10-methylenetetrahydrofolate reductase
    is( $f->start , 2491 );
    is( $f->end , 3423 );
    cmp_ok( $f->strand, '>', 0 );
    is( ($f->get_tag_values('PID'))[0],'24217063' );
    is( ($f->get_tag_values('Gene'))[0], 'metF' );
    is( ($f->get_tag_values('Synonym'))[0], 'LB002' );
    ok( ! $f->has_tag('Code') );
    is( ($f->get_tag_values('COG'))[0],'COG0685E' );
    is( ($f->get_tag_values('Product'))[0], '5,10-methylenetetrahydrofolate reductase' );   
  }
}
is($fcount , 367);
