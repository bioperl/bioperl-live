# -*-Perl-*-
# $Id$
use strict;
use vars qw($DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
use constant NUMTESTS => 33;
my $error;
BEGIN {
  eval { require Test; };
  if( $@ ) {
	  use lib 't';
  }
  $error = 0;
  use Test;
  plan tests => NUMTESTS;
  unless( eval "require Graph; require Bio::FeatureIO; 1;" ) {
	  warn("Graph not installed.  Bio::FeatureIO is not installed.\n");
	  $error = 1;
	  for ( 1..NUMTESTS ) {
		  skip("Graph is not installed. Bio::FeatureIO::gff cannot be run",1);
	  }
  }
}

if( $error ==  1 ) {
	exit(0);
}
END {
	foreach ( $Test::ntest..NUMTESTS) {
		skip('Cannot complete FeatureIO tests',1);
	}
}

use Bio::Root::IO;
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
ok($fcount == 15);

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
    ok( $f->start == 2491 );
    ok( $f->end == 3423 );
    ok( $f->strand > 0 );
    ok( ($f->get_tag_values('PID'))[0] eq '24217063' );
    ok( ($f->get_tag_values('Gene'))[0] eq 'metF' );
    ok( ($f->get_tag_values('Synonym'))[0] eq 'LB002' );
    ok( not $f->has_tag('Code') );
    ok( ($f->get_tag_values('COG'))[0] eq 'COG0685E' );
    ok( ($f->get_tag_values('Product'))[0] eq '5,10-methylenetetrahydrofolate reductase' );   
  }
}
ok($fcount == 367);
