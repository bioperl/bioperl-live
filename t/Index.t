# -*-Perl-*-
# $Id$

use strict;
my $exit;
BEGIN {
   eval { require Test; };
   use vars qw($NUMTESTS);
   $NUMTESTS = 41;
   if ( $@ ) {
      use lib 't';
   }
   use Test;
   eval { require Bio::Index::Fasta;
	  require Bio::Index::Qual;
	  require Bio::Index::SwissPfam;
	  require Bio::Index::EMBL;
	  require Bio::Index::GenBank;
	  require Bio::Index::Swissprot;
	  require DB_File;
	  require Storable;
	  require File::Temp;
       };
   if ( $@ ) {
      warn("Module DB_File or Storable or File::Temp not installed - skipping tests");
      $exit = 1;
   }
   plan tests => $NUMTESTS;
}

exit(0) if $exit;

use Bio::Root::IO;
use Bio::DB::InMemoryCache;

eval { require Bio::DB::FileCache };

use vars qw ($dir);

($Bio::Root::IO::FILESPECLOADED && File::Spec->can('curdir') &&
($dir = File::Spec->curdir) ) || ($dir = `pwd`) || ($dir = '.');
chomp( $dir );

my $ind = Bio::Index::Fasta->new(-filename => 'Wibbl',
				 -write_flag => 1,
				 -verbose => 0);

$ind->make_index(Bio::Root::IO->catfile($dir,"t","data","multifa.seq"));
$ind->make_index(Bio::Root::IO->catfile($dir,"t","data","seqs.fas"));

ok ( -e "Wibbl" || -e "Wibbl.pag" );
my $seq = $ind->fetch('HSEARLOBE');
ok($seq->length,321);
$seq = $ind->fetch('HSMETOO');
ok($seq->length,134);
$seq = $ind->fetch('MMWHISK');
ok($seq->length,62);
$seq = $ind->fetch('gi|238775|bbs|65126');
ok($seq->length,70);

my $stream = $ind->get_PrimarySeq_stream();
$seq = $stream->next_seq;
ok ($seq->isa('Bio::PrimarySeqI'));

$ind = Bio::Index::Fasta->new(-filename => 'multifa_index',
			      -write_flag => 1,
			      -verbose => 0);
$ind->make_index(Bio::Root::IO->catfile($dir,"t","data","multifa.seq.qual"));

ok ( -e "multifa_index" );

$ind = Bio::Index::Qual->new(-filename => 'multifa_qual_index',
			     -write_flag => 1,
			     -verbose => 0);
$ind->make_index(Bio::Root::IO->catfile($dir,"t","data","multifa.seq.qual"));

ok ( -e "multifa_qual_index" );

ok ( defined($seq) && $seq->isa('Bio::SeqI'));
$seq = $ind->fetch('HSEARLOBE');
ok($seq->length,321);
$seq = $ind->fetch('HSMETOO');
ok($seq->length,134);
$seq = $ind->fetch('MMWHISK');
ok($seq->length,62);

$ind = Bio::Index::SwissPfam->new(-filename   => 'Wibbl2',
				  -write_flag =>1);
$ind->make_index(Bio::Root::IO->catfile($dir,"t","data","swisspfam.data"));

ok ( -e "Wibbl2" || -e "Wibbl2.pag" );

$ind = Bio::Index::EMBL->new(-filename   => 'Wibbl3',
			     -write_flag =>1);
$ind->make_index(Bio::Root::IO->catfile($dir,"t","data","test.embl"));
ok ( -e "Wibbl3" || -e "Wibbl3.pag" );
ok ($ind->fetch('AL031232')->length, 4870);

$ind = Bio::Index::Swissprot->new(-filename   => 'Wibbl4',
				  -write_flag =>1);
$ind->make_index(Bio::Root::IO->catfile($dir,"t","data","roa1.swiss"));
ok ( -e "Wibbl4" || -e "Wibbl4.pag" );
ok ($ind->fetch('P09651')->display_id(), 'ROA1_HUMAN');

my $gb_ind = Bio::Index::GenBank->new(-filename   => 'Wibbl5',
				      -write_flag =>1,
				      -verbose    => 0);
$gb_ind->make_index(Bio::Root::IO->catfile($dir,"t","data","roa1.genbank"));
ok ( -e "Wibbl5" || -e "Wibbl5.pag" );
$seq =$gb_ind->fetch('AI129902');
ok ($seq->length, 37);
ok ($seq->species->binomial, 'Homo sapiens');

my $cache = Bio::DB::InMemoryCache->new( -seqdb => $gb_ind );

ok ( $cache->get_Seq_by_id('AI129902') );

if (Bio::DB::FileCache->can('new')) {
   $cache = Bio::DB::FileCache->new(-seqdb => $gb_ind,
				    -keep  => 1,
				    -file  => 'filecache.idx');
   my $seq = $cache->get_Seq_by_id('AI129902');
   ok ( $seq);
   ok ( $seq->length, 37);
   ok ( lc($seq->seq()), 'ctccgcgccaactccccccaccccccccccacacccc');

   my ( $f1 ) = $seq->get_SeqFeatures();
   ok ( ($f1->each_tag_value('sex'))[0], 'female');
   ok ( ($f1->each_tag_value('lab_host'))[0], 'DH10B');
   my $species = $seq->species;
   ok( $species );
   ok( $species->binomial, 'Homo sapiens');
   ok( $species->species(), 'sapiens');
   ok( $species->genus(), 'Homo');
   ok ($species->common_name(), 'human');

   $cache = undef;
   $cache = Bio::DB::FileCache->new(-seqdb => $gb_ind,
				    -keep  => 0,
				    -file  => 'filecache.idx');
   $seq = $cache->get_Seq_by_id('AI129902');
   ok ( $seq);
   ok ( $seq->length, 37);
   ok ( lc($seq->seq()), 'ctccgcgccaactccccccaccccccccccacacccc');

   ( $f1 ) = $seq->get_SeqFeatures();
   ok ( ($f1->each_tag_value('sex'))[0], 'female');
   ok ( ($f1->each_tag_value('lab_host'))[0], 'DH10B');
   $species = $seq->species;
   ok( $species );
   ok( $species->binomial, 'Homo sapiens');
   ok( $species->species(), 'sapiens');
   ok( $species->genus(), 'Homo');
   ok ($species->common_name(), 'human');
} else {
   skip('Bio::DB::FileCache not loaded because one or more of Storable, DB_File or File::Temp not installed',1);
}

foreach my $root ( qw( Wibbl Wibbl2 Wibbl3 Wibbl4 Wibbl5 multifa_index multifa_qual_index ) ) {
   if( -e $root ) { unlink $root;}
   if( -e "$root.pag") { unlink "$root.pag";}
   if( -e "$root.dir") { unlink "$root.dir";}
}
