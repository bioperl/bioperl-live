# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
   use lib '.';
   use Bio::Root::Test;
   
   test_begin(-tests => 73,
              -requires_modules => [qw(DB_File
                                       Storable
                                       Fcntl)]);
   
   use_ok('Bio::Index::Fasta');
   use_ok('Bio::Index::Qual');
   use_ok('Bio::Index::SwissPfam');
   use_ok('Bio::Index::EMBL');
   use_ok('Bio::Index::GenBank');
   use_ok('Bio::Index::Stockholm');
   use_ok('Bio::Index::Swissprot');
   use_ok('Bio::Index::Hmmer');
   use_ok('Bio::DB::InMemoryCache');
   use_ok('Bio::DB::InMemoryCache');
}

my $ind = Bio::Index::Fasta->new(-filename => 'Wibbl',
                                 -write_flag => 1,
                                 -verbose => 0);
$ind->make_index(test_input_file('multifa.seq'));
$ind->make_index(test_input_file('seqs.fas'));

ok ( -e "Wibbl" || -e "Wibbl.pag" );
my $seq = $ind->fetch('HSEARLOBE');
is($seq->length,321);
is($seq->primary_id(),'HSEARLOBE');
$seq = $ind->fetch('HSMETOO');
is($seq->length,134);
is($seq->primary_id(),'HSMETOO');
$seq = $ind->fetch('MMWHISK');
is($seq->length,62);
is($seq->primary_id(),'MMWHISK');
$seq = $ind->fetch('gi|238775|bbs|65126');
is($seq->length,70);

my $stream = $ind->get_PrimarySeq_stream();
$seq = $stream->next_seq;
isa_ok $seq, 'Bio::PrimarySeqI';

$ind = Bio::Index::Fasta->new(-filename => 'multifa_index',
                              -write_flag => 1,
                              -verbose => 0);
$ind->make_index(test_input_file('multifa.seq.qual'));

ok ( -e "multifa_index" );

$ind = Bio::Index::Qual->new(-filename => 'multifa_qual_index',
                             -write_flag => 1,
                             -verbose => 0);
$ind->make_index(test_input_file('multifa.seq.qual'));

ok ( -e "multifa_qual_index" );

ok ( defined($seq) );
isa_ok $seq, 'Bio::SeqI';
$seq = $ind->fetch('HSEARLOBE');
is($seq->length,321);
is($seq->primary_id(),'HSEARLOBE');
$seq = $ind->fetch('HSMETOO');
is($seq->length,134);
is($seq->primary_id(),'HSMETOO');
$seq = $ind->fetch('MMWHISK');
is($seq->length,62);
is($seq->primary_id(),'MMWHISK');
$seq = $ind->fetch('NONEXISTENT_SEQ');
ok(! defined $seq);

$ind = Bio::Index::SwissPfam->new(-filename => 'Wibbl2',
                                  -write_flag =>1);
$ind->make_index(test_input_file('swisspfam.data'));

ok ( -e "Wibbl2" || -e "Wibbl2.pag" );

$ind = Bio::Index::EMBL->new(-filename   => 'Wibbl3',
                             -write_flag =>1);
$ind->make_index(test_input_file('test.embl'));
ok ( -e "Wibbl3" || -e "Wibbl3.pag" );
is ($ind->fetch('AL031232')->length, 4870);

$ind = Bio::Index::Swissprot->new(-filename   => 'Wibbl4',
                                  -write_flag => 1);
$ind->make_index(test_input_file('roa1.swiss'));
ok ( -e "Wibbl4" || -e "Wibbl4.pag" );
$seq = $ind->fetch('ROA1_HUMAN');
is ($seq->display_id(), 'ROA1_HUMAN');
$seq = $ind->fetch('P09651');
is ($seq->display_id(), 'ROA1_HUMAN');

# test id_parser
$ind = Bio::Index::Swissprot->new(-filename   => 'Wibbl4',
                                  -write_flag => 1);
$ind->id_parser(\&get_id);
$ind->make_index(test_input_file('roa1.swiss'));
ok ( -e "Wibbl4" || -e "Wibbl4.pag" );
$seq = $ind->fetch('X12671');
is ($seq->length,371);


my $gb_ind = Bio::Index::GenBank->new(-filename => 'Wibbl5',
                                      -write_flag =>1,
                                      -verbose    => 0);
$gb_ind->make_index(test_input_file('roa1.genbank'));
ok ( -e "Wibbl5" || -e "Wibbl5.pag" );
$seq = $gb_ind->fetch('AI129902');
is ($seq->length, 37);
is ($seq->species->binomial, 'Homo sapiens');
$seq = $gb_ind->fetch(3598416);
is ($seq->seq,"CTCCGCGCCAACTCCCCCCACCCCCCCCCCACACCCC");

my $cache = Bio::DB::InMemoryCache->new( -seqdb => $gb_ind );

ok ( $cache->get_Seq_by_id('AI129902') );

SKIP: {
   test_skip(-tests => 22, -requires_module => 'Bio::DB::FileCache');

   $cache = Bio::DB::FileCache->new(-seqdb => $gb_ind,
                                    -keep  => 1,
                                    -file  => 'filecache.idx');
   # problem:
   my $seq = $cache->get_Seq_by_id('AI129902');
   ok ( $seq);
   is ( $seq->length, 37);
   is ( lc($seq->seq()), 'ctccgcgccaactccccccaccccccccccacacccc');

   my ( $f1 ) = $seq->get_SeqFeatures();
   is ( ($f1->get_tag_values('sex'))[0], 'female');
   is ( ($f1->get_tag_values('lab_host'))[0], 'DH10B');
   my $species = $seq->species;
   ok( $species );
   is( $species->binomial, 'Homo sapiens');
   is( $species->species(), 'sapiens');
   is( $species->genus(), 'Homo');
   # changes in GenBank file SOURCE line
   # this is now the abbreviated name
   ok defined($species->name('abbreviated'));
   is ($species->name('abbreviated')->[0], 'human');

   $cache = undef;
   $cache = Bio::DB::FileCache->new(-seqdb => $gb_ind,
                                    -keep  => 0,
                                    -file  => 'filecache.idx');
   $seq = $cache->get_Seq_by_id('AI129902');
   ok ( $seq);
   is ( $seq->length, 37);
   is ( lc($seq->seq()), 'ctccgcgccaactccccccaccccccccccacacccc');

   ( $f1 ) = $seq->get_SeqFeatures();
   is ( ($f1->get_tag_values('sex'))[0], 'female');
   is ( ($f1->get_tag_values('lab_host'))[0], 'DH10B');
   $species = $seq->species;
   ok( $species );
   is( $species->binomial, 'Homo sapiens');
   is( $species->species(), 'sapiens');
   is( $species->genus(), 'Homo');
   # changes in GenBank file SOURCE line
   # this is now the abbreviated name
   ok defined($species->name('abbreviated'));
   is ($species->name('abbreviated')->[0], 'human');
}

# test id_parser
$gb_ind = Bio::Index::GenBank->new(-filename => 'Wibbl5',
                                   -write_flag =>1,
                                   -verbose    => 0);
$gb_ind->id_parser(\&get_id);
$gb_ind->make_index(test_input_file('roa1.genbank'));
ok ( -e "Wibbl5" || -e "Wibbl5.pag" );
$seq = $gb_ind->fetch('alpha D-globin');
is ($seq->length,141);

# test Stockholm
my $st_ind = Bio::Index::Stockholm->new(-filename => 'Wibbl6',
                                        -write_flag => 1,
                                        -verbose    => 0);
isa_ok $st_ind, 'Bio::Index::Stockholm';
$st_ind->make_index(test_input_file('testaln.stockholm'));
ok ( -e "Wibbl6" );
my $aln = $st_ind->fetch_aln('PF00244');
isa_ok($aln,'Bio::SimpleAlign');

# test Hmmer
my $hmmer_ind = Bio::Index::Hmmer->new(-filename => 'Wibbl7',
                                       -write_flag => 1,
                                       -verbose    => 0);
isa_ok $hmmer_ind, 'Bio::Index::Hmmer';
$hmmer_ind->make_index(test_input_file('hmmpfam_multiresult.out'));
ok ( -e "Wibbl7" );
my $hmm_result = $hmmer_ind->fetch_report('lcl|gi|340783625|Plus1');
is ($hmm_result->query_description, 'megaplasmid, complete sequence [UNKNOWN]');




sub get_id {
   my $line = shift;
   return $1 if ($line =~ /product="([^"]+)"/);
   return $1 if ($line =~ /^DR\s+EMBL;\s+([^;]+)/);
}

END {
   cleanup();
}

sub cleanup {
   for my $root ( qw( Wibbl Wibbl2 Wibbl3 Wibbl4 Wibbl5 Wibbl6 Wibbl7
                      multifa_index multifa_qual_index ) ) {
      unlink $root if( -e $root );
      unlink "$root.pag" if( -e "$root.pag");
      unlink "$root.dir" if( -e "$root.dir");
   }
}
