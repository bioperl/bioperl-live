# -*-Perl-*-
# $Id$

use strict;
my $exit;
BEGIN {     
    eval { require Test; };
    use vars qw($NUMTESTS);
    $NUMTESTS = 35;
    if( $@ ) {
	use lib 't';
    }
    use Test;
    $exit = 0;
    eval { require Bio::Index::Fasta;
           require Bio::Index::Qual;
	   require Bio::Index::SwissPfam;
	   require Bio::Index::EMBL;
	   require Bio::Index::GenBank;
	   require Bio::Index::Swissprot;
	   require DB_File;	   
       };
    if( $@ ) {
	$exit = 1;
    }
    plan tests => $NUMTESTS;
}
END { 
    foreach ( $Test::ntest..$NUMTESTS) {
	skip('DB_File not loaded because one or more of Storable, DB_File or File::Temp not installed',1);

    } 
    foreach my $root ( qw( Wibbl Wibbl2 Wibbl3 Wibbl4 Wibbl5 
			   ) ) {
	if( -e "$root" ) { unlink $root;}
	if( -e "$root.pag") { unlink "$root.pag";}
	if( -e "$root.dir") { unlink "$root.dir";}
   }
 }

exit(0) if $exit;

use Bio::Root::IO;
use Bio::DB::InMemoryCache;

eval { require Bio::DB::FileCache };

use vars qw ($dir);

($Bio::Root::IO::FILESPECLOADED && File::Spec->can('cwd') && 
 ($dir = File::Spec->cwd) ) ||
    ($dir = `pwd`) || ($dir = '.');

chomp( $dir );
{
    my $ind = Bio::Index::Fasta->new(-filename => 'Wibbl', 
				     -write_flag => 1,
				     -verbose => 0);
    $ind->make_index(Bio::Root::IO->catfile($dir,"t","data","multifa.seq"));
    $ind->make_index(Bio::Root::IO->catfile($dir,"t","data","seqs.fas"));
    $ind->make_index(Bio::Root::IO->catfile($dir,"t","data","multi_1.fa"));
}

ok ( -e "Wibbl" || -e "Wibbl.pag" );

{
    my %t_seq = (
        HSEARLOBE               => 321,
        HSMETOO                 => 134,
        MMWHISK                 => 62,
        'gi|238775|bbs|65126'   => 70,
    );

    my $ind = Bio::Index::Abstract->new(-FILENAME => 'Wibbl');

    my $ok_3 = 1;
    while (my($name, $length) = each %t_seq) {
        my $seq = $ind->fetch($name);
        if( defined($seq) and $seq->isa('Bio::SeqI') ) {
            my $r_length = $seq->length;
	    unless ($r_length == $length) {
                warn "$name - retrieved length '$r_length' doesn't match known length '$length'\n";
                $ok_3 = 0;
            }
        } else {
            warn "Didn't get sequence '$name' from index\n";
            $ok_3 = 0;
        }
    }
    ok $ok_3;

    my $stream = $ind->get_PrimarySeq_stream();
    my $ok_4 = 1;
    while( my $seq2 = $stream->next_seq ) {
	unless ($seq2->isa('Bio::PrimarySeqI')) {
	    $ok_4 = 0;
	    last; # no point continuing...
	}
    }
    ok $ok_4;
}

{
    my $ind = Bio::Index::Fasta->new(-filename => 'multifa_index',
				     -write_flag => 1,
				     -verbose => 0);
    $ind->make_index(Bio::Root::IO->catfile($dir,"t","data","multifa.seq.qual"));
}

ok ( -e "multifa_index" );

{
    my $ind = Bio::Index::Qual->new(-filename => 'multifa_qual_index',
                                    -write_flag => 1,
                                    -verbose => 0);
    $ind->make_index(Bio::Root::IO->catfile($dir,"t","data","multifa.seq.qual"));
}

ok ( -e "multifa_qual_index" );

{
    my %t_seq = (
        HSEARLOBE               => 321,
        HSMETOO                 => 134,
        MMWHISK                 => 62,
    );

    my $ind_f = Bio::Index::Fasta->new(-FILENAME => 'multifa_index');
    my $ind_q = Bio::Index::Qual->new(-FILENAME => 'multifa_qual_index');

    my $ok_3 = 1;
    my $ok_4 = 1;
    while (my($name, $length) = each %t_seq) {
        my $seq = $ind_f->fetch($name);
	my $qual = $ind_q->fetch($name);
        if( defined($seq) and $seq->isa('Bio::SeqI') ) {
            my $r_length = $seq->length;
	    unless ($r_length == $length) {
                warn "$name - retrieved seq length '$r_length' doesn't match known length '$length'\n";
                $ok_3 = 0;
            }
        } else {
            warn "Didn't get sequence '$name' from fasta index\n";
            $ok_3 = 0;
        }
	ok $ok_3;
	if( defined($qual) and $qual->isa('Bio::Seq::QualI') ) {
            my $r_length = $qual->length;
	    unless ($r_length == $length) {
                warn "$name - retrieved qual length '$r_length' doesn't match known length '$length'\n";
                $ok_4 = 0;
            }
        } else {
            warn "Didn't get sequence '$name' from qual index\n";
            $ok_4 = 0;
        }

	ok $ok_4;
    }
}

{
    my $ind = Bio::Index::SwissPfam->new(-filename=>'Wibbl2', 
					 -write_flag=>1);
    $ind->make_index(Bio::Root::IO->catfile($dir,"t","data","swisspfam.data"));
    ok ( -e "Wibbl2" || -e "Wibbl2.pag" );
}

{
    my $ind = Bio::Index::EMBL->new(-filename=>'Wibbl3', 
				    -write_flag=>1);
    $ind->make_index(Bio::Root::IO->catfile($dir,"t","data","test.embl"));
    ok ( -e "Wibbl3" || -e "Wibbl3.pag" );
    ok $ind->fetch('AL031232')->length, 4870;
}

{
    my $ind = Bio::Index::Swissprot->new(-filename=>'Wibbl4', 
				    -write_flag=>1);
    $ind->make_index(Bio::Root::IO->catfile($dir,"t","data","roa1.swiss"));
    ok ( -e "Wibbl4" || -e "Wibbl4.pag" );
    ok ($ind->fetch('P09651')->display_id(), 'ROA1_HUMAN');
}

my $gb_ind;
{
    $gb_ind = Bio::Index::GenBank->new(-filename=>'Wibbl5', 
				       -write_flag=>1, 
				       -verbose => 0);
    $gb_ind->make_index(Bio::Root::IO->catfile($dir,"t","data","roa1.genbank"));
    ok ( -e "Wibbl5" || -e "Wibbl5.pag" );
    my $seq =$gb_ind->fetch('AI129902'); 
    ok ($seq->length, 37);
    ok ($seq->species->binomial, 'Homo sapiens');
}

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
