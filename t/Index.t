use strict;
BEGIN {     
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    plan tests => 11;
}

use Bio::Root::IO;
use Bio::Index::Fasta;
use Bio::Index::SwissPfam;
use Bio::Index::EMBL;
use Bio::Index::GenBank;
use Bio::Index::Swissprot;
use Bio::DB::InMemoryCache;
use vars qw ($dir);

($Bio::Root::IO::FILESPECLOADED && File::Spec->can('cwd') && ($dir = File::Spec->cwd) ) ||
    ($dir = `pwd`) || ($dir = '.');
 
END {  unlink qw( Wibbl Wibbl2 Wibbl3 Wibbl4 Wibbl5); }

chomp( $dir );
{
    my $ind = Bio::Index::Fasta->new(-filename => 'Wibbl', 
				     -write_flag => 1,
				     -verbose => 0);
    $ind->make_index(Bio::Root::IO->catfile($dir,"t","data","multifa.seq"));
    $ind->make_index(Bio::Root::IO->catfile($dir,"t","data","seqs.fas"));
    $ind->make_index(Bio::Root::IO->catfile($dir,"t","data","multi_1.fa"));
}

ok ( -e "Wibbl" );

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
    while( my $seq2 = $stream->next_primary_seq ) {
	unless ($seq2->isa('Bio::PrimarySeqI')) {
	    $ok_4 = 0;
	    last; # no point continuing...
	}
    }
    ok $ok_4;
}

{
    my $ind = Bio::Index::SwissPfam->new(-filename=>'Wibbl2', 
					 -write_flag=>1);
    $ind->make_index(Bio::Root::IO->catfile($dir,"t","data","swisspfam.data"));
    ok ( -e "Wibbl2" );
}

{
    my $ind = Bio::Index::EMBL->new(-filename=>'Wibbl3', 
				    -write_flag=>1);
    $ind->make_index(Bio::Root::IO->catfile($dir,"t","data","test.embl"));
    ok ( -e "Wibbl3" );
    ok $ind->fetch('AL031232')->length, 4870;
}

{
    my $ind = Bio::Index::Swissprot->new(-filename=>'Wibbl4', 
				    -write_flag=>1);
    $ind->make_index(Bio::Root::IO->catfile($dir,"t","data","roa1.swiss"));
    ok ( -e "Wibbl4" );
    ok ($ind->fetch('P09651')->display_id(), 'ROA1_HUMAN');
}

my $gb_ind;
{
    $gb_ind = Bio::Index::GenBank->new(-filename=>'Wibbl5', 
				       -write_flag=>1, 
				       -verbose => 0);
    $gb_ind->make_index(Bio::Root::IO->catfile($dir,"t","data","roa1.genbank"));
    ok ( -e "Wibbl5" );
    ok ($gb_ind->fetch('AI129902')->length, 37);
}

my $cache = Bio::DB::InMemoryCache->new( -seqdb => $gb_ind );

ok ( $cache->get_Seq_by_id('AI129902') );




