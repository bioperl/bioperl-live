#-*-perl-*-
# $Id$

use strict;

BEGIN {
    use Bio::Root::Test;
    test_begin(
        -tests            => 61,
        -requires_modules => [qw(XML::LibXML XML::LibXML::Reader XML::Writer)]
    );

    use_ok('Bio::SeqIO');
}
use_ok('Bio::PrimarySeq');

my $verbose = test_debug();

SKIP: {

    # XML library version checks
    if ( 1000 * $] < 5008 ) {
        skip( "Reader interface only supported in Perl >= 5.8", 96 );
    }
    elsif ( XML::LibXML::LIBXML_VERSION() <= 20620 ) {
        skip( "Reader not supported for libxml2 <= 2.6.20", 96 );
    }

    if ($verbose) {
        diag( "libxml version: ", XML::LibXML::LIBXML_VERSION() );
    }

    # checks that your module is there and loads ok
    use_ok('Bio::SeqIO::seqxml');

    # read data
    ok(
        my $seq_stream = Bio::SeqIO->new(
            -file    => test_input_file("seqxml.xml"),
            -format  => 'seqxml',
            -verbose => $verbose,
        ),
        'stream ok',
    );

    # check metadata
    is( $seq_stream->seqXMLversion, '0.3',     'seqXML version' );
    is( $seq_stream->source,        'Ensembl', 'source' );
    is( $seq_stream->sourceVersion, '56',      'source version' );

    # now get and check the sequence entry itself
    my $seq_obj = $seq_stream->next_seq;
    isa_ok( $seq_obj, 'Bio::Seq' );
    is( $seq_obj->display_id, 'ENST00000308775',           'display id' );
    is( $seq_obj->primary_id, 'ENST00000308775',           'primary id' );
    is( $seq_obj->desc,       'dystroglycan 1',            'description' );
    is( $seq_obj->seq,        'AAGGC----UGAUGUC.....ACAU', 'sequence' );
    is( $seq_obj->length,     25,                          'length' );

    my ($source) = $seq_obj->get_Annotations('source');
    if ($source) { is($source->value, 'Ensembl', 'entry source'); }

    # species
    isa_ok( $seq_obj->species, 'Bio::Species', 'species' );
    is( $seq_obj->species->node_name,    'Homo sapiens', 'species name' );
    is( $seq_obj->species->ncbi_taxid, '9606',         'NCBI tax id' );

    # alternative IDs
    my @dblinks = $seq_obj->get_Annotations('dblink');
    my $dblink  = shift @dblinks;
    isa_ok( $dblink, 'Bio::Annotation::DBLink' );
    is( $dblink->database,   'RefSeq',   'dblink source' );
    is( $dblink->primary_id, 'NM_004393', 'dblink ID' );

    # properties
    my @annotations = $seq_obj->get_Annotations();
    foreach my $annot_obj (@annotations) {
        next if ( $annot_obj->tagname eq 'dblink' );
        next if ( $annot_obj->tagname eq 'source' );        
        isa_ok( $annot_obj, 'Bio::Annotation::SimpleValue' );
        if ( $annot_obj->tagname eq 'has_splice_variants' ) {
            is( $annot_obj->value, undef, 'boolean property' );
        }
        elsif ( $annot_obj->tagname eq 'prediction_method' ) {
            is( $annot_obj->value, 'manual curation', 'property with value' );
        }
    }

    # write data
    my $outfile = test_output_file();
    ok(
        my $seq_writer = Bio::SeqIO->new(
            -file          => ">$outfile",
            -format        => 'seqxml',
            -verbose       => $verbose,
            -source        => 'Ensembl',
            -sourceVersion => '56',
            -seqXMLversion => '0.3',
        ),
        'writer ok',
    );
    $seq_writer->flush;    # to make sure output is written to file
    ok( -s $outfile, 'outfile is created' );

    # check metadata
    is( $seq_writer->seqXMLversion, '0.3',     'seqXML version' );
    is( $seq_writer->source,        'Ensembl', 'source' );
    is( $seq_writer->sourceVersion, '56',      'source version' );
    is( $seq_writer->schemaLocation, 'http://www.seqxml.org/0.3/seqxml.xsd', 'schemaLocation' );

    # write one sequence entry to file
    $seq_writer->write_seq($seq_obj);
    $seq_writer->close;
    if ( $verbose > 0 ) {
        diag("writing first seqXML outfile");
        diag(`cat $outfile`);
    }

    # verify written data by roundtripping it
    {
        my $new_in = Bio::SeqIO->new(
            -file   => $outfile,
            -format => 'seqxml'
        );

        my $new_seqobj = $new_in->next_seq;
        isa_ok( $new_seqobj, 'Bio::Seq' );
        is( $new_seqobj->display_id, 'ENST00000308775', 'display id' );
        is( $new_seqobj->primary_id, 'ENST00000308775', 'primary id' );
        is( $new_seqobj->desc,       'dystroglycan 1',  'description' );
        is( $new_seqobj->seq, 'AAGGC----UGAUGUC.....ACAU', 'sequence' );
        is( $new_seqobj->length, 25, 'length' );

        my ($new_source) = $new_seqobj->get_Annotations('source');
        if ($new_source) { is($new_source->value, 'Ensembl', 'entry source'); }


        # species
        isa_ok( $new_seqobj->species, 'Bio::Species', 'species' );
        is( $new_seqobj->species->node_name,    'Homo sapiens', 'species name' );
        is( $new_seqobj->species->ncbi_taxid, '9606',         'NCBI tax id' );

        # alternative IDs
        my @dblinks = $new_seqobj->get_Annotations('dblink');
        my $dblink  = shift @dblinks;
        isa_ok( $dblink, 'Bio::Annotation::DBLink' );
        is( $dblink->database,   'RefSeq',   'dblink source' );
        is( $dblink->primary_id, 'NM_004393', 'dblink ID' );

        # properties
        my @annotations = $new_seqobj->get_Annotations();
        foreach my $annot_obj (@annotations) {
            next if ( $annot_obj->tagname eq 'dblink' );
            next if ( $annot_obj->tagname eq 'source' );
            isa_ok( $annot_obj, 'Bio::Annotation::SimpleValue' );
            if ( $annot_obj->tagname eq 'has_splice_variants' ) {
                is( $annot_obj->value, undef, 'boolean property' );
            }
            elsif ( $annot_obj->tagname eq 'prediction_method' ) {
                is(
                    $annot_obj->value,
                    'manual curation',
                    'property with value'
                );
            }
        }
    }

    # write data from a Seq object created from a fasta file
    {
        # forcing a Bio::Seq to be created
        # due to SeqIO::fasta creating a PrimarySeq by default
        # as of r16838
        my $factory = Bio::Seq::SeqFactory->new(-type => 'Bio::Seq');
        
        my $seq_stream = Bio::SeqIO->new(
            -file   => test_input_file("test.fasta"),
            -format => 'fasta',
            -seqfactory => $factory,
        );

        my $outfile = test_output_file();
        my $writer  = Bio::SeqIO->new(
            -file   => ">$outfile",
            -format => 'seqxml'
        );
        $writer->flush;
        ok( -s $outfile, 'outfile is created' );

        while ( my $seq_obj = $seq_stream->next_seq ) {
            $writer->write_seq($seq_obj);
        }
        $writer->close;
        if ( $verbose > 0 ) {
            diag(`cat $outfile`);
        }

        # now read that newly made seqxml back in
        my $in = Bio::SeqIO->new(
            -file   => $outfile,
            -format => 'seqxml'
        );

        # check header
        is( $in->seqXMLversion, '0.3', 'seqXML version' );
        is( $in->source,        undef, 'source' );
        is( $in->sourceVersion, undef, 'source version' );

        # check first sequence entry
        my $seqxml_obj = $in->next_seq;
        is( $seqxml_obj->display_id, 'roa1_drome', 'display id' );
        is( $seqxml_obj->primary_id, 'roa1_drome', 'primary id' );
        is( $seqxml_obj->desc, 'Rea guano receptor type III >> 0.1',
            'description' );
        is(
            $seqxml_obj->seq,
'MVNSNQNQNGNSNGHDDDFPQDSITEPEHMRKLFIGGLDYRTTDENLKAHEKWGNIVDVVVMKDPRTKRSRGFGFITYSHSSMIDEAQKSRPHKIDGRVEPKRAVPRQDIDSPNAGATVKKLFVGALKDDHDEQSIRDYFQHFGNIVDNIVIDKETGKKRGFAFVEFDDYDPVDKVVLQKQHQLNGKMVDVKKALPKNDQQGGGGGRGGPGGRAGGNRGNMGGGNYGNQNGGGNWNNGGNNWGNNRGNDNWGNNSFGGGGGGGGGYGGGNNSWGNNNPWDNGNGGGNFGGGGNNWNGGNDFGGYQQNYGGGPQRGGGNFNNNRMQPYQGGGGFKAGGGNQGNYGNNQGFNNGGNNRRY',
            'sequence'
        );
        is( $seqxml_obj->length, 358, 'length' );

        # check second sequence entry
        my $seqxml_obj2 = $in->next_seq;
        is( $seqxml_obj2->display_id, 'roa2_drome',       'display id' );
        is( $seqxml_obj2->primary_id, 'roa2_drome',       'primary id' );
        is( $seqxml_obj2->desc,       'Rea guano ligand', 'description' );
        is(
            $seqxml_obj2->seq,
'MVNSNQNQNGNSNGHDDDFPQDSITEPEHMRKLFIGGLDYRTTDENLKAHEKWGNIVDVVVMKDPTSTSTSTSTSTSTSTSTMIDEAQKSRPHKIDGRVEPKRAVPRQDIDSPNAGATVKKLFVGALKDDHDEQSIRDYFQHLLLLLLLDLLLLDLLLLDLLLFVEFDDYDPVDKVVLQKQHQLNGKMVDVKKALPKNDQQGGGGGRGGPGGRAGGNRGNMGGGNYGNQNGGGNWNNGGNNWGNNRGNDNWGNNSFGGGGGGGGGYGGGNNSWGNNNPWDNGNGGGNFGGGGNNWNGGNDFGGYQQNYGGGPQRGGGNFNNNRMQPYQGGGGFKAGGGNQGNYGNNQGFNNGGNNRRY',
            'sequence'
        );
        is( $seqxml_obj2->length, 358, 'length' );

    }

}
