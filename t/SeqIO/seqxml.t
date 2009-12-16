#-*-perl-*-
# $Id$

use strict;

BEGIN {
    use Bio::Root::Test;
    test_begin(
        -tests            => 41,
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

    my $seq_obj = $seq_stream->next_seq;
    isa_ok( $seq_obj, 'Bio::Seq' );
    is( $seq_obj->display_id, 'ENSG00000173402',           'display id' );
    is( $seq_obj->primary_id, 'ENSG00000173402',           'primary id' );
    is( $seq_obj->desc,       'dystroglycan 1',            'description' );
    is( $seq_obj->seq,        'AAGGC----UGAUGUC.....ACAU', 'sequence' );
    is( $seq_obj->length,     25,                          'length' );

    # species
    isa_ok( $seq_obj->species, 'Bio::Species', 'species' );
    is( $seq_obj->species->species,    'Homo sapiens', 'species name' );
    is( $seq_obj->species->ncbi_taxid, '9606',         'NCBI tax id' );

    # alternative IDs
    my @dblinks = $seq_obj->get_Annotations('dblink');
    my $dblink  = shift @dblinks;
    isa_ok( $dblink, 'Bio::Annotation::DBLink' );
    is( $dblink->database,   'GenBank',   'dblink source' );
    is( $dblink->primary_id, 'NM_004393', 'dblink ID' );

    # properties
    my @annotations = $seq_obj->get_Annotations();
    foreach my $annot_obj (@annotations) {
        next if ( $annot_obj->tagname eq 'dblink' );
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
            -seqXMLversion => '0.1',
        ),
        'writer ok',
    );
    $seq_writer->flush;    # to make sure output is written to file
    ok( -s $outfile, 'outfile is created' );
    if ( $verbose > 0 ) {
        diag(`cat $outfile`);
    }

    # check metadata
    is( $seq_writer->seqXMLversion, '0.1',     'seqXML version' );
    is( $seq_writer->source,        'Ensembl', 'source' );
    is( $seq_writer->sourceVersion, '56',      'source version' );

    # write one sequence entry to file
    $seq_writer->write_seq($seq_obj);
    $seq_writer->close;

    # verify written data by roundtripping it
    {
        my $new_in = Bio::SeqIO->new(-file   => $outfile,
                                     -format => 'seqxml');

        my $new_seqobj = $new_in->next_seq;
        isa_ok( $new_seqobj, 'Bio::Seq' );
        is( $new_seqobj->display_id, 'ENSG00000173402',           'display id' );
        is( $new_seqobj->primary_id, 'ENSG00000173402',           'primary id' );
        is( $new_seqobj->desc,       'dystroglycan 1',            'description' );
        is( $new_seqobj->seq,        'AAGGC----UGAUGUC.....ACAU', 'sequence' );
        is( $new_seqobj->length,     25,                          'length' );

        # species
        isa_ok( $new_seqobj->species, 'Bio::Species', 'species' );
        is( $new_seqobj->species->species,    'Homo sapiens', 'species name' );
        is( $new_seqobj->species->ncbi_taxid, '9606',         'NCBI tax id' );

        # alternative IDs
        my @dblinks = $new_seqobj->get_Annotations('dblink');
        my $dblink  = shift @dblinks;
        isa_ok( $dblink, 'Bio::Annotation::DBLink' );
        is( $dblink->database,   'GenBank',   'dblink source' );
        is( $dblink->primary_id, 'NM_004393', 'dblink ID' );

        # properties
        my @annotations = $new_seqobj->get_Annotations();
        foreach my $annot_obj (@annotations) {
            next if ( $annot_obj->tagname eq 'dblink' );
            isa_ok( $annot_obj, 'Bio::Annotation::SimpleValue' );
            if ( $annot_obj->tagname eq 'has_splice_variants' ) {
                is( $annot_obj->value, undef, 'boolean property' );
            }
            elsif ( $annot_obj->tagname eq 'prediction_method' ) {
                is( $annot_obj->value, 'manual curation', 'property with value' );
            }
        }
    }

}
