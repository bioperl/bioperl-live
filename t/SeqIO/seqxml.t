#-*-perl-*-
# $Id$

use strict;

BEGIN {
    use Bio::Root::Test;
    test_begin(
        -tests            => 20,
        -requires_modules => [qw(XML::LibXML XML::LibXML::Reader)]
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
  	diag("libxml version: ", XML::LibXML::LIBXML_VERSION()); 
    }

    use_ok('Bio::SeqIO::seqxml');    # checks that your module is there and loads ok
    
    # read data
    ok(
        my $seq_stream = Bio::SeqIO->new(
            -file   => test_input_file("seqxml.xml"),
            -format => 'seqxml'
        ),
        'stream ok',
    );
    
    my $seq_obj = $seq_stream->next_seq;
    isa_ok($seq_obj, 'Bio::Seq');
    is($seq_obj->display_id, 'ENSG00000173402', 'display id');
    is($seq_obj->primary_id, 'ENSG00000173402', 'primary id');
    is($seq_obj->desc, 'dystroglycan 1', 'description');
    is($seq_obj->seq, 'AAGGC----UGAUGUC.....ACAU', 'sequence');
    is($seq_obj->length, 25, 'length');
    
    # species
    isa_ok($seq_obj->species, 'Bio::Species', 'species');
    is($seq_obj->species->species, 'Homo sapiens', 'species name');
    is($seq_obj->species->ncbi_taxid, '9606', 'NCBI tax id');

    # alternative IDs
    my @dblinks = $seq_obj->get_Annotations('dblink');
    my $dblink = shift @dblinks;
    isa_ok($dblink, 'Bio::Annotation::DBLink');
    is($dblink->database, 'GenBank', 'dblink source');
    is($dblink->primary_id, 'NM_004393', 'dblink ID');
    
    # properties
    my $expected = {
        'has_splice_variants' => undef,
        'prediction_method'   => 'manual curation',
    };
    my @annotations = $seq_obj->get_Annotations();
    foreach my $annot_obj (@annotations) {
        next if ($annot_obj->tagname eq 'dblink');
        isa_ok($annot_obj, 'Bio::Annotation::SimpleValue');
        if ($annot_obj->tagname eq 'has_splice_variants') {
            is($annot_obj->value, undef, 'boolean property');
        }
        elsif ($annot_obj->tagname eq 'prediction_method') {
            is($annot_obj->value, 'manual curation', 'property with value');
        }
    }
}
