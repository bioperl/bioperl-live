# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 175);
    
    use_ok('Bio::Phenotype::OMIM::OMIMparser');
}


my $omim_parser = Bio::Phenotype::OMIM::OMIMparser->new( -genemap  => test_input_file('omim_genemap_test'),
                                                         -omimtext => test_input_file('omim_text_test'));

isa_ok( $omim_parser, "Bio::Phenotype::OMIM::OMIMparser");

my $omim_entry = $omim_parser->next_phenotype();

isa_ok($omim_entry, "Bio::Phenotype::OMIM::OMIMentry");

is( $omim_entry->MIM_number(), "100500" );
is( $omim_entry->title(), "*100500 title" );
is( $omim_entry->alternative_titles_and_symbols(), ";;title1;;\ntitle2;;\ntitle3" );
is( $omim_entry->more_than_two_genes(), 0 );
is( $omim_entry->is_separate(), 1 );
is( $omim_entry->description(), undef); # "DESCRIPTION1\nDESCRIPTION2" );
is( $omim_entry->mapping_method(), "M method 1" );
is( $omim_entry->gene_status(), "C" );
is( $omim_entry->comment(), "comment1" );
is( $omim_entry->edited(), undef); # "ed1\ned2\ned3" );
is( $omim_entry->created(), undef); # "cd1\ncd2\ncd3" );
is( $omim_entry->contributors, undef); # "cn1\ncn2\ncn3" );
is( $omim_entry->additional_references(), "sa" );
is( ref($omim_entry->clinical_symptoms()), 'HASH' );
is( $omim_entry->species()->binomial(), "Homo sapiens" );


my $mini_mim = $omim_entry->miniMIM();

isa_ok( $mini_mim,"Bio::Phenotype::OMIM::MiniMIMentry" );
is( $mini_mim->description(), "Mini MIM text" );
is( $mini_mim->created(), "Mini MIM - cd" );
is( $mini_mim->contributors(), "Mini MIM - cn" );
is( $mini_mim->edited(), "Mini MIM - ed" );


my @corrs      = $omim_entry->each_Correlate();

is( $corrs[ 0 ]->name(), "mousecorrelate1" );
is( $corrs[ 0 ]->type(), "OMIM mouse correlate" );
is( $corrs[ 0 ]->species()->binomial(), "Mus musculus" );


my @cps        = $omim_entry->each_CytoPosition();

is( $cps[ 0 ]->value(), "1pter-p36.14" );


my @gss        = $omim_entry->each_gene_symbol();

is( $gss[ 0 ], "gene-symbol1" );


my @refs       = $omim_entry->each_Reference();

is( $refs[ 0 ]->authors(), "Author11, A. A.; Author12, A. A." );
is( $refs[ 0 ]->title(), "Title 1." );
is( $refs[ 0 ]->location(), "Am. J. Med. Genet1. 11 11-111 \(1981\)" );

is( $refs[ 1 ]->authors(), "Author21, A. A.; Author22, A. A." );
is( $refs[ 1 ]->title(), "Title 2." );
is( $refs[ 1 ]->location(), "Am. J. Med. Genet2. 12 22-222 \(1982\)" );

is( $refs[ 2 ]->authors(), "Author31, A. A.; Author32, A. A." );
is( $refs[ 2 ]->title(), "Title 3." );
is( $refs[ 2 ]->location(), "Am. J. Med. Genet3. 13 33-333 \(1983\)" );

is( $refs[ 3 ]->authors(), "" );
is( $refs[ 3 ]->title(), "other reference undef format" );
is( $refs[ 3 ]->location(), "" );



my @avs        = $omim_entry->each_AllelicVariant();

is( $avs[ 0 ]->number(), ".0001" );
is( $avs[ 0 ]->title(), "ALCOHOL INTOLERANCE, ACUTE" );
is( $avs[ 0 ]->symbol(), "ALDH2" );
is( $avs[ 0 ]->description(), "AV1-text" );
is( $avs[ 0 ]->aa_ori(), "GLU" );
is( $avs[ 0 ]->aa_mut(), "LYS" );
is( $avs[ 0 ]->position(), "487" );
is( $avs[ 0 ]->additional_mutations(), "" );


is( $avs[ 1 ]->number(), ".0002" );
is( $avs[ 1 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
is( $avs[ 1 ]->symbol(), "CHRNA1" );
is( $avs[ 1 ]->description(), "AV2-text" );
is( $avs[ 1 ]->aa_ori(), "VAL" );
is( $avs[ 1 ]->aa_mut(), "MET" );
is( $avs[ 1 ]->position(), "156" );
is( $avs[ 1 ]->additional_mutations(), "" );


is( $avs[ 2 ]->number(), ".0003" );
is( $avs[ 2 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
is( $avs[ 2 ]->symbol(), "CHRNE" );
is( $avs[ 2 ]->description(), "AV2-text a\nAV2-text b" );
is( $avs[ 2 ]->aa_ori(), "ARG" );
is( $avs[ 2 ]->aa_mut(), "LEU" );
is( $avs[ 2 ]->position(), "147" );
is( $avs[ 2 ]->additional_mutations(), "" );


is( $avs[ 3 ]->number(), ".0004" );
is( $avs[ 3 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
is( $avs[ 3 ]->symbol(), "CHRNE" );
is( $avs[ 3 ]->description(), "Sieb et al. (2000) found that a brother and sister with congenital\nmyasthenic syndrome (601462) were compound heterozygotes for a deletion\nof 911T and a splicing mutation (IVS4+1G-A; 100725.0007)." );
is( $avs[ 3 ]->aa_ori(), "" );
is( $avs[ 3 ]->aa_mut(), "" );
is( $avs[ 3 ]->position(), "" );
is( $avs[ 3 ]->additional_mutations(), "1-BP DEL, 911T" );


is( $avs[ 4 ]->number(), ".0005" );
is( $avs[ 4 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
is( $avs[ 4 ]->symbol(), "CHRNE" );
is( $avs[ 4 ]->description(), "See 100725.0006 and Sieb et al. (2000)." );
is( $avs[ 4 ]->aa_ori(), "" );
is( $avs[ 4 ]->aa_mut(), "" );
is( $avs[ 4 ]->position(), "" );
is( $avs[ 4 ]->additional_mutations(), "IVS4DS, G-A, +1" );



is( $avs[ 5 ]->number(), ".0006" );
is( $avs[ 5 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
is( $avs[ 5 ]->symbol(), "CHRNE" );
is( $avs[ 5 ]->description(), "AV6-text" );
is( $avs[ 5 ]->aa_ori(), "" );
is( $avs[ 5 ]->aa_mut(), "" );
is( $avs[ 5 ]->position(), "" );
is( $avs[ 5 ]->additional_mutations(), "1-BP DEL, 1030C" );





my $omim_entry2 = $omim_parser->next_phenotype();


isa_ok( $omim_entry2, "Bio::Phenotype::OMIM::OMIMentry" );

is( $omim_entry2->MIM_number(), "100501" );
is( $omim_entry2->title(), "#100501 second entry" );
is( $omim_entry2->alternative_titles_and_symbols(), ";;title1;;\ntitle2;;\ntitle3" );
is( $omim_entry2->more_than_two_genes(), 1 );
is( $omim_entry2->is_separate(), 0 );
is( $omim_entry2->description(), undef); # "DESCRIPTION1\nDESCRIPTION2" );
is( $omim_entry2->mapping_method(), "M method 2" );
is( $omim_entry2->gene_status(), "C" );
is( $omim_entry2->comment(), "comment2" );
is( $omim_entry2->edited(), undef); # "ed1\ned2\ned3" );
is( $omim_entry2->created(), undef); # "cd1\ncd2\ncd3" );
is( $omim_entry2->contributors(), undef); # "cn1\ncn2\ncn3" );
is( $omim_entry2->additional_references(), "sa" );

my $cs = $omim_entry2->clinical_symptoms();
is( ref($cs), 'HASH' );
is( $omim_entry2->species()->binomial(), "Homo sapiens" );


$mini_mim   = $omim_entry2->miniMIM();

isa_ok( $mini_mim, "Bio::Phenotype::OMIM::MiniMIMentry" );
is( $mini_mim->description(), "Mini MIM text" );
is( $mini_mim->created(), "Mini MIM - cd" );
is( $mini_mim->contributors(), "Mini MIM - cn" );
is( $mini_mim->edited(), "Mini MIM - ed" );


@corrs      = $omim_entry2->each_Correlate();

is( $corrs[ 0 ]->name(), "mousecorrelate2" );
is( $corrs[ 0 ]->type(), "OMIM mouse correlate" );
is( $corrs[ 0 ]->species()->binomial(), "Mus musculus" );


@cps        = $omim_entry2->each_CytoPosition();

is( $cps[ 0 ]->value(), "1pter-p36.15" );


@gss        = $omim_entry2->each_gene_symbol();

is( $gss[ 0 ], "gene-symbol2" );


@refs       = $omim_entry2->each_Reference();

is( $refs[ 0 ]->authors(), "Author11, A. A.; Author12, A. A." );
is( $refs[ 0 ]->title(), "Title 1." );
is( $refs[ 0 ]->location(), "Am. J. Med. Genet1. 11 11-111 \(1981\)" );

is( $refs[ 1 ]->authors(), "Author21, A. A.; Author22, A. A." );
is( $refs[ 1 ]->title(), "Title 2." );
is( $refs[ 1 ]->location(), "Am. J. Med. Genet2. 12 22-222 \(1982\)" );

is( $refs[ 2 ]->authors(), "Author31, A. A.; Author32, A. A." );
is( $refs[ 2 ]->title(), "Title 3." );
is( $refs[ 2 ]->location(), "Am. J. Med. Genet3. 13 33-333 \(1983\)" );

is( $refs[ 3 ]->authors(), "" );
is( $refs[ 3 ]->title(), "other reference undef format" );
is( $refs[ 3 ]->location(), "" );



@avs        = $omim_entry2->each_AllelicVariant();

is( $avs[ 0 ]->number(), ".0001" );
is( $avs[ 0 ]->title(), "ALCOHOL INTOLERANCE, ACUTE" );
is( $avs[ 0 ]->symbol(), "ALDH2" );
is( $avs[ 0 ]->description(), "AV1-text" );
is( $avs[ 0 ]->aa_ori(), "GLU" );
is( $avs[ 0 ]->aa_mut(), "LYS" );
is( $avs[ 0 ]->position(), "487" );
is( $avs[ 0 ]->additional_mutations(), "" );


is( $avs[ 1 ]->number(), ".0002" );
is( $avs[ 1 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
is( $avs[ 1 ]->symbol(), "CHRNA1" );
is( $avs[ 1 ]->description(), "AV2-text" );
is( $avs[ 1 ]->aa_ori(), "VAL" );
is( $avs[ 1 ]->aa_mut(), "MET" );
is( $avs[ 1 ]->position(), "156" );
is( $avs[ 1 ]->additional_mutations(), "" );


is( $avs[ 2 ]->number(), ".0003" );
is( $avs[ 2 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
is( $avs[ 2 ]->symbol(), "CHRNE" );
is( $avs[ 2 ]->description(), "AV2-text a\nAV2-text b" );
is( $avs[ 2 ]->aa_ori(), "ARG" );
is( $avs[ 2 ]->aa_mut(), "LEU" );
is( $avs[ 2 ]->position(), "147" );
is( $avs[ 2 ]->additional_mutations(), "" );


is( $avs[ 3 ]->number(), ".0004" );
is( $avs[ 3 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
is( $avs[ 3 ]->symbol(), "CHRNE" );
is( $avs[ 3 ]->description(), "Sieb et al. (2000) found that a brother and sister with congenital\nmyasthenic syndrome (601462) were compound heterozygotes for a deletion\nof 911T and a splicing mutation (IVS4+1G-A; 100725.0007)." );
is( $avs[ 3 ]->aa_ori(), "" );
is( $avs[ 3 ]->aa_mut(), "" );
is( $avs[ 3 ]->position(), "" );
is( $avs[ 3 ]->additional_mutations(), "1-BP DEL, 911T" );


is( $avs[ 4 ]->number(), ".0005" );
is( $avs[ 4 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
is( $avs[ 4 ]->symbol(), "CHRNE" );
is( $avs[ 4 ]->description(), "See 100725.0006 and Sieb et al. (2000)." );
is( $avs[ 4 ]->aa_ori(), "" );
is( $avs[ 4 ]->aa_mut(), "" );
is( $avs[ 4 ]->position(), "" );
is( $avs[ 4 ]->additional_mutations(), "IVS4DS, G-A, +1" );



is( $avs[ 5 ]->number(), ".0006" );
is( $avs[ 5 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
is( $avs[ 5 ]->symbol(), "CHRNE" );
is( $avs[ 5 ]->description(), "AV6-text" );
is( $avs[ 5 ]->aa_ori(), "" );
is( $avs[ 5 ]->aa_mut(), "" );
is( $avs[ 5 ]->position(), "" );
is( $avs[ 5 ]->additional_mutations(), "1-BP DEL, 1030C" );


# catch missing linebreak
throws_ok { my $omim_parser = Bio::Phenotype::OMIM::OMIMparser->new(
    -genemap  => test_input_file('omim_genemap_test_linebreak'),
    -omimtext => test_input_file('omim_text_test'));
} qr/linebreak/, 'missing linebreak caught';
