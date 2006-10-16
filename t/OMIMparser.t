# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## # $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan tests => 173;
}

use File::Spec;
use Bio::Phenotype::OMIM::OMIMparser;


my $omim_parser = Bio::Phenotype::OMIM::OMIMparser->new( -genemap  => File::Spec->catfile(qw(t data omim_genemap_test)),
                                                         -omimtext => File::Spec->catfile(qw(t data omim_text_test)) );


ok( $omim_parser->isa( "Bio::Phenotype::OMIM::OMIMparser" ) );

my $omim_entry = $omim_parser->next_phenotype();


ok( $omim_entry->isa( "Bio::Phenotype::OMIM::OMIMentry" ) );

ok( $omim_entry->MIM_number(), "100500" );
ok( $omim_entry->title(), "*100500 title" );
ok( $omim_entry->alternative_titles_and_symbols(), ";;title1;;\ntitle2;;\ntitle3" );
ok( $omim_entry->more_than_two_genes(), 0 );
ok( $omim_entry->is_separate(), 1 );
ok( $omim_entry->description(), undef); # "DESCRIPTION1\nDESCRIPTION2" );
ok( $omim_entry->mapping_method(), "M method 1" );
ok( $omim_entry->gene_status(), "C" );
ok( $omim_entry->comment(), "comment1" );
ok( $omim_entry->edited(), undef); # "ed1\ned2\ned3" );
ok( $omim_entry->created(), undef); # "cd1\ncd2\ncd3" );
ok( $omim_entry->contributors, undef); # "cn1\ncn2\ncn3" );
ok( $omim_entry->additional_references(), "sa" );
ok( ref($omim_entry->clinical_symptoms()), 'HASH' );
ok( $omim_entry->species()->binomial(), "Homo sapiens" );


my $mini_mim = $omim_entry->miniMIM();

ok( $mini_mim->isa( "Bio::Phenotype::OMIM::MiniMIMentry" ) );
ok( $mini_mim->description(), "Mini MIM text" );
ok( $mini_mim->created(), "Mini MIM - cd" );
ok( $mini_mim->contributors(), "Mini MIM - cn" );
ok( $mini_mim->edited(), "Mini MIM - ed" );


my @corrs      = $omim_entry->each_Correlate();

ok( $corrs[ 0 ]->name(), "mousecorrelate1" );
ok( $corrs[ 0 ]->type(), "OMIM mouse correlate" );
ok( $corrs[ 0 ]->species()->binomial(), "Mus musculus" );


my @cps        = $omim_entry->each_CytoPosition();

ok( $cps[ 0 ]->value(), "1pter-p36.14" );


my @gss        = $omim_entry->each_gene_symbol();

ok( $gss[ 0 ], "gene-symbol1" );


my @refs       = $omim_entry->each_Reference();

ok( $refs[ 0 ]->authors(), "Author11, A. A.; Author12, A. A." );
ok( $refs[ 0 ]->title(), "Title 1." );
ok( $refs[ 0 ]->location(), "Am. J. Med. Genet1. 11 11-111 \(1981\)" );

ok( $refs[ 1 ]->authors(), "Author21, A. A.; Author22, A. A." );
ok( $refs[ 1 ]->title(), "Title 2." );
ok( $refs[ 1 ]->location(), "Am. J. Med. Genet2. 12 22-222 \(1982\)" );

ok( $refs[ 2 ]->authors(), "Author31, A. A.; Author32, A. A." );
ok( $refs[ 2 ]->title(), "Title 3." );
ok( $refs[ 2 ]->location(), "Am. J. Med. Genet3. 13 33-333 \(1983\)" );

ok( $refs[ 3 ]->authors(), "" );
ok( $refs[ 3 ]->title(), "other reference undef format" );
ok( $refs[ 3 ]->location(), "" );



my @avs        = $omim_entry->each_AllelicVariant();

ok( $avs[ 0 ]->number(), ".0001" );
ok( $avs[ 0 ]->title(), "ALCOHOL INTOLERANCE, ACUTE" );
ok( $avs[ 0 ]->symbol(), "ALDH2" );
ok( $avs[ 0 ]->description(), "AV1-text" );
ok( $avs[ 0 ]->aa_ori(), "GLU" );
ok( $avs[ 0 ]->aa_mut(), "LYS" );
ok( $avs[ 0 ]->position(), "487" );
ok( $avs[ 0 ]->additional_mutations(), "" );


ok( $avs[ 1 ]->number(), ".0002" );
ok( $avs[ 1 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
ok( $avs[ 1 ]->symbol(), "CHRNA1" );
ok( $avs[ 1 ]->description(), "AV2-text" );
ok( $avs[ 1 ]->aa_ori(), "VAL" );
ok( $avs[ 1 ]->aa_mut(), "MET" );
ok( $avs[ 1 ]->position(), "156" );
ok( $avs[ 1 ]->additional_mutations(), "" );


ok( $avs[ 2 ]->number(), ".0003" );
ok( $avs[ 2 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
ok( $avs[ 2 ]->symbol(), "CHRNE" );
ok( $avs[ 2 ]->description(), "AV2-text a\nAV2-text b" );
ok( $avs[ 2 ]->aa_ori(), "ARG" );
ok( $avs[ 2 ]->aa_mut(), "LEU" );
ok( $avs[ 2 ]->position(), "147" );
ok( $avs[ 2 ]->additional_mutations(), "" );


ok( $avs[ 3 ]->number(), ".0004" );
ok( $avs[ 3 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
ok( $avs[ 3 ]->symbol(), "CHRNE" );
ok( $avs[ 3 ]->description(), "Sieb et al. (2000) found that a brother and sister with congenital\nmyasthenic syndrome (601462) were compound heterozygotes for a deletion\nof 911T and a splicing mutation (IVS4+1G-A; 100725.0007)." );
ok( $avs[ 3 ]->aa_ori(), "" );
ok( $avs[ 3 ]->aa_mut(), "" );
ok( $avs[ 3 ]->position(), "" );
ok( $avs[ 3 ]->additional_mutations(), "1-BP DEL, 911T" );


ok( $avs[ 4 ]->number(), ".0005" );
ok( $avs[ 4 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
ok( $avs[ 4 ]->symbol(), "CHRNE" );
ok( $avs[ 4 ]->description(), "See 100725.0006 and Sieb et al. (2000)." );
ok( $avs[ 4 ]->aa_ori(), "" );
ok( $avs[ 4 ]->aa_mut(), "" );
ok( $avs[ 4 ]->position(), "" );
ok( $avs[ 4 ]->additional_mutations(), "IVS4DS, G-A, +1" );



ok( $avs[ 5 ]->number(), ".0006" );
ok( $avs[ 5 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
ok( $avs[ 5 ]->symbol(), "CHRNE" );
ok( $avs[ 5 ]->description(), "AV6-text" );
ok( $avs[ 5 ]->aa_ori(), "" );
ok( $avs[ 5 ]->aa_mut(), "" );
ok( $avs[ 5 ]->position(), "" );
ok( $avs[ 5 ]->additional_mutations(), "1-BP DEL, 1030C" );





my $omim_entry2 = $omim_parser->next_phenotype();


ok( $omim_entry2->isa( "Bio::Phenotype::OMIM::OMIMentry" ) );

ok( $omim_entry2->MIM_number(), "100501" );
ok( $omim_entry2->title(), "#100501 second entry" );
ok( $omim_entry2->alternative_titles_and_symbols(), ";;title1;;\ntitle2;;\ntitle3" );
ok( $omim_entry2->more_than_two_genes(), 1 );
ok( $omim_entry2->is_separate(), 0 );
ok( $omim_entry2->description(), undef); # "DESCRIPTION1\nDESCRIPTION2" );
ok( $omim_entry2->mapping_method(), "M method 2" );
ok( $omim_entry2->gene_status(), "C" );
ok( $omim_entry2->comment(), "comment2" );
ok( $omim_entry2->edited(), undef); # "ed1\ned2\ned3" );
ok( $omim_entry2->created(), undef); # "cd1\ncd2\ncd3" );
ok( $omim_entry2->contributors(), undef); # "cn1\ncn2\ncn3" );
ok( $omim_entry2->additional_references(), "sa" );

my $cs = $omim_entry2->clinical_symptoms();
ok( ref($cs), 'HASH' );
ok( $omim_entry2->species()->binomial(), "Homo sapiens" );


$mini_mim   = $omim_entry2->miniMIM();

ok( $mini_mim->isa( "Bio::Phenotype::OMIM::MiniMIMentry" ) );
ok( $mini_mim->description(), "Mini MIM text" );
ok( $mini_mim->created(), "Mini MIM - cd" );
ok( $mini_mim->contributors(), "Mini MIM - cn" );
ok( $mini_mim->edited(), "Mini MIM - ed" );


@corrs      = $omim_entry2->each_Correlate();

ok( $corrs[ 0 ]->name(), "mousecorrelate2" );
ok( $corrs[ 0 ]->type(), "OMIM mouse correlate" );
ok( $corrs[ 0 ]->species()->binomial(), "Mus musculus" );


@cps        = $omim_entry2->each_CytoPosition();

ok( $cps[ 0 ]->value(), "1pter-p36.15" );


@gss        = $omim_entry2->each_gene_symbol();

ok( $gss[ 0 ], "gene-symbol2" );


@refs       = $omim_entry2->each_Reference();

ok( $refs[ 0 ]->authors(), "Author11, A. A.; Author12, A. A." );
ok( $refs[ 0 ]->title(), "Title 1." );
ok( $refs[ 0 ]->location(), "Am. J. Med. Genet1. 11 11-111 \(1981\)" );

ok( $refs[ 1 ]->authors(), "Author21, A. A.; Author22, A. A." );
ok( $refs[ 1 ]->title(), "Title 2." );
ok( $refs[ 1 ]->location(), "Am. J. Med. Genet2. 12 22-222 \(1982\)" );

ok( $refs[ 2 ]->authors(), "Author31, A. A.; Author32, A. A." );
ok( $refs[ 2 ]->title(), "Title 3." );
ok( $refs[ 2 ]->location(), "Am. J. Med. Genet3. 13 33-333 \(1983\)" );

ok( $refs[ 3 ]->authors(), "" );
ok( $refs[ 3 ]->title(), "other reference undef format" );
ok( $refs[ 3 ]->location(), "" );



@avs        = $omim_entry2->each_AllelicVariant();

ok( $avs[ 0 ]->number(), ".0001" );
ok( $avs[ 0 ]->title(), "ALCOHOL INTOLERANCE, ACUTE" );
ok( $avs[ 0 ]->symbol(), "ALDH2" );
ok( $avs[ 0 ]->description(), "AV1-text" );
ok( $avs[ 0 ]->aa_ori(), "GLU" );
ok( $avs[ 0 ]->aa_mut(), "LYS" );
ok( $avs[ 0 ]->position(), "487" );
ok( $avs[ 0 ]->additional_mutations(), "" );


ok( $avs[ 1 ]->number(), ".0002" );
ok( $avs[ 1 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
ok( $avs[ 1 ]->symbol(), "CHRNA1" );
ok( $avs[ 1 ]->description(), "AV2-text" );
ok( $avs[ 1 ]->aa_ori(), "VAL" );
ok( $avs[ 1 ]->aa_mut(), "MET" );
ok( $avs[ 1 ]->position(), "156" );
ok( $avs[ 1 ]->additional_mutations(), "" );


ok( $avs[ 2 ]->number(), ".0003" );
ok( $avs[ 2 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
ok( $avs[ 2 ]->symbol(), "CHRNE" );
ok( $avs[ 2 ]->description(), "AV2-text a\nAV2-text b" );
ok( $avs[ 2 ]->aa_ori(), "ARG" );
ok( $avs[ 2 ]->aa_mut(), "LEU" );
ok( $avs[ 2 ]->position(), "147" );
ok( $avs[ 2 ]->additional_mutations(), "" );


ok( $avs[ 3 ]->number(), ".0004" );
ok( $avs[ 3 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
ok( $avs[ 3 ]->symbol(), "CHRNE" );
ok( $avs[ 3 ]->description(), "Sieb et al. (2000) found that a brother and sister with congenital\nmyasthenic syndrome (601462) were compound heterozygotes for a deletion\nof 911T and a splicing mutation (IVS4+1G-A; 100725.0007)." );
ok( $avs[ 3 ]->aa_ori(), "" );
ok( $avs[ 3 ]->aa_mut(), "" );
ok( $avs[ 3 ]->position(), "" );
ok( $avs[ 3 ]->additional_mutations(), "1-BP DEL, 911T" );


ok( $avs[ 4 ]->number(), ".0005" );
ok( $avs[ 4 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
ok( $avs[ 4 ]->symbol(), "CHRNE" );
ok( $avs[ 4 ]->description(), "See 100725.0006 and Sieb et al. (2000)." );
ok( $avs[ 4 ]->aa_ori(), "" );
ok( $avs[ 4 ]->aa_mut(), "" );
ok( $avs[ 4 ]->position(), "" );
ok( $avs[ 4 ]->additional_mutations(), "IVS4DS, G-A, +1" );



ok( $avs[ 5 ]->number(), ".0006" );
ok( $avs[ 5 ]->title(), "MYASTHENIC SYNDROME, SLOW-CHANNEL CONGENITAL" );
ok( $avs[ 5 ]->symbol(), "CHRNE" );
ok( $avs[ 5 ]->description(), "AV6-text" );
ok( $avs[ 5 ]->aa_ori(), "" );
ok( $avs[ 5 ]->aa_mut(), "" );
ok( $avs[ 5 ]->position(), "" );
ok( $avs[ 5 ]->additional_mutations(), "1-BP DEL, 1030C" );









