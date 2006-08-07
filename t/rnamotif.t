# -*-Perl-*-

## Bioperl Test Harness Script for Modules

## $Id$



# Before `make install' is performed this script should be runnable with

# `make test'. After `make install' it should work as `perl test.t'



my $error = 0;



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

    plan tests => 24; 

}



if( $error == 1 ) {

    exit(0);

}



use Bio::Tools::RNAMotif;

use Bio::Root::IO;

my $verbose = $ENV{'BIOPERLDEBUG'} || 0;



my $parser = Bio::Tools::RNAMotif->new(

			-verbose => $verbose,

            -file => Bio::Root::IO->catfile('t','data','trna.strict.rnamotif'),

            -motiftag => 'tRNA_gene',

            -desctag => 'tRNA');



my @genes;

while( my $gene = $parser->next_prediction ) {

    push @genes, $gene;

}

#tests 1-12 

ok($genes[1]->display_name, 'tRNA');

ok($genes[12]->seq_id, 'M33910');

ok($genes[6]->primary_tag, 'tRNA_gene');

ok($genes[22]->start, 464);

ok($genes[8]->end, 585);

ok($genes[9]->strand, 1);

ok($genes[90]->get_Annotations('sequence'), 'cggatt ta ttg ggcg taa a gggct cgtaggc ggctc gtcgcgtccggtgtgaaagtc catc gcttaac ggtg gatctg cgcc');

ok($genes[84]->get_Annotations('descfile'), 'trna.strict.descr');

ok($genes[4]->get_Annotations('descline'),'gi|173683|gb|M10671|ACSTRW Avian oncornavirus Trp-tRNA');

ok($genes[26]->get_Annotations('secstructure'), 'h5 ss h5 ss h3 ss h5 ss h3 ss h5 ss h3 h3 ss');

ok($genes[4]->score, '0.000');

ok($genes[4]->source_tag, 'RNAMotif');



@genes=();



my $parser2 = Bio::Tools::RNAMotif->new(

			-verbose => $verbose,

            -file => Bio::Root::IO->catfile('t','data','sprintf.rnamotif'),

            -motiftag => 'term',

            -desctag => 'stem_loop');



while( my $gene = $parser2->next_prediction ) {

    push @genes, $gene;

}



#tests 13-24

ok($genes[1]->display_name, 'stem_loop');

ok($genes[12]->seq_id, 'M82700');

ok($genes[6]->primary_tag, 'term');

ok($genes[22]->start, 141);

ok($genes[8]->end, 154);

ok($genes[9]->strand, -1);

ok($genes[90]->get_Annotations('sequence'), 'ggggaag cttg cttcccc');

ok($genes[84]->get_Annotations('descfile'), 'sprintf.descr');

ok($genes[4]->get_Annotations('descline'),'gi|173741|gb|M83548|AQF16SRRN Aquifex pyrophilus 16S ribosomal RNA (16S rRNA)');

ok($genes[26]->get_Annotations('secstructure'), 'h5 ss h3');

ok($genes[4]->score, '-12.100,5,gaaa');

ok($genes[4]->source_tag, 'RNAMotif');

