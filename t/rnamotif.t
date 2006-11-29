# -*-Perl-*-
# Bioperl Test Script for RNA Motif Modules
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test::More; };
    if( $@ ) {
    use lib 't/lib';
    }

    use Test::More;
	
	eval {
		require Bio::Tools::RNAMotif;
	};
	if ($@) {
		plan skip_all => 'Bio::Tools::RNAMotif failed to load, DB_File probably not installed. This means that the module is not usable. Skipping tests';
	}
	else {
		plan tests => 72;
	}
}

use Bio::Tools::ERPIN;
use Bio::Root::IO;
my $verbose = $ENV{'BIOPERLDEBUG'} || 0;

### RNAMotif.pm tests ###

my $parser = new Bio::Tools::RNAMotif(
		-verbose => $verbose,
        -file => Bio::Root::IO->catfile('t','data','trna.strict.rnamotif'),
        -motiftag => 'tRNA_gene',
        -desctag => 'tRNA');


my @genes;
while( my $gene = $parser->next_prediction ) {
    push @genes, $gene;
}
#tests 1-12 
is($genes[1]->display_name, 'tRNA','RNAMotif::display_name()');
is($genes[12]->seq_id, 'M33910','RNAMotif::seq_id()');
is($genes[6]->primary_tag, 'tRNA_gene','RNAMotif::primary_tag()');
is($genes[22]->start, 464,'RNAMotif::start()');
is($genes[8]->end, 585,'RNAMotif::end()');
is($genes[9]->strand, 1,'RNAMotif::strand()');
is($genes[90]->get_Annotations('sequence'),
   'cggatt ta ttg ggcg taa a gggct cgtaggc ggctc'.
   ' gtcgcgtccggtgtgaaagtc catc gcttaac ggtg gatctg cgcc',
   "RNAMotif::get_Annotations('sequence')");
is($genes[84]->get_Annotations('descfile'), 'trna.strict.descr',
   "RNAMotif::get_Annotations('descfile')");
is($genes[4]->get_Annotations('descline'),
   'gi|173683|gb|M10671|ACSTRW Avian oncornavirus Trp-tRNA',
   "RNAMotif::get_Annotations('descline')");
is($genes[26]->get_Annotations('secstructure'),
   'h5 ss h5 ss h3 ss h5 ss h3 ss h5 ss h3 h3 ss',
   "RNAMotif::get_Annotations('secstructure')");
is($genes[4]->score, '0.000','RNAMotif::score()');
is($genes[4]->source_tag, 'RNAMotif','RNAMotif::source_tag()');

@genes=();

$parser = Bio::Tools::RNAMotif->new(
			-verbose => $verbose,
            -file => Bio::Root::IO->catfile('t','data','sprintf.rnamotif'),
            -motiftag => 'term',
            -desctag => 'stem_loop');

while( my $gene = $parser->next_prediction ) {
    push @genes, $gene;
}

#tests 13-24
is($genes[1]->display_name, 'stem_loop','RNAMotif::display_name()');
is($genes[12]->seq_id, 'M82700','RNAMotif::seq_id()');
is($genes[6]->primary_tag, 'term','RNAMotif::primary_tag()');
is($genes[22]->start, 141,'RNAMotif::start()');
is($genes[8]->end, 154,'RNAMotif::end()');
is($genes[9]->strand, -1,'RNAMotif::strand()');
is($genes[90]->get_Annotations('sequence'), 'ggggaag cttg cttcccc',
   "RNAMotif::get_Annotations('sequence')");
is($genes[84]->get_Annotations('descfile'), 'sprintf.descr',
   "RNAMotif::get_Annotations('descfile')");
is($genes[4]->get_Annotations('descline'),
   'gi|173741|gb|M83548|AQF16SRRN Aquifex pyrophilus 16S ribosomal RNA (16S rRNA)',
   "RNAMotif::get_Annotations('descline')");
is($genes[26]->get_Annotations('secstructure'), 'h5 ss h3',
   "RNAMotif::get_Annotations('secstructure')");
is($genes[4]->score, '-12.100,5,gaaa','RNAMotif::score()');
is($genes[4]->source_tag, 'RNAMotif','RNAMotif::source_tag()');

### ERPIN.pm tests ###

@genes = ();

my @erpinstats = (
['30260185','5181155','5181183',1,'CTTT.aacc--.CAACC.CCGTGA.GGTTG.a.GAAG',0,
 'gi|30260185|gb|AE016879.1| Bacillus anthracis str. Ames, complete genome',
 0,'1.68e-05'],
['30260185','3709092','3709121',-1,'CTTT.taatt-.CAGTC.CTGTGA.GACCG.g.AAAG',0,
 'gi|30260185|gb|AE016879.1| Bacillus anthracis str. Ames, complete genome',
 0,'5.61e-05'],
['30260185','3710524','3710553',-1,'TTTT.aaatg-.TAGTC.CTGTGA.GGCTG.c.CAAA',0,
 'gi|30260185|gb|AE016879.1| Bacillus anthracis str. Ames, complete genome',
 0,'1.31e-04'],
['30260185','3711223','3711251',-1,'CTTT.aaca--.CAGCC.CCGTGA.GGTTG.a.GAAG',0,
 'gi|30260185|gb|AE016879.1| Bacillus anthracis str. Ames, complete genome',
 0,'4.44e-06']
);

$parser = Bio::Tools::ERPIN->new(
			-verbose => $verbose,
            -file => Bio::Root::IO->catfile('t','data','testfile.erpin'),
            -motiftag => 'protein_bind',
			-desctag =>  'pyrR_BL');

while( my $gene = $parser->next_prediction ) {
	my @stats = @{ shift @erpinstats };
	is($gene->display_name, 'pyrR_BL','ERPIN::display_name()');
	is($gene->seq_id, shift @stats,'ERPIN::seq_id()');
	is($gene->primary_tag, 'protein_bind','ERPIN::primary_tag()');
	is($gene->start, shift @stats,'ERPIN::start()');
	is($gene->end, shift @stats,'ERPIN::end()');
	is($gene->strand, shift @stats,'ERPIN::strand()');
	is($gene->get_Annotations('sequence'), shift @stats,
	   "ERPIN::get_Annotations('sequence')");
	is($gene->get_Annotations('descfile'), shift @stats,
	   "ERPIN::get_Annotations('descfile')");
	is($gene->get_Annotations('descline'), shift @stats,
	   "ERPIN::get_Annotations('descline')");
	is($gene->get_Annotations('secstructure'), shift @stats,
	   "ERPIN::get_Annotations('secstructure')");
	is($gene->score, shift @stats,'ERPIN::score()');
	is($gene->source_tag, 'ERPIN','ERPIN::source_tag()');
}

### Infernal.pm tests ###
### FASTR.pm tests ###