# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 116,
               -requires_module => 'Bio::Tools::RNAMotif');
	
    use_ok('Bio::Tools::ERPIN');
    use_ok('Bio::Tools::Infernal');
}

my $verbose = test_debug();

### RNAMotif.pm tests ###

my $parser = Bio::Tools::RNAMotif->new(
		-verbose => $verbose,
        -file => test_input_file('trna.strict.rnamotif'),
        -motiftag => 'tRNA_gene',
        -desctag => 'tRNA');

my @genes;
while( my $gene = $parser->next_prediction ) {
    push @genes, $gene;
}

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
            -file => test_input_file('sprintf.rnamotif'),
            -motiftag => 'term',
            -desctag => 'stem_loop');

while( my $gene = $parser->next_prediction ) {
    push @genes, $gene;
}

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
            -file => test_input_file('testfile.erpin'),
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

my @stats = (
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

$parser = Bio::Tools::Infernal->new(
            -file => test_input_file('test.infernal'),
            -motiftag => 'misc_binding',
            -desctag => 'Purine riboswitch',
            -cm    => 'Purine',
            -rfam  =>  'RF00167',
            -minscore => 20);

my $gene = $parser->next_prediction;
# get query (model) data
is($gene->display_name, 'Purine riboswitch','Infernal::display_name()');
is($gene->seq_id, 'RF00167','Infernal::seq_id()');
is($gene->primary_tag, 'misc_binding','Infernal::primary_tag()');
is($gene->source_tag, 'Infernal 0.71','Infernal::source_tag()');
is($gene->start, '1','Infernal::start()');
is($gene->end, '102','Infernal::end()');
is($gene->strand, '0','Infernal::strand()');
is($gene->score, '78.40','Infernal::strand()');
# get hit data
$gene->invert;
is($gene->display_name, 'Purine riboswitch','Infernal::display_name()');
is($gene->seq_id, '2239287','Infernal::seq_id()');
is($gene->primary_tag, 'misc_binding','Infernal::primary_tag()');
is($gene->source_tag, 'Infernal 0.71','Infernal::source_tag()');
is($gene->start, '15589','Infernal::start()');
is($gene->end, '15691','Infernal::end()');
is($gene->strand, '1','Infernal::strand()');
is($gene->score, '78.40','Infernal::strand()');

is($gene->get_Annotations('model'),
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGACUAcG.aGuaAauauuaaauauuu',
   "Infernal::get_Annotations('model')");
is($gene->get_Annotations('midline'),
   ' A+ A+A+ AAAA A   :CUC:UAUAAU: :GGGAAUAUGGCCC: :AGUUUCUACC:GGCAACCGUAAAUUGCC:GACUA:G AG: AA + ++  +++++',
   "Infernal::get_Annotations('midline')");
is($gene->get_Annotations('hit'),
   'CAUGAAAUCAAAACACGACCUCAUAUAAUCUUGGGAAUAUGGCCCAUAAGUUUCUACCCGGCAACCGUAAAUUGCCGGACUAUGcAGGGAAGUGAUCGAUAAA',
   "Infernal::get_Annotations('hit')");
is($gene->get_Annotations('secstructure'),
   ':::::::::::::::::((((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,)))).))))::::::::::::::',
   "Infernal::get_Annotations('secstructure')");
is($gene->get_Annotations('seq_name'),
   'gi|2239287|gb|U51115.1|BSU51115',
   "Infernal::get_Annotations('seq_name')");

$gene = $parser->next_prediction;
# get query (model) data
is($gene->display_name, 'Purine riboswitch','Infernal::display_name()');
is($gene->seq_id, 'RF00167','Infernal::seq_id()');
is($gene->primary_tag, 'misc_binding','Infernal::primary_tag()');
is($gene->source_tag, 'Infernal 0.71','Infernal::source_tag()');
is($gene->start, '1','Infernal::start()');
is($gene->end, '102','Infernal::end()');
is($gene->strand, '0','Infernal::strand()');
is($gene->score, '81.29','Infernal::strand()');

$gene->invert; # switch to get hit data
is($gene->display_name, 'Purine riboswitch','Infernal::display_name()');
is($gene->seq_id, '2239287','Infernal::seq_id()');
is($gene->primary_tag, 'misc_binding','Infernal::primary_tag()');
is($gene->source_tag, 'Infernal 0.71','Infernal::source_tag()');
is($gene->start, '11655','Infernal::start()');
is($gene->end, '11756','Infernal::end()');
is($gene->strand, '1','Infernal::strand()');
is($gene->score, '81.29','Infernal::strand()');

is($gene->get_Annotations('model'),
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGACUAcGaGuaAauauuaaauauuu',
   "Infernal::get_Annotations('model')");
is($gene->get_Annotations('midline'),
   'A AAAU AAA+AA A+   : CGUAUAAU::CG:GAAUAUGGC:CG::AGU UCUACCA:GC ACCGUAAAU GC:UGACUACG :   AU+U +++  UUU',
   "Infernal::get_Annotations('midline')");
is($gene->get_Annotations('hit'),
   'AGAAAUCAAAUAAGAUGAAUUCGUAUAAUCGCGGGAAUAUGGCUCGCAAGUCUCUACCAAGCUACCGUAAAUGGCUUGACUACGUAAACAUUUCUUUCGUUU',
   "Infernal::get_Annotations('hit')");
is($gene->get_Annotations('secstructure'),
   ':::::::::::::::::((((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,))))))))::::::::::::::',
   "Infernal::get_Annotations('secstructure')");
is($gene->get_Annotations('seq_name'),
   'gi|2239287|gb|U51115.1|BSU51115',
   "Infernal::get_Annotations('seq_name')");
