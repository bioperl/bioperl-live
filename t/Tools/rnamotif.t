# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 0);
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


my $val;
is($genes[1]->display_name, 'tRNA','RNAMotif::display_name()');
is($genes[12]->seq_id, 'M33910','RNAMotif::seq_id()');
is($genes[6]->primary_tag, 'tRNA_gene','RNAMotif::primary_tag()');
is($genes[22]->start, 464,'RNAMotif::start()');
is($genes[8]->end, 585,'RNAMotif::end()');
is($genes[9]->strand, 1,'RNAMotif::strand()');
is(($genes[90]->get_tag_values('sequence'))[0],
   'cggatt ta ttg ggcg taa a gggct cgtaggc ggctc'.
   ' gtcgcgtccggtgtgaaagtc catc gcttaac ggtg gatctg cgcc',
   "RNAMotif::get_tag_values('sequence')");
is(($genes[90]->get_tag_values('descfile'))[0], 'trna.strict.descr',
   "RNAMotif::get_tag_values('descfile')");
is(($genes[4]->get_tag_values('descline'))[0],
   'gi|173683|gb|M10671|ACSTRW Avian oncornavirus Trp-tRNA',
   "RNAMotif::get_tag_values('descline')");
is(($genes[26]->get_tag_values('secstructure'))[0],
   'h5 ss h5 ss h3 ss h5 ss h3 ss h5 ss h3 h3 ss',
   "RNAMotif::get_tag_values('secstructure')");
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
is(($genes[90]->get_tag_values('sequence'))[0], 'ggggaag cttg cttcccc',
   "RNAMotif::get_tag_values('sequence')");
is(($genes[84]->get_tag_values('descfile'))[0], 'sprintf.descr',
   "RNAMotif::get_tag_values('descfile')");
is(($genes[4]->get_tag_values('descline'))[0],
   'gi|173741|gb|M83548|AQF16SRRN Aquifex pyrophilus 16S ribosomal RNA (16S rRNA)',
   "RNAMotif::get_tag_values('descline')");
is(($genes[26]->get_tag_values('secstructure'))[0], 'h5 ss h3',
   "RNAMotif::get_Annotations('secstructure')");
is($genes[4]->score, '-12.100,5,gaaa','RNAMotif::score()');
is($genes[4]->source_tag, 'RNAMotif','RNAMotif::source_tag()');

### ERPIN.pm tests ###

@genes = ();

my @erpinstats = (
['30260185','5181155','5181183',1,'CTTT.aacc--.CAACC.CCGTGA.GGTTG.a.GAAG',
 'gi|30260185|gb|AE016879.1| Bacillus anthracis str. Ames, complete genome',
 '1.68e-05'],
['30260185','3709092','3709121',-1,'CTTT.taatt-.CAGTC.CTGTGA.GACCG.g.AAAG',
 'gi|30260185|gb|AE016879.1| Bacillus anthracis str. Ames, complete genome',
 '5.61e-05'],
['30260185','3710524','3710553',-1,'TTTT.aaatg-.TAGTC.CTGTGA.GGCTG.c.CAAA',
 'gi|30260185|gb|AE016879.1| Bacillus anthracis str. Ames, complete genome',
 '1.31e-04'],
['30260185','3711223','3711251',-1,'CTTT.aaca--.CAGCC.CCGTGA.GGTTG.a.GAAG',
 'gi|30260185|gb|AE016879.1| Bacillus anthracis str. Ames, complete genome',
 '4.44e-06']
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
	($val) = $gene->get_tag_values('sequence');
	is($val, shift @stats, "ERPIN::get_tag_values('sequence')");
	eval {($val) = $gene->get_tag_values('descfile')}; # no descfile
	like($@, qr(asking for tag value that does not exist descfile));
	($val) = $gene->get_tag_values('descline'); 
	is($val, shift @stats, "ERPIN::get_tag_values('descline')");
	eval {($val) = $gene->get_tag_values('secstructure')}; # no secstructure 
	like($@, qr(asking for tag value that does not exist secstructure));
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
($val) = $gene->get_tag_values('model');
is($val,
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGACUAcG.aGuaAauauuaaauauuu',
   "Infernal::get_tag_values('model')");
($val) = $gene->get_tag_values('midline');
is($val,
   ' A+ A+A+ AAAA A   :CUC:UAUAAU: :GGGAAUAUGGCCC: :AGUUUCUACC:GGCAACCGUAAAUUGCC:GACUA:G AG: AA + ++  +++++',
   "Infernal::get_tag_values('midline')");
($val) = $gene->get_tag_values('hit');
is($val,
   'CAUGAAAUCAAAACACGACCUCAUAUAAUCUUGGGAAUAUGGCCCAUAAGUUUCUACCCGGCAACCGUAAAUUGCCGGACUAUGcAGGGAAGUGAUCGAUAAA',
   "Infernal::get_tag_values('hit')");
($val) = $gene->get_tag_values('secstructure');
is($val,
   ':::::::::::::::::((((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,)))).))))::::::::::::::',
   "Infernal::get_tag_values('secstructure')");
($val) = $gene->get_tag_values('seq_name'),
is($val, 'gi|2239287|gb|U51115.1|BSU51115',
   "Infernal::get_tag_values('seq_name')");

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
($val) = $gene->get_tag_values('model');
is($val,
   'aAaaauaaAaaaaaaaauaCuCgUAUAaucucgggAAUAUGGcccgagaGUuUCUACCaGgcaaCCGUAAAuugcCuGACUAcGaGuaAauauuaaauauuu',
   "Infernal::get_tag_values('model')");
($val) = $gene->get_tag_values('midline');
is($val,
   'A AAAU AAA+AA A+   : CGUAUAAU::CG:GAAUAUGGC:CG::AGU UCUACCA:GC ACCGUAAAU GC:UGACUACG :   AU+U +++  UUU',
   "Infernal::get_tag_values('midline')");
($val) = $gene->get_tag_values('hit');
is($val,
   'AGAAAUCAAAUAAGAUGAAUUCGUAUAAUCGCGGGAAUAUGGCUCGCAAGUCUCUACCAAGCUACCGUAAAUGGCUUGACUACGUAAACAUUUCUUUCGUUU',
   "Infernal::get_tag_values('hit')");
($val) = $gene->get_tag_values('secstructure');
is($val,
   ':::::::::::::::::((((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,))))))))::::::::::::::',
   "Infernal::get_tag_values('secstructure')");
($val) = $gene->get_tag_values('seq_name');
is($val,
   'gi|2239287|gb|U51115.1|BSU51115',
   "Infernal::get_tag_values('seq_name')");
