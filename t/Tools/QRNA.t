# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 30);
	
	use_ok('Bio::Tools::QRNA');
}

my $inputfilename= test_input_file('ecoli-trna-qrna.out');
ok my $parser = Bio::Tools::QRNA->new(-file => $inputfilename);

my $rnacount = 0;
while( my $f = $parser->next_feature ) {
    if( $f->primary_tag eq 'RNA' ) { # winning model is primary tag
	if( ! $rnacount ) { # 1st time through let's test
	    is($f->feature1->start,4);
	    is($f->feature1->end,  70);
	    is($f->score, 22.147);
	    is($f->feature1->seq_id,'DA0780-1-');
	    
	    is($f->feature2->start, 4);
	    is($f->feature2->end,  70);
	    is($f->feature2->seq_id, 'ECOLI-3979754-');
	    is(($f->get_tag_values('alignment_len'))[0], 70);
	    is(($f->get_tag_values('alignment_pid'))[0], '72.86');
	    is(($f->get_tag_values('COD_score'))[0], '16.954');
	    is(($f->get_tag_values('COD_logoddspost'))[0], '-4.365');
	    is(($f->get_tag_values('OTH_score'))[0], '21.319');
	    is(($f->get_tag_values('OTH_logoddspost'))[0], '0.000');
	}
	$rnacount++;
    }
}
is($rnacount, 21);
$inputfilename= test_input_file('qrna-relloc.out');
$parser = Bio::Tools::QRNA->new(-file => $inputfilename);

my $qrna = $parser->next_feature;
is($qrna->primary_tag, 'COD');
is($qrna->source_tag, 'qrna');
is($qrna->feature1->seq_id, 'Contig1');
is($qrna->feature2->seq_id, 'chr5.pseudo');
is($qrna->feature1->start, 24732);
is($qrna->feature1->end, 24881);

is($qrna->feature2->start, 527251);
is($qrna->feature2->end, 527400);

is($parser->seq_file,'tst.out');
is($parser->RNA_model, '/mix_tied_linux.cfg');
is($parser->PAM_model, 'BLOSUM62 scaled by 1.000');
is($parser->program_name, 'qrna');
is($parser->program_version, '1.2b');
is($parser->program_date, 'Tue Dec 18 15:04:38 CST 2001');
