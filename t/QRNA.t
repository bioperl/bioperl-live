# -*-Perl-*-
## Bioperl Test Harness Script for Modules
# $Id$

use strict;
BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    use vars qw($NTESTS);
    $NTESTS = 29;
    plan tests => $NTESTS;
}
use Bio::Tools::QRNA;
use Bio::Root::IO;

my $inputfilename= Bio::Root::IO->catfile("t","data","ecoli-trna-qrna.out");
my $parser = new Bio::Tools::QRNA(-file => $inputfilename);
ok($parser);
my $rnacount = 0;
while( my $f = $parser->next_feature ) {
    if( $f->primary_tag eq 'RNA' ) { # winning model is primary tag
	if( ! $rnacount ) { # 1st time through let's test
	    ok($f->feature1->start,4);
	    ok($f->feature1->end,  70);
	    ok($f->score, 22.147);
	    ok($f->feature1->seq_id,'DA0780-1-');
	    
	    ok($f->feature2->start, 4);
	    ok($f->feature2->end,  70);
	    ok($f->feature2->seq_id, 'ECOLI-3979754-');
	    ok(($f->get_tag_values('alignment_len'))[0], 70);
	    ok(($f->get_tag_values('alignment_pid'))[0], '72.86');
	    ok(($f->get_tag_values('COD_score'))[0], '16.954');
	    ok(($f->get_tag_values('COD_logoddspost'))[0], '-4.365');
	    ok(($f->get_tag_values('OTH_score'))[0], '21.319');
	    ok(($f->get_tag_values('OTH_logoddspost'))[0], '0.000');
	}
	$rnacount++;
    }
}
ok($rnacount, 21);
$inputfilename= Bio::Root::IO->catfile("t","data","qrna-relloc.out");
$parser = new Bio::Tools::QRNA(-file => $inputfilename);

my $qrna = $parser->next_feature;
ok($qrna->primary_tag, 'COD');
ok($qrna->source_tag, 'qrna');
ok($qrna->feature1->seq_id, 'Contig1');
ok($qrna->feature2->seq_id, 'chr5.pseudo');
ok($qrna->feature1->start, 24732);
ok($qrna->feature1->end, 24881);

ok($qrna->feature2->start, 527251);
ok($qrna->feature2->end, 527400);

ok($parser->seq_file,'tst.out');
ok($parser->RNA_model, '/mix_tied_linux.cfg');
ok($parser->PAM_model, 'BLOSUM62 scaled by 1.000');
ok($parser->program_name, 'qrna');
ok($parser->program_version, '1.2b');
ok($parser->program_date, 'Tue Dec 18 15:04:38 CST 2001');

