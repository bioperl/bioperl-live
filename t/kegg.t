# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

use strict;

BEGIN {
	eval { require Test; };
	if ( $@ ) {
		use lib 't';
	}
	use Test;
	plan tests => 12;
}

use Bio::SeqIO;

ok(1);

my $verbose = $ENV{'BIOPERLDEBUG'} || 0;

my $kegg = Bio::SeqIO->new(-format => 'kegg',
									-verbose => $verbose,
                           -file => Bio::Root::IO->catfile
                           ("t","data","AHCYL1.kegg"));
ok($kegg);
$kegg = $kegg->next_seq();
ok($kegg);
ok($kegg->accession, '10768');
ok($kegg->display_id, 'AHCYL1');
ok($kegg->alphabet, 'dna');
ok($kegg->seq);
ok($kegg->translate->seq);

ok(($kegg->annotation->get_Annotations('description'))[0]->text,
   'S-adenosylhomocysteine hydrolase-like 1 [EC:3.3.1.1]');

ok(($kegg->annotation->get_Annotations('pathway'))[0]->text,
   'Metabolism; Amino Acid Metabolism; Methionine metabolism');

ok( (grep {$_->database eq 'KO'}
     $kegg->annotation->get_Annotations('dblink'))[0]->comment, 
    'adenosylhomocysteinase' );

ok( (grep {$_->database eq 'PATH'} 
     $kegg->annotation->get_Annotations('dblink'))[0]->primary_id,
    'hsa00271' );
