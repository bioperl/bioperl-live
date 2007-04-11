# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

use strict;

BEGIN {
	eval { require Test::More; };
	if ( $@ ) {
		use lib 't/lib';
	}
	use Test::More;
	plan tests => 13;
	use_ok('Bio::SeqIO');
}

my $verbose = $ENV{'BIOPERLDEBUG'} || 0;

my $io = Bio::SeqIO->new(-format => 'kegg',
								 -verbose => $verbose,
								 -file => Bio::Root::IO->catfile
								 ("t","data","AHCYL1.kegg"));
ok($io);
my $kegg = $io->next_seq();
ok($kegg);
is($kegg->accession, '10768');
is($kegg->display_id, 'AHCYL1');
is($kegg->alphabet, 'dna');
is($kegg->seq,'atgtcgatgcctgacgcgatgccgctgcccggggtcggggaggagctgaagcaggccaaggagatcgaggacgccgagaagtactccttcatggccaccgtcaccaaggcgcccaagaagcaaatccagtttgctgatgacatgcaggagttcaccaaattccccaccaaaactggccgaagatctttgtctcgctcgatctcacagtcctccactgacagctacagttcagctgcatcctacacagatagctctgatgatgaggtttctccccgagagaagcagcaaaccaactccaagggcagcagcaatttctgtgtgaagaacatcaagcaggcagaatttggacgccgggagattgagattgcagagcaagacatgtctgctctgatttcactcaggaaacgtgctcagggggagaagcccttggctggtgctaaaatagtgggctgtacacacatcacagcccagacagcggtgttgattgagacactctgtgccctgggggctcagtgccgctggtctgcttgtaacatctactcaactcagaatgaagtagctgcagcactggctgaggctggagttgcagtgttcgcttggaagggcgagtcagaagatgacttctggtggtgtattgaccgctgtgtgaacatggatgggtggcaggccaacatgatcctggatgatgggggagacttaacccactgggtttataagaagtatccaaacgtgtttaagaagatccgaggcattgtggaagagagcgtgactggtgttcacaggctgtatcagctctccaaagctgggaagctctgtgttccggccatgaacgtcaatgattctgttaccaaacagaagtttgataacttgtactgctgccgagaatccattttggatggcctgaagaggaccacagatgtgatgtttggtgggaaacaagtggtggtgtgtggctatggtgaggtaggcaagggctgctgtgctgctctcaaagctcttggagcaattgtctacattaccgaaatcgaccccatctgtgctctgcaggcctgcatggatgggttcagggtggtaaagctaaatgaagtcatccggcaagtcgatgtcgtaataacttgcacaggaaataagaatgtagtgacacgggagcacttggatcgcatgaaaaacagttgtatcgtatgcaatatgggccactccaacacagaaatcgatgtgaccagcctccgcactccggagctgacgtgggagcgagtacgttctcaggtggaccatgtcatctggccagatggcaaacgagttgtcctcctggcagagggtcgtctactcaatttgagctgctccacagttcccacctttgttctgtccatcacagccacaacacaggctttggcactgatagaactctataatgcacccgaggggcgatacaagcaggatgtgtacttgcttcctaagaaaatggatgaatacgttgccagcttgcatctgccatcatttgatgcccaccttacagagctgacagatgaccaagcaaaatatctgggactcaacaaaaatgggccattcaaacctaattattacagatactaa');
ok($kegg->translate->seq);

is(($kegg->annotation->get_Annotations('description'))[0]->text,
   'S-adenosylhomocysteine hydrolase-like 1 [EC:3.3.1.1]');

is(($kegg->annotation->get_Annotations('pathway'))[0]->text,
   'Metabolism; Amino Acid Metabolism; Methionine metabolism');

is( (grep {$_->database eq 'KO'}
     $kegg->annotation->get_Annotations('dblink'))[0]->comment, 
    'adenosylhomocysteinase' );

is( (grep {$_->database eq 'PATH'} 
     $kegg->annotation->get_Annotations('dblink'))[0]->primary_id,
    'hsa00271' );

is( ($kegg->annotation->get_Annotations('aa_seq'))[0]->text,
'MSMPDAMPLPGVGEELKQAKEIEDAEKYSFMATVTKAPKKQIQFADDMQEFTKFPTKTGRRSLSRSISQSSTDSYSSAASYTDSSDDEVSPREKQQTNSKGSSNFCVKNIKQAEFGRREIEIAEQDMSALISLRKRAQGEKPLAGAKIVGCTHITAQTAVLIETLCALGAQCRWSACNIYSTQNEVAAALAEAGVAVFAWKGESEDDFWWCIDRCVNMDGWQANMILDDGGDLTHWVYKKYPNVFKKIRGIVEESVTGVHRLYQLSKAGKLCVPAMNVNDSVTKQKFDNLYCCRESILDGLKRTTDVMFGGKQVVVCGYGEVGKGCCAALKALGAIVYITEIDPICALQACMDGFRVVKLNEVIRQVDVVITCTGNKNVVTREHLDRMKNSCIVCNMGHSNTEIDVTSLRTPELTWERVRSQVDHVIWPDGKRVVLLAEGRLLNLSCSTVPTFVLSITATTQALALIELYNAPEGRYKQDVYLLPKKMDEYVASLHLPSFDAHLTELTDDQAKYLGLNKNGPFKPNYYRY');
