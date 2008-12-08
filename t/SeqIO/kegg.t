# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 16);
	
	use_ok('Bio::SeqIO::kegg');
}

my $verbose = test_debug();

my $io = Bio::SeqIO->new(-format => 'kegg',
								 -verbose => $verbose,
								 -file => test_input_file('AHCYL1.kegg'));
isa_ok($io, 'Bio::SeqIO');
my $kegg = $io->next_seq();
ok($kegg);
isa_ok($kegg, 'Bio::Seq::RichSeq');
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

is( ($kegg->annotation->get_Annotations('position'))[0]->text, '1:join(26813..26932,50794..50905,52974..53117,54122..54222,54657..
            54759,56523..56617,57185..57291,58104..58220,58427..58490,59255..
            59343,59706..59776,60133..60227,60312..60410,60811..60879,61308..
            61386,62491..62611,63434..63440)');

is($kegg->species->common_name , 'H.sapiens');