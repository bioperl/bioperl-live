# -*-Perl-*- Test Harness script for Bioperl
# $Id: SeqWords.t 11525 2007-06-27 10:16:38Z sendu $

use strict;

BEGIN { 
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 5);
	
	use_ok('Bio::SeqIO');
	use_ok('Bio::Tools::SeqWords');
}

my $str = Bio::SeqIO->new(-file=> test_input_file('multifa.seq'), '-format' => 'Fasta');
my $seqobj= $str->next_seq();
ok defined $seqobj, 'new Bio::Root::IO object';

my $words = Bio::Tools::SeqWords->new('-seq' => $seqobj);
my $hash = $words->count_words(6);
ok (defined $words, 'new Bio::Tools::SeqWords object');
ok (defined $hash, 'count_words');
