# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 5);
	
	use_ok('Bio::SeqIO::fastq');
}

my $DEBUG = test_debug();

# bug 2335
my $in_qual  = Bio::SeqIO->new('-file' => test_input_file('bug2335.fastq'),
							   '-format' => 'fastq',
							  );
isa_ok($in_qual, 'Bio::SeqIO');

my $qual = $in_qual->next_seq();
isa_ok($qual, 'Bio::Seq::Quality');

my @quals = @{$qual->qual()};

is(@quals, 111, 'number of qual values');

my $qualslice = join(',',@quals[0..10]);
is($qualslice, '31,23,32,23,31,22,27,28,32,24,25', 'qual slice');
