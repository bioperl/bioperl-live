# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 21);
	
	use_ok('Bio::SeqIO');
	use_ok('Bio::Seq::PrimaryQual');
}

my $DEBUG = test_debug();

my $in_qual  = Bio::SeqIO->new('-file' => test_input_file('qualfile.qual'),
			       '-format' => 'qual');
ok($in_qual);

my @quals;

my $first = 1;
while ( my $qual = $in_qual->next_seq() ) {
		# ::dumpValue($qual);
	isa_ok($qual, 'Bio::Seq::PrimaryQual');
    @quals = @{$qual->qual()};
    if( $DEBUG ) {
	warn($qual->id()."\n");
	
	warn("(".scalar(@quals).") quality values.\n");
    }
    if( $first ) { 
		is(@quals, 484);
    }
    $first = 0;
}

# in October 2004, Carlos Mauricio La Rota posted a problem with descriptions
# this routine is to test that

@quals = 10..20;
# this one has a forced header
my $seq = Bio::Seq::PrimaryQual->new(
                    -qual =>   \@quals,
                    -header   =>   "Hank is a good cat. I gave him a bath yesterday.");
my $out = Bio::SeqIO->new(-file  =>   ">".test_output_file(),
                         -format   =>   'qual');
# yes, that works
is $seq->header, 'Hank is a good cat. I gave him a bath yesterday.';
@quals = @{$seq->qual()};
is scalar(@quals), 11;
ok $out->write_seq($seq);
$seq->header('');
is $seq->header, '';
$seq->id('Hank1');
is $seq->id, 'Hank1';
# yes, that works
ok $out->write_seq($seq);

# bug 2335

my $in_qual  = Bio::SeqIO->new('-file' => test_input_file('bug2335.fastq'),
			       '-format' => 'fastq');
ok($in_qual);

my @quals;

my $qual = $in_qual->next_seq();
isa_ok($qual, 'Bio::Seq::Quality');
@quals = @{$qual->qual()};
if( $first ) { 
	is(@quals, 111);
}
my $qualslice = join(',',@quals[0..10]);
is($qualslice, '31,23,32,23,31,22,27,28,32,24,25');
