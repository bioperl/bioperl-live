# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 9);
	
	use_ok('Bio::SeqIO');
}

my $DEBUG = test_debug();

print("Checking to see if SeqWithQuality objects can be created from a file...\n") if ($DEBUG);
my $in_phd  = Bio::SeqIO->new('-file' => test_input_file('phredfile.phd'),
			      '-format'  => 'phd',
			      '-verbose' => $DEBUG);
isa_ok($in_phd,'Bio::SeqIO::phd');

my @phreds;
my $phd = $in_phd->next_seq();
is($phd->{comments}->{'QUALITY_LEVELS'}, 99, "Did you get the 'QUALITY_LEVELS' comment?");
isa_ok($phd,"Bio::Seq::Quality");

my $position = 6;

if( $DEBUG ) {
    print("What is the base at position $position (using subseq)?\n");
    print($phd->subseq($position,$position)."\n");
    print("What is the base at position $position (using baseat)?\n");
    print($phd->baseat($position)."\n");
    print("What is the quality at $position? (using subqual)\n");
}
my @qualsretr = @{$phd->subqual($position,$position)};
if( $DEBUG ) {
    print($qualsretr[0]."\n");
    print("What is the quality at $position? (using qualat)\n");
    print($phd->qualat($position)."\n");
}

print("OK. Now testing write_phd...\n") if($DEBUG);

my $outfile = test_output_file();
my $out_phd = Bio::SeqIO->new(-file => ">$outfile",
			      '-format' => 'phd');
isa_ok($out_phd,"Bio::SeqIO::phd");

$out_phd->write_seq(	-SeqWithQuality		=>	$phd,
			-CHROMAT_FILE		=>	$phd->id(),
			-ABI_THUMBPRINT		=>	"",
			-PHRED_VERSION		=>	"",
			-CALL_METHOD		=>	"",
			-QUALITY_LEVELS		=>	"",
			-TIME			=>	"",
			-TRACE_ARRAY_MIN_INDEX	=>	"",
			-TRACE_ARRAY_MAX_INDEX	=>	"",
			-CHEM	=>	"",
			-DYE	=>	""	
			);
ok( -s $outfile);

# Bug 2120

my @qual = q(9 9 12 12 8 8 9 8 8 8 9);
my @trace = q(113 121 130 145 153 169 177 203 210 218 234);

$in_phd  = Bio::SeqIO->new('-file' => test_input_file('bug2120.phd'),
			      '-format'  => 'phd',
			      '-verbose' => $DEBUG);

my $seq = $in_phd->next_seq;
is($seq->subseq(10,20),'gggggccttat','$seq->subseq()');
my @seq_qual =$seq->subqual_text(10,20);
is_deeply(\@seq_qual,\@qual,'$seq->subqual_tex()');
my @seq_trace = $seq->subtrace_text(10,20);
is_deeply(\@seq_trace,\@trace,'$seq->subqual_tex()');
