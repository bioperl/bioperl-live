# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#

use strict;
use vars qw($DEBUG);
use Bio::Root::IO;
$DEBUG = $ENV{'BIOPERLDEBUG'};

BEGIN {
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
    eval { require Test::More; };
    if( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
    plan tests => 9;
}
END  {
    unlink qw(write_phd.phd);
}

print("Checking if the Bio::SeqIO::phd module could be used, even though it shouldn't be directly used...\n") if ( $DEBUG);
        # test 1
use_ok('Bio::SeqIO::phd');

print("Checking to see if SeqWithQuality objects can be created from a file...\n") if ($DEBUG);
my $in_phd  = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile("t","data",
								"phredfile.phd"),
			      '-format'  => 'phd',
			      '-verbose' => $DEBUG || 0);
isa_ok($in_phd,'Bio::SeqIO::phd');

my @phreds;
print("I saw these in qualfile.qual:\n") if($DEBUG);
my $phd = $in_phd->next_seq();
print("Did you get the 'QUALITY_LEVELS' comment?\n") if ($DEBUG);
is($phd->{comments}->{'QUALITY_LEVELS'}, 99);
print("Checking to see if this is the right reference...\n") if( $DEBUG);
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

my $out_phd = Bio::SeqIO->new(-file => ">write_phd.phd",
			      '-format' => 'phd');
print("Did it return the right reference?\n") if($DEBUG);
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
ok( -e "write_phd.phd");

# Bug 2120

my @qual = q(9 9 12 12 8 8 9 8 8 8 9);
my @trace = q(113 121 130 145 153 169 177 203 210 218 234);

$in_phd  = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile("t","data",
								"bug2120.phd"),
			      '-format'  => 'phd',
			      '-verbose' => $DEBUG || 0);

my $seq = $in_phd->next_seq;
is($seq->subseq(10,20),'gggggccttat','$seq->subseq()');
my @seq_qual =$seq->subqual_text(10,20);
is_deeply(\@seq_qual,\@qual,'$seq->subqual_tex()');
my @seq_trace = $seq->subtrace_text(10,20);
is_deeply(\@seq_trace,\@trace,'$seq->subqual_tex()');
