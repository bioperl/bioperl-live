# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#


use strict;
use vars qw($DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'};

BEGIN {
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan tests => 5;
}
END  {
    unlink qw(write_phd.phd);
}

print("Checking if the Bio::SeqIO::phd module could be used, even though it shouldn't be directly use'd...\n") if ( $DEBUG);
        # test 1
use Bio::SeqIO::phd;
ok(1);

print("Checking to see if SeqWithQuality objects can be created from a file...\n") if ($DEBUG);
my $in_phd  = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile("t","data",
								"phredfile.phd"),
			      '-format'  => 'phd',
			      '-verbose' => $DEBUG || 0);
ok(1);

my @phreds;
print("I saw these in qualfile.qual:\n") if($DEBUG);
my $phd = $in_phd->next_phd();
print("Checking to see if this is the right reference...\n") if( $DEBUG);
ok(ref($phd) eq "Bio::Seq::SeqWithQuality");

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
ok(ref($out_phd) eq "Bio::SeqIO::phd");

$out_phd->write_phd(	-SeqWithQuality		=>	$phd,
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
