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
    plan tests => 12;
}

END {
    unlink qw(write_qual.qual );
}

print ("\$DEBUG is this: $DEBUG\n");

print("Checking if the Bio::SeqIO::Qual module could be used, even though it shouldn't be directly use'd...\n") if ( $DEBUG );
        # test 1
use Bio::SeqIO::qual;
ok(1);

print("Checking to see if PrimaryQual.pm can be used...\n") if ( $DEBUG );
use Bio::Seq::PrimaryQual;
ok(1);

print("Checking to see if PrimaryQual objects can be created from a file...\n") if ( $DEBUG );
my $in_qual  = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile("t","data",
								 "qualfile.qual"),
			       '-format' => 'qual');
ok(1);

my @quals;
print("I saw these in qualfile.qual:\n") if $DEBUG;
my $first = 1;
while ( my $qual = $in_qual->next_seq() ) {
		# ::dumpValue($qual);

    ok(1);
    @quals = @{$qual->qual()};
    if( $DEBUG ) {
	print($qual->id()."\n");
	
	print("(".scalar(@quals).") quality values.\n");
    }
    if( $first ) { 
	ok(@quals, 484);
    }
    $first = 0;
}

# in October 2004, Carlos Mauricio La Rota posted a problem with descriptions
# this routine is to test that

my @quals = 10..20;
# this one has a forced header
my $seq = new Bio::Seq::PrimaryQual(
                    -qual =>   \@quals,
                    -header   =>   "Hank is a good cat. I gave him a bath yesterday.");
my $out = new Bio::SeqIO(     -fh  =>   \*STDOUT,
                         -format   =>   'qual');
# yes, that works
# $out->write_seq($seq);
$seq->header('');
$seq->id('Hank1');
# yes, that works
# $out->write_seq($seq);

