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
    eval { require Test::More; };
    if( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
    plan tests => 18;
	use_ok('Bio::SeqIO::qual');
	use_ok('Bio::Seq::PrimaryQual');
}

END {
    unlink qw(write_qual1.qual write_qual2.qual);
}

warn("Checking to see if PrimaryQual objects can be created from a file...\n") if ( $DEBUG );
my $in_qual  = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile("t","data",
								 "qualfile.qual"),
			       '-format' => 'qual');
ok($in_qual);

my @quals;
warn("I saw these in qualfile.qual:\n") if $DEBUG;
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
my $seq = new Bio::Seq::PrimaryQual(
                    -qual =>   \@quals,
                    -header   =>   "Hank is a good cat. I gave him a bath yesterday.");
my $out = new Bio::SeqIO(-file  =>   '>write_qual2.qual',
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

