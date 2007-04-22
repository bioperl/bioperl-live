# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
## $Id$

use strict;
use vars qw($NUMTESTS $DEBUG $ERROR);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

BEGIN {
    # to handle systems with no installed Test::More module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test::More; };
    $ERROR = 0;
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;

    $NUMTESTS = 11;
    plan tests => $NUMTESTS;
}

use_ok('Bio::Tools::SiRNA');
use_ok('Bio::Seq');
use_ok('Bio::SeqIO');

my $input = Bio::SeqIO->new( -file 	=> File::Spec->catfile(qw(t data NM_002254.gb)),
			     -format 	=> 'Genbank' );
my $seq = $input->next_seq;

isa_ok( $input, 'Bio::SeqIO' ) ;


#object creation
my $sirna = Bio::Tools::SiRNA->new( -target 	=> $seq,
                                     );
isa_ok( $sirna, 'Bio::Tools::SiRNA' ) ;

# first test - cds only
my @pairs = $sirna->design;
is ( scalar(@pairs), 65, "CDS only: got ". scalar(@pairs) );


# next test - include 3prime utr
my @feats = $seq->remove_SeqFeatures;
foreach my $feat (@feats) {
    $seq->add_SeqFeature($feat) unless
	($feat->primary_tag eq 'Target' or $feat->isa('Bio::SeqFeature::SiRNA::Pair'));
}
ok( $sirna->include_3pr(1) ) ;
@pairs = $sirna->design;
print "With 3p UTR: got ",scalar(@pairs),"\n" if $DEBUG;
is( scalar(@pairs), 140 );


#third test - naked sequence
my $newseq = Bio::Seq->new( -seq => $seq->seq);
isa_ok($newseq, 'Bio::Seq') ;

ok( $sirna->target($newseq) );
@pairs = $sirna->design;
print "Bare sequence: got ",scalar(@pairs),"\n" if $DEBUG;
is ( scalar(@pairs),  142 ) ;
