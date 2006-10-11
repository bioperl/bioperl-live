# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
## $Id$

use strict;
use vars qw($NUMTESTS $DEBUG $ERROR);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    $ERROR = 0;
    if( $@ ) {
	use lib 't';
    }
    use Test;

    $NUMTESTS = 7;
    plan tests => $NUMTESTS;
}

use Bio::Tools::SiRNA;
use Bio::Seq;
use Bio::SeqIO;


# modules compile
ok 1;

my $input = Bio::SeqIO->new( -file 	=> File::Spec->catfile(qw(t data NM_002254.gb)),
			     -format 	=> 'Genbank' );
my $seq = $input->next_seq;

#object creation
ok my $sirna = Bio::Tools::SiRNA->new( -target 	=> $seq,
                                     );

# first test - cds only
my @pairs = $sirna->design;
ok scalar(@pairs), 65, "CDS only: got ". scalar(@pairs) ;


# next test - include 3prime utr
my @feats = $seq->remove_SeqFeatures;
foreach my $feat (@feats) {
    $seq->add_SeqFeature($feat) unless
	($feat->primary_tag eq 'Target' or $feat->isa('Bio::SeqFeature::SiRNA::Pair'));
}
ok $sirna->include_3pr(1);
@pairs = $sirna->design;
print "With 3p UTR: got ",scalar(@pairs),"\n" if $DEBUG;
ok  scalar(@pairs), 140;


#third test - naked sequence
my $newseq = Bio::Seq->new( -seq => $seq->seq);
ok $sirna->target($newseq);
@pairs = $sirna->design;
print "Bare sequence: got ",scalar(@pairs),"\n" if $DEBUG;
ok scalar(@pairs),  142;
