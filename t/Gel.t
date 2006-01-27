# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

use strict;

BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    plan test => 7 }

use Bio::PrimarySeq;
use Bio::Restriction::Analysis;
use Bio::Tools::Gel;
ok(1);



my $seq1 = Bio::PrimarySeq->new(-id=>'groundhog day',
                                -seq=>'AAAAAAAAAGAATTCTTTTTTTTTTTTTTGAATTCGGGGGGGGGGGGGGGGGGGG');


my $ra=Bio::Restriction::Analysis->new(-seq=>$seq1);
ok my @cuts = $ra->fragments('EcoRI'), 3;


ok my $gel = Bio::Tools::Gel->new(-seq=>\@cuts,-dilate=>10);
ok my %bands = $gel->bands;
my @bands = (26, 27, 30);
my $c = 0;
foreach my $band (sort {$b <=> $a} keys %bands){
    #print $band,"\t",  sprintf("%.1f", $bands{$band}), "\n";
    ok $bands[$c],  sprintf("%.0f", $bands{$band});
    $c++;
}

