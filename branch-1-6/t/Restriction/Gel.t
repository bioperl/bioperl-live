# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 9);
	
    use_ok('Bio::PrimarySeq');
    use_ok('Bio::Restriction::Analysis');
    use_ok('Bio::Tools::Gel');
}

my $seq1 = Bio::PrimarySeq->new(-id=>'groundhog day',
                                -seq=>'AAAAAAAAAGAATTCTTTTTTTTTTTTTTGAATTCGGGGGGGGGGGGGGGGGGGG');

my $ra=Bio::Restriction::Analysis->new(-seq=>$seq1);
is my @cuts = $ra->fragments('EcoRI'), 3;

ok my $gel = Bio::Tools::Gel->new(-seq=>\@cuts,-dilate=>10);
ok my %bands = $gel->bands;
my @bands = (26, 27, 30);
my $c = 0;
foreach my $band (sort {$b <=> $a} keys %bands){
    #print $band,"\t",  sprintf("%.1f", $bands{$band}), "\n";
    is $bands[$c],  sprintf("%.0f", $bands{$band});
    $c++;
}
