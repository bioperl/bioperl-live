# -*-Perl-*-
## Bioperl Test Harness Script for Modules


use strict;
BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    use vars qw($NTESTS);
    $NTESTS = 9;
    plan tests => $NTESTS;
}
use Bio::Tools::Run::RepeatMasker;
use Bio::SeqIO;

END {
    for ( $Test::ntest..$NTESTS ) {
        skip("RepeatMasker program not found. Skipping. (Be sure you have the phrap package )",1);
    }
}
my @params=("mam" => 1,"noint"=>1);
my $fact = Bio::Tools::Run::RepeatMasker->new(@params);

ok ($fact->mam, 1);
ok ($fact->noint,1);
my $inputfilename= Bio::Root::IO->catfile("t","data","repeatmasker.fa");

my $in  = Bio::SeqIO->new(-file => "$inputfilename" , '-format' => 'fasta');
my $seq = $in->next_seq();
my @feats = $fact->mask($seq);
ok ($feats[0]->feature1->start, 1337);
ok ($feats[0]->feature1->end, 1407);
ok ($feats[0]->feature1->strand, 1);
ok ($feats[1]->feature1->start, 1712);
ok ($feats[1]->feature1->end, 2225);
my ($rn) = $feats[0]->each_tag_value('repeat_name');
my ($rc) = $feats[0]->each_tag_value('repeat_class');
ok ($rn,'(TTAGGG)n');
ok ($rc,'Simple_repeat');





