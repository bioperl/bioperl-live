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
    $NTESTS = 15;
    plan tests => $NTESTS;
}
use Bio::Tools::Run::Alignment::DBA;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Root::IO;

END {
    for ( $Test::ntest..$NTESTS ) {
        skip("dba program not found. Skipping. (Be sure you have the wise package > 2.2.0)",1);
    }
}

ok(1);
my $verbose = -1;
my @params = ('matchA' => 0.75, 'matchB' => '0.55','dymem'=>'linear');
my  $factory = Bio::Tools::Run::Alignment::DBA->new(@params);
ok $factory->isa('Bio::Tools::Run::Alignment::DBA');
my $bequiet = 1;
$factory->quiet($bequiet);  # Suppress clustal messages to terminal


#test with one file with 2 sequences
my $inputfilename_1a = Bio::Root::IO->catfile("t","data","dba1a.fa");
my $inputfilename_1b = Bio::Root::IO->catfile("t","data","dba1b.fa");
my $inputfilename2 = Bio::Root::IO->catfile("t","data","dba2.fa");
my $aln;
my $dba_present = $factory->exists_dba();
unless ($dba_present) {
    warn("Clustalw program not found. Skipping tests $Test::ntest to $NTESTS.\n");
    exit 0;
}
my @hsps = $factory->align($inputfilename2);
ok($hsps[0]->isa("Bio::Search::HSP::GenericHSP"));
ok($hsps[0]->query->start,4);
ok($hsps[0]->query->end,209);
ok($hsps[0]->gaps,6);


#test with 2 files of 1 sequence each
my @files = ($inputfilename_1a,$inputfilename_1b);
@hsps = $factory->align(\@files);
ok($hsps[0]->query->start,3);
ok($hsps[0]->query->end,88);
ok($hsps[0]->gaps,0);
ok($hsps[1]->hit->start,90);
ok($hsps[1]->hit->end,195);
ok($hsps[1]->gaps,0);

#test with an array of 2 PrimarySeqI objects

my $str = Bio::SeqIO->new(-file=> Bio::Root::IO->catfile("t","data","dba2.fa"), '-format' => 'Fasta');
my @seq_array =();

while ( my $seq = $str->next_seq() ) {
  push (@seq_array, $seq) ;
}
@hsps = $factory->align(\@seq_array);
ok($hsps[0]->query->start,4);
ok($hsps[0]->query->end,209);
ok($hsps[0]->gaps,6);




