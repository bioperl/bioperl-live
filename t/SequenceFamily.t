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
    $NTESTS = 17;
    plan tests => $NTESTS;
}

use Bio::Root::IO;
use Bio::SeqIO;
use Bio::Cluster::SequenceFamily;


my $seqio= new Bio::SeqIO('-format' => 'swiss',
                           '-file'   => Bio::Root::IO->catfile('t','data','sequencefamily.dat'));
my @mem;
while(my $seq = $seqio->next_seq){
    push @mem, $seq;
}
my $family = Bio::Cluster::SequenceFamily->new(-family_id=>"Family_1",
                                       -description=>"SomeFamily",
                                       -annotation_score=>"100",
                                       -family_score=>"50",
                                       -version=>"1.0",
                                       -members=>\@mem);
ok $family->description, "SomeFamily";
ok $family->annotation_score,100;
ok $family->size, 5;
ok $family->family_id,"Family_1";
ok $family->version, "1.0";

$family->add_members($mem[0]);
$family->add_members($mem[1]);
ok $family->size, 7;
ok $family->cluster_score, "50";
ok $family->family_score, "50";

my @members = $family->get_members(-ncbi_taxid=>9606);

foreach my $mem(@members){
    ok $mem->species->ncbi_taxid, 9606;
}

@members = $family->get_members(-binomial=>"Homo sapiens");

foreach my $mem(@members){
    ok $mem->species->binomial, "Homo sapiens";
}


$family->flush_members();

ok $family->size, 0;








