# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests           => 17,
               -requires_module => 'Data::Stag');

    use_ok('Bio::SeqIO');
    use_ok('Bio::Cluster::SequenceFamily');
}

my $seqio= Bio::SeqIO->new('-format' => 'swiss',
                           '-file'   => test_input_file('sequencefamily.dat'));
my @mem;
while(my $seq = $seqio->next_seq){
    push @mem, $seq;
}
my $family = Bio::Cluster::SequenceFamily->new(
    -family_id=>"Family_1",
    -description=>"SomeFamily",
    -annotation_score=>"100",
    -family_score=>"50",
    -version=>"1.0",
    -members=>\@mem,
);
is $family->description, "SomeFamily";
is $family->annotation_score,100;
is $family->size, 5;
is $family->family_id,"Family_1";
is $family->version, "1.0";

$family->add_members($mem[0]);
$family->add_members($mem[1]);
is $family->size, 7;
is $family->cluster_score, "50";
is $family->family_score, "50";

my @members = $family->get_members(-ncbi_taxid=>9606);

foreach my $mem(@members){
    is $mem->species->ncbi_taxid, 9606;
}

@members = $family->get_members(-binomial=>"Homo sapiens");

foreach my $mem(@members){
    is $mem->species->binomial, "Homo sapiens";
}


$family->flush_members();

is $family->size, 0;
