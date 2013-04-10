# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    test_begin(-tests => 137);
    use_ok('Bio::SeqIO');
}

my $verbosity = test_debug;

my ($seqio, $seq); # predeclare variables for strict

# default mode
ok $seqio = Bio::SeqIO->new('-file' => test_input_file('test.genbank'), 
                            '-format' => 'GenBank');
$seqio->verbose($verbosity);

my $numseqs = 0;
my @loci = qw(U63596 U63595 M37762 NT_010368 L26462);
my @numfeas = (3,1,6,3,26);

while ($seq = $seqio->next_seq) {
    is $seq->accession_number, $loci[$numseqs++];
    ok $seq->annotation->get_Annotations;
    is scalar($seq->top_SeqFeatures), $numfeas[$numseqs-1];
    ok $seq->species->binomial;
    ok $seq->seq;
    ok $seq->desc;
    ok $seq->id;
}
is $numseqs, 5;

# minimalistic mode
$seqio = Bio::SeqIO->new('-file' => test_input_file('test.genbank'), 
                         '-format' => 'GenBank');
ok $seqio;
$seqio->verbose($verbosity);
ok my $seqbuilder = $seqio->sequence_builder;
isa_ok $seqbuilder, "Bio::Factory::ObjectBuilderI";
$seqbuilder->want_none;
$seqbuilder->add_wanted_slot('display_id','accession_number','desc');

$numseqs = 0;

while ($seq = $seqio->next_seq) {
    is $seq->accession_number, $loci[$numseqs++];
    is scalar(grep { ! ($_->tagname eq "keyword" ||
                        $_->tagname eq "date_changed" ||
                        $_->tagname eq "secondary_accession"); }
               $seq->annotation->get_Annotations), 0;
    if ($numseqs <= 3) {
        is scalar($seq->top_SeqFeatures), 0;
    }
    else {
        is scalar($seq->top_SeqFeatures), $numfeas[$numseqs-1];
    }
    is $seq->species, undef;
    is $seq->seq, undef;
    ok $seq->desc;
    ok $seq->id;
    # switch on features for the last 2 seqs
    $seqbuilder->add_wanted_slot('features') if $numseqs == 3;
}
is $numseqs, 5;

# everything but no sequence, and no features
$seqio = Bio::SeqIO->new('-file' => test_input_file('test.genbank'), 
                         '-format' => 'GenBank');
ok $seqio;
$seqio->verbose($verbosity);
$seqbuilder = $seqio->sequence_builder;
# want-all is default
$seqbuilder->add_unwanted_slot('seq','features');

$numseqs = 0;

while ($seq = $seqio->next_seq) {
    is $seq->accession_number, $loci[$numseqs++];
    ok scalar($seq->annotation->get_Annotations);
    if ($numseqs <= 3) {
        is scalar($seq->top_SeqFeatures), 0;
    }
    else {
        is scalar($seq->top_SeqFeatures), $numfeas[$numseqs-1];
    }
    ok $seq->species->binomial;
    is $seq->seq, undef;
    ok $seq->desc;
    ok $seq->id;
    # switch on features for the last 2 seqs
    if ($numseqs == 3) {
        $seqbuilder->add_unwanted_slot(
            grep { $_ ne 'features'; } $seqbuilder->remove_unwanted_slots
        );
    }
}
is $numseqs, 5;

# skip sequences less than 100bp or accession like 'NT_*'
$seqio = Bio::SeqIO->new('-file' => test_input_file('test.genbank'), 
                         '-format' => 'GenBank');
ok $seqio;
$seqio->verbose($verbosity);
$seqbuilder = $seqio->sequence_builder;
# we could have as well combined the two conditions into one, but we want to
# test the implicit AND here
$seqbuilder->add_object_condition(sub {
    my $h = shift;
    return 0 if($h->{'-length'} < 100);
    return 1;
});
$seqbuilder->add_object_condition(sub {
    my $h = shift;
    return 0 if($h->{'-display_id'} =~ /^NT_/);
    return 1;
});

$numseqs = 0;
my $i = 0;

while ($seq = $seqio->next_seq) {
    $numseqs++;
    is $seq->accession_number, $loci[$i];
    ok scalar($seq->annotation->get_Annotations);
    is scalar($seq->top_SeqFeatures), $numfeas[$i];
    ok $seq->species->binomial;
    ok $seq->seq;
    ok $seq->desc;
    ok $seq->id;
    $i += 2;
}
is $numseqs, 3;
