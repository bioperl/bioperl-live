# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

use strict;
use vars qw($DEBUG $TESTCOUNT);
BEGIN {     
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    $TESTCOUNT = 101;
    plan tests => $TESTCOUNT;
}

use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::IO;

ok(1);

my $verbosity = -1;   # Set to -1 for release version, so warnings aren't printed

my ($seqio,$seq); # predeclare variables for strict

# default mode
$seqio = Bio::SeqIO->new('-file'=> Bio::Root::IO->catfile(
						   "t","data","test.genbank"), 
			 '-format' => 'GenBank');
ok $seqio;
$seqio->verbose($verbosity);

my $numseqs = 0;
my @loci = qw(U63596 U63595 M37762 NT_010368 L26462);
my @numfeas = (3,1,6,3,26);

while($seq = $seqio->next_seq()) {
    ok ($seq->accession_number, $loci[$numseqs++]);
    ok ($seq->annotation->get_Annotations());
    ok (scalar($seq->top_SeqFeatures), $numfeas[$numseqs-1]);
    ok ($seq->species->binomial);
    ok ($seq->seq);
}
ok ($numseqs, 5);

# minimalistic mode
$seqio = Bio::SeqIO->new('-file'=> Bio::Root::IO->catfile(
						   "t","data","test.genbank"), 
			 '-format' => 'GenBank');
ok $seqio;
$seqio->verbose($verbosity);
my $seqbuilder = $seqio->sequence_builder();
ok $seqbuilder;
ok $seqbuilder->isa("Bio::Factory::ObjectBuilderI");
$seqbuilder->want_none();
$seqbuilder->add_wanted_slot('display_id','accession_number','desc');

$numseqs = 0;

while($seq = $seqio->next_seq()) {
    ok ($seq->accession_number, $loci[$numseqs++]);
    ok (scalar($seq->annotation->get_Annotations()), 0);
    if($numseqs <= 3) {
	ok (scalar($seq->top_SeqFeatures), 0);
    } else {
	ok (scalar($seq->top_SeqFeatures), $numfeas[$numseqs-1]);
    }
    ok ($seq->species, undef);
    ok ($seq->seq, undef);
    # switch on features for the last 2 seqs
    $seqbuilder->add_wanted_slot('features') if $numseqs == 3;
}
ok ($numseqs, 5);

# everything but no sequence, and no features
$seqio = Bio::SeqIO->new('-file'=> Bio::Root::IO->catfile(
						   "t","data","test.genbank"), 
			 '-format' => 'GenBank');
ok $seqio;
$seqio->verbose($verbosity);
$seqbuilder = $seqio->sequence_builder();
# want-all is default
$seqbuilder->add_unwanted_slot('seq','features');

$numseqs = 0;

while($seq = $seqio->next_seq()) {
    ok ($seq->accession_number, $loci[$numseqs++]);
    ok scalar($seq->annotation->get_Annotations());
    if($numseqs <= 3) {
	ok (scalar($seq->top_SeqFeatures), 0);
    } else {
	ok (scalar($seq->top_SeqFeatures), $numfeas[$numseqs-1]);
    }
    ok $seq->species->binomial;
    ok ($seq->seq, undef);
    # switch on features for the last 2 seqs
    if($numseqs == 3) {
	$seqbuilder->add_unwanted_slot(
	     grep { $_ ne 'features'; } $seqbuilder->remove_unwanted_slots());
    }
}
ok ($numseqs, 5);

# skip sequences less than 100bp or accession like 'NT_*'
$seqio = Bio::SeqIO->new('-file'=> Bio::Root::IO->catfile(
						   "t","data","test.genbank"), 
			 '-format' => 'GenBank');
ok $seqio;
$seqio->verbose($verbosity);
$seqbuilder = $seqio->sequence_builder();
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

while($seq = $seqio->next_seq()) {
    $numseqs++;
    ok ($seq->accession_number, $loci[$i]);
    ok scalar($seq->annotation->get_Annotations());
    ok (scalar($seq->top_SeqFeatures), $numfeas[$i]);
    ok $seq->species->binomial;
    ok $seq->seq;
    $i += 2;
}
ok ($numseqs, 3);

