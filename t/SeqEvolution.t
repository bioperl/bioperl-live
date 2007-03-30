# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test::More.pm is located)
    # as a fallback
    eval { require Test::More; };
    if( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
    use vars qw($NUMTESTS $DEBUG);

    $NUMTESTS = 39;

    $DEBUG   = $ENV{'BIOPERLDEBUG'} || 0;
    plan tests => $NUMTESTS;
    use_ok('Bio::SeqEvolution::Factory');
    use_ok('Bio::PrimarySeq');
}

#
# point mutations
#

use Data::Dumper;

ok my $evolve = Bio::SeqEvolution::Factory->new (-verbose => 1);
ok ref $evolve eq 'Bio::SeqEvolution::DNAPoint', 'ref';
#print Dumper $evolve;

ok $evolve = Bio::SeqEvolution::DNAPoint->new (-verbose => 2);
ok ref $evolve eq 'Bio::SeqEvolution::DNAPoint', 'ref';
#print Dumper $evolve;

ok $evolve = Bio::SeqEvolution::Factory->new (-verbose => 3,
                                             -type => 'Bio::SeqEvolution::DNAPoint');
ok ref $evolve  eq 'Bio::SeqEvolution::DNAPoint', 'ref';
#print Dumper $evolve;


my $seq = new Bio::PrimarySeq(-id=>'test',
                              -seq=>'aaacccgggtta'
                             );

#ok $evolve = Bio::SeqEvolution::DNAPoint->new (-verbose => 0,
ok $evolve = Bio::SeqEvolution::Factory->new (-verbose => 0,
                                              -foo => 'bar',
                                              -rate => 5,
                                              -seq => $seq
                                             );

ok $evolve->seq->id, 'test';
#print Dumper $evolve;
ok $evolve->rate == 5,               'get rate()';
ok $evolve->rate(2) == 2,            'get and set rate()';


# next_seq
# each_seq(3);

ok $evolve->identity(80) == 80, 'identity()';
is $evolve->identity(undef), undef, 'identity()';

is $evolve->pam, undef, 'pam()';
ok $evolve->pam(80) == 80, 'pam()';

is $evolve->mutation_count, undef, 'mutation_count()';
ok $evolve->mutation_count(5) == 5, , 'mutation_count()';



ok $evolve->seq_type eq 'Bio::PrimarySeq', 'seq_type()';
ok $evolve->seq_type('Bio::Seq') eq 'Bio::Seq', 'seq_type()';

ok my $newseq = $evolve->next_seq, 'next_seq()';

#print Dumper $evolve;
#print Dumper $newseq;

foreach ( $evolve->each_mutation('string')) {
    #print Dumper $_;
    ok $_;
}
ok $evolve->each_mutation == 5, 'each_mutation()';

ok $evolve = Bio::SeqEvolution::Factory->new (-verbose => 0,
                                              -rate => 5,
                                              -seq => $seq,
                                              -identity => 50     ###
                                             );
ok $newseq = $evolve->next_seq;
ok $evolve->get_alignment_identity <= 50, "get_alignment_identity()";

ok $evolve = Bio::SeqEvolution::Factory->new (-verbose => 0,
                                              -rate => 5,
                                              -seq => $seq,
                                              -pam => 50     ###
                                             );
ok $newseq = $evolve->next_seq;
ok $evolve->get_mutation_counter == 6, 'get_mutation_counter()';
$newseq = $evolve->next_seq;
ok $evolve->get_sequence_counter == 2, 'get_sequence_counter()';
ok $evolve->reset_sequence_counter(), 'reset_sequence_counter()';
ok $evolve->get_sequence_counter == 0, 'get_sequence_counter() == 0';

ok $evolve = Bio::SeqEvolution::Factory->new (-verbose => 0,
                                              -rate => 5,
                                              -seq => $seq,
                                              -pam => 50
                                             );
#print Dumper $evolve;
ok $newseq = $evolve->next_seq;
# ok $evolve->get_alignment_identity <= 50, "get_alignment_identity()";
ok ref $evolve->get_alignmet eq 'Bio::SimpleAlign', 'get_alignment()';
#
#print "-----------------------------------------\n";
#print $evolve->{'_ori_locatableseq'}->seq, "\n";
#print $evolve->{'_mut_locatableseq'}->seq, "\n";
#print $evolve->{'_align'}->overall_percentage_identity, "\n";
#print $evolve->get_mutation_counter, "\n";
#print "-----------------------------------------\n";
#






exit;

#
# indels
#

#use Bio::SeqEvolution::DNAIndel;
#ok $evolve = Bio::SeqEvolution::DNAIndel->new (-verbose => 0,
#                                               -mutation_count => 3,
#                                               -rate => 1,
#                                               -seq => $seq
#                                             );
#
#
#
#ok $newseq = $evolve->next_seq;
#ok $evolve->get_mutation_counter, 3;
##print Dumper $evolve;
#
#foreach ( $evolve->each_mutation) {
##    print Dumper $_;
##    print $_->sysname, "\n";
#    ok $_->sysname;
#}
#
#
#ok $evolve = Bio::SeqEvolution::DNAIndel->new (-verbose => 0,
#                                               -duplication => 1,
##                                               -identity =>50,
#                                               -mutation_count => 3,
#                                               -rate => 1,
#                                               -seq => $seq
#                                             );
##$evolve->{'_mut_string'} = $evolve->{'_ori_string'};
#
##ok $newseq = $evolve->mutate(12);
#ok $newseq = $evolve->next_seq;
##print Dumper $evolve;
##print $evolve->{'_ori_locatableseq'}->seq, "\n";
##print $evolve->{'_mut_locatableseq'}->seq, "\n";
#
#print "-----------------------------------------\n";
#print $evolve->{'_ori_locatableseq'}->seq, "\n";
#print $evolve->{'_mut_locatableseq'}->seq, "\n";
#print $evolve->{'_align'}->overall_percentage_identity, "\n";
#print "-----------------------------------------\n";
#
#
##
## mixer, simulation
##
#
#
#ok my $evolve2 = Bio::SeqEvolution::Factory->new (-verbose => 0,
#                                                  -rate => 2,
#                                                  -seq => $seq,
#                                                  -set_mutated_seq => $newseq,
##                                                  -identity => 50
#                                                  -mutation_count => 1     ###
#                                             );
#ok $evolve2->seq_type('Bio::LocatableSeq');
#
#print "=-----------------------------------------\n";
#print $evolve2->{'_ori_locatableseq'}->seq, "\n";
#print $evolve2->{'_mut_locatableseq'}->seq, "\n";
#print $evolve2->{'_align'}->overall_percentage_identity, "\n";
#print "-----------------------------------------\n";
#
#ok $newseq = $evolve2->next_seq;
#
#
#print "-----------------------------------------\n";
#print $evolve2->{'_ori_locatableseq'}->seq, "\n";
#print $evolve2->{'_mut_locatableseq'}->seq, "\n";
#print $evolve2->{'_align'}->overall_percentage_identity, "\n";
#print "-----------------------------------------\n";
#
