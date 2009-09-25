# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 16,
               -requires_modules => [qw(IO::String
                                        LWP::UserAgent
                                        HTTP::Request::Common)]);
    
	use_ok('Bio::DB::RefSeq');
    use_ok('Bio::DB::GenBank');
    use_ok('Bio::DB::EMBL');
}

my $verbose = test_debug() || -1;

my ($db,$seq,$db2,$seq2,$seqio);
# get a single seq

$seq = $seqio = undef;

#test redirection from GenBank and EMBL
#GenBank
ok $db = Bio::DB::GenBank->new('-verbose'=> $verbose, -redirect_refseq => 1);
#EMBL
ok $db2 = Bio::DB::EMBL->new('-verbose'=> $verbose, -redirect_refseq => 1);     

eval {
    $seq = $db->get_Seq_by_acc('NT_006732');
    $seq2 = $db2->get_Seq_by_acc('NT_006732');
};
ok $@;

SKIP: {
    test_skip(-tests => 10, -requires_networking => 1);
    
    eval {
        ok($seq = $db->get_Seq_by_acc('NM_006732'));
        is($seq->length, 3776);
        ok $seq2 = $db2->get_Seq_by_acc('NM_006732');
        is($seq2->length, 3776);
    };
    skip "Warning: Couldn't connect to RefSeq with Bio::DB::RefSeq.pm!", 10 if $@;
    
    eval { 
        ok defined($db = Bio::DB::RefSeq->new(-verbose=>$verbose)); 
        ok(defined($seq = $db->get_Seq_by_acc('NM_006732')));
        is( $seq->length, 3776);
        ok defined ($db->request_format('fasta'));
        ok(defined($seq = $db->get_Seq_by_acc('NM_006732')));
        is( $seq->length, 3776); 
    };
    skip "Warning: Couldn't connect to RefSeq with Bio::DB::RefSeq.pm!", 6 if $@;
}
