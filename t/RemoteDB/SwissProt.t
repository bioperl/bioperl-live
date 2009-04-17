# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 19,
			   -requires_modules => [qw(IO::String
									    LWP::UserAgent
										HTTP::Request::Common)],
			   -requires_networking => 1);
	
	use_ok('Bio::DB::SwissProt');
}

ok my $gb = Bio::DB::SwissProt->new(-retrievaltype =>'pipeline',
                                 -delay => 0);

my %expected_lengths = (
                        'NDP_MOUSE' => 131, 
                        'NDP_HUMAN' => 133, 
                        'BOLA_HAEIN'=> 103, 
                        'YNB3_YEAST'=> 125, 
                        'O39869'    => 56,  
                        'DEGP_CHLTR'=> 497, 
                        );

my ($seq, $seqio);

SKIP: {
    eval {$seq = $gb->get_Seq_by_id('YNB3_YEAST');};
    skip "Couldn't connect to SwissProt with Bio::DB::Swiss.pm. Skipping those tests", 14 if $@;
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
    is $seq->division, 'YEAST';
    
    eval {$seq = $gb->get_Seq_by_acc('P43780');};
    skip "Couldn't connect to SwissProt with Bio::DB::Swiss.pm. Skipping those tests", 12 if $@;
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
    eval {$seq = $gb->get_Seq_by_acc('O39869');};
    skip "Couldn't connect to SwissProt with Bio::DB::Swiss.pm. Skipping those tests", 11 if $@;
    is $seq->length, $expected_lengths{$seq->accession_number}, $seq->accession_number;
    is $seq->accession_number, 'O39869';
    is $seq->division, '9PICO';
    
    # test for bug #958
    eval {$seq = $gb->get_Seq_by_id('P18584');};
    skip "Couldn't connect to SwissProt with Bio::DB::Swiss.pm. Skipping those tests", 8 if $@;
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
    is $seq->display_id, 'DEGP_CHLTR';
    is $seq->division, 'CHLTR';

    ok $gb = Bio::DB::SwissProt->new('-retrievaltype' => 'tempfile', '-delay' => 0);
    eval {$seqio = $gb->get_Stream_by_id(['NDP_MOUSE', 'NDP_HUMAN']);};
    skip "Couldn't connect to SwissProt with Bio::DB::Swiss.pm. Skipping those tests", 4 if $@;
    undef $gb; # testing to see if we can remove gb
    ok $seq = $seqio->next_seq();
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
    ok $seq = $seqio->next_seq();
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
}

# test idtracker() method

ok $gb = Bio::DB::SwissProt->new(-retrievaltype =>'pipeline',
                                 -delay => 0);

SKIP: {
    my $newid;
    # check old ID
    eval {$newid = $gb->idtracker('myod_pig');};
    skip("Couldn't connect to SwissProt with Bio::DB::Swiss.pm. Skipping those tests", 2) if $@;
    is($newid, 'MYOD1_PIG');
    # check ID that is current
    eval {$newid = $gb->idtracker('YNB3_YEAST');};
    skip("Couldn't connect to SwissProt with Bio::DB::Swiss.pm. Skipping those tests", 1) if $@;
    is($newid, 'YNB3_YEAST');    
}

1;
