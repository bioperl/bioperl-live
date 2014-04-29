# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests               => 23,
               -requires_modules    => [qw(IO::String
                                           LWP::UserAgent
                                           HTTP::Request::Common)],
               -requires_networking => 1);

    use_ok('Bio::DB::SwissProt');
}

ok my $gb = Bio::DB::SwissProt->new(-retrievaltype => 'pipeline',
                                    -delay         => 0);

my %expected_lengths = (
                        'NDP_MOUSE'   => 131,
                        'NDP_HUMAN'   => 133,
                        'BOLA_HAEIN'  => 103,
                        'YNB3_YEAST'  => 125,
                        'O39869'      => 56,
                        'DEGP_CHLTR'  => 497,
                        'DEGPL_CHLTR' => 497
                        );

my ($seq, $seqio);

SKIP: {
    eval {$seq = $gb->get_Seq_by_id('YNB3_YEAST');};
    skip "Couldn't connect to SwissProt with Bio::DB::SwissProt.pm. Skipping those tests", 14 if $@;
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
    is $seq->division, 'YEAST';

    eval {$seq = $gb->get_Seq_by_acc('P43780');};
    skip "Couldn't connect to SwissProt with Bio::DB::SwissProt.pm. Skipping those tests", 12 if $@;
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
    eval {$seq = $gb->get_Seq_by_acc('O39869');};
    skip "Couldn't connect to SwissProt with Bio::DB::SwissProt.pm. Skipping those tests", 11 if $@;
    is $seq->length, $expected_lengths{$seq->accession_number}, $seq->accession_number;
    is $seq->accession_number, 'O39869';
    is $seq->division, '9PICO';

    # test for bug #958
    eval {$seq = $gb->get_Seq_by_id('P18584');};
    skip "Couldn't connect to SwissProt with Bio::DB::SwissProt.pm. Skipping those tests", 8 if $@;
    ok exists $expected_lengths{$seq->display_id}, 'P18584';
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
    is $seq->division, 'CHLTR';

    ok $gb = Bio::DB::SwissProt->new('-retrievaltype' => 'tempfile', '-delay' => 0);
    eval {$seqio = $gb->get_Stream_by_id(['NDP_MOUSE', 'NDP_HUMAN']);};
    skip "Couldn't connect to SwissProt with Bio::DB::SwissProt.pm. Skipping those tests", 4 if $@;
    undef $gb; # testing to see if we can remove gb
    ok $seq = $seqio->next_seq();
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
    ok $seq = $seqio->next_seq();
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
}

# test idtracker() method
ok $gb = Bio::DB::SwissProt->new(-retrievaltype => 'pipeline',
                                 -delay         => 0,
                                 -verbose       => 2);

SKIP: {
    my $map;
    # check old ID
    eval {$map = $gb->id_mapper(-from => 'ACC+ID',
                                -to   => 'ACC',
                                -ids  => [qw(MYOD1_PIG PYRC_YEAST)]
                                )};
    skip("Problem with idtracker(), skipping these tests: $@", 6) if $@;

    cmp_ok(@{$map->{MYOD1_PIG}}, '>=', 1);
    is($map->{MYOD1_PIG}[0], 'P49811');
    cmp_ok(@{$map->{PYRC_YEAST}}, '>=', 1);
    is($map->{PYRC_YEAST}[0], 'P20051');

    eval {$map = $gb->id_mapper(-from => 'ACC+ID',
                                -to   => 'EMBL',
                                -ids  => [qw(PYRC_YEAST)]
                                )};
    skip("Problem with idtracker(), skipping these tests: $@", 2) if $@;

    cmp_ok(@{$map->{PYRC_YEAST}}, '>=', 2);
    like($map->{PYRC_YEAST}[0], qr/^[A-Z0-9]/);
}

1;
