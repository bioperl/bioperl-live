# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 18,
			   -requires_modules => [qw(IO::String
									    LWP::UserAgent
										HTTP::Request::Common)],
			   -requires_networking => 1);
	
	use_ok('Bio::DB::Query::GenBank');
	use_ok('Bio::DB::GenBank');
}

my %expected_lengths = (
                        'MUSIGHBA1' => 408,  
                        'AF303112'  => 1611, 
                        'AF041456'  => 1156, 
                        'AY080910'  => 798,  
                        'AY080909'  => 1042, 
                        'AF155220'  => 1172, 
                        'AF442768'  => 2547, 
                        );

my ($gb, $seq, $seqio, $seqin, $query);

# test query facility
ok $query = Bio::DB::Query::GenBank->new('-db'      => 'nucleotide',
                                         '-query'   => 'Onchocerca volvulus[Organism]',
                                         '-mindate' => '2002/1/1',
                                         '-maxdate' => '2002/12/31'), 'Bio::DB::Query::GenBank';
SKIP: {
    cmp_ok $query->count, '>', 0;
    my @ids = $query->ids;
    cmp_ok @ids, '>', 0;
    is @ids, $query->count;
    ok $gb = Bio::DB::GenBank->new('-delay' => 0);
    eval {$seqio = $gb->get_Stream_by_query($query);};
    skip "Couldn't connect to complete GenBank query tests. Skipping those tests", 5 if $@;
    my $done = 0;
    while (my $s = $seqio->next_seq) {
        is $s->length, $expected_lengths{$s->display_id}, $s->display_id;
        undef $gb; # test the case where the db is gone, 
        # but the pipeline should remain until seqio goes away
        $done++;
    }
    skip('No seqs returned', 5) if !$done;
    is $done, 1;
}

$seq = $seqio = undef;

# test query facility (again)
ok $query = Bio::DB::Query::GenBank->new('-db'  => 'nucleotide',
                                         '-ids' => [qw(J00522 AF303112 2981014)]);
SKIP: {
    cmp_ok $query->count, '>', 0;
    my @ids = $query->ids;
    cmp_ok @ids, '>', 0;
    is @ids, $query->count;
    $gb = Bio::DB::GenBank->new('-delay' => 0);
    eval {$seqio = $gb->get_Stream_by_query($query);};
    skip "Couldn't connect to complete GenBank query tests. Skipping those tests: $@", 4 if $@;
    my $done = 0;
    while (my $s = $seqio->next_seq) {
        is $s->length, $expected_lengths{$s->display_id}, $s->display_id;
        $done++;
    }
    skip('No seqs returned', 4) if !$done;
    is $done, 3;
    $seqio->close(); # the key to preventing errors during make test, no idea why
}

$seq = $seqio = undef;

# and yet again, for bug 2133
$query = Bio::DB::Query::GenBank->new('-query'  => 'AF303112',
                                      '-ids' => [qw(J00522 AF303112 2981014)]);
is $query->query, 'J00522[PACC]|AF303112[PACC]|2981014[UID]';
