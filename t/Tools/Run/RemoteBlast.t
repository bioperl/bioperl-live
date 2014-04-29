# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests               => 21,
               -requires_modules    => [qw(IO::String
                                           LWP
                                           LWP::UserAgent)],
               -requires_networking => 1);

    use_ok('Bio::Tools::Run::RemoteBlast');
}

my $prog  = 'blastp';
my $db    = 'swissprot';
my $e_val = '1e-10';
my $v     = test_debug();

SKIP: {
    my $remote_blast
        = Bio::Tools::Run::RemoteBlast->new(-verbose => $v,
                                            -prog    => $prog,
                                            -data    => $db,
                                            -expect  => $e_val,
                                            );
    $remote_blast->submit_parameter('ENTREZ_QUERY', 'Escherichia coli[ORGN]');

    my $inputfilename = test_input_file('ecolitst.fa');
    ok( -e $inputfilename);
    ok(1, 'Text BLAST');

    my $attempt = 1;
    my $status;
    eval{
        $status = $remote_blast->submit_blast($inputfilename);
    };

    ok($status,'BLAST text output submitted');
    skip("Error accessing remote BLAST interface: $@", 5) if $@;

    cmp_ok($remote_blast->get_rtoe, '>', 0, "Estimated retrieval time");
    diag("Time remaining: ".$remote_blast->get_rtoe) if $v;

    my @rids = $remote_blast->each_rid;
    is(@rids, 1, 'should only be one RID');
    skip("Wrong number of RIDs: ".scalar(@rids), 3) if @rids != 1;

    diag("Retrieving $rids[0]...\n") if $v;
    my $max_attempts = 20;
    my $rc;
    while (defined($rc = $remote_blast->retrieve_blast($rids[0]))) {
        if ( !ref($rc) ) {
            if ( $rc < 0 ) {
                skip("need a better solution for when 'Server failed to return any data'",3);
            }
            sleep 30;
            diag("Retrieval attempt: $attempt\n") if $v;
            $attempt++ < $max_attempts ? next : last;
        }
        else {
            last;
        }
    }

    if ($rc) {
        ok(1,'retrieve_blast succeeded');
        $remote_blast->remove_rid($rids[0]);
        my $count = 0;
        isa_ok($rc, 'Bio::SearchIO');
        while (my $result = $rc->next_result) {
            while ( my $hit = $result->next_hit ) {
                $count++;
                next unless ( $v > 0);
                print "sbjct name is ", $hit->name, "\n";
                while ( my $hsp = $hit->next_hsp ) {
                    print "score is ", $hsp->bits, "\n";
                }
            }
        }
        is($count, 3, 'HSPs returned');
    }
    elsif ($attempt > $max_attempts) {
        # have a test fail here (there should not be repeated failed attempts to
        # get reports)
        ok(0,'Exceeded maximum attempts on server to retrieve report');
        diag("Timeout, did not return report after ".($attempt - 1)." attempts");
        skip('Remote server timeout problems', 2);
    }
    else {
        # have a test fail here (whatever is returned should be eval as true and
        # be a SearchIO)
        ok(0,"Other problem on remote server, no report returned: $rc");
        skip('Possible remote server problems', 2);
    }
}

SKIP: {
    # test blasttable
    my $remote_blast
        = Bio::Tools::Run::RemoteBlast->new(-verbose    => $v,
                                            -prog       => $prog,
                                            -data       => $db,
                                            -readmethod => 'blasttable',
                                            -expect     => $e_val,
                                            );
    $remote_blast->submit_parameter('ENTREZ_QUERY', 'Escherichia coli[ORGN]');
    $remote_blast->retrieve_parameter('ALIGNMENT_VIEW', 'Tabular');

    my $inputfilename = test_input_file('ecolitst.fa');
    my $attempt = 1;

    my $status;
    eval{
        $status = $remote_blast->submit_blast($inputfilename);
    };

    ok($status,'Tabular BLAST submitted');
    skip("Error accessing remote BLAST interface: $@", 5) if $@;

    cmp_ok($remote_blast->get_rtoe, '>', 0, "Estimated retrieval time");
    diag("Time remaining: ".$remote_blast->get_rtoe) if $v;

    my @rids = $remote_blast->each_rid;
    is(@rids, 1, 'should only be one RID');
    skip("Wrong number of RIDs: ".scalar(@rids), 3) if @rids != 1;

    diag("waiting [$rids[0]]...") if $v;
    my $max_attempts = 20;
    my $rc;
    while (defined($rc = $remote_blast->retrieve_blast($rids[0]))) {
        if ( !ref($rc) ) {
            if ( $rc < 0 ) {
                skip("need a better solution for when 'Server failed to return any data'",3);
            }
            sleep 30;
            diag("Retrieval attempt: $attempt") if $v;
            $attempt++ < $max_attempts ? next : last;
        }
        else {
            last;
        }
    }

    if ($rc) {
        ok(1,'retrieve_blast succeeded');
        $remote_blast->remove_rid($rids[0]);
        my $count = 0;
        isa_ok($rc, 'Bio::SearchIO');
        while (my $result = $rc->next_result) {
            while ( my $hit = $result->next_hit ) {
                $count++;
                next unless ( $v > 0);
                print "sbjct name is ", $hit->name, "\n";
                while ( my $hsp = $hit->next_hsp ) {
                    print "score is ", $hsp->bits, "\n";
                }
            }
        }
        is($count, 3, 'HSPs returned');
    }
    elsif ($attempt > $max_attempts) {
        # have a test fail here (there should not be repeated failed attempts to
        # get reports)
        ok(0,'Exceeded maximum attempts on server to retrieve report');
        diag("Timeout, did not return report after ".($attempt - 1)." attempts");
        skip('Remote server timeout problems', 2);
    }
    else {
        # have a test fail here (whatever is returned should be eval as true and
        # be a SearchIO)
        ok(0,"Other problem on remote server, no report returned: $rc");
        skip('Possible remote server problems', 2);
    }
}

SKIP: {
    test_skip(-tests           => 6,
              -requires_module => 'Bio::SearchIO::blastxml');

    my $remote_blast
        = Bio::Tools::Run::RemoteBlast->new(-prog       => $prog,
                                            -data       => $db,
                                            -readmethod => 'xml',
                                            -expect     => $e_val,
                                            );
    $remote_blast->submit_parameter('ENTREZ_QUERY',
                                    'Escherichia coli[ORGN]');
    $remote_blast->retrieve_parameter('FORMAT_TYPE', 'XML');

    my $inputfilename = test_input_file('ecolitst.fa');
    my $attempt = 1;

    my $status;
    eval{
        $status = $remote_blast->submit_blast($inputfilename);
    };

    ok($status,'XML BLAST submitted');
    skip("Error accessing remote BLAST interface: $@", 5) if $@;

    cmp_ok($remote_blast->get_rtoe, '>', 0, "Estimated retrieval time");
    diag("Time remaining: ".$remote_blast->get_rtoe) if $v;

    my @rids = $remote_blast->each_rid;
    is(@rids, 1, 'should only be one RID');
    skip("Wrong number of RIDs: ".scalar(@rids), 3) if @rids != 1;

    diag("waiting [$rids[0]]...") if $v;
    # This blast usually takes 2x time compared to the other two
    my $max_attempts = 30;
    my $rc;
    while (defined($rc = $remote_blast->retrieve_blast($rids[0]))) {
        if ( !ref($rc) ) {
            if ( $rc < 0 ) {
                skip("need a better solution for when 'Server failed to return any data'",3);
            }
            sleep 30;
            diag("Retrieval attempt: $attempt") if $v;
            $attempt++ < $max_attempts ? next : last;
        }
        else {
            last;
        }
    }

    if ($rc) {
        ok(1,'retrieve_blast succeeded');
        $remote_blast->remove_rid($rids[0]);
        my $count = 0;
        isa_ok($rc, 'Bio::SearchIO');
        while (my $result = $rc->next_result) {
            while ( my $hit = $result->next_hit ) {
                $count++;
            }
        }
        is($count, 3, 'HSPs returned');
    }
    elsif ($attempt > $max_attempts) {
        # have a test fail here (there should not be repeated failed attempts to
        # get reports)
        ok(0,'Exceeded maximum attempts on server to retrieve report');
        diag("Timeout, did not return report after ".($attempt - 1)." attempts");
        skip('Remote server timeout problems', 2);
    }
    else {
        # have a test fail here (whatever is returned should be eval as true and
        # be a SearchIO)
        ok(0,"Other problem on remote server, no report returned: $rc");
        skip('Possible remote server problems', 2);
    }
}
