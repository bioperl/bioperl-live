# -*-Perl-*- Test Harness script for Bioperl
# $Id$

# malcolm_cook@stowers.org: this test is in a separate file from
# RemoteBlast.t (on which it is modelled) since there is some sort of
# side-effecting between the multiple remote blasts that is causing
# this test to fail, if it comes last, or the other test to fail, if
# this one comes first.  THIS IS A BUG EITHER IN REMOTE BLAST OR MY
# UNDERSTANDING, i.e. of how to initialize it.

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests               => 7,
               -requires_modules    => [qw(IO::String
                                           LWP
                                           LWP::UserAgent)],
               -requires_networking => 1);

    use_ok('Bio::Tools::Run::RemoteBlast');
}

my $v = test_debug();
my $inputfilename = test_input_file('ecolitst.fa');
ok( -e $inputfilename);

my $remote_rpsblast
    = Bio::Tools::Run::RemoteBlast->new(-verbose    => test_debug(),
                                        -prog       => 'blastp',
                                        -data       => 'cdsearch/cdd',
                                        -readmethod => 'blasttable',
                                        -expect     => '1e-10',
                                        );
$remote_rpsblast->retrieve_parameter('ALIGNMENT_VIEW', 'Tabular');

# This is the key to getting job run using rpsblast:
$Bio::Tools::Run::RemoteBlast::HEADER{'SERVICE'} = 'rpsblast';

my $attempt = 1;

SKIP: {
    my $status;
    eval{
        $status = $remote_rpsblast->submit_blast($inputfilename);
    };

    ok($status,'rpsblast blasttable submitted');
    skip("Error accessing remote BLAST interface: $@", 4) if $@;

    my @rids = $remote_rpsblast->each_rid;
    is(@rids, 1, 'should only be one RID');
    skip("Wrong number of RIDs: ".scalar(@rids), 3) if @rids != 1;

    diag("Retrieving $rids[0]...\n") if $v;
    my $rc;
    while (defined($rc = $remote_rpsblast->retrieve_blast($rids[0]))) {
        if ( !ref($rc) ) {
            if ( $rc < 0 ) {
                skip("need a better solution for when 'Server failed to return any data'",2);
            }
            sleep 5;
            diag("Retrieval attempt: $attempt\n") if $v;
            $attempt++ < 10 ? next : last;
        }
        else {
            last;
        }
    }
    $remote_rpsblast->remove_rid($rids[0]);

    if ($rc) {
        ok(1,'retrieve_blast succeeded');
        $remote_rpsblast->remove_rid($rids[0]);
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
        cmp_ok($count, '>=', 45, 'HSPs returned');
    }
    elsif ($attempt > 10) {
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

# To be a good citizen, we should restore the default NCBI service
# ('plain') for future tests
$Bio::Tools::Run::RemoteBlast::HEADER{'SERVICE'} = 'plain';
