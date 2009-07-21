# -*-Perl-*- Test Harness script for Bioperl
# 	$Id$	

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

    test_begin(-tests => 3,
               -requires_modules => [qw(IO::String LWP LWP::UserAgent)],
               -requires_networking => 1);

    use_ok('Bio::Tools::Run::RemoteBlast');
}

my $v = test_debug();
my $inputfilename = test_input_file('ecolitst.fa');
ok( -e $inputfilename);	

my $remote_rpsblast = Bio::Tools::Run::RemoteBlast->new
  ('-verbose'    => test_debug(),
   '-prog'       => 'blastp',
   '-data'       => 'cdsearch/cdd',
   '-readmethod' => 'blasttable',
   '-expect'     => '1e-10',
  );

$remote_rpsblast->retrieve_parameter('ALIGNMENT_VIEW', 'Tabular');

# This is the key to getting job run using rpsblast:
$Bio::Tools::Run::RemoteBlast::HEADER{'SERVICE'} = 'rpsblast'; 

ok($remote_rpsblast->submit_blast($inputfilename),'rpsblast blasttable submitted');

print STDERR "waiting..." if( $v > 0 );
while ( my @rids = $remote_rpsblast->each_rid ) {
  foreach my $rid ( @rids ) {
    my $rc = $remote_rpsblast->retrieve_blast($rid);
    if ( !ref($rc) ) {
      if ( $rc < 0 ) {
	die "need a better solution for when 'Server failed to return any data'";
      }
      $remote_rpsblast->remove_rid($rid);
      print STDERR "." if ( $v > 0 );
      sleep 5;
    } else {
      ok(1,'retrieve_blast succeeded');
      $remote_rpsblast->remove_rid($rid);
      my $count = 0;
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
      is($count, 44, 'correct result count'); # of course, this could change whenever CDD changes
    }
  }
}

# To be a good citizen, we should restore the default NCBI service
# ('plain') for future tests
$Bio::Tools::Run::RemoteBlast::HEADER{'SERVICE'} = 'plain';
