# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 16,
               -requires_modules => [qw(IO::String LWP LWP::UserAgent)],
               -requires_networking => 1);
    
    use_ok('Bio::Tools::Run::RemoteBlast');
}

my $prog = 'blastp';
my $db   = 'swissprot';
my $e_val= '1e-10';
my $v = test_debug();
my $remote_blast = Bio::Tools::Run::RemoteBlast->new('-verbose' => $v,
													'-prog' => $prog,
													'-data' => $db,
													'-expect' => $e_val,
							  );
$remote_blast->submit_parameter('ENTREZ_QUERY', 
										  'Escherichia coli[ORGN]');
my $inputfilename = test_input_file('ecolitst.fa');
ok( -e $inputfilename);	

ok(1, 'Text BLAST');

my $r = $remote_blast->submit_blast($inputfilename);
ok($r);
print STDERR "waiting..." if( $v > 0 );
while ( my @rids = $remote_blast->each_rid ) {
	foreach my $rid ( @rids ) {
		my $rc = $remote_blast->retrieve_blast($rid);
		if( !ref($rc) ) {
			if( $rc < 0 ) { 		
				$remote_blast->remove_rid($rid);
				# need a better solution for when 'Server failed to return any data'
			}
			print STDERR "." if ( $v > 0 );
			sleep 5;
		} else { 
			ok(1);
			$remote_blast->remove_rid($rid);
			my $result = $rc->next_result;
			like($result->database_name, qr/swissprot/i);
			my $count = 0;
			while( my $hit = $result->next_hit ) {		
				$count++;
				next unless ( $v > 0);
				print "sbjct name is ", $hit->name, "\n";
				while( my $hsp = $hit->next_hsp ) {
					print "score is ", $hsp->score, "\n";
				} 
			}
			is($count, 3);
		}
	}
}

# test blasttable

ok(1, 'Tabular BLAST');

my $remote_blast2 = Bio::Tools::Run::RemoteBlast->new
  ('-verbose'    => $v,
	'-prog'       => $prog,
	'-data'       => $db,
	'-readmethod' => 'blasttable',
	'-expect'     => $e_val,
  );
$remote_blast2->submit_parameter('ENTREZ_QUERY', 'Escherichia coli[ORGN]');

$remote_blast2->retrieve_parameter('ALIGNMENT_VIEW', 'Tabular');

$inputfilename = test_input_file('ecolitst.fa');

$r = $remote_blast2->submit_blast($inputfilename);
ok($r);
print STDERR "waiting..." if( $v > 0 );
while ( my @rids = $remote_blast2->each_rid ) {
	foreach my $rid ( @rids ) {
		my $rc = $remote_blast2->retrieve_blast($rid);
		if( !ref($rc) ) {
			if( $rc < 0 ) { 		
				$remote_blast2->remove_rid($rid);
				# need a better solution for when 'Server failed to return any data'
			}
			print STDERR "." if ( $v > 0 );
			sleep 5;
		} else { 
			ok(1);
			$remote_blast2->remove_rid($rid);
			my $count = 0;
			while (my $result = $rc->next_result) {
				while( my $hit = $result->next_hit ) {		
					$count++;
					next unless ( $v > 0);
					print "sbjct name is ", $hit->name, "\n";
					while( my $hsp = $hit->next_hsp ) {
						print "score is ", $hsp->score, "\n";
					} 
				}
			}
			is($count, 3);
		}
	}
}

SKIP: {
    test_skip(-tests => 5, -requires_module => 'Bio::SearchIO::blastxml');
	
	my $remote_blastxml = Bio::Tools::Run::RemoteBlast->new('-prog' => $prog,
		'-data'       => $db,
		'-readmethod' => 'xml',
		'-expect'     => $e_val,
	);
	$remote_blastxml->submit_parameter('ENTREZ_QUERY', 
									'Escherichia coli[ORGN]');
	
	$remote_blastxml->retrieve_parameter('FORMAT_TYPE', 'XML');
	$inputfilename = test_input_file('ecolitst.fa');
	
	ok(1, 'XML BLAST');
	
	$r = $remote_blastxml->submit_blast($inputfilename);
	ok($r);
	print STDERR "waiting..." if( $v > 0 );
	while ( my @rids = $remote_blastxml->each_rid ) {
		foreach my $rid ( @rids ) {
			my $rc = $remote_blastxml->retrieve_blast($rid);
			if( !ref($rc) ) {
				if( $rc < 0 ) { 		
					$remote_blastxml->remove_rid($rid);
					# need a better solution for when 'Server failed to return any data'
				}
				print STDERR "." if ( $v > 0 );
				sleep 5;
			} else { 
				ok(1);
				$remote_blastxml->remove_rid($rid);
				my $result = $rc->next_result;
				like($result->database_name, qr/swissprot/i);
				my $count = 0;
				while( my $hit = $result->next_hit ) {		
					$count++;
					next unless ( $v > 0);
					print "sbjct name is ", $hit->name, "\n";
					while( my $hsp = $hit->next_hsp ) {
						print "score is ", $hsp->score, "\n";
					} 
				}
				is($count, 3);
			}
		}
	}
}
