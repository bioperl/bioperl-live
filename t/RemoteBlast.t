# -*-Perl-*- for my emacs
# $Id$

use strict;
use vars qw($NUMTESTS $DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
my $error;
BEGIN { 
	eval { require Test::More; };
	if( $@ ) {
		use lib 't/lib';
	}
	use Test::More;
	$NUMTESTS = 21; 
	plan tests => $NUMTESTS;

}

my $actually_submit = $DEBUG > 0;

require_ok('Bio::Tools::Run::RemoteBlast');
require_ok('Bio::SeqIO');
require_ok('Bio::AlignIO');
require_ok('Bio::Seq');
require_ok('Bio::Root::IO');
use_ok('Env');


SKIP: {
	my $prog = 'blastp';
	my $db   = 'swissprot';
	my $e_val= '1e-10';
	my $v = $DEBUG > 1;
	my  $remote_blast = Bio::Tools::Run::RemoteBlast->new('-verbose' => $v,
														'-prog' => $prog,
														'-data' => $db,
														'-expect' => $e_val,
								  );
	$remote_blast->submit_parameter('ENTREZ_QUERY', 
											  'Escherichia coli[ORGN]');
	my $inputfilename = Bio::Root::IO->catfile("t","data","ecolitst.fa");
	ok( -e $inputfilename);	
	eval { 	require IO::String;
			require LWP;
			require LWP::UserAgent;
		 };
	if( $@ ) {
		skip("IO::String or LWP::UserAgent not installed. This means the Bio::Tools::Run::RemoteBlast is not usable. Skipping tests", 14);
	}
	
	ok(1, 'Text BLAST');
	
	if( $actually_submit == 0 ) {
		print STDERR "Skipping submitting remote BLAST to avoid Time-out\n" if( $DEBUG );
		skip('Skip to avoid timeout',1);
	} else {
		my $r = $remote_blast->submit_blast($inputfilename);
		ok($r);
		print STDERR "waiting..." if( $v > 0 );
		while ( my @rids = $remote_blast->each_rid ) {
			foreach my $rid ( @rids ) {
				my $rc = $remote_blast->retrieve_blast($rid);
				if( !ref($rc) ) {
					if( $rc < 0 ) { 		
						$remote_blast->remove_rid($rid);
						ok(0);
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

	$inputfilename = Bio::Root::IO->catfile("t","data","ecolitst.fa");

	my $r = $remote_blast2->submit_blast($inputfilename);
	ok($r);
	print STDERR "waiting..." if( $v > 0 );
	while ( my @rids = $remote_blast2->each_rid ) {
		foreach my $rid ( @rids ) {
			my $rc = $remote_blast2->retrieve_blast($rid);
			if( !ref($rc) ) {
				if( $rc < 0 ) { 		
					$remote_blast2->remove_rid($rid);
					ok(0);
				}
				print STDERR "." if ( $v > 0 );
				sleep 10;
			} else { 
				ok(1);
				$remote_blast2->remove_rid($rid);
				my $result = $rc->next_result;
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
	
	my $remote_blastxml = Bio::Tools::Run::RemoteBlast->new
	  ('-prog'       => $prog,
		'-data'       => $db,
		'-readmethod' => 'xml',
		'-expect'     => $e_val,
	  );
	$remote_blast->submit_parameter('ENTREZ_QUERY', 
											  'Escherichia coli[ORGN]');
	
	$remote_blastxml->retrieve_parameter('FORMAT_TYPE', 'XML');
	$inputfilename = Bio::Root::IO->catfile("t","data","ecolitst.fa");
	
	eval {require Bio::SearchIO::blastxml;};
	if ($@) {
		skip('Skip blastxml tests probably because XML::SAX not installed',5);
	}
	
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
					ok(0);
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

1;