# -*-Perl-*- for my emacs
# $Id$

use strict;
use vars qw($NUMTESTS $DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
my $error;
BEGIN { 
	eval { require Test; };
	if( $@ ) {
		use lib 't';
	}
	use Test;
	$error = 0;
	$NUMTESTS = 13; 
	plan tests => $NUMTESTS;
	eval { require IO::String;
			 require LWP;
			 require LWP::UserAgent;
			 1;
		 };
	if( $@ ) {
		print STDERR "IO::String,LWP,LWP::UserAgent not installed. This means cannot query remote webserver. Skipping tests.\n";
		for( $Test::ntest..$NUMTESTS ) {
			skip("IO::String or LWP::UserAgent not installed. This means the Bio::Tools::Run::RemoteBlast is not usable. Skipping tests",1);
		}
		$error = 1; 
	}
}

END {     
	for ( $Test::ntest..$NUMTESTS ) {
		skip("Unable to complete RemoteBlast tests",1);
	}
}

my $actually_submit = $DEBUG > 0;

if( $error ==  1 ) {
    exit(0);
}

require Bio::Tools::Run::RemoteBlast;

require Bio::SeqIO;
require Bio::AlignIO;
require Bio::Seq;
require Bio::Root::IO;
use Env;

ok(1);

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

if( $actually_submit == 0 ) {
	print STDERR "Skipping submitting remote BLAST to avoid Time-out\n" if( $DEBUG );
	foreach( $Test::ntest..$NUMTESTS) { 
		skip('Skip to avoid timeout',1);
	}
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
				ok($result->database_name, qr/swissprot/i);
				my $count = 0;
				while( my $hit = $result->next_hit ) {		
					$count++;
					next unless ( $v > 0);
					print "sbjct name is ", $hit->name, "\n";
					while( my $hsp = $hit->next_hsp ) {
						print "score is ", $hsp->score, "\n";
					} 
				}
				ok($count, 3);
			}
		}
	}
}

# test blasttable

my $remote_blast2 = Bio::Tools::Run::RemoteBlast->new
  ('-verbose'    => $v,
	'-prog'       => $prog,
	'-data'       => $db,
	'-readmethod' => 'blasttable',
	'-expect'     => $e_val,     
  );
$remote_blast2->submit_parameter('ENTREZ_QUERY', 
											'Escherichia coli[ORGN]');

$remote_blast2->retrieve_parameter('ALIGNMENT_VIEW', 'Tabular');

$inputfilename = Bio::Root::IO->catfile("t","data","ecolitst.fa");

if( $actually_submit == 0 ) {
	print STDERR "Skipping submitting remote BLAST to avoid Time-out\n" if( $DEBUG );
	foreach( $Test::ntest..$NUMTESTS) { 
		skip('Skip to avoid timeout',1);
	}
} else {
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
				sleep 5;
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
				ok($count, 3);
			}
		}
	}
}


my $remote_blastxml = Bio::Tools::Run::RemoteBlast->new
  ('-verbose'    => $v,
	'-prog'       => $prog,
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
	foreach( $Test::ntest..$NUMTESTS) { 
		skip('Skip blastxml tests probably because XML::SAX not installed',1);
	}
}
elsif( $actually_submit == 0 ) {
	print STDERR "Skipping submitting remote BLAST to avoid Time-out\n" if( $DEBUG );
	foreach( $Test::ntest..$NUMTESTS) { 
		skip('Skip to avoid timeout',1);
	}
}
else {
	my $r = $remote_blastxml->submit_blast($inputfilename);
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
				ok($result->database_name, qr/swissprot/i);
				my $count = 0;
				while( my $hit = $result->next_hit ) {		
					$count++;
					next unless ( $v > 0);
					print "sbjct name is ", $hit->name, "\n";
					while( my $hsp = $hit->next_hsp ) {
						print "score is ", $hsp->score, "\n";
					} 
				}
				ok($count, 3);
			}
		}
	}
}
