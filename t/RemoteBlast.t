# -*-Perl-*- for my emacs

use strict;
use vars qw($NUMTESTS);
my $error;
BEGIN { 
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    $error = 0;
    $NUMTESTS = 6; 
    plan tests => 6;
    eval { require 'IO/String.pm' };
    if( $@ ) {
	print STDERR "IO::String not installed. This means the Bio::DB::* modules are not usable. Skipping tests.\n";
	for( 1..$NUMTESTS ) {
	    skip(1,"IO::String not installed. This means the Bio::DB::* modules are not usable. Skipping tests");
	}
	$error = 1; 
    }
}

my $actually_submit = 0;

if( $error ==  1 ) {
    exit(0);
}

require Bio::Tools::Run::RemoteBlast;

require Bio::SeqIO;
require Bio::AlignIO;
require Bio::Seq;
require Bio::Root::IO;

ok(1);

my $prog = 'blastp';
my $db   = 'ecoli';
my $e_val= '1e-10';
my $v = 0;
my  $remote_blast = Bio::Tools::Run::RemoteBlast->new('-verbose' => $v,
						      '-prog' => $prog,
						      '-data' => $db,
						      '-expect' => $e_val,
						      );
my $inputfilename = Bio::Root::IO->catfile("t","data","ecolitst.fa");
ok( -e $inputfilename);

if( $actually_submit == 0 ) {
    print STDERR "Skipping submitting remote BLAST to avoid Time-out\n";
    ok(1);
    ok(1);
    ok(1);
    ok(1);
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
		ok($rc->database, qr/E. coli/);
		my $count = 0;
		while( my $sbjct = $rc->nextSbjct ) {		
		$count++;
		next unless ( $v > 0);
		print "sbjct name is ", $sbjct->name, "\n";
		while( my $hsp = $sbjct->nextHSP ) {
		    print "score is ", $hsp->score, "\n";
		} 
	    }
		ok($count, 3);
	    }
	}
    }
}

