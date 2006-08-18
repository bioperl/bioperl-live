# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use lib '..','.','./blib/lib';
use vars qw($NUMTESTS $DEBUG $error $asn_error);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

BEGIN { 
   # to handle systems with no installed Test module
   # we include the t dir (where a copy of Test.pm is located)
   # as a fallback
   eval { require Test; };
   $error = 0;
   if ( $@ ) {
      use lib 't';
   }
   use Test;
   $NUMTESTS = 121;
   plan tests => $NUMTESTS;

   eval { require IO::String;
			 require LWP::UserAgent;
			 require HTTP::Request::Common; };
   if ( $@ ) {
      print STDERR "IO::String or LWP::UserAgent or HTTP::Request not installed. This means the Bio::DB::* modules are not usable. Skipping tests.\n";
      for ( 1..$NUMTESTS ) {
			skip("IO::String, LWP::UserAgent,or HTTP::Request not installed",1);
      }
      $error = 1;
   }
	eval { require Bio::ASN1::EntrezGene; };
	if ($@) {
      print STDERR "Bio::ASN1::EntrezGene not found, Bio::DB::EntrezGene module is not usable. Skipping tests.\n";
		$asn_error = 1;
	}
}

END { 
   foreach ( $Test::ntest..$NUMTESTS) {
      skip('Unable to run all of the DB tests',1);
   }
}

exit (0) if ($error);

require Bio::DB::GenBank;
require Bio::DB::GenPept;
require Bio::DB::SwissProt;
require Bio::DB::EntrezGene;
require Bio::DB::GDB;
require Bio::DB::MeSH;

my $testnum;
my $verbose = 0;

my %expected_lengths = ('NDP_MOUSE' => 131,
								  'NDP_HUMAN' => 133,
								  'MUSIGHBA1' => 408,
								  'AF303112'  => 1611,
								  'J00522'    => 408,
								  'AF303112'  => 1611,
								  'AF303112.1' => 1611,
								  '2981014'   => 1156,
								  'AF041456'  => 1156,
								  'AY080910'  => 798,
								  'AY080909'  => 1042,
								  'AF155220'  => 1172,
								  '405830'    => 1743,
								  'CELRABGDI' => 1743,
								  '195055'    => 136,
								  'AAD15290'  => 136,
								  'AAC06201'  => 353,
								  'P43780'    => 103,
								  'BOLA_HAEIN'=> 103,
								  'YNB3_YEAST'=> 125,
								  'O39869'    => 56,
								  'P18584'    => 497,
								  'DEGP_CHLTR'=> 497,
								  'AF442768'  => 2547,
								  'P31383'    => 635,
                                  'CH402638'  => 5041,
							  );

if( $DEBUG ) {
	#
	# Bio::DB::GenBank
	#
	my ($gb,$seq,$seqio,$query);
	# get a single seq
	eval {
		ok defined ( $gb = new Bio::DB::GenBank('-verbose'=>$verbose,'-delay'=>0) );
		$seq = $gb->get_Seq_by_id('MUSIGHBA1');
		$seq ? ok 1 : exit 0;
		ok( $seq->length, $expected_lengths{$seq->display_id});
		ok( defined ($seq = $gb->get_Seq_by_acc('AF303112')));
		ok( $seq->length, $expected_lengths{$seq->display_id});
		ok( defined ($seq = $gb->get_Seq_by_version('AF303112.1')));
		ok( $seq->length, $expected_lengths{$seq->display_id});
		ok( defined ($seq = $gb->get_Seq_by_gi('405830')));
		ok( $seq->length, $expected_lengths{$seq->display_id});
	};
	if ($@) {
		if( $DEBUG ) { 
			warn "Warning: Couldn't connect to Genbank with Bio::DB::GenBank.pm!\nError: $@\nDo you have network access? Skipping all other tests";
		}
		foreach ( $Test::ntest..$NUMTESTS ) { skip('no network access',1); }
		exit(0);
	}
	$seq = $seqio = undef;

	eval {
		ok( defined($seqio = $gb->get_Stream_by_id([ qw(J00522 AF303112 
																		2981014)])));
		while( my $s = $seqio->next_seq ) {
			ok( $s->length, $expected_lengths{$s->display_id});
		}
	};

	if ($@) {
		if( $DEBUG ) { warn "Batch access test failed.\nError: $@\n"; }
		foreach ( $Test::ntest..$NUMTESTS ) { skip('no network access',1); }
		exit(0);
	}
	$seq = $seqio = undef;

	# test the temporary file creation and fasta
	eval {
		ok defined ( $gb = new Bio::DB::GenBank('-verbose' =>$verbose,
															 '-format' => 'fasta',
															 '-retrievaltype' => 'tempfile',
															 '-delay' => 0) );
		ok( defined ($seq = $gb->get_Seq_by_id('MUSIGHBA1')));
		# last part of id holds the key
		ok($seq->length, $expected_lengths{(split(/\|/,$seq->display_id))[-1]});
		$seq = $gb->get_Seq_by_acc('AF303112');
		ok( defined $seq);
		# last part of id holds the key
		ok($seq->length, $expected_lengths{(split(/\|/,$seq->display_id))[-1]});
		# batch mode requires genbank format
		$gb->request_format("gb");
		ok(defined($seqio = $gb->get_Stream_by_id([ qw(J00522 AF303112 
																	  2981014)])));
		while( my $s = $seqio->next_seq ) {
			ok( $s->length, $expected_lengths{$s->display_id});
			undef $gb;		# test the case where the db is gone, 
			# but a temp file should remain until seqio goes away. 
		}
	};

	if ($@) {
		if( $DEBUG ) {
			warn "Warning: Couldn't connect to complete GenBank tests with a tempfile with Bio::DB::GenBank.pm!\n $@\n";
		}
		foreach ( $Test::ntest..$NUMTESTS ) { 
			skip('could not connect to Genbank',1); 
		}
	}
	$seq = $seqio = undef;

	# test pipeline creation
	eval {
		ok defined ( $gb = new Bio::DB::GenBank('-verbose' =>$verbose,
															 '-retrievaltype' => 'pipeline',
															 '-delay'  => 0,
						) );
		ok( defined ($seq = $gb->get_Seq_by_id('MUSIGHBA1')));
		ok($seq->length, $expected_lengths{$seq->display_id});
		$seq = $gb->get_Seq_by_acc('AF303112');
		ok( defined $seq);
		ok($seq->length, $expected_lengths{$seq->display_id});
		ok(defined($seqio = $gb->get_Stream_by_id([ qw(J00522 AF303112 
																	  2981014)])));
		while( my $s = $seqio->next_seq ) {
			ok( $s->length, $expected_lengths{$s->display_id});
			undef $gb;		# test the case where the db is gone, 
			# but the pipeline should remain until seqio goes away
		}
	};

	if ($@) {
		if( $DEBUG ) {
			warn "Warning: Couldn't connect to complete GenBank tests with a pipeline with Bio::DB::GenBank.pm!\n $@\n";
		}
		foreach ( $Test::ntest..$NUMTESTS ) { 
			skip('could not connect to Genbank',1); 
		}
	}
	$seq = $seqio = undef;

	# test query facility
	eval {
		ok defined ( $query = Bio::DB::Query::GenBank->new('-verbose' => $verbose,
																			'-db'      => 'nucleotide',
																			'-query'   => 'Onchocerca volvulus[Organism]',
																			'-mindate' => '2002/1/1',
																			'-maxdate' => '2002/12/31'));
		ok $query->count > 0;
		my @ids = $query->ids;
		ok @ids > 0;
		ok @ids == $query->count;
		ok defined ($gb = Bio::DB::GenBank->new('-verbose' =>$verbose,
															 '-delay'  => 0,
															));
		ok defined ($seqio = $gb->get_Stream_by_query($query));
		while( my $s = $seqio->next_seq ) {
			ok( $s->length, $expected_lengths{$s->display_id});
			undef $gb;		# test the case where the db is gone, 
			# but the pipeline should remain until seqio goes away
		}
	};

	if ($@) {
		if( $DEBUG ) {
			warn "Warning: Couldn't connect to complete GenBank query tests!\n $@\n";
		}
		foreach ( $Test::ntest..$NUMTESTS ) { 
			skip('could not connect to Genbank',1); 
		}
	}
	$seq = $seqio = undef;

	# test query facility
	eval {
		ok defined ( $query = Bio::DB::Query::GenBank->new('-verbose' => $verbose,
																			'-db'      => 'nucleotide',
																			'-ids'     => [qw(J00522
																									AF303112
																									2981014)]));
		ok $query->count > 0;
		my @ids = $query->ids;
		ok @ids > 0;
		ok @ids == $query->count;
		ok defined ($gb = Bio::DB::GenBank->new('-verbose' =>$verbose,
															 '-delay'  => 0,
															));
		ok defined ($seqio = $gb->get_Stream_by_query($query));
		while( my $s = $seqio->next_seq ) {
			ok( $s->length, $expected_lengths{$s->display_id});
		}
		$seqio->close(); # the key to preventing errors during make test, no idea why
	};

	if ($@) {
		if( $DEBUG ) {
			warn "Warning: Couldn't connect to complete GenBank query tests!\n $@\n";
		}
		foreach ( $Test::ntest..$NUMTESTS ) { 
			skip('could not connect to Genbank',1); 
		}
	}
	$seq = $seqio = undef;

	# bug 1405
	$gb = Bio::DB::GenBank->new(-format     => 'Fasta',
										 -seq_start  => 2,
										 -seq_stop   => 7 );

	$seq = $gb->get_Seq_by_acc("A11111");
	ok($seq->length,6);
	# complexity tests
	$gb = Bio::DB::GenBank->new(-format     => 'Fasta',
										 -complexity  => 0);
   
	my $seqin = $gb->get_Stream_by_acc("5");
	my @result = (1136, 'dna', 342, 'protein');
	while ($seq = $seqin->next_seq) {
	  ok($seq->length,shift(@result));
	  ok($seq->alphabet,shift(@result));
	}
    ok @result, 0;
	#
	# Real batch retrieval using epost/efetch 
	# these tests may change if integrated further into Bio::DB::Gen*
	# Currently only useful for retrieving GI's via get_seq_stream
	#
   $gb = Bio::DB::GenBank->new();
   
   $seqin = $gb->get_seq_stream(-uids 	=> [4887706 ,431229, 147460],
								-mode	=> 'batch');
   @result = ('M59757', 12611 ,'X76083', 3140, 'J01670', 1593);
   while ($seq = $seqin->next_seq) {
	  ok($seq->accession,shift(@result));
	  ok($seq->length,shift(@result));
   }
   ok @result, 0;
   
   $gb = Bio::DB::GenPept->new();
   
   $seqin = $gb->get_seq_stream(-uids 	=> [2981015, 1621261, 195055],
								-mode	=> 'batch');
   @result = ('AAC06201', 353, 'CAB02640', 193, 'AAD15290', 136);
   while ($seq = $seqin->next_seq) {
	  ok($seq->accession,shift(@result));
	  ok($seq->length,shift(@result));
   }
   ok @result, 0;
   
	#
	# Bio::DB::GenPept
	#
	eval { 
		ok defined($gb = new Bio::DB::GenPept('-verbose'=>$verbose,'-delay'=>0)); 
		ok( defined($seq = $gb->get_Seq_by_id('195055')));
		ok( $seq->length, $expected_lengths{$seq->display_id});
		$seq = $gb->get_Seq_by_acc('AAC06201');
		ok(defined $seq);
		ok( $seq->length, $expected_lengths{$seq->display_id});
		$seqio = $gb->get_Stream_by_id([ qw(AAC06201 195055)]);
		ok( defined $seqio);
		while( my $s = $seqio->next_seq ) {
			ok( $s->length, $expected_lengths{$s->display_id});
		}

		# swissprot genpept parsing   
		ok( defined($seq = $gb->get_Seq_by_acc('2AAA_YEAST') ));
		ok($seq->length, $expected_lengths{$seq->display_id}, 
			$expected_lengths{$seq->display_id});

		# test dbsource stuff
		# small chance this might change but hopefully not
		my @annot = $seq->annotation->get_Annotations('dblink');
		ok(scalar @annot, 31);	# changed 16-Aug-06
		ok($annot[0]->database, 'swissprot');
		ok($annot[0]->primary_id, '2AAA_YEAST');
		ok( ($seq->annotation->get_Annotations('swissprot_dates'))[0]->value, 'Jul 1, 1993');
	};

	if ($@) {
		if( $DEBUG ) { 
			warn "Warning: Couldn't connect to Genbank with Bio::DB::GenPept.pm!\n$@";
		}
		foreach( $Test::ntest..$NUMTESTS ) { 
			skip('could not connect with GenPept',1); 
		}
		exit(0);
	}
	$seq = $seqio = undef;

	#
	# Bio::DB::SwissProt
	#
	eval {
		ok defined($gb = new Bio::DB::SwissProt(-verbose => $verbose,
															 -retrievaltype =>'pipeline',
															 -delay => 0));
		ok(defined($seq = $gb->get_Seq_by_id('YNB3_YEAST')));
		ok( $seq->length, $expected_lengths{$seq->display_id});
		ok($seq->division, 'YEAST');

		ok(defined($seq = $gb->get_Seq_by_acc('P43780')));
		ok( $seq->length, $expected_lengths{$seq->display_id});
		ok( defined( $seq = $gb->get_Seq_by_acc('O39869')));
		ok( $seq->length, $expected_lengths{$seq->accession_number});
		ok($seq->accession_number, 'O39869');
		ok($seq->division, '9PICO');

		# test for bug #958
		$seq = $gb->get_Seq_by_id('P18584');
		ok( defined $seq );
		ok( $seq->length, $expected_lengths{$seq->display_id});

		#skip($seq->primary_id =~ /^Bio::Seq/, $seq->primary_id, 'DEGP');
		ok( $seq->display_id, 'DEGP_CHLTR');
		ok( $seq->division, 'CHLTR');

		ok( defined($gb = new Bio::DB::SwissProt('-verbose'=>$verbose, 
															  '-retrievaltype' => 'tempfile',
															  '-delay' => 0,
															 )));
		ok(defined($seqio = $gb->get_Stream_by_id(['NDP_MOUSE', 'NDP_HUMAN'])));
		undef $gb;		# testing to see if we can remove gb
		ok( defined($seq = $seqio->next_seq()));
		ok( $seq->length, $expected_lengths{$seq->display_id});
		ok( defined($seq = $seqio->next_seq()));
		ok( $seq->length, $expected_lengths{$seq->display_id});
	};

	if ($@) {
		if( $DEBUG ) { 
			print STDERR "Warning: Couldn't connect to SwissProt with Bio::DB::Swiss.pm!\n$@";
		}
		foreach ( $Test::ntest..$NUMTESTS) { 
			skip('could not connect to swissprot',1);
		}
		exit(0);
	}
	$seq = $seqio = undef;
   # test contig retrieval
   eval {
       ok defined ( $gb = new Bio::DB::GenBank('-verbose' => $verbose,
                                               '-delay'  => 0,
                                               '-format' => 'gbwithparts'
                       ) );
       ok( defined ($seq = $gb->get_Seq_by_id('CH402638')));
       ok($seq->length, $expected_lengths{$seq->display_id});
	   # now to check that postprocess_data in NCBIHelper catches CONTIG...
	   ok defined ( $gb = new Bio::DB::GenBank('-verbose' => -1,
                                               '-delay'  => 0,
                                               '-format' => 'gb'
                       ) );
       ok( defined ($seq = $gb->get_Seq_by_id('CH402638')));
       ok($seq->length, $expected_lengths{$seq->display_id});
       };
   if ($@) {
       if( $DEBUG ) { 
           print STDERR "Couldn't connect to GenBank with Bio::DB::GenBank.pm!\n$@";
       }
       foreach ( $Test::ntest..$NUMTESTS) { 
           skip('Could not connect to GenBank',1);
       }
       exit(0);
   }
   	$seq = $seqio = undef;

	#
	# Bio::DB::GDB
	#
	my $gdb = new Bio::DB::GDB;  
	my $info = $gdb->get_info(-type => 'marker', 
									  -id => 'D1S243'); 
	ok $info->{gdbid},'GDB:188393';

	#
	# Bio::DB::EntrezGene
	#
	unless ($asn_error) {
		eval {
			ok( defined($gb = Bio::DB::EntrezGene->new(-verbose => $verbose, 
																	 -retrievaltype => 'tempfile',
																	 -delay => 0 )));
			ok(defined($seqio = $gb->get_Stream_by_id([2,3064])));
			$seq = $seqio->next_seq;
			ok( $seq->display_id, "A2M");
			ok($seq->accession_number, 2);
			$seq = $seqio->next_seq;
			ok( $seq->display_id, "HD");
			ok($seq->accession_number, 3064);
			$seq = $gb->get_Seq_by_id(6099);
			ok( $seq->display_id, "RP");
			ok($seq->accession_number, 6099);
		};
	}

	if ($@) {
		if( $DEBUG ) { 
			print STDERR "Warning: Couldn't connect to Entrez with Bio::DB::EntrezGene!\n$@";
		}
		foreach ( $Test::ntest..$NUMTESTS) { 
			skip('Could not connect to Entrez Gene',1);
		}
		exit(0);
	}
	$seq = $seqio = undef;

	#
	# Bio::DB::MeSH
	#
	eval {
		ok my $mesh = new Bio::DB::MeSH(-verbose => $verbose);
		ok my $t = $mesh->get_exact_term('Dietary Fats');
		ok $t->each_twig(), 2;

		$t = $mesh->get_exact_term("Sinus Thrombosis, Intracranial");
		ok $t->description, "Thrombus formation in an intracranial venous sinus, including the superior sagittal, cavernous, lateral, and petrous sinuses. Etiologies include thrombosis due to infection,  DEHYDRATION, coagulation disorders (see  THROMBOPHILIA), and  CRANIOCEREBRAL TRAUMA.";
		ok $t->id, "D012851";
	};
	if ($@) {
		if( 1 ) { 
			print STDERR "Warning: Couldn't connect to MeSH with Bio::DB::MeSH!\n$@";
		}
		foreach ( $Test::ntest..$NUMTESTS) { 
			skip('Could not connect to MeSH',1);
		}
		exit(0);
	}
	$seq = $seqio = undef;


} else {
	for ( $Test::ntest..$NUMTESTS) {
		skip("Skipping tests which require remote servers - set env variable BIOPERLDEBUG to test",1);
	}
}
