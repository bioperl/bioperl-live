##-*-Perl-*-
# $Id$
# test for Bio::DB::Registry, Flat::BinarySearch, and Flat::BDB
# written by Brian Osborne

use strict;
use vars qw($NUMTESTS $old_search_path $no_DB_File $no_LWP $DEBUG);

use File::Spec;

BEGIN {
	$DEBUG = $ENV{BIOPERLDEBUG} || 0;
	eval { require Test; };
	if( $@ ) {
		use lib 't';
	}
	use Test;
	$no_DB_File = 0;
	$old_search_path = $ENV{OBDA_SEARCH_PATH} if defined $ENV{OBDA_SEARCH_PATH};
	$ENV{OBDA_SEARCH_PATH} = 't/data/registry/flat;t/data/registry/bdb';
	eval { require DB_File;
			 require BerkeleyDB;
		 };
	if ($@) {
		$ENV{OBDA_SEARCH_PATH} = 't/data/registry/flat';
		$no_DB_File = 1;
	}
	eval { require LWP::UserAgent;
			 require HTTP::Request::Common; };
	if( $@ ) {
		$no_LWP =1;
	}
	plan tests => ($NUMTESTS=13);
}

require Bio::DB::Registry;
use Bio::DB::Flat;
use Bio::Root::IO;

my $tmpdir = File::Spec->catfile(qw(t tmp));
mkdir($tmpdir,0777);
ok (-d $tmpdir);

my $flat = Bio::DB::Flat->new(-directory  => $tmpdir,
			      -dbname     => 'testflat',
			      -format     => 'fasta',
			      -index      => 'binarysearch',
                              -write_flag => 1 );
my $entries = $flat->build_index(File::Spec->catfile(qw(t data cysprot.fa)));
ok $entries == 7;

if ($no_DB_File){
    for (1..2) {
	skip("DB_File or BerkeleyDB not found, skipping DB_File tests",1);
    }
} else {
    my $bdb = Bio::DB::Flat->new(-directory  => $tmpdir,
				 -dbname     => 'testbdb',
				 -format     => 'fasta',
				 -index      => 'bdb',
				 -write_flag => 1 );
    ok defined($bdb);
    $entries = $bdb->build_index(File::Spec->catfile(qw(t data cysprot.fa)));
    ok $entries == 7;
}

if( $no_LWP ) {
    for (1..9 ) {
	skip("No LWP::UserAgent or HTTP::Request::Common modules installed skipping tests",1);
    }
} else {
    my $registry = Bio::DB::Registry->new;
    ok defined($registry);
    my @available_services = $registry->services;
    
    ok grep /testflat/,@available_services;
    my $db = $registry->get_database('testflat');
    ok defined($db);
    my $seq = $db->get_Seq_by_id("ALEU_HORVU");
    ok defined($seq);
    my $sequence = $seq->seq;
    ok $sequence eq "MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLGALGRTRHALRFARFAVRYGKSYESAAEVRRRFRIFSESLEEVRSTNRKGLPYRLGINRFSDMSWEEFQATRLGAAQTCSATLAGNHLMRDAAALPETKDWREDGIVSPVKNQAHCGSCWTFSTTGALEAAYTQATGKNISLSEQQLVDCAGGFNNFGCNGGLPSQAFEYIKYNGGIDTEESYPYKGVNGVCHYKAENAAVQVLDSVNITLNAEDELKNAVGLVRPVSVAFQVIDGFRQYKSGVYTSDHCGTTPDDVNHAVLAVGYGVENGVPYWLIKNSWGADWGDNGYFKMEMGKNMCAIATCASYPVVAA";

    if ($no_DB_File) {
	for (1..4) {
	    skip("DB_File or BerkeleyDB not found, skipping DB_File tests",1);
	}
    } else {
	ok grep /testbdb/,@available_services;
	$db = $registry->get_database('testbdb');
	ok defined($db);
	$seq = $db->get_Seq_by_id("ALEU_HORVU");
	ok defined($seq);
	$sequence = $seq->seq;
	ok $sequence eq "MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLGALGRTRHALRFARFAVRYGKSYESAAEVRRRFRIFSESLEEVRSTNRKGLPYRLGINRFSDMSWEEFQATRLGAAQTCSATLAGNHLMRDAAALPETKDWREDGIVSPVKNQAHCGSCWTFSTTGALEAAYTQATGKNISLSEQQLVDCAGGFNNFGCNGGLPSQAFEYIKYNGGIDTEESYPYKGVNGVCHYKAENAAVQVLDSVNITLNAEDELKNAVGLVRPVSVAFQVIDGFRQYKSGVYTSDHCGTTPDDVNHAVLAVGYGVENGVPYWLIKNSWGADWGDNGYFKMEMGKNMCAIATCASYPVVAA";
    }
}
END {
	cleanup();
}

sub cleanup {
	eval {
		Bio::Root::IO->rmtree($tmpdir) if (-d $tmpdir);
	};
}

__END__
