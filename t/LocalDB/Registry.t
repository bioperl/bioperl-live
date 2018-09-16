# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 14);
	
	$ENV{OBDA_SEARCH_PATH} = 't/data/registry/flat;t/data/registry/bdb';
	
	use_ok('Bio::DB::Registry');
	use_ok('Bio::DB::Flat');
}

# we need a temp directory t/tmp since t/tmp is specified in the registry files
my $tmpdir = File::Spec->catfile(qw(t tmp));
mkdir($tmpdir,0777);

SKIP: {
	skip "unable to create temp dir '$tmpdir', skipping tests", 12 unless -d $tmpdir;
	
	my $flat = Bio::DB::Flat->new(-directory  => $tmpdir,
					  -dbname     => 'testflat',
					  -format     => 'fasta',
					  -index      => 'binarysearch',
								  -write_flag => 1 );
	my $entries = $flat->build_index(test_input_file('cysprot.fa'));
	is $entries, 7;
	
	SKIP: {
		test_skip(-tests => 2, -requires_modules => [qw(DB_File)]);
		
		my $bdb = Bio::DB::Flat->new(-directory  => $tmpdir,
					 -dbname     => 'testbdb',
					 -format     => 'fasta',
					 -index      => 'bdb',
					 -write_flag => 1 );
		ok defined($bdb);
		$entries = $bdb->build_index(test_input_file('cysprot.fa'));
		is $entries, 7;
	}
	
	SKIP: {
		test_skip(-tests => 9,
                  -requires_modules => [qw(LWP::UserAgent HTTP::Request::Common)],
                  -requires_networking => 1);
		
		my $registry = Bio::DB::Registry->new();
		ok defined($registry);
		my @available_services = $registry->services;
		
		ok grep /testflat/,@available_services;
		my $db = $registry->get_database('testflat');
		ok defined($db);
		my $seq = $db->get_Seq_by_id("ALEU_HORVU");
		ok defined($seq);
		my $sequence = $seq->seq;
		is $sequence, "MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLGALGRTRHALRFARFAVRYGKSYESAAEVRRRFRIFSESLEEVRSTNRKGLPYRLGINRFSDMSWEEFQATRLGAAQTCSATLAGNHLMRDAAALPETKDWREDGIVSPVKNQAHCGSCWTFSTTGALEAAYTQATGKNISLSEQQLVDCAGGFNNFGCNGGLPSQAFEYIKYNGGIDTEESYPYKGVNGVCHYKAENAAVQVLDSVNITLNAEDELKNAVGLVRPVSVAFQVIDGFRQYKSGVYTSDHCGTTPDDVNHAVLAVGYGVENGVPYWLIKNSWGADWGDNGYFKMEMGKNMCAIATCASYPVVAA";
	
		SKIP: {
			test_skip(-tests => 4, -requires_modules => [qw(DB_File)]);
			
			ok grep /testbdb/,@available_services;
			$db = $registry->get_database('testbdb');
			ok defined($db);
			$seq = $db->get_Seq_by_id("ALEU_HORVU");
			ok defined($seq);
			$sequence = $seq->seq;
			is $sequence, "MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLGALGRTRHALRFARFAVRYGKSYESAAEVRRRFRIFSESLEEVRSTNRKGLPYRLGINRFSDMSWEEFQATRLGAAQTCSATLAGNHLMRDAAALPETKDWREDGIVSPVKNQAHCGSCWTFSTTGALEAAYTQATGKNISLSEQQLVDCAGGFNNFGCNGGLPSQAFEYIKYNGGIDTEESYPYKGVNGVCHYKAENAAVQVLDSVNITLNAEDELKNAVGLVRPVSVAFQVIDGFRQYKSGVYTSDHCGTTPDDVNHAVLAVGYGVENGVPYWLIKNSWGADWGDNGYFKMEMGKNMCAIATCASYPVVAA";
		}
	}
}

END {
	File::Path::rmtree($tmpdir) if ($tmpdir && (-d $tmpdir));
}
