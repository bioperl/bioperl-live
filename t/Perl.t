# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 31,
			   -requires_module => 'IO::String');
	
	use_ok('Bio::Perl');
}

# Bio::Perl isn't OO so we don't see Bio::Perl->new() here

my ($seq_object,$filename,$blast_report,@seq_object_array);

# will guess file format from extension
$filename = test_input_file('cysprot1.fa');
ok ($seq_object = read_sequence($filename));
isa_ok $seq_object, 'Bio::SeqI';

# forces genbank format
$filename = test_input_file('AF165282.gb');
ok  ($seq_object = read_sequence($filename,'genbank'));
isa_ok $seq_object, 'Bio::SeqI';

# reads an array of sequences
$filename = test_input_file('amino.fa');
is (@seq_object_array = read_all_sequences($filename,'fasta'), 2);
isa_ok $seq_object_array[0], 'Bio::SeqI';
isa_ok $seq_object_array[1], 'Bio::SeqI';

$filename = test_output_file();
ok write_sequence(">$filename",'genbank',$seq_object);
ok ($seq_object = new_sequence("ATTGGTTTGGGGACCCAATTTGTGTGTTATATGTA","myname","AL12232"));
isa_ok $seq_object, 'Bio::SeqI';

my $trans;

ok ($trans = translate($seq_object));

isa_ok $trans, 'Bio::SeqI';

ok ($trans = translate("ATTGGTTTGGGGACCCAATTTGTGTGTTATATGTA"));

isa_ok $trans, 'Bio::PrimarySeqI';

ok ($trans = translate_as_string($seq_object));

is $trans, 'IGLGTQFVCYM';

$trans = '';

ok ($trans = translate_as_string("ATTGGTTTGGGGACCCAATTTGTGTGTTATATGTA"));

is $trans, 'IGLGTQFVCYM';

# we need to keep tests that depend on net connection at the end
# these now run only with BIOPERLDEBUG set

SKIP: {
	test_skip(-tests => 12, -requires_networking => 1);
	
	# swissprot
	SKIP: {
		eval {
			$seq_object = get_sequence('swissprot',"ROA1_HUMAN");
		};
		if ($@) {
			skip("problem connecting to SwissProt:$@",2);
		} else {
			ok $seq_object;
			isa_ok $seq_object, 'Bio::SeqI';
		}
	}
	
    # embl
	SKIP: {
		eval {
			$seq_object = get_sequence('embl',"BUM");
		};
		if ($@) {
			skip("problem connecting to EMBL:$@",2);
		} else {
			ok $seq_object;
			isa_ok $seq_object, 'Bio::SeqI';
		}
	}

	# genbank	
	SKIP: {
		eval {
			$seq_object = get_sequence('genbank',"AI129902");
		};
		if ($@) {
			skip("problem connecting to GenBank:$@",2);
		} else {
			ok $seq_object;
			isa_ok $seq_object, 'Bio::SeqI';
		}
	}

    # refseq
	SKIP: {
		eval {
			$seq_object = get_sequence('genbank',"NM_006732");
		};
		if( $@ ) {
			skip("problem connecting to RefSeq:$@",2);
		} else {
			ok $seq_object;
			isa_ok $seq_object, 'Bio::SeqI';
		}
	}
	
	# genpept
	SKIP: {
		eval {
			$seq_object = get_sequence('genpept',"AAC06201");
		};
		if ($@) {
			skip("problem connecting to RefSeq:$@",2);
		} else {
			ok $seq_object;
			isa_ok $seq_object, 'Bio::SeqI';
		}
	}
    
    # blast
    SKIP: {
		eval {
			$blast_report = blast_sequence($seq_object, 0);
		};
		if ($@) {
			skip("problem connecting to NCBI BLAST:$@",2);
		} else {
			ok $blast_report;
			isa_ok $blast_report, 'Bio::Search::Result::ResultI';
		}
    }
}
