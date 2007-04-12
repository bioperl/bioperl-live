# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

my $error;

BEGIN { 
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test::More; };
	if( $@ ) {
		use lib 't/lib';
	}
	use Test::More;

	$NUMTESTS = 30;
	plan tests => $NUMTESTS;
	eval { require IO::String };
	if( $@ ) {
		$error = "IO::String not installed. This means the Bio::DB::* modules are not usable. Skipping some tests.";
		$DEBUG = 0;
	} elsif (!$DEBUG) {
		$error = "BIOPERLDEBUG must be set to 1 for running tests requiring network access. Skipping some tests.";
	}
	use_ok('Bio::Perl');
	use_ok('File::Spec');
}

END {
	# clean up after oneself
	unlink (  'Perltmp' );
}

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my ($seq_object,$filename,@seq_object_array);

# will guess file format from extension
$filename = File::Spec->catfile(qw(t data cysprot1.fa));
ok ($seq_object = read_sequence($filename));
isa_ok $seq_object, 'Bio::SeqI';

# forces genbank format
$filename = File::Spec->catfile(qw(t data AF165282.gb));
ok  ($seq_object = read_sequence($filename,'genbank'));
isa_ok $seq_object, 'Bio::SeqI';

# reads an array of sequences
$filename = File::Spec->catfile(qw(t data amino.fa));
is (@seq_object_array = read_all_sequences($filename,'fasta'), 2);
isa_ok $seq_object_array[0], 'Bio::SeqI';
isa_ok $seq_object_array[1], 'Bio::SeqI';

$filename = 'Perltmp';
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
	skip($error, 10) unless $DEBUG;

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
}
