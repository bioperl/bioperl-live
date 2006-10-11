# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $DEBUG $BIODBTESTS);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

my $error;

BEGIN { 
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test; };
	$error = 0;
	if( $@ ) {
		use lib 't';
	}
	use Test;

	$NUMTESTS = 14;
	$BIODBTESTS = 5;
	plan tests => $NUMTESTS;
	eval { require IO::String };
	if( $@ ) {
		print STDERR "IO::String not installed. This means the Bio::DB::* modules are not usable. Skipping some tests.\n";
		for( 1..$BIODBTESTS ) {
			skip("IO::String not installed. This means the Bio::DB::* modules are not usable. Skipping some tests",1);
		}
		$error = 1;
	}
}

END {
	# clean up after oneself
	unlink (  'Perltmp' );
	for ( $Test::ntest..$NUMTESTS ) {
		skip("Unable to run database access tests",1);
	}
}

use Bio::Perl;
use File::Spec;

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my ($seq_object,$filename,@seq_object_array);



# will guess file format from extension
$filename = File::Spec->catfile(qw(t data cysprot1.fa));
ok ($seq_object = read_sequence($filename)); 
# forces genbank format
$filename = File::Spec->catfile(qw(t data AF165282.gb));
ok  ($seq_object = read_sequence($filename,'genbank')); 
# reads an array of sequences
$filename = File::Spec->catfile(qw(t data amino.fa));
ok (@seq_object_array = read_all_sequences($filename,'fasta'), 2); 
$filename = 'Perltmp';
ok write_sequence(">$filename",'genbank',$seq_object);
ok ($seq_object = new_sequence("ATTGGTTTGGGGACCCAATTTGTGTGTTATATGTA","myname","AL12232"));

my $trans;

ok ($trans = translate($seq_object));

ok ($trans = translate("ATTGGTTTGGGGACCCAATTTGTGTGTTATATGTA"));

ok ($trans = translate_as_string($seq_object));

ok ($trans = translate_as_string("ATTGGTTTGGGGACCCAATTTGTGTGTTATATGTA"));


# we need to keep tests that depend on net connection at the end

unless ( $error ) {
    # swissprot
    eval {
	ok ($seq_object = get_sequence('swissprot',"ROA1_HUMAN"));
    };
    if ($@) {
	if($DEBUG) {
	    warn "Warning: Couldn't connect to SWISS-PROT! Do you have network access?\n";
        }
	exit 0;
    }

    # embl
    eval {
	ok ($seq_object = get_sequence('embl',"BUM"));
    };
    if ($@) {
	if($DEBUG ) {
	    warn "Warning: Couldn't connect to EMBL! Do you have network access?\n";
	}
        exit 0;
    }

    # genbank
    eval {
	ok ($seq_object = get_sequence('genbank',"AI129902"));
    };
    if ($@) {
	if($DEBUG) {
	    warn "Warning: Couldn't connect to GenBank! Do you have network access?\n";
	}
        exit 0;
    }

    # refseq
    eval {
	ok ($seq_object = get_sequence('genbank',"NM_006732"));
    };
    if ($@) {
	if( $DEBUG ) {
	    warn "Warning: Couldn't connect to RefSeq! Do you have network access?\n";
	}
        exit 0;
    }

        # genbank
    eval {
	ok ($seq_object = get_sequence('genpept',"AAC06201"));
    };
    if ($@) {
	if($DEBUG) {
	    warn "Warning: Couldn't connect to GenPept! Do you have network access?\n";
	}
        exit 0;
    }

}
