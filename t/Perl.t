# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $BIODBTESTS);

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

    $NUMTESTS = 9;
    $BIODBTESTS = 4;
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
}

use Bio::Perl qw( read_sequence 
		  read_all_sequences 
		  write_sequence 
		  new_sequence 
		  get_sequence );

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my ($seq_object,$filename,@seq_object_array);


unless ( $error ) {
    # swissprot
    eval {
	ok ($seq_object = get_sequence('swissprot',"ROA1_HUMAN"));
    };
    if ($@) {
	warn "Warning: Couldn't connect to SWISS-PROT! Do you have network access?\n";
        skip(1,1, 'no network access');
    }

    # embl
    eval {
	ok ($seq_object = get_sequence('embl',"BUM"));
    };
    if ($@) {
	warn "Warning: Couldn't connect to EMBL! Do you have network access?\n";
        skip(1,1, 'no network access');
    }

    # genbank
    eval {
	ok ($seq_object = get_sequence('genbank',"AI129902"));
    };
    if ($@) {
	warn "Warning: Couldn't connect to GenBank! Do you have network access?\n";
        skip(1,1, 'no network access');
    }

    # refseq
    eval {
	ok ($seq_object = get_sequence('genbank',"NM_006732"));
    };
    if ($@) {
	warn "Warning: Couldn't connect to RefSeq! Do you have network access?\n";
        skip(1,1, 'no network access');
    }
}




# will guess file format from extension
$filename = 't/data/cysprot1.fa';
ok ($seq_object = read_sequence($filename)); 
# forces genbank format
$filename = 't/data/AF165282.gb';
ok  ($seq_object = read_sequence($filename,'genbank')); 
# reads an array of sequences
$filename = 't/data/amino.fa';
ok (@seq_object_array = read_all_sequences($filename,'fasta'), 2); 
$filename = 'Perltmp';
ok write_sequence(">$filename",'genbank',$seq_object);
ok ($seq_object = new_sequence("ATTGGTTTGGGGACCCAATTTGTGTGTTATATGTA","myname","AL12232"));

