# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS);

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

    $NUMTESTS = 7;
    plan tests => $NUMTESTS;
    eval { require 'IO/String.pm' };
    if( $@ ) {
	print STDERR "IO::String not installed. This means the Bio::DB::* modules are not usable. Skipping tests.\n";
	for( 1..$NUMTESTS ) {
	    skip(1,"IO::String not installed. This means the Bio::DB::* modules are not usable. Skipping tests");
	}
       $error = 1; 
    }
}

if( $error ==  1 ) {
    exit(0);
}

my $testnum;
my $verbose = 0;

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 


#First of all we need to create an flat db
use Bio::DB::Flat;
use Cwd;
my $cd = cwd();
my $tmpdir = $cd."/t/tmp";
&maketmpdir();
my $db = Bio::DB::Flat->new(-directory  => $tmpdir,
                            -index => 'bdb',
			    -format => 'fasta',
			    -verbose => 1,
                 	    -write_flag => 1
                            );
ok($db);
my $dir = $cd."/t/data/AAC12660.fa";
my $result = $db->build_index(glob($dir));
ok($result);

#Now let's get the sequence out again
my $seq = $db->get_Seq_by_id('AAC12660');
ok($seq);
ok($seq->length,504);

&cleanup();
&maketmpdir();
$db = Bio::DB::Flat->new(-directory  => $tmpdir,
                         -index => 'bdb',
                         -format => 'embl',
                         -verbose => 1,
                         -write_flag => 1
                            );

$dir= $cd."/t/data/factor7.embl";
$result = $db->build_index(glob($dir));
ok($result);
$seq = $db->get_Seq_by_id('HSCFVII');
ok($seq);
ok($seq->length,12850);
&cleanup();
#&maketmpdir();
#my $db = Bio::DB::Flat->new(-directory  => $tmpdir,
#                            -index => 'flat',
#                            -format => 'fasta',
#                            -verbose => 1,
#                            -write_flag => 1
#                            );


sub maketmpdir {
     system("mkdir $tmpdir");
}
sub cleanup {
     system("rm -rf $tmpdir");
} 
