# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$NUMTESTS $DEBUG $ERROR);
use lib '../';
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    $ERROR = 0;
    if( $@ ) {
	use lib 't';
    }
    use Test;

    $NUMTESTS = 14;
    plan tests => $NUMTESTS;

    eval {
	require Class::AutoClass;
	require Clone; 
    }; 
    if( $@ ) {
        warn("Class::AutoClass or Clone not installed. This means that the module is not usable. Skipping tests");
	$ERROR = 1;
    }
}



exit 0 if $ERROR ==  1;
require Bio::Graph::ProteinGraph;
ok 1;

my $verbose = 0;
$verbose = 1 if $DEBUG;


ok 1;


