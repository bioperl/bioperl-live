# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#

#####
#
# this test script tests working within the chromat_dir,phd_dir,edit_dir structure
# it also tests the ability of Trim.pm to do its thing
#
#####


use ExtUtils::testlib;
use strict;
require 'dumpvar.pl';

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan tests => 7;
}
my $DEBUG = 0;
print("------------- 05test_trim.t -------------\n") ;
print("Checking if the Bio::Tools::Alignment::Consed module could be used...\n") if($DEBUG);
use Bio::Tools::Alignment::Consed;
use Bio::Root::IO;
ok(1);

print("Checking if the Bio::Tools::Alignment::Trim module could be used...\n") if $DEBUG;
use Bio::Tools::Alignment::Trim;

ok(1);

	# scope some variables
my($o_consed,@singlets,@singletons,@pairs,@doublets,@multiplets,$invoker);

	# instantiate a new object

$o_consed = Bio::Tools::Alignment::Consed->new(	-acefile	=>Bio::Root::IO->catfile("t","data","consed_project","edit_dir","test_project.fasta.screen.ace.1"));
print("Checking if a new CSM::Consed object was created...\n") if( $DEBUG);
ok defined $o_consed;

	# set the verbosity to a valid value (0)
my $verbosity = $o_consed->set_verbose(0);

	#
print("Checking if the new object is a reference to a Bio::Tools::Alignment::Consed object...\n") if($DEBUG);
	# test 3
ok ref($o_consed),'Bio::Tools::Alignment::Consed';

print("Checking if singlets can be successfully set...\n") if ($DEBUG);
	# test 4
$invoker = $o_consed->set_singlets("verbosely");
ok $invoker != 1;

print("Checking if singlets quality can be set...\n") if ($DEBUG);
ok !($o_consed->set_singlet_quality());

print("Checking if singlet and singleton qualities can be set...\n") if( $DEBUG);
ok !($o_consed->set_trim_points_singlets_and_singletons());


