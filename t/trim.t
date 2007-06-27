# -*-Perl-*- Test Harness script for Bioperl
# $Id$

#####
#
# this test script tests working within the chromat_dir,phd_dir,edit_dir structure
# it also tests the ability of Trim.pm to do its thing
#
#####

use strict;

BEGIN {
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 6);
	
	use_ok('Bio::Tools::Alignment::Consed');
}

# scope some variables
my($o_consed,@singlets,@singletons,@pairs,@doublets,@multiplets,$invoker);

# instantiate a new object
$o_consed = Bio::Tools::Alignment::Consed->new(	-acefile	=>test_input_file("consed_project","edit_dir","test_project.fasta.screen.ace.1"));
ok ( $o_consed && defined $o_consed, 'create new CSM::Consed object' ) ;

# set the verbosity to a valid value (0)
my $verbosity = $o_consed->set_verbose(0);

# test 3
isa_ok($o_consed,'Bio::Tools::Alignment::Consed');

# test 4
$invoker = $o_consed->set_singlets("verbosely");
ok ( $invoker != 1, 'singlets can be successfully set');

ok ( ! $o_consed->set_singlet_quality(), 'singlet quality can be set' );
ok ( ! $o_consed->set_trim_points_singlets_and_singletons(), 'singleton/singlet qualities can be set' );

# TODO? huh? are these tests really valid, is Trim even being tested here?
