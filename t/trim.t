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

use strict;

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    
    eval { 
	    require Test::More ; 
	} ;

    if( $@ ) {
        use lib 't/lib' ;
    }
    
    use Test::More tests => 8 ;

}

my $DEBUG = $ENV{'BIOPERLDEBUG'} || 0 ;

use_ok('Bio::Tools::Alignment::Consed', 'Bio::Tools::Alignment::Consed can be used');
use_ok('Bio::Root::IO', 'Bio::Root::IO can be used');
use_ok('Bio::Tools::Alignment::Trim', 'Bio::Tools::Alignment::Trim can be used');

# scope some variables
my($o_consed,@singlets,@singletons,@pairs,@doublets,@multiplets,$invoker);

# instantiate a new object
$o_consed = Bio::Tools::Alignment::Consed->new(	-acefile	=>Bio::Root::IO->catfile("t","data","consed_project","edit_dir","test_project.fasta.screen.ace.1"));
ok ( $o_consed && defined $o_consed, 'create new CSM::Consed object' ) ;

# set the verbosity to a valid value (0)
my $verbosity = $o_consed->set_verbose(0);

# test 3
isa_ok($o_consed,'Bio::Tools::Alignment::Consed', 'creating new Bio::Tools::Alignment::Consed object') ;

# test 4
$invoker = $o_consed->set_singlets("verbosely");
ok ( $invoker != 1, 'singlets can be successfully set' ) ;

ok ( ! $o_consed->set_singlet_quality(), 'singlet quality can be set' );
ok ( ! $o_consed->set_trim_points_singlets_and_singletons(), 'singleton/singlet qualities can be set' );


