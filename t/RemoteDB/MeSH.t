# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 5,
			   -requires_modules => [qw(IO::String
									    LWP::UserAgent
										HTTP::Request::Common)],
			   -requires_networking => 1);
	
    use_ok('Bio::DB::MeSH');
}

#
# Bio::DB::MeSH
#
ok my $mesh = Bio::DB::MeSH->new();
SKIP: {
    my $t;
    eval {$t = $mesh->get_exact_term('Dietary Fats');};
    skip "Couldn't connect to MeSH with Bio::DB::MeSH. Skipping those tests", 3 if $@;
    is $t->each_twig(), 2;
    eval {$t = $mesh->get_exact_term("Sinus Thrombosis, Intracranial");};
    skip "Couldn't connect to MeSH with Bio::DB::MeSH. Skipping those tests", 2 if $@;
    is $t->description, "Thrombus formation in an intracranial venous sinus, including the superior sagittal, cavernous, lateral, and petrous sinuses. Etiologies include thrombosis due to infection,  DEHYDRATION, coagulation disorders (see  THROMBOPHILIA), and  CRANIOCEREBRAL TRAUMA.";
    is $t->id, "D012851";
}
