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
    is $t->description, "Formation or presence of a blood clot (THROMBUS) in the CRANIAL SINUSES, large endothelium-lined venous channels situated within the SKULL. Intracranial sinuses, also called cranial venous sinuses, include the superior sagittal, cavernous, lateral, petrous sinuses, and many others. Cranial sinus thrombosis can lead to severe HEADACHE; SEIZURE; and other neurological defects.";
    is $t->id, "D012851";
}
