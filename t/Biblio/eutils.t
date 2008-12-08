# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 4,
			   -requires_modules => [qw(IO::String LWP::Simple)],
			   -requires_networking => 1);
	
	use_ok('Bio::Biblio');
}

## End of black magic.

my $db;

my $verbose = test_debug();
SKIP: {
	ok ($db = Bio::Biblio->new(-access => 'eutils',
					   -verbose=>$verbose));
	eval { 
		ok(defined($db->find('"Day A"[AU] AND ("Database Management Systems"[MH] OR "Databases, Genetic"[MH] OR "Software"[MH] OR "Software Design"[MH])')));
	};
	
	if ($@) {
		skip("Warning: Couldn't connect to Eutils server!\n$@\n",1);
	}
	
	# these aren't exactly the most stringent of tests...
	my $ct = 0;
	while(my $xml = $db->get_next) {
		$ct++
	}
	# bullet-proof this, though it really needs more stringent tests...
	cmp_ok($ct, '>=', 4)
}
