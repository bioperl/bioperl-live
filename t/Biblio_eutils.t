# This is -*-Perl-*- code
# $Id$

use strict;

BEGIN { 
    use lib 't/lib';
	use BioperlTest;
	
	test_begin(-tests => 6,
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
	
	while(my $xml = $db->get_next) {
		ok(1);
	}
}
