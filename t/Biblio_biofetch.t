# This is -*-Perl-*- code
# $Id$

use strict;

BEGIN { 
	use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 13,
			   -requires_modules => [qw(IO::String
									    LWP::Simple
										HTTP::Request::Common)]);
	
	use_ok('Bio::Biblio');
	use_ok('Bio::Biblio::IO');
}

my ($db,$ref,$refio);
# get a single ref

my $verbose =  test_debug();

$ref = $refio = undef;

SKIP: {
	# check BioFetch access method
	ok ($db = Bio::Biblio->new(-access => 'biofetch',
							   -verbose => $verbose));
	eval { 
		$ref = $db->get_by_id('10592273');
	};
	
	if ($@) {
		skip( "Warning: Couldn't connect to BioFetch server with Bio::DB::Biblio::biofetch! $@", 10); 
	}
	ok(defined($ref)); 
	is $ref->identifier, '10592273';
	$ref = $refio = undef;
	
	ok defined($db = Bio::Biblio->new(-access => 'biofetch',
									  -verbose => $verbose)); 
	
	my $ids = ['10592273', '9613206'];
	eval {
		$refio = $db->get_all($ids);
	};
	
	if ($@) {
		skip("Batch access test failed. Error: $@", 7);
	}
	
	ok(defined($refio));
	is($refio->next_bibref->identifier, '9613206');
	is($refio->next_bibref->identifier, '10592273');
	
	ok defined($db = Bio::Biblio->new(-access => 'biofetch',
									  -verbose => $verbose)); 
	eval {
		$refio = $db->get_Stream_by_id(['10592273', '9613206']);
	};
	
	if ($@) {    
		skip("Batch access test failed. Error: $@", 3);
	}
	
	ok(defined($refio));	
	is($refio->next_bibref->identifier, '9613206');
	is($refio->next_bibref->identifier, '10592273');
}
