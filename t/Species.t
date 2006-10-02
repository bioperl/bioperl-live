# -*-Perl-*-
## $Id$

use strict;
use vars qw($NUMTESTS $DEBUG);

BEGIN {
	$NUMTESTS = 20;
	$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
	
	eval {require Test::More;};
	if ($@) {
		use lib 't';
	}
	use Test::More;
    
	plan tests => $NUMTESTS;
}

use_ok('Bio::Species');
use_ok('Bio::DB::Taxonomy');

ok my $sps = Bio::Species->new();
$sps->classification(qw( sapiens Homo Hominidae
			 Catarrhini Primates Eutheria Mammalia Vertebrata
			 Chordata Metazoa Eukaryota));

is $sps->binomial, 'Homo sapiens';

ok $sps->sub_species('sapiensis');
is $sps->binomial, 'Homo sapiens';
is $sps->binomial('FULL'), 'Homo sapiens sapiensis';
is $sps->sub_species, 'sapiensis';

$sps->classification(qw( sapiens Homo Hominidae
			 Catarrhini Primates Eutheria Mammalia Vertebrata
			 Chordata Metazoa Eukaryota));
is $sps->binomial, 'Homo sapiens';


# test cmd line initializtion
ok my $species = new Bio::Species( -classification => 
				[ qw( sapiens Homo Hominidae
				      Catarrhini Primates Eutheria 
				      Mammalia Vertebrata
				      Chordata Metazoa Eukaryota) ] );
is $species->binomial, 'Homo sapiens';
is $species->species, 'sapiens';
is $species->genus, 'Homo';


# A Bio::Species isa Bio::Taxon, so test some things from there briefly
is $species->scientific_name, 'sapiens';
is $species->rank, 'species';

# We can make a species object from just an id an db handle
SKIP: {
    skip "Skipping tests which require network access, set BIOPERLDEBUG=1 to test", 5 unless $DEBUG;
    $species = new Bio::Species(-id => 51351);
    my $taxdb = new Bio::DB::Taxonomy(-source => 'entrez');
    eval {$species->db_handle($taxdb);};
    skip "Unable to connect to entrez database; no network or server busy?", 5 if $@;
    is $species->binomial, 'Brassica rapa subsp.';
    is $species->binomial('FULL'), 'Brassica rapa subsp. pekinensis';
    is $species->genus, 'Brassica';
    is $species->species, 'rapa subsp.';
    is $species->sub_species, 'pekinensis';
}
