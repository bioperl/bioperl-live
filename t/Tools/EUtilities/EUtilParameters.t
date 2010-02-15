# -*-Perl-*- Test Harness script for Bioperl
# $Id: esearch.t 15112 2008-12-08 18:12:38Z sendu $
#

use strict;
use warnings;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 13,
			   -requires_modules =>
               [qw(URI HTTP::Request)]);
	
    use_ok('Bio::Tools::EUtilities::EUtilParameters');
}

my @ids = qw(6679096 31543332 134288853 483581 20805941 187951953 169158074
123228044 148676374 114326469 148707003 187952787 123233807 148694865 148694864
148694863 148694861 148694862 8705244 8568086);

my %params = (-eutil => 'efetch',
             -db => 'nucleotide',
             -id => \@ids,
             -email => 'me@foo.bar',
             -retmode => 'xml');

my $pobj = Bio::Tools::EUtilities::EUtilParameters->new(%params);

# initial 'primed' state
is($pobj->parameters_changed, 1);

my $request = $pobj->to_request; # 'exhaust' state 
isa_ok($request, 'HTTP::Request');
is($request->url, 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'.
   'db=nucleotide&retmode=xml&id=6679096%2C31543332%2C134288853%2C483581%2C'.
   '20805941%2C187951953%2C169158074%2C123228044%2C148676374%2C114326469%2C'.
   '148707003%2C187952787%2C123233807%2C148694865%2C148694864%2C148694863%2C'.
   '148694861%2C148694862%2C8705244%2C8568086&tool=BioPerl&email=me%40foo.bar');
is($pobj->to_string(), 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/'.
   'efetch.fcgi?db=nucleotide&retmode=xml&id=6679096%2C31543332%2C134288853%2C'.
   '483581%2C20805941%2C187951953%2C169158074%2C123228044%2C148676374%2C'.
   '114326469%2C148707003%2C187952787%2C123233807%2C148694865%2C148694864%2C'.
   '148694863%2C148694861%2C148694862%2C8705244%2C8568086'.
   '&tool=BioPerl&email=me%40foo.bar');
is($pobj->parameters_changed, 0);

# state won't change if the same parameters are passed
$pobj->set_parameters(%params);
is($pobj->parameters_changed, 0);
$pobj->retmode('xml');
is($pobj->parameters_changed, 0);

# reprime state with new value
$pobj->retmode('text');
is($pobj->parameters_changed, 1);

is(join(',',$pobj->available_parameters('epost')),
   'db,retmode,id,tool,email,WebEnv,query_key', 'available_parameters');
is(join(',',$pobj->available_parameters('efetch')),
   'db,retmode,id,retmax,retstart,rettype,strand,seq_start,seq_stop,complexity,report,tool,email,WebEnv,query_key', 'available_parameters');

my %data = $pobj->get_parameters;
is_deeply($data{id}, $params{-id}, 'get_parameters');
is($data{email}, $params{-email}, 'get_parameters');
