# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;

    plan tests => 24 }

use Bio::Biblio::Person;
ok(1);
use Bio::Biblio::Organisation;
ok(1);
use Bio::Biblio::Service;
ok(1);

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my($obj);
ok $obj = new Bio::Biblio::Organisation (-name => 'EBI');
ok $obj->set_name, 'EBI';
ok $obj->get_name('EMBL-EBI'), 'EMBL-EBI';

ok $obj = new Bio::Biblio::Service (-name => 'EBI');
ok $obj->set_name, 'EBI';
ok $obj->get_name('EMBL-EBI'), 'EMBL-EBI';

ok $obj = new Bio::Biblio::Person (-name => 'Lehvaslaiho');
ok $obj->get_lastname, 'Lehvaslaiho';
ok $obj->set_lastname('Lehva'), 'Lehva';
ok $obj->set_lastname(), 'Lehva';
ok $obj->get_lastname('Lehvaslaiho'), 'Lehvaslaiho';

ok $obj->set_firstname('Heikki'), 'Heikki';
ok $obj->get_firstname(), 'Heikki';

ok $obj->set_middlename('O'), 'O';
ok $obj->get_middlename(), 'O';

ok $obj->set_email('heikki@ebi.a.uk'), 'heikki@ebi.a.uk';
ok $obj->get_email(), 'heikki@ebi.a.uk';

ok $obj->set_address('Genome Campus, Cambridge, UK'), 'Genome Campus, Cambridge, UK';
ok $obj->get_address(), 'Genome Campus, Cambridge, UK';

ok $obj->set_affiliation('Dr'), 'Dr';
ok $obj->get_affiliation(), 'Dr';
