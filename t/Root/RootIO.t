##############################################
# tests http retrieval
##############################################

use strict;
use warnings;

BEGIN {
    use Bio::Root::Test;
    test_begin(-tests => 3,
               -requires_modules => ["LWP::Protocol::https"],
               -requires_networking => 1);
    use_ok 'Bio::Root::IO';
}

my $TESTURL = 'https://bioperl.org';

my $rio;

ok $rio = Bio::Root::IO->new(-url=>$TESTURL), 'default -url method';
lives_ok {$rio = Bio::Root::IO->new(-url=>$TESTURL)};
