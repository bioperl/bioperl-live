##############################################
# tests http retrieval
##############################################

use strict;
use warnings;

BEGIN {
    use Bio::Root::Test;
    test_begin(-tests => 3,
	       -requires_networking => 1);
    use_ok 'Bio::Root::IO';
}

my $TESTURL = 'http://www.google.com/index.html';

my $rio;

ok $rio = Bio::Root::IO->new(-url=>$TESTURL), 'default -url method';
lives_ok {$rio = Bio::Root::IO->new(-url=>$TESTURL)};
