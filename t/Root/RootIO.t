##############################################
# tests http retrieval
##############################################

use strict;
use warnings;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    test_begin(-tests => 2,
	       -requires_networking => 1);
    use_ok 'Bio::Root::IO';
}

my $TESTURL = 'http://www.google.com/index.html';

my $rio = Bio::Root::IO->new();

ok $rio = Bio::Root::IO->new(-url=>$TESTURL), 'default -url method';
lives_ok {$rio = Bio::Root::IO->new(-url=>$TESTURL)};
