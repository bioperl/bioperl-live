##############################################
# tests http retrieval
##############################################

use strict;
use warnings;
use Test::More;
use Test::Exception;

use Bio::Root::IO;

my $TESTURL = 'http://www.google.com/index.html';

my $rio = Bio::Root::IO->new();

ok $rio = Bio::Root::IO->new(-url=>$TESTURL), 'default -url method';
lives_ok {$rio = Bio::Root::IO->new(-url=>$TESTURL)};

done_testing;
