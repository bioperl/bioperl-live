use strict;
use warnings;
use File::Spec::Functions qw/catfile/;
use Bio::Root::Test;
use FindBin;
test_begin(-tests => 1);

use Bio::SearchIO;

my $searchio = Bio::SearchIO->new(
    -format => 'psl',
    -file   => catfile($FindBin::Bin, 'headerless.psl'),
);

lives_ok { my $result = $searchio->next_result };
