# -*-Perl-*- Test Harness script for Bioperl

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    use Bio::SeqFeature::Generic;
    use Bio::Location::Split;

    test_begin(-tests => 17);
}

my $DEBUG = test_debug();

my $orig = Bio::SeqFeature::Generic->new(
                                       -start => 40,
                                       -end => 80,
                                       -strand => 1,
                                       -primary => 'exon',
                                       -source  => 'internal',
                                       -tag => {
                                           silly => 20,
                                           new => 1
                                           }
                                       );

# ----------
# Verify simple attributes work and are independent of each other
ok(my $clone = $orig->clone(),                'clone()');
ok($clone->start(140),                        'start() clone set');
is($clone->start(),     140,                  'start() clone get');
is($orig->start(),      40,                   'start() original unchanged');

# ----------
# Verify that arguments passed into clone() are applied to the cloned object
# and that the attributes are still independent.
ok($clone = $orig->clone(-start => 150, -end => 157),   'clone() with arguments');
is($orig->start(),       40,                  'start() orig get');
is($orig->end(),         80,                  'end() orig get');
is($clone->start(),     150,                  'start() clone get');
is($clone->end(),       157,                  'end() clone get');
ok($clone->start(140),                        'start() clone set');
is($clone->start(),     140,                  'start() clone get');
is($orig->start(),       40,                  'start() original unchanged');

# ----------
# Verify that object attributes can be cloned, and are independent after cloning
my $splitlocation = Bio::Location::Split->new();
$splitlocation->add_sub_Location(Bio::Location::Simple->new(
   -start=>1, -end=>30, -strand=>1
));
$splitlocation->add_sub_Location(Bio::Location::Simple->new(
   -start=>50, -end=>61, -strand=>1
));
ok($orig->location($splitlocation),                      'location() Bio::Location::Split');
ok($clone = $orig->clone(),                           'clone()');
ok(($clone->location->sub_Location())[1]->start(51),     'start() clone set');
is(($clone->location->sub_Location())[1]->start,     51, 'start() clone get');
is(($orig->location->sub_Location())[1]->start,      50, 'start() original unchanged');


