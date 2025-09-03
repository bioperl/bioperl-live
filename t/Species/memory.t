#!/usr/bin/env perl

use strict;
use warnings;

use Bio::Root::Test;
use Bio::Species;

eval { require Test::Memory::Cycle; 1; };
my $CYCLE = $@ ? 0 : 1;
eval { require Test::Weaken; 1; };
my $WEAKEN = $@ ? 0 : 1;

test_begin(-tests => 6);

SKIP: {
    skip("Test::Memory::Cycle not installed", 3) if !$CYCLE;

    # Intentional circular ref (should leak)
    my ($a, $b); $a = \$b; $b = \$a;
    ok(Test::Memory::Cycle::memory_cycle_exists($a), "Intentional circular reference leaks");

    # Bio::Species should not leak
    my $species = Bio::Species->new(
        -classification => [qw(sapiens Homo)],
        -common_name    => 'human'
    );
    ok(!Test::Memory::Cycle::memory_cycle_exists($species), "Bio::Species does not leak");

    # GitHub issue #81 regression
    ok(!Test::Memory::Cycle::memory_cycle_exists(Bio::Species->new(-classification => ['A'])), "Regression test for #81");
}

SKIP: {
    skip("Test::Weaken not installed", 3) if !$WEAKEN;

    # Deliberate leak
    ok(Test::Weaken::leaks({
        constructor => sub { my ($a, $b); $a = \$b; $b = \$a; }
    }), "Deliberate circular ref detected by Test::Weaken");

    # Bio::Species should not leak
    ok(!Test::Weaken::leaks({
        constructor => sub {
            Bio::Species->new(
                -classification => [qw(sapiens Homo)],
                -common_name    => 'human'
            )
        }
    }), "Bio::Species passes Test::Weaken");

    # Regression test for #81
    ok(!Test::Weaken::leaks({
        constructor => sub { Bio::Species->new(-classification => ['A']) }
    }), "Regression check for #81 via Test::Weaken");
}

