#!/usr/bin/env perl
use strict; use warnings;
use Test::More;
use File::Spec;
use File::Basename;

eval "use Test::Pod 1.44";
plan skip_all => "Test::Pod 1.44 required for testing POD" if $@;
plan skip_all => "Perl 5.18 or greater required for testing POD" if $] < 5.018;

my $bindir=File::Spec->catfile(dirname(__FILE__), '../bin');
my $blib=File::Spec->catfile(dirname(__FILE__), '../blib');

my @poddirs = qw( ../blib ../bin );
all_pod_files_ok($blib, $bindir);
