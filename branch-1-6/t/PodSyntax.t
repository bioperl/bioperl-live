# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use Test::More;
	eval 'use Test::Pod 1.00';
	plan (skip_all => 'Test::Pod 1.00 required for testing POD' ) if $@;
}

# check pod is syntactically correct
all_pod_files_ok( all_pod_files(qw(Bio scripts examples maintenance)) )
