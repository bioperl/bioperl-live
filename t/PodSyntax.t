# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use Test::More;
	eval 'use Test::Pod 1.00';
	plan (skip_all => 'Test::Pod 1.00 required for testing POD' ) if $@;
}

# check pod is syntactically correct
# Files in Bio::Root come from another git project. Those use an extended POD
# format that Test::Pod doesn't think is valid.
my @pod_files = grep(!m{^Bio/Root},
		     all_pod_files(qw(Bio scripts examples maintenance)));

all_pod_files_ok( @pod_files )
