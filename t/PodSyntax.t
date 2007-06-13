use strict;

BEGIN {
	use Test::More;
	eval 'use Test::Pod 1.00';
	plan (skip_all => 'Test::Pod 1.00 required for testing POD' ) if $@;
  
	use vars qw($NTESTS $error);
	$error = 0;
}


# check pod is syntactically correct
all_pod_files_ok( all_pod_files('.') )
