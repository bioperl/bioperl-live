# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(
		#-tests => 10,
		-requires_modules => [qw(REST::Client)]
    );

    use_ok('Bio::DB::NextProt');
}
