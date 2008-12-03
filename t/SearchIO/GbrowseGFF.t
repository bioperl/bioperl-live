# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 3);
	
    use_ok('Bio::SearchIO');
}

my $in = Bio::SearchIO->new(-format => 'blast',
			    -file   => test_input_file('brassica_ATH.WUBLASTN'));
my $out = Bio::SearchIO->new(-output_format  => 'GbrowseGFF',
			    -prefix => 'Sequence',
			    -output_cigar   => 1,
			    -output_signif  => 1,
			    -file           => ">".test_output_file());
ok($out);
while( my $r = $in->next_result ) {
    ok($out->write_result($r));
}
