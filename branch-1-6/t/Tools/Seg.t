# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 15);
	
	use_ok('Bio::Tools::Seg');
}

my ($infile, $parser) ;

$infile = test_input_file('seg.out');
ok ($parser = Bio::Tools::Seg->new(-file=>$infile), 'parser defined') ;

my @feat;
while ( my $feat = $parser->next_result ) {
	push @feat, $feat;
}

is scalar(@feat), 3;

# seq 0
#>LBL_0012(32-46) complexity=2.47 (12/2.20/2.50)
#gdggwtfegwggppe

# seq 1
#>LBL_0012(66-80) complexity=2.31 (12/2.20/2.50)
#kfssrasakavakks

# seq 2
#>LBL_0012(123-138) complexity=2.31 (12/2.20/2.50)
#svivsqsqgvvkgvgv

my $raa_testdata = [
	[ 'LBL_0012', 32, 46, 2.47   ],
	[ 'LBL_0012', 66, 80, 2.31   ],
	[ 'LBL_0012', 123, 138, 2.31 ],
] ;

for (0..( scalar(@feat)-1 )) {
	is ( $feat[$_]->seq_id, $raa_testdata->[$_]->[0], "seq id for seq $_ identified" ) ;
	is ( $feat[$_]->start,  $raa_testdata->[$_]->[1], "start for seq $_ identified"  ) ;
	is ( $feat[$_]->end,    $raa_testdata->[$_]->[2], "end for seq $_ identified"    ) ;
	is ( $feat[$_]->score,  $raa_testdata->[$_]->[3], "score for seq $_ identified"  ) ;
}
