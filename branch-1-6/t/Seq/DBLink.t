# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 131);
	
    use_ok('Bio::SeqIO');
}

my @genbank_files = qw/BK000016-tpa.gbk ay116458.gb ay149291.gb NC_006346.gb ay007676.gb dq519393.gb P35527.gb/;

# bug 2152
#------------------------------------
my $verbose = test_debug();
for my $in ( @genbank_files ) {
    my $infile = test_input_file($in);
    my $seqio = Bio::SeqIO->new(
                            -format  =>'genbank',
                            -verbose => $verbose,
                            -file    => $infile,
    );
    my $seq = $seqio->next_seq;
    my @values = $seq->annotation->get_Annotations('dblink');
    foreach my $value (@values) {
        my $output = $value->display_text;
        ok(defined $output, '"'.$output . '"');          # check value is not empty
        ok(index($output,'::') < 0       , 'no double colon'  );  # these tests seems silly
        ok( substr($output,-1) ne ':'    , 'no trailing colon');  #    but all have been known to occur
        ok(index($output,'  ') < 0       , 'no double space'  );  #
        my @parts = split(/:/,$output, 2);
        ok( scalar(@parts) == 2, 'dblink value is splittable');
    }
}
#------------------------------------
