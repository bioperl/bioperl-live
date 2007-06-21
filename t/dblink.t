# This is -*-Perl-*- code
# $Id$
# Original code donated by Erikjan

use strict;

BEGIN {
    use lib 't/lib';
	use BioperlTest;
	
	test_begin(-tests => 134);
	
    use_ok('Bio::Annotation::DBLink');
    use_ok('Bio::SeqIO::genbank');
    use_ok('Bio::SeqIO::swiss');
    use_ok('Bio::Root::IO');
}

my @genbank_files = qw/BK000016-tpa.gbk ay116458.gb ay149291.gb NC_006346.gb ay007676.gb dq519393.gb P35527.gb/;

# bug 2152
#------------------------------------
my $verbose = 1;
for my $in ( @genbank_files ) {
    my $infile = Bio::Root::IO->catfile("t","data",$in);
    my $seqio = Bio::SeqIO->new(
                            -format  =>'genbank',
                            -verbose => $verbose,
                            -file    => $infile,
    );
    my $seq = $seqio->next_seq;
    my @values = $seq->annotation->get_Annotations('dblink');
    foreach my $value (@values) {
        ok(defined $value, '"'.$value . '"');          # check value is not empty
        ok(index($value,'::') < 0       , 'no double colon'  );  # these tests seems silly
        ok( substr($value,-1) ne ':'    , 'no trailing colon');  #    but all have been known to occur
        ok(index($value,'  ') < 0       , 'no double space'  );  #
        my @parts = split(/:/,$value, 2);
        ok( scalar(@parts) == 2, 'dblink value is splittable');
    }
}
#------------------------------------
