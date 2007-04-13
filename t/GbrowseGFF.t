#-*-Perl-*-
## Bioperl Test Harness Script for Modules

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;

BEGIN {     
    eval { require Test::More; };
    if( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
    plan tests => 5;
    use_ok('Bio::SearchIO');
    use_ok('Bio::SearchIO::Writer::GbrowseGFF');
    use_ok('Bio::Root::IO');
}

END {
    unlink(Bio::Root::IO->catfile(qw(t data gbrowsegff.out)) );
}

my $in = Bio::SearchIO->new(-format => 'blast',
			    -file   => Bio::Root::IO->catfile(
				 qw(t data brassica_ATH.WUBLASTN)));
my $out = new Bio::SearchIO(-output_format  => 'GbrowseGFF',
			    -prefix => 'Sequence',
			    -output_cigar   => 1,
			    -output_signif  => 1,
			    -file           => ">".Bio::Root::IO->catfile
			    (qw(t data gbrowsegff.out) ));
ok($out);
while( my $r = $in->next_result ) {
    ok($out->write_result($r));
}
