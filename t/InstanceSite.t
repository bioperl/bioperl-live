#-*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

use strict;

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;

    plan tests => 6;
	use_ok('Bio::Matrix::PSM::InstanceSite');
}

my %params=(-seq=>'TATAAT',-id=>"TATAbox1", -accession_number=>'ENSG00000122304', -mid=>'TB1',
            -desc=>'TATA box, experimentally verified in PRM1 gene',-relpos=>-35, -start=>1965);

ok my $instance=new  Bio::Matrix::PSM::InstanceSite(%params);
is $instance->seq, 'TATAAT';
is $instance->subseq(1,3),'TAT';
is $instance->accession_number, 'ENSG00000122304';
is $instance->end(1999), 1999;

