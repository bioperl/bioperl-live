# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 6);
	
	use_ok('Bio::Matrix::PSM::InstanceSite');
}

my %params=(-seq=>'TATAAT',-id=>"TATAbox1", -accession_number=>'ENSG00000122304', -mid=>'TB1',
            -desc=>'TATA box, experimentally verified in PRM1 gene',-relpos=>-35, -start=>1965);

ok my $instance = Bio::Matrix::PSM::InstanceSite->new(%params);
is $instance->seq, 'TATAAT';
is $instance->subseq(1,3),'TAT';
is $instance->accession_number, 'ENSG00000122304';
is $instance->end(1970), 1970;
