use strict;

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;

    plan tests => 8;
}

use Bio::Matrix::PSM::InstanceSite;
ok(1);

  my %params=(-seq=>'TATAAT',-id=>"TATAbox1", -accession_number=>'ENSG00000122304', -mid=>'TB1',
              -desc=>'TATA box, experimentally verified in PRM1 gene',-relpos=>-35, -start=>1965);
my $instance=new  Bio::Matrix::PSM::InstanceSite(%params);
ok $instance;

my $x=$instance->seq;
ok $x, 'TATAAT';

my $x=$instance->subseq(1,3);
ok $x,'TAT';

my $accession=$instance->accession_number;
ok $accession;

ok $instance->end(1999);

