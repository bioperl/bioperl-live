# -*-Perl-*-
## Bioperl Test Harness Script for Modules


use strict;
BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    use vars qw($NTESTS);
    $NTESTS = 10;
    plan tests => $NTESTS;
}

use Bio::Tools::Run::TribeMCL;
use Bio::SearchIO;

END {
    for ( $Test::ntest..$NTESTS ) {
        skip("dba program not found. Skipping. (Be sure you have the wise package > 2.2.0)",1);
    }
}

my $blast_out = Bio::Root::IO->catfile("t","data","TribeMCL.bls");

#do from raw blast output
my @params=('blastfile'=>$blast_out,I=>'3.0');
my $fact = Bio::Tools::Run::TribeMCL->new(@params);

ok $fact->isa('Bio::Tools::Run::TribeMCL');
ok ($fact->blastfile,$blast_out);

my $bequiet =1 ;
$fact->quiet($bequiet);

my $fam = $fact->run;
ok ($fam->[0]->[0], 'ENSANGP00000008485');
ok ($fam->[1]->[0], 'COE1_MOUSE');
ok ($fam->[2]->[0], 'ENSANGP00000019582');

#do from searchio

my $sio = Bio::SearchIO->new(-format=>'blast',
                                                -file=>$blast_out);
my @params=('searchio'=>$sio,I=>'3.0');
my $fact = Bio::Tools::Run::TribeMCL->new(@params);
ok $fact->isa('Bio::Tools::Run::TribeMCL');
ok ($fact->searchio->isa("Bio::SearchIO"));
my $bequiet =1 ;
$fact->quiet($bequiet);

my $fam = $fact->run;
ok ($fam->[0]->[0], 'ENSANGP00000008485');
ok ($fam->[1]->[0], 'COE1_MOUSE');
ok ($fam->[2]->[0], 'ENSANGP00000019582');


