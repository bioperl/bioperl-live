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
    $NTESTS = 6;
    plan tests => $NTESTS;
}
use Bio::Tools::Run::Eponine;
use Bio::SeqIO;

END {
    for ( $Test::ntest..$NTESTS ) {
        skip("Eponine program not found. Skipping. (Be sure you have the Eponine-scanpackage )",1);
    }
}
my $inputfilename= Bio::Root::IO->catfile("t","data","eponine.fa");
my $fact = Bio::Tools::Run::Eponine->new("threshold" => 0.999,"EPOJAR" =>"/home/tania/progs/eponine-scan.jar", "seq" =>$inputfilename);
#my $fact = Bio::Tools::Run::Eponine->new("threshold" => 0.998,"java"=>'',"EPOJAR" =>"/home/tania/progs/eponine-scan.jar", "seq" =>$inputfilename);
#my $fact = Bio::Tools::Run::Eponine->new("threshold" => 0.998,"java"=>'/usr/java/jre1.3.1.02/bin/java',"EPOJAR" =>"/home/tania/progs/eponine-scan.jar", "seq" =>$inputfilename);
ok ($fact->threshold, 0.999);
#ok ($fact->java,1);

my @feats = $fact->run_eponine();
ok ($feats[0]->start, 69);
ok ($feats[0]->end, 69);
ok ($feats[0]->strand, 1);
ok ($feats[1]->start,178 );
ok ($feats[1]->end, 180);





