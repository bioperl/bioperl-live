# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

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

use Bio::Tools::Run::Phylo::Phylip::ProtDist;
use Bio::Tools::Run::Alignment::Clustalw; 
END {     
    for ( $Test::ntest..$NTESTS ) {
	skip("ProtDist not found. Skipping.",1);
    }
}

ok(1);
my $verbose = -1;
my @params = ('idlength'=>30,'model'=>'cat','gencode'=>'U','category'=>'H','probchange'=>'0.4','trans'=>'5','freq'=>'0.25,0.5,0.125,0.125');
my $dist_factory = Bio::Tools::Run::Phylo::Phylip::ProtDist->new(@params);

ok $dist_factory->isa('Bio::Tools::Run::Phylo::Phylip::ProtDist');

my $model = 'CAT';
$dist_factory->model($model);

my $new_model= $dist_factory->model();
ok $new_model , 'CAT', " couldn't set factory parameter";

my $gencode = 'M';
$dist_factory->gencode($gencode);

my $new_gencode= $dist_factory->gencode();
ok $new_gencode, 'M', " couldn't set factory parameter";


my $category= "H";
$dist_factory->category($category);

my $new_category= $dist_factory->category();
ok $new_category, "H", " couldn't set factory parameter";

my $probchange= 0.4;
$dist_factory->probchange($probchange);

my $new_probchange= $dist_factory->probchange();
ok $new_probchange, 0.4, " couldn't set factory parameter";

my $trans= 5;
$dist_factory->trans($trans);

my $new_trans= $dist_factory->trans();
ok $new_trans, 5, " couldn't set factory parameter";

my $freq= "0.25,0.5,0.125,0.125";
$dist_factory->freq($freq);

my $new_freq= $dist_factory->freq();
ok $new_freq, "0.25,0.5,0.125,0.125", " couldn't set factory parameter";

my $bequiet = 1;
$dist_factory->quiet($bequiet);  # Suppress protpars messages to terminal 

my $inputfilename = Bio::Root::IO->catfile("t","data","protpars.phy");
my $matrix;
my $protdist_present = $dist_factory->exists_protdist();
unless ($protdist_present) {
    warn("protdist program not found. Skipping tests $Test::ntest to $NTESTS.\n");    
    exit 0;
}

$matrix = $dist_factory->create_distance_matrix($inputfilename);
ok ($matrix->{'ENSP000003'}{'SINFRUP001'},0.73389,"failed creating distance matrix");

my $inputfilename = Bio::Root::IO->catfile("t","data","cysprot.fa");
my @params = ('ktuple' => 2, 'matrix' => 'BLOSUM', 
              -verbose => $verbose);
my  $align_factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);
my $aln = $align_factory->align($inputfilename);
$matrix = $dist_factory->create_distance_matrix($aln);
my $hello;

ok ($matrix->{'ALEU_HORVU'}{'CATL_HUMAN'},3.07376,"failed creating distance matrix");
