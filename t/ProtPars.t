# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

use vars qw($DEBUG );
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

use strict;
BEGIN {
    eval { require Test; };
    if( $@ ) { 
	use lib 't';
    }
    use Test;
    use vars qw($NTESTS);
    $NTESTS = 7;
    plan tests => $NTESTS;
}

use Bio::Tools::Run::Phylo::Phylip::ProtPars;
use Bio::Tools::Run::Alignment::Clustalw; 
END {     
    for ( $Test::ntest..$NTESTS ) {
	skip("Protpars not found. Skipping.",1);
    }
}

ok(1);
my $verbose = -1;
my @params = ('threshold'=>10,'jumble'=>'17,10',outgroup=>2,'idlength'=>10);
my $tree_factory = Bio::Tools::Run::Phylo::Phylip::ProtPars->new(@params);

ok $tree_factory->isa('Bio::Tools::Run::Phylo::Phylip::ProtPars');

my $threshold = 5;
$tree_factory->threshold($threshold);

my $new_threshold= $tree_factory->threshold();
ok $new_threshold, 5, " couldn't set factory parameter";

my $outgroup = 3;
$tree_factory->outgroup($outgroup);

my $new_outgroup= $tree_factory->outgroup();
ok $new_outgroup, 3, " couldn't set factory parameter";


my $jumble = "7,5";
$tree_factory->jumble($jumble);

my $new_jumble= $tree_factory->jumble();
ok $new_jumble, "7,5", " couldn't set factory parameter";

my $bequiet = 1;
$tree_factory->quiet($bequiet);  # Suppress protpars messages to terminal 

my $inputfilename = Bio::Root::IO->catfile("t","data","protpars.phy");
my $tree;
my $protpars_present = $tree_factory->exists_protpars();
unless ($protpars_present) {
    if($DEBUG ){ 
	warn("protpars program not found. Skipping tests $Test::ntest to $NTESTS.\n");    
    }
    exit 0;
}


$tree = $tree_factory->create_tree($inputfilename);

# have to sort the since there is polytomy here.
my @nodes = sort { defined $a->id && defined $b->id && $a->id cmp $b->id } $tree->get_nodes();
ok ($nodes[3]->id, 'SINFRUP002', 
    "failed creating tree by protpars");

$inputfilename = Bio::Root::IO->catfile("t","data","cysprot.fa");
@params = ('ktuple' => 2, 'matrix' => 'BLOSUM', 
	   -verbose => $verbose);
my  $align_factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);
my $clustal_present = $align_factory->exists_clustal();

unless ($clustal_present) {
    warn("Clustalw program not found. Skipping tests $Test::ntest to $NTESTS.\n") if( $DEBUG);    
    exit 0;
}

my $aln = $align_factory->align($inputfilename);

$tree = $tree_factory->create_tree($aln);

@nodes = sort { defined $a->id && defined $b->id && $a->id cmp $b->id } $tree->get_nodes();
ok ($nodes[6]->id, 'CYS1_DICDI', 
    "failed creating tree by protpars");
