# -*-Perl-*-
#Some simple tests for meme and transfac parsers

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

    plan tests => 42;
}
use Bio::Matrix::PSM::IO;

ok(1);
#Let's try meme here
my $psmIO =  new Bio::Matrix::PSM::IO(-format=>'meme', 
				      -file=>Bio::Root::IO->catfile(qw(t data meme.dat)));
ok $psmIO;

my @inputfile=grep(/datafile/i,$psmIO->unstructured);
ok @inputfile;

my $release=$psmIO->release;
ok $release;

my @ids=$psmIO->hid;
ok @ids,4;

my %weights=$psmIO->weight;
ok %weights;

my %seq = $psmIO->seq;
ok %seq,'0';#Meme doesn't have seq

ok $psmIO->version,'3.0';

my $psm = $psmIO->next_psm;
ok $psm;

my $matrix=$psm->matrix;
ok $matrix;
my $psm2=$psm;
$psm2->matrix($matrix);
ok $psm,$psm2;

my %psm_header=$psm->header;
ok $psm_header{IC},38.1;
ok $psm_header{sites},4;
ok $psm_header{width},25;
ok $psm_header{e_val},'1.2e-002';


#Quick check if returned object works
my $IUPAC=$psm->IUPAC;
ok $IUPAC,'CAGAAAAATWVAATYCCCACCHCCC';
ok $IUPAC,$psm2->IUPAC;
ok $IUPAC,$matrix->IUPAC;

my $instances=$psm->instances;
ok $instances;

foreach my $instance (@{$instances}) {
  my $id=$instance->primary_id;
  last if (ok $id);
}

ok $psm->header('e_val');
#Meme parser should be OK if tests passed


#Now we are going to try transfac

$psmIO =  new Bio::Matrix::PSM::IO(-format=>'transfac', 
				   -file=> Bio::Root::IO->catfile(qw(t data transfac.dat)));
ok $psmIO;

my $version=$psmIO->version;
ok !$version;

ok $psmIO->release, '6.4--2002-12-02';

@ids     = $psmIO->hid;
ok @ids, '1';

$psm     = $psmIO->next_psm;
ok $psm;

%weights = $psmIO->weight;
ok !$weights{''};

%seq     = $psmIO->seq;
ok scalar keys %seq, 0;

#Quick check if returned object works
$IUPAC   = $psm->IUPAC;
ok $IUPAC,'NNNNNNNNNNNN';

#Now we are going to try mast
$psmIO =  new Bio::Matrix::PSM::IO(-format=>'mast', 
				   -file=>Bio::Root::IO->catfile(qw(t data mast.dat)));
ok $psmIO;

@inputfile = grep(/datafile/i,$psmIO->unstructured);
ok !@inputfile;

ok( $psmIO->release, '2002/04/02 0:11:59');

@ids     = $psmIO->hid;
ok @ids,4;

%weights = $psmIO->weight;
ok !%weights; #Mast doesn't have weights

ok %seq    = $psmIO->seq;

foreach my $id ($psmIO->hid) {
    ok $seq{$id};
}
ok $psm=$psmIO->next_psm;

my %instances=$psmIO->instances;
ok %instances;

ok $psmIO->version, '3.0';


