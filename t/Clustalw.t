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
    
    plan tests => 9;
}

use Bio::Tools::Run::Alignment::Clustalw; 
use Bio::SimpleAlign; 
use Bio::AlignIO; 
use Bio::SeqIO; 
use Bio::Root::IO;

ok(1);

my @params = ('ktuple' => 2, 'matrix' => 'BLOSUM');
my  $factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);

ok $factory->isa('Bio::Tools::Run::Alignment::Clustalw');

my $ktuple = 3;
$factory->ktuple($ktuple);

my $new_ktuple = $factory->ktuple();
ok $new_ktuple, 3, " couldn't set factory parameter";


my $what_matrix = $factory->matrix();
ok $what_matrix, 'BLOSUM', "couldn't get factory parameter";

my $bequiet = 1;
$factory->quiet($bequiet);  # Suppress clustal messages to terminal

my $inputfilename = Bio::Root::IO->catfile("t","data","cysprot.fa");
my $aln;

my $clustal_present = Bio::Tools::Run::Alignment::Clustalw->exists_clustal();
unless ($clustal_present) {
	warn "Clustalw program not found. Skipping tests 5 to 9.\n";
    	skip(1,1);
	skip(1,1);
	skip(1,1);
	skip(1,1);
	skip(1,1);
	exit 0;
}
$aln = $factory->align($inputfilename);

ok ($aln->{'_order'}->{'0'}, 'CATH_HUMAN/1-335', 
    "failed clustalw alignment using input file");

my $str = Bio::SeqIO->new(-file=> Bio::Root::IO->catfile("t","data","cysprot.fa"), 
			  '-format' => 'Fasta');
my @seq_array =();

while ( my $seq = $str->next_seq() ) {
	push (@seq_array, $seq) ;
    }

$aln = $factory->align(\@seq_array);
	
ok ($aln->{'_order'}->{'0'}, 'CATH_HUMAN/1-335', 
    "failed clustalw alignment using BioSeq array ");
	
my $profile1 = Bio::Root::IO->catfile("t","data","cysprot1a.msf");
my $profile2 = Bio::Root::IO->catfile("t","data","cysprot1b.msf");
$aln = $factory->profile_align($profile1,$profile2);

ok( $aln->{'_order'}->{'1'}, 'CATH_HUMAN/1-335', 
    " failed clustalw profile alignment using input file" );

my $str1 = Bio::AlignIO->new(-file=> Bio::Root::IO->catfile("t","data","cysprot1a.msf"));
my $aln1 = $str1->next_aln();
my $str2 = Bio::AlignIO->new(-file=> Bio::Root::IO->catfile("t","data","cysprot1b.msf"));
my $aln2 = $str2->next_aln();

$aln = $factory->profile_align($aln1,$aln2);
ok($aln->{'_order'}->{'1'}, 'CATH_HUMAN/1-335', 
   "failed clustalw profile alignment using SimpleAlign input ");

$str1 = Bio::AlignIO->new(-file=> Bio::Root::IO->catfile("t","data","cysprot1a.msf"));
$aln1 = $str1->next_aln();
$str2 = Bio::SeqIO->new(-file=> Bio::Root::IO->catfile("t","data","cysprot1b.fa"));
my $seq = $str2->next_seq();
$aln = $factory->profile_align($aln1,$seq);

ok ($aln->{'_order'}->{'1'},  'CATH_HUMAN/1-335', 
    "failed adding new sequence to alignment");
