# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

use strict;
use vars qw($NUMTESTS);
BEGIN { 
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;

    $NUMTESTS = 15; 
    plan tests => $NUMTESTS; 
}

END { unlink qw(cysprot.dnd cysprot1a.dnd) }

use Bio::Tools::Run::Alignment::TCoffee;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Root::IO;

END {     
    for ( $Test::ntest..$NUMTESTS ) {
	skip("TCoffee program not found. Skipping.\n",1);
    }
}

ok(1);


my @params = ('ktuple' => 2, 'matrix' => 'BLOSUM');
my  $factory = Bio::Tools::Run::Alignment::TCoffee->new(@params);

ok ($factory =~ /Bio::Tools::Run::Alignment::TCoffee/);


my $ktuple = 3;
$factory->ktuple($ktuple);

my $new_ktuple = $factory->ktuple();
ok $new_ktuple, 3, " couldn't set factory parameter";


my $what_matrix = $factory->matrix();
ok $what_matrix, 'BLOSUM', "couldn't get factory parameter";

my $bequiet = 1;
$factory->quiet($bequiet);  # Suppress tcoffee messages to terminal

my $inputfilename = Bio::Root::IO->catfile("t","data","cysprot.fa");
my $aln;


my $coffee_present = Bio::Tools::Run::Alignment::TCoffee->exists_tcoffee();
unless ($coffee_present) {
    warn "tcoffee program not found. Skipping tests $Test::ntest to $NUMTESTS.\n";
    exit(0);
}
my $version = $factory->version;
ok ($version >= 1.22, 1, "Code tested only on t_coffee versions > 1.22" );
$aln = $factory->align($inputfilename);
ok $aln->no_sequences, 7;

my $str = Bio::SeqIO->new('-file' => 
			Bio::Root::IO->catfile("t","data","cysprot.fa"), 
			  '-format' => 'Fasta');
my @seq_array =();

while ( my $seq = $str->next_seq() ) {
    push (@seq_array, $seq) ;
}

my $seq_array_ref = \@seq_array;

$aln = $factory->align($seq_array_ref);
ok $aln->no_sequences, 7;
	
my $profile1 = Bio::Root::IO->catfile("t","data","cysprot1a.msf");
my $profile2 = Bio::Root::IO->catfile("t","data","cysprot1b.msf");
$aln = $factory->profile_align($profile1,$profile2);
ok $aln->no_sequences, 7;


my $str1 = Bio::AlignIO->new(-file=> Bio::Root::IO->catfile("t","data","cysprot1a.msf"));
my $aln1 = $str1->next_aln();
ok $aln1->no_sequences, 3;

my $str2 = Bio::AlignIO->new(-file=> Bio::Root::IO->catfile("t","data","cysprot1b.msf"));
my $aln2 = $str2->next_aln();
ok $aln2->no_sequences, 4;

$aln = $factory->profile_align($aln1,$aln2);
ok $aln->no_sequences, 7;


$str1 = Bio::AlignIO->new(-file=> Bio::Root::IO->catfile("t","data","cysprot1a.msf"));
$aln1 = $str1->next_aln();
$str2 = Bio::SeqIO->new(-file=> Bio::Root::IO->catfile("t","data","cysprot1b.fa"));
my $seq = $str2->next_seq();

ok $aln1->no_sequences, 3;
ok( int($aln1->average_percentage_identity), 39);
$aln = $factory->profile_align($aln1,$seq);
ok( $aln->no_sequences, 4);
if( $version <= 1.22 ) {
    ok( int($aln->overall_percentage_identity), 18);    
    ok( int($aln->average_percentage_identity), 44);
} else {
    ok( int($aln->overall_percentage_identity), 21);
    ok( int($aln->average_percentage_identity), 39);    
}
