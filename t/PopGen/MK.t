# -*-Perl-*- Test Harness script for Bioperl
# $Id$

# This will outline many tests for the population genetics
# objects in the Bio::PopGen namespace

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 46);
    use_ok('Bio::AlignIO');
    use_ok('Bio::PopGen::Statistics');
    use_ok('Bio::PopGen::Utilities');
    
}
my $verbose = test_debug();


# - McDonald Kreitman tests -

my $stats = Bio::PopGen::Statistics->new(-verbose => $verbose);
isa_ok($stats, 'Bio::PopGen::Statistics');
my $alnio = Bio::AlignIO->new(-format => 'fasta',
			    -file   => test_input_file('CG2865.fasaln'));
my $aln = $alnio->next_aln;
isa_ok($aln,'Bio::SimpleAlign');
my $population = Bio::PopGen::Utilities->aln_to_population(-alignment => $aln,
							   -site_model => 'codon',
							   -include_monomorphic => 1);
isa_ok($population,'Bio::PopGen::Population');

my @marker_names = $population->get_marker_names;
my @inds = $population->get_Individuals;

is(scalar @marker_names, 434, 'Marker Names');
is(scalar @inds, 6,'Number of Inds');

my (@ingroup_seqs,@outgroup_seqs1, @outgroup_seqs2);

for my $ind ( $population->get_Individuals ) {
    my $id = $ind->unique_id;
    # do we allow ingroup to be a list as well?
    my $pushed = 0;
    if( $id =~ /sim/ ) {
	push @ingroup_seqs, $ind;
	$pushed++;
    }
    if( ! $pushed ) {
	if( $id =~ /mel/ ) {
	    push @outgroup_seqs1, $ind;
	    $pushed++;
	} elsif( $id =~ /yak/ ) {
	    push @outgroup_seqs2, $ind;
	    $pushed++;
	}
    }
    #    if( ! $pushed ) {
    #	warn("sequence $id was not grouped, ignoring...\n");
    #	    push @ingroup_seqs, $ind;
    #}   
}
is(scalar @ingroup_seqs, 4, 'number of ingroup sequences');
is(scalar @outgroup_seqs1, 1, 'number of outgroup1 sequences');
is(scalar @outgroup_seqs2, 1, 'number of outgroup2 sequences');

my $polarized = 0;

my @counts = $stats->mcdonald_kreitman(-ingroup   => \@ingroup_seqs,
				       -outgroup  => \@outgroup_seqs1,
				       -polarized => $polarized);
is($counts[0], 0, 'NSpoly');
is($counts[1], 1, 'NSfixed');
is($counts[2], 3, 'Spoly');
is($counts[3], 7, 'Sfixed');
my $mk;
 SKIP: {
     test_skip(-tests => 1, 
	       -requires_module => 'Text::NSP::Measures::2D::Fisher2::twotailed');
     skip "Some problem with Bio::PopGen::Statistics::has_twotailed", 1 unless $Bio::PopGen::Statistics::has_twotailed;
     
     $mk = $stats->mcdonald_kreitman_counts(@counts);
     is($mk, 1, 'McDonald Kreitman');
 }
@counts = $stats->mcdonald_kreitman(-ingroup   => \@ingroup_seqs,
				    -outgroup  => \@outgroup_seqs2,
				    -polarized => $polarized);
is($counts[0], 0, 'NSpoly');
is($counts[1], 6, 'NSfixed');
is($counts[2], 3, 'Spoly');
is($counts[3], 16, 'Sfixed');

SKIP: {
	test_skip(-tests => 1, 
			  -requires_module => 'Text::NSP::Measures::2D::Fisher2::twotailed');
	skip "Some problem with Bio::PopGen::Statistics::has_twotailed", 1 unless $Bio::PopGen::Statistics::has_twotailed;
	
    $mk = $stats->mcdonald_kreitman_counts(@counts);
    is(sprintf("%.2f",$mk), 0.55, 'McDonald Kreitman');
}

@counts = $stats->mcdonald_kreitman(-ingroup  => \@ingroup_seqs,
				    -outgroup => [@outgroup_seqs1,
						  @outgroup_seqs2],
				    -polarized=> 1);
is($counts[0], 0, 'NSpoly');
is($counts[1], 1, 'NSfixed');
is($counts[2], 3, 'Spoly');
is($counts[3], 1, 'Sfixed');


# test 2nd aln file
$alnio = Bio::AlignIO->new(-format => 'fasta',
			   -file   => test_input_file('CG11099.fasaln'));
$aln = $alnio->next_aln;
isa_ok($aln,'Bio::SimpleAlign');
$population = Bio::PopGen::Utilities->aln_to_population(-alignment => $aln,
							-site_model => 'codon',
							-include_monomorphic => 1);
isa_ok($population,'Bio::PopGen::Population');

@marker_names = $population->get_marker_names;
@inds = $population->get_Individuals;

is(scalar @marker_names, 378, 'Marker Names');
is(scalar @inds, 7,'Number of Inds');

@ingroup_seqs = ();
@outgroup_seqs1 = ();
@outgroup_seqs2 = ();

for my $ind ( $population->get_Individuals ) {
    my $id = $ind->unique_id;
    # do we allow ingroup to be a list as well?
    my $pushed = 0;
    if( $id =~ /sim/ ) {
	push @ingroup_seqs, $ind;
	$pushed++;
    }
    if( ! $pushed ) {
	if( $id =~ /mel/ ) {
	    push @outgroup_seqs1, $ind;
	    $pushed++;
	} elsif( $id =~ /yak/ ) {
	    push @outgroup_seqs2, $ind;
	    $pushed++;
	}
    }
    #    if( ! $pushed ) {
    #	warn("sequence $id was not grouped, ignoring...\n");
    #	    push @ingroup_seqs, $ind;
    #}   
}
is(scalar @ingroup_seqs, 5, 'number of ingroup sequences');
is(scalar @outgroup_seqs1, 1, 'number of outgroup1 sequences');
is(scalar @outgroup_seqs2, 1, 'number of outgroup2 sequences');

$polarized = 0;

@counts = $stats->mcdonald_kreitman(-ingroup   => \@ingroup_seqs,
				    -outgroup  => \@outgroup_seqs1,
				    -polarized => $polarized);
is($counts[0], 9, 'NSpoly');
is($counts[1], 1, 'NSfixed');
is($counts[2], 26, 'Spoly');
is($counts[3], 17, 'Sfixed');


SKIP: {
	test_skip(-tests => 1, 
			  -requires_module => 'Text::NSP::Measures::2D::Fisher2::twotailed');
	skip "Some problem with Bio::PopGen::Statistics::has_twotailed", 1 unless $Bio::PopGen::Statistics::has_twotailed;
	
    $mk = $stats->mcdonald_kreitman_counts(@counts);
    is(sprintf("%.2f",$mk), 0.14, 'McDonald Kreitman');
}

@counts = $stats->mcdonald_kreitman(-ingroup  => \@ingroup_seqs,
				    -outgroup => \@outgroup_seqs2,
				    -polarized=> $polarized);
is($counts[0], 9, 'NSpoly');
is($counts[1], 10, 'NSfixed');
is($counts[2], 26, 'Spoly');
is($counts[3], 42, 'Sfixed');

SKIP: {
	test_skip(-tests => 1, 
			  -requires_module => 'Text::NSP::Measures::2D::Fisher2::twotailed');
	skip "Some problem with Bio::PopGen::Statistics::has_twotailed", 1 unless $Bio::PopGen::Statistics::has_twotailed;
	
    $mk = $stats->mcdonald_kreitman_counts(@counts);
    is(sprintf("%.2f",$mk), '0.60', 'McDonald Kreitman');
}

@counts = $stats->mcdonald_kreitman(-ingroup  => \@ingroup_seqs,
				    -outgroup => [@outgroup_seqs1,
						  @outgroup_seqs2],
				    -polarized=> 1);
is($counts[0], 6, 'NSpoly');
is($counts[1], 0, 'NSfixed');
is($counts[2], 17, 'Spoly');
is($counts[3], 1, 'Sfixed');
