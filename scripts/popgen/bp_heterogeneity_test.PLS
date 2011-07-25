#!/usr/bin/perl
# -*-Perl-*- (for my emacs)

use strict;
use warnings;

=head1 NAME 

bp_heterogeneity_test - a test for distinguishing between selection and population expansion.

=head1 SYNOPSIS

heterogenetity_test -mut_1/--mutsyn synonymous_mut_count -mut_2/--mutnon nonsyn_mut_count -s/--smaplesize sample_size [-i/--iterations iterations] [-o/--observed observed_D] [-v/--verbose] [--silent ] [-m/--method tajimaD or fuD] [--precision]

=head2 DESCRIPTION

This is an implementation of the Heterogenetity test as described in
Hahn MW, Rausher MD, and Cunningham CW. 2002. Genetics 161(1):11-20. 

=head2 OPTIONS

 Options in brackets above are optional

 -s or --samplesize samplesize 
 -mut_1 or --mutsyn synonymous mutation count 
 -mut_2 or --mutnon nonsynonmous mutation count 
 -i or --iterations number of iterations 
 -o or --observed   observed D 
 -m or --method     tajimaD or fuD  for Tajima's D or Fu and Li's D
 -v or --verbose    print out extra verbose messages
 --silent           Be extra quiet
 --precision        Level of precision - specify the number of digits 
                   (default 4)

=head2 AUTHOR Matthew Hahn E<lt>matthew.hahn-at-duke.eduE<gt>

For more information contact:

Matthew Hahn, E<lt>matthew.hahn-at-duke.eduE<gt>
Jason Stajich E<lt>jason-at-bioperl-dot-orgE<gt>

=cut

use Getopt::Long;
use Bio::PopGen::Simulation::Coalescent;
use Bio::PopGen::Statistics;
use Bio::PopGen::Individual;
use Bio::PopGen::Genotype;

my $sample_size = 4;
my $mut_count_1 = 10; # synonymous
my $mut_count_2 = 20; # non-synonymous
my $iterations = 1;
my $verbose = 0;
my $observedD = undef;
my $method = 'fuD';
my $help = 0;
my $precision = '4'; # Let's make the random precision between
                     # 0->1 to 1000th digits

GetOptions( 
	    's|samplesize|samp_size:i' => \$sample_size,
	    'mut_1|mutsyn:i'           => \$mut_count_1,
	    'mut_2|mutnon:i'           => \$mut_count_2, 
	    'i|iterations:i'           => \$iterations,
	    'o|obsered|observedD:f'    => \$observedD, 
	    'v|verbose'                => \$verbose,
	    'm|method:s'               => \$method,
	    'h|help'                   => \$help,
	    'silent'                   => sub { $verbose = -1; },
	    'p|precision:i'            => \$precision,
	    );

if( $help ) {
    system("perldoc",$0);
    exit(0);
}

if( $method ne 'fuD' and $method ne 'tajimaD' ) {
    die("available methods are [fu and li's D] (fuD) and Tajima's D (tajimaD)");
}
my @D_distribution;  
printf("sample size is %d iteration count = %d\n", $sample_size, 
       $iterations);

my $sim = new Bio::PopGen::Simulation::Coalescent
    (-sample_size => $sample_size);

for(my $iter = 0; $iter < $iterations; $iter++ ) {
    my $tree = $sim->next_tree($sample_size);
    my $f1 = 0;
    if( $mut_count_1 > 0 ) {
	$sim->add_Mutations($tree,$mut_count_1,$precision);

	my @leaves = $tree->get_leaf_nodes;
	# the outgroup is just an individual with the ancestral state 
	# (no mutations)
	my $outgroup = new Bio::PopGen::Individual();
	foreach my $m ( $leaves[0]->get_marker_names ) {
	    $outgroup->add_Genotype(Bio::PopGen::Genotype->new
				    (-marker_name=> $m,
				     -alleles    => [ 0 ]));
	}
	if( $method eq 'fuD' ) {
	    $f1 = Bio::PopGen::Statistics->fu_and_li_D(\@leaves,[$outgroup]);
	} elsif( $method eq 'tajimaD' ) {
	    $f1 = Bio::PopGen::Statistics->tajima_D(\@leaves);
	}
	print "(mutation count = $mut_count_1) D=$f1\n" 
	    if( $verbose >= 0);
    }
    
    my $f2 = 0;
    if( $mut_count_2 > 0 ) {
	$sim->add_Mutations($tree,$mut_count_2,$precision);
	my @leaves = $tree->get_leaf_nodes;
        # the outgroup is just an individual with the ancestral state 
	# (no mutations)
	my $outgroup = new Bio::PopGen::Individual();
	foreach my $m ( $leaves[0]->get_marker_names ) {
	    $outgroup->add_Genotype(Bio::PopGen::Genotype->new
				    (-marker_name=> $m,
				     -alleles    => [ 0 ]));
	}
	if( $method eq 'fuD' ) {
	    $f2 = Bio::PopGen::Statistics->fu_and_li_D(\@leaves,[$outgroup]);
	} elsif( $method eq 'tajimaD' ) {
	    $f2 = Bio::PopGen::Statistics->tajima_D(\@leaves);
	}
	print "(mutation count = $mut_count_2) D=$f2\n" if( $verbose >= 0);

    }
    my $deltaD = ( $f1 - $f2 );
    push @D_distribution, $deltaD;
    if( $iter % 10 == 0 && $iter > 0 ) { 
	print STDERR "iter = $iter\n"; 
    }
}

if( defined $observedD && $iterations > 1 ) { 
    my @sortedD = sort { $a <=> $b } @D_distribution;
    my $i;
    for($i = 0; $i < scalar @sortedD; $i++ ) {
	if( $sortedD[$i] > $observedD ) { 
	    last;
	}
    }
    
    printf( "index %d value=%.4f out of %d total (obs=%.4f)\n", 
	    $i, $sortedD[$i], scalar @sortedD, $observedD);
}

