# -*-Perl-*- Test Harness script for Bioperl
# $Id$

# This will outline many tests for the population genetics
# objects in the Bio::PopGen namespace

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 105);

    use_ok('IO::String');
    use_ok('Bio::PopGen::Individual');
    use_ok('Bio::PopGen::Genotype');
    use_ok('Bio::PopGen::Population');
    use_ok('Bio::PopGen::IO');
    use_ok('Bio::PopGen::PopStats');
    use_ok('Bio::AlignIO');
    use_ok('Bio::PopGen::Statistics');
    use_ok('Bio::PopGen::Utilities');
}

# test Fu and Li's D using data from the paper

is(sprintf("%.3f",
	   Bio::PopGen::Statistics->fu_and_li_D_counts(24, 18, 9)),
   -1.529);

is(sprintf("%.3f",
	   Bio::PopGen::Statistics->fu_and_li_D_star_counts(24, 18, 10)),
   -1.558);

is(sprintf("%.3f",
	   Bio::PopGen::Statistics->fu_and_li_F_counts(24, 3.16, 18, 9)),
   -1.735);

is(sprintf("%.2f",
	   Bio::PopGen::Statistics->fu_and_li_F_star_counts(24, 3.16, 18, 10)),
   -1.71);

my $FILE1 = test_output_file();

my @individuals = ( Bio::PopGen::Individual->new(-unique_id => '10a'));
ok($individuals[0]);

my @genotypes = ( Bio::PopGen::Genotype->new(-marker_name    => 'Mkr1',
					    -individual_id  => '10a',
					    -alleles => [ qw(A a)]),
		  Bio::PopGen::Genotype->new(-marker_name    => 'Mkr2',
					    -individual_id  => '10a',
					    -alleles => [ qw(B B)]),
		  Bio::PopGen::Genotype->new(-marker_name    => 'Mkr3',
					    -individual_id  => '10a',
					    -alleles => [ qw(A a)]));
is(($genotypes[1]->get_Alleles)[0], 'B');

$individuals[0]->add_Genotype(@genotypes);
is($individuals[0]->get_Genotypes,3);
is($individuals[0]->get_Genotypes(-marker => 'Mkr3')->get_Alleles(),2);
my @alleles = $individuals[0]->get_Genotypes(-marker => 'Mkr2')->get_Alleles();
is($alleles[0], 'B');

					     
my $population = Bio::PopGen::Population->new(-name        => 'TestPop1',
					     -source      => 'testjasondata',
					     -description => 'throw away example',
					     -individuals => \@individuals);

is(scalar ($population->get_Individuals()), 1);
is($population->name, 'TestPop1');
is($population->source, 'testjasondata');
is($population->description, 'throw away example');

my @genotypes2 = ( Bio::PopGen::Genotype->new(-marker_name   => 'Mkr1',
					     -individual_id => '11',
					     -alleles       => [ qw(A A)]),
		   Bio::PopGen::Genotype->new(-marker_name   => 'Mkr2',
					     -individual_id => '11',
					     -alleles       => [ qw(B B)]),
		   Bio::PopGen::Genotype->new(-marker_name   => 'Mkr3',
					     -individual_id => '11',
					     -alleles       => [ qw(a a)]),
		   Bio::PopGen::Genotype->new(-marker_name   => 'Mkr4',
					     -individual_id => '11',
					     -alleles       => [ qw(C C)])
		   );
push @individuals, Bio::PopGen::Individual->new(-genotypes   => \@genotypes2,
					       -unique_id   => '11');
$population->add_Individual($individuals[1]);

is(scalar ($population->get_Individuals()), 2);
my ($found_ind) = $population->get_Individuals(-unique_id => '10a');
is($found_ind->unique_id, '10a');
is(scalar($population->get_Individuals(-marker => 'Mkr4')) , 1);
is(scalar($population->get_Individuals(-marker => 'Mkr3')) , 2);

my @g = $population->get_Genotypes(-marker => 'Mkr4');

is($g[0]->individual_id, '11');
is(($g[0]->get_Alleles())[0], 'C');

my $marker = $population->get_Marker('Mkr3');
ok($marker);

is($marker->marker_coverage, 2);
@alleles = $marker->get_Alleles;
is(@alleles,2);
my %af = $marker->get_Allele_Frequencies();
is($af{'a'}, 0.75);
is($af{'A'}, 0.25);

$population->remove_Individuals('10a');
$marker = $population->get_Marker('Mkr3');
is($marker->marker_coverage, 1);
%af = $marker->get_Allele_Frequencies();

is($af{'a'}, 1);
is($af{'A'}, undef);

# Read in data from a file
my $io = Bio::PopGen::IO->new(-format => 'csv',
			     -file   => test_input_file('popgen_saureus.dat'));

my @inds;
while( my $ind = $io->next_individual ) {
    push @inds, $ind;
}

my @mrsainds = grep { $_->unique_id =~ /^MRSA/ } @inds;
my @mssainds = grep { $_->unique_id =~ /^MSSA/ } @inds;
my @envinds = grep { $_->unique_id =~ /^NC/ } @inds;

is(scalar @mrsainds, 9);
is(scalar @mssainds, 10);
is(scalar @envinds, 5);

my $mrsapop = Bio::PopGen::Population->new(-name        => 'MRSA',
					  -description => 'Resistant S.aureus',
					  -individuals => \@mrsainds);

my $mssapop = Bio::PopGen::Population->new(-name        => 'MSSA',
					  -description =>'Suceptible S.aureus',
					  -individuals => \@mssainds);

my $envpop = Bio::PopGen::Population->new(-name        => 'NC',
					 -description => 'WT isolates',
					  -individuals => \@envinds);

my $stats = Bio::PopGen::PopStats->new(-haploid => 1);
my $fst = $stats->Fst([$mrsapop,$mssapop],[qw(AFLP1)]);
# We're going to check the values against other programs first
is(sprintf("%.3f",$fst),0.077,'mrsa,mssa aflp1'); 
  
$fst = $stats->Fst([$envpop,$mssapop,$mrsapop],[qw(AFLP1 )]);
is(sprintf("%.3f",$fst),0.035,'all pops, aflp1'); 

$fst = $stats->Fst([$mrsapop,$envpop],[qw(AFLP1 AFLP2)]);
is(sprintf("%.3f",$fst),0.046,'mrsa,envpop aflp1,aflp2');

# Read in data from a file
$io = Bio::PopGen::IO->new(-format => 'csv',
			  -file   => test_input_file('popgen_saureus.multidat'));

@inds = ();
while( my $ind = $io->next_individual ) {
    push @inds, $ind;
}

@mrsainds = grep { $_->unique_id =~ /^MRSA/ } @inds;
@mssainds = grep { $_->unique_id =~ /^MSSA/ } @inds;
@envinds = grep { $_->unique_id =~ /^NC/ } @inds;

is(scalar @mrsainds, 7);
is(scalar @mssainds, 10);
is(scalar @envinds, 5);

$mrsapop = Bio::PopGen::Population->new(-name        => 'MRSA',
				       -description => 'Resistant S.aureus',
				       -individuals => \@mrsainds);

$mssapop = Bio::PopGen::Population->new(-name        => 'MSSA',
				       -description =>'Suceptible S.aureus',
				       -individuals => \@mssainds);

$envpop = Bio::PopGen::Population->new(-name        => 'NC',
				      -description => 'WT isolates',
				      -individuals => \@envinds);

$stats = Bio::PopGen::PopStats->new(-haploid => 1);
my @all_bands = map { 'B' . $_ } 1..20;
my @mkr1     = map { 'B' . $_ } 1..13;
my @mkr2     = map { 'B' . $_ } 14..20;

# still wrong ? seems to work now? --sendubala
$fst = $stats->Fst([$mrsapop,$mssapop],[@all_bands ]);
#TODO: {
#    local $TODO = 'Possible bug with Fst() output';
is(sprintf("%.3f",$fst),'-0.001','mssa,mrsa all_bands'); # We're going to check the values against other programs first
#}
$fst = $stats->Fst([$envpop,$mssapop],[ @mkr1 ]);
is(sprintf("%.3f",$fst),0.023,'env,mssa mkr1'); # We're going to check the values against other programs first

$fst = $stats->Fst([$envpop,$mssapop,$mrsapop],[ @all_bands ]);
is(sprintf("%.3f",$fst),0.071,'env,mssa,mrsa all bands'); # We're going to check the values against other programs first

$fst = $stats->Fst([$envpop,$mssapop,$mrsapop],[ @mkr2 ]);
is(sprintf("%.3f",$fst),0.076, 'env,mssa,mrsa mkr2'); # We're going to check the values against other programs first

$fst = $stats->Fst([$mrsapop,$envpop],[@all_bands ]);
is(sprintf("%.3f",$fst),0.241,'mrsa,nc all_bands'); # We're going to check the values against other programs first

# test overall allele freq setting for a population

my $poptst1 = Bio::PopGen::Population->new(-name => 'tst1');
my $poptst2 = Bio::PopGen::Population->new(-name => 'tst2');

$poptst1->set_Allele_Frequency(-frequencies => 
			       { 'marker1' => { 'a' => '0.20',
						'A' => '0.80' },
				 'marker2' => { 'A' => '0.10',
						'B' => '0.20',
						'C' => '0.70' }
			     });

my $mk1 = $poptst1->get_Marker('marker1');
my %f1 = $mk1->get_Allele_Frequencies;
is($f1{'a'}, '0.20');
is($f1{'A'}, '0.80');
my $mk2 = $poptst1->get_Marker('marker2');
my %f2 = $mk2->get_Allele_Frequencies;
is($f2{'C'}, '0.70');

$poptst2->set_Allele_Frequency(-name      => 'marker1',
			       -allele    => 'A',
			       -frequency => '0.60');
$poptst2->set_Allele_Frequency(-name      => 'marker1',
			       -allele    => 'a',
			       -frequency => '0.40');

#TODO: {
#    local $TODO = 'Fst not calculated yet for just allele freqs';
#    ok 0;
#    #$fst = $stats->Fst([$poptst1,$poptst2],[qw(marker1 marker2) ]);
#}

$io = Bio::PopGen::IO->new(-format => 'csv',
			  -file   => ">$FILE1");

$io->write_individual(@inds);
$io->close();
ok( -s $FILE1);
$io = Bio::PopGen::IO->new(-format => 'csv',
			  -file   => ">$FILE1");

$io->write_population(($mssapop,$mrsapop));
$io->close();
ok( -s $FILE1);

$io = Bio::PopGen::IO->new(-format => 'prettybase',
			  -file   => ">$FILE1");

$io->write_individual(@inds);
$io->close();
ok( -s $FILE1);

$io = Bio::PopGen::IO->new(-format => 'prettybase',
			  -file   => ">$FILE1");

$io->write_population(($mssapop,$mrsapop));
$io->close();
ok( -s $FILE1);


# Let's do PopGen::Statistics tests here

$io = Bio::PopGen::IO->new(-format          => 'prettybase',
			  -no_header       => 1,
			  -file            => test_input_file('popstats.prettybase'));
my (@ingroup,@outgroup);
my $sitecount;
while( my $ind = $io->next_individual ) {
    if($ind->unique_id =~ /out/) {
	push @outgroup, $ind;
    } else { 
	push @ingroup, $ind;
	$sitecount = scalar $ind->get_marker_names() unless defined $sitecount;
    }
}
$stats = Bio::PopGen::Statistics->new();

# Real data and values courtesy M.Hahn and DNASP

is($stats->pi(\@ingroup),2);
is(Bio::PopGen::Statistics->pi(\@ingroup,$sitecount),0.4);

is(Bio::PopGen::Statistics->theta(\@ingroup),1.92);
is(Bio::PopGen::Statistics->theta(\@ingroup,$sitecount),0.384);

# Test with a population object
my $ingroup  = Bio::PopGen::Population->new(-individuals => \@ingroup);
my $outgroup = Bio::PopGen::Population->new(-individuals => \@outgroup);

is($stats->pi($ingroup),2);
is(Bio::PopGen::Statistics->pi($ingroup,$sitecount),0.4);

is(Bio::PopGen::Statistics->theta($ingroup),1.92);
is(Bio::PopGen::Statistics->theta($ingroup,$sitecount),0.384);

my $haploidpop = $ingroup->haploid_population;
is(sprintf("%.5f",Bio::PopGen::Statistics->tajima_D($haploidpop)), 0.27345);

# to fix
is(sprintf("%.5f",Bio::PopGen::Statistics->tajima_D(\@ingroup)),0.27345);
is(sprintf("%.5f",Bio::PopGen::Statistics->tajima_D($ingroup)),0.27345);

is(sprintf("%.5f",Bio::PopGen::Statistics->fu_and_li_D_star(\@ingroup)),
   0.27345);
is(sprintf("%.5f",Bio::PopGen::Statistics->fu_and_li_D_star($ingroup)),
   0.27345);

is(sprintf("%.5f",Bio::PopGen::Statistics->fu_and_li_F_star(\@ingroup)),
     0.27834);
is(sprintf("%.5f",Bio::PopGen::Statistics->fu_and_li_F_star($ingroup)),
   0.27834);

is((Bio::PopGen::Statistics->derived_mutations(\@ingroup,\@outgroup))[0], 1);
is((Bio::PopGen::Statistics->derived_mutations($ingroup,\@outgroup))[0], 1);
is((Bio::PopGen::Statistics->derived_mutations(\@ingroup,$outgroup))[0], 1);
is((Bio::PopGen::Statistics->derived_mutations($ingroup,$outgroup))[0], 1);

# expect to have 1 external mutation
@ingroup = $haploidpop->get_Individuals;

is(sprintf("%.5f",Bio::PopGen::Statistics->fu_and_li_D(\@ingroup,1)),0.75653);
is(sprintf("%.5f",Bio::PopGen::Statistics->fu_and_li_D($ingroup,1)),0.75653);

is(sprintf("%.5f",Bio::PopGen::Statistics->fu_and_li_D(\@ingroup,
						       \@outgroup)),0.75653);
is(sprintf("%.5f",Bio::PopGen::Statistics->fu_and_li_D($ingroup,
						       \@outgroup)),0.75653);
is(sprintf("%.5f",Bio::PopGen::Statistics->fu_and_li_D($ingroup,
						       $outgroup)),0.75653);
is(sprintf("%.5f",Bio::PopGen::Statistics->fu_and_li_D(\@ingroup,
						       $outgroup)),0.75653);

is(sprintf("%.5f",Bio::PopGen::Statistics->fu_and_li_F(\@ingroup,1)),
     0.77499);
is(sprintf("%.5f",Bio::PopGen::Statistics->fu_and_li_F($haploidpop,1)),0.77499);
is(sprintf("%.5f",Bio::PopGen::Statistics->fu_and_li_F($ingroup,
						       \@outgroup)),0.77499);
is( sprintf("%.5f",Bio::PopGen::Statistics->fu_and_li_F($ingroup,
							 $outgroup)),0.77499);


# Test composite LD

$io = Bio::PopGen::IO->new(-format => 'prettybase',
			  -file   => test_input_file('compLD_test.prettybase'));

my $pop = $io->next_population;

my %LD = $stats->composite_LD($pop);

is($LD{'01'}->{'02'}->[1], 10);
is($LD{'01'}->{'03'}->[1], 0);
is($LD{'02'}->{'03'}->[1], 0);

# Test composite LD

$io = Bio::PopGen::IO->new(-format => 'prettybase',
			  -file   => test_input_file('compLD_missingtest.prettybase'));

$pop = $io->next_population;

%LD = $stats->composite_LD($pop);

is(sprintf("%.4f",$LD{'ProC9198EA'}->{'ProcR2973EA'}->[0]), -0.0375);
is(sprintf("%.2f",$LD{'ProC9198EA'}->{'ProcR2973EA'}->[1]), 2.56);



# build a population from an alignment

my $alnin = Bio::AlignIO->new(-format => 'clustalw',
			      -file   => test_input_file('T7.aln'));
my $aln = $alnin->next_aln;
$population = Bio::PopGen::Utilities->aln_to_population(-alignment => $aln);
is($population->get_number_individuals,9);
#warn($aln->match_line,"\n");
my $matchline = $aln->match_line;
is( $population->get_marker_names, $matchline =~ tr/ //);
for my $name ( $population->get_marker_names ) {
    my $marker = $population->get_Marker($name); 
#    warn("$name ",join(" ",$marker->get_Alleles()),"\n");
}


# test Rich's phase and hap parsers

$io = Bio::PopGen::IO->new(-format   => 'hapmap',
			  -verbose  => 1,
			  -no_header=> 1,
			  -starting_column => 10,
			  -file     => test_input_file('example.hap'));

# Some IO might support reading in a population at a time

my @population;
while( my $ind = $io->next_individual ) {
    push @population, $ind;
}
is(@population, 90);
is($population[3]->unique_id, 'NA06994');
is($population[3]->get_Genotypes, 34);
$population = Bio::PopGen::Population->new(-individuals => \@population);

is(sprintf("%.3f",$stats->pi($population)),12.266);
# if forced haploid population is called within pi
# need to decide about that...
# is(sprintf("%.3f",$stats->pi($population)),12.266);

is(sprintf("%.3f",$stats->theta($population)),5.548);
#TODO: {
#    local $TODO = 'May be TJd inconsistency, need to recalculate';
    is(sprintf("%.3f",$stats->tajima_D($population)),'2.926');
    is(sprintf("%.3f",$stats->tajima_D($population->haploid_population)),3.468);
#}

# test converting from hapmap to phase
my $string;
my $out = IO::String->new($string);
Bio::PopGen::IO->new(-fh => $out, -format => 'phase')->write_individual($population[0]);
is($string, "1
34
P rs1022827 rs1323262 rs1359058 rs1475202 rs1570473 rs2025307 rs2025308 rs2296049 rs2296050 rs2296054 rs3739586 rs4400444 rs4584192 rs4740849 rs4740850 rs4740851 rs4740866 rs4742215 rs4742219 rs4742220 rs4742222 rs4742223 rs4742225 rs4742227 rs4742236 rs4742292 rs732118 rs732119 rs745876 rs745877 rs881684 rs912174 rs912175 rs962817
SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
#NA06985
G C A G G A A A C T C T A C A A T G C A A C G G A A A T A T A A C C
G G T G G A A G G T C T G G G T T T T G G T T G G A A T A T A C G T
");


$io = Bio::PopGen::IO->new(-format => 'phase',
			   -file   => test_input_file('example.phase'));

# Some IO might support reading in a population at a time

@population = ();
while( my $ind = $io->next_individual ) {
    push @population, $ind;
}
is(@population, 3);

# test writing in phase format
$population = Bio::PopGen::Population->new(-individuals => \@population);
$string = '';
$out = IO::String->new($string);
Bio::PopGen::IO->new(-fh => $out, -format => 'phase')->write_population($population);
is($string, "3
5
P 1313 1500 2023 300 5635
SSSMM
#1
1 0 1 12 3
0 1 0 11 3
#2
1 1 1 12 2
0 0 0 12 3
#3
? 0 0 -1 2
? 1 1 -1 13
");

# test diploid data

# bug 2492
{
    my $in = Bio::PopGen::IO->new(-format=>"csv", -fh=>\*DATA);
    my $pop = $in->next_population;
    is(sprintf("%.3f",$stats->pi($pop)),0.833,'Pi on 3-allele data');
    is(sprintf("%.3f",$stats->theta($pop)),0.545,'Theta on 3-allele data');
}
__DATA__
SAMPLE,Site-1
seq_1,G
seq_2,C
seq_3,T
seq_4,G
