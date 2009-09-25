# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 24,
			   -requires_module => 'IO::String');
	
	use_ok('Bio::LiveSeq::Mutator');
	use_ok('Bio::LiveSeq::IO::BioPerl');
	use_ok('Bio::Variation::IO');
}

$a = Bio::LiveSeq::Mutator->new();
ok $a;

is $a->numbering, 'coding';
ok $a->numbering('coding 1');
is $a->numbering, 'coding 1';

my $mt = Bio::LiveSeq::Mutation->new();
ok $mt->seq('g');
$mt->pos(100);
ok ($a->add_Mutation($mt));
my @each = $a->each_Mutation;
is( (scalar @each), 1 );
my $mt_b = pop @each;
is($mt_b->seq, 'g');
my $filename=test_input_file('ar.embl');
my $loader=Bio::LiveSeq::IO::BioPerl->load('-file' => "$filename");
my $gene_name='AR'; # was G6PD

my $gene=$loader->gene2liveseq('-gene_name' => $gene_name);
ok($gene);
ok $a->gene($gene);

my $results = $a->change_gene();
ok($results);

# bug 1701 - mutations on intron/exon boundaries where codon is split 

$loader = Bio::LiveSeq::IO::BioPerl->load( -db   => 'EMBL',
                                -file => test_input_file('ssp160.embl.1')
                        );
# move across intron/exon boundaries, check expected mutations
my @positions = (3128..3129,3188..3189);
my @bases = (qw(C C T T));
my @expected = (qw(T683T T684P T684I T684T));
my $ct = 0;

for my $pos (@positions) {
    # reset gene
    my $gene = $loader->gene2liveseq( -gene_name => 'ssp160');
    my $mutation = Bio::LiveSeq::Mutation->new( -seq => $bases[$ct],
                                                -pos => $pos,
                          );
    my $mutate = Bio::LiveSeq::Mutator->new( -gene      => $gene,
                                             -numbering => 'entry',
                           );
    
    $mutate->add_Mutation( $mutation );

    my $results = $mutate->change_gene();
    
    ok(defined($results));
    is($expected[$ct], $results->trivname);
    $ct++;
}

my $s;
my $io = IO::String->new($s);
my $out = Bio::Variation::IO->new('-fh'   => $io,
                  '-format' => 'flat'
                  );
ok($out->write($results));
ok ($s=~/DNA/ && $s=~/RNA/ && $s=~/AA/);
