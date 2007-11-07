# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 36);
	
	use_ok('Bio::Seq');
	use_ok('Bio::SeqIO');
	use_ok('Bio::SeqFeature::Slim');
}

# predeclare variables for strict
my ($feat,$str,$feat2,$pair,$comp_obj1,$comp_obj2,@sft); 

my $DEBUG = test_debug();

$feat = Bio::SeqFeature::Slim->new(-start => 40,
				   -end   => 80,
				   -strand => 1,
				   -primary => 'exon',
				   -source  => 'internal',
				   -tag => {
					'silly' => 20,
					'new' => 1
				   });

is $feat->start, 40, 'start of feature location';
is $feat->end, 80, 'end of feature location';
is $feat->primary_tag, 'exon', 'primary tag';
is $feat->source_tag, 'internal', 'source tag';
$str = $feat->gff_string() || ""; # placate -w
is $feat->length, 41, 'length';
is $feat->strand, 1, 'strand';

my @tags = $feat->get_all_tags;
is @tags, 2, 'number of tags';
is $feat->has_tag('silly'), 1, 'has silly tag';
my ($v) = $feat->get_tag_values('silly');
is $v, 20, 'tag value';

my $gf = $feat->create_seqfeature_generic;
is ref($gf) , 'Bio::SeqFeature::Generic', 'Create generic SF';

my $seq = Bio::Seq->new(-seq => 'ACGT'x100,
			-id  => 'generic');
$seq->add_SeqFeature($feat);

is $feat->entire_seq->seq, $seq->seq, 'Entire sequence object';
is $feat->seq->length, $feat->length, 'Feature length';

my $geneid = 'gene001';
my $mrnaid = 'mRNA001';

my $mRNA = Bio::SeqFeature::Slim->new(-start => 20,
				      -end   => 70,
				      -strand=> 1,
				      -score => 0.70,
				      -primary=> 'mRNA',
				      -source => 'Curated',
				      -display_name => 'BTB_0001',
				      -id           => $mrnaid,
				      -parent       => $geneid,
    );

is $mRNA->score, 0.70, 'Score';
is $mRNA->display_name, 'BTB_0001', 'Display name';
is $mRNA->parent_id, $geneid, 'gene parent id';
is $mRNA->primary_id, $mrnaid, 'mRNA id';

my $c = 0;
my @exons = ([20,29,1],[40,70,1]);
for my $subf ( @exons ) {
    my $f = Bio::SeqFeature::Slim->new(-start => $subf->[0],
				       -end   => $subf->[1],
				       -strand=> $subf->[2],
				       -primary=> 'CDS',
				       -source => 'Curated',
				       -parent => $mrnaid,
				       -id     => sprintf('cds%03d',$c++),
	);
    $mRNA->add_SeqFeature($f);
}
is $mRNA->start, 20, 'mRNA start';
is $mRNA->end, 70, 'mRNA end';
is $mRNA->length, 51, '2 exon mRNA length';

my $i = 0;
for my $cds ( $mRNA->get_SeqFeatures ) {
    is $cds->primary_tag, 'CDS', 'primary tag of cds feature';
    is $cds->start, $exons[$i]->[0], 'exon start';
    is $cds->end, $exons[$i]->[1], 'exon end';
    is $cds->strand, $exons[$i]->[2], 'strand';
    is $cds->parent_id, $mrnaid, 'cds parent id';
    $i++;
}
is $i,2, '2 exons seen';
 
$mRNA->add_SeqFeature(Bio::SeqFeature::Slim->new(
			  -start => 80,
			  -end   => 88,
			  -strand=> 1,
			  -primary=> 'CDS',
			  -source => 'Curated',
			  -parent => $mrnaid,
			  -id     => sprintf('cds%03d',$c++),
		      ),'EXPAND');
is $mRNA->end, 88, '3 exon mRNA end';
is $mRNA->start, 20, '3 exon mRNA start';
is $mRNA->length, 69, '3 exon mRNA length';
