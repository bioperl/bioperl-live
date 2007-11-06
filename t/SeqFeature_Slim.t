# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 7);
	
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
