

use strict;

BEGIN { 
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    plan tests => 12;
}

use Bio::LiveSeq::Mutator;
use Bio::LiveSeq::IO::BioPerl;
use Bio::LiveSeq::Gene;
use Bio::Root::IO;


$a = Bio::LiveSeq::Mutator->new();
ok $a;

ok $a->numbering, 'coding';
ok $a->numbering('coding 1');
ok $a->numbering, 'coding 1';

use Bio::LiveSeq::Mutation;
my $mt = new Bio::LiveSeq::Mutation;
ok $mt->seq('g');
$mt->pos(100);
ok ($a->add_Mutation($mt));
my @each = $a->each_Mutation;
ok( (scalar @each), 1 );
my $mt_b = pop @each;
ok($mt_b->seq, 'g');
my $filename=Bio::Root::IO->catfile("t","ar.embl");
my $loader=Bio::LiveSeq::IO::BioPerl->load('-file' => "$filename");
my $gene_name='AR'; # was G6PD

my $gene=$loader->gene2liveseq('-gene_name' => $gene_name, 
			       '-getswissprotinfo' => 0);
ok($gene);
ok $a->gene($gene);

my $results = $a->change_gene();

ok($results);
use Bio::Variation::IO;
if ($results) {    
    my $out = Bio::Variation::IO->new( '-format' => 'flat');
    ok($out->write($results));    
}
