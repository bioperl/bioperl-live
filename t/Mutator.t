# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

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
    plan tests => 12;
}

use Bio::LiveSeq::Mutator;
use Bio::LiveSeq::IO::BioPerl;
use Bio::LiveSeq::Gene;


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
#print STDERR  ref $mt_b, "\n";
ok($mt_b->seq, 'g');
#my $filename='/home/heikki/src/bioperl-live/t/g6pd.embl';
my $filename='t/ar.embl';
my $loader=Bio::LiveSeq::IO::BioPerl->load('-file' => "$filename");
my $gene_name='AR'; # was G6PD

my $gene=$loader->gene2liveseq('-gene_name' => $gene_name, 
			       '-getswissprotinfo' => 0);
ok($gene);
#print STDERR "Gene: ",$gene->name,"\n";
ok $a->gene($gene);

my $results = $a->change_gene();

ok($results);
use Bio::Variation::IO;
if ($results) {    
    my $out = Bio::Variation::IO->new( '-format' => 'flat');
    ok($out->write($results));    
}
