# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..12\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use lib '../';
use Bio::LiveSeq::Mutator;
use Bio::LiveSeq::IO::BioPerl;
use Bio::LiveSeq::Gene;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

$a = Bio::LiveSeq::Mutator->new();
test 2, $a;

test 3, $a->numbering eq 'coding';
test 4, $a->numbering('coding 1');
test 5, $a->numbering eq 'coding 1';

use Bio::LiveSeq::Mutation;
$mt = new Bio::LiveSeq::Mutation;
test 6, $mt->seq('g');
$mt->pos(100);
test 7, ($a->add_Mutation($mt));
@each = $a->each_Mutation;
test 8, (scalar @each == 1 );
$mt_b = pop @each;
#print STDERR  ref $mt_b, "\n";
test 9, ($mt_b->seq == 'g');
#my $filename='/home/heikki/src/bioperl-live/t/g6pd.embl';
my $filename='t/ar.embl';
my $loader=Bio::LiveSeq::IO::BioPerl->load(-file => "$filename");
my $gene_name='G6PD';
my $gene_name='AR';
test 10, my $gene=$loader->gene2liveseq(-gene_name => $gene_name, 
			       -getswissprotinfo => 0);
#print STDERR "Gene: ",$gene->name,"\n";
test 11, $a->gene($gene);

test 12, $results = $a->change_gene();

print "\n";
use Bio::Variation::IO;
if ($results) {    
    my $out = Bio::Variation::IO->new( '-format' => 'flat');
    $out->write($results);
}

#constuctor test
#$b = Bio::LiveSeq::Mutation->new(-seq=>'AC',
#				     -seqori => 'GG'
#				     -pos => 5,
#				     -len => 2
#				  );
#test 8,  $b;




