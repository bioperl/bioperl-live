# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

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

    plan tests => 90;
}

use Bio::Symbol::Alphabet;
use Bio::Symbol::Symbol;
use Bio::Symbol::DNAAlphabet;
use Bio::Symbol::ProteinAlphabet;

my $A = new Bio::Symbol::Symbol(-token => 'A' );
my $U = new Bio::Symbol::Symbol(-token => 'U' );
my $G = new Bio::Symbol::Symbol(-token => 'G' );
my $T = new Bio::Symbol::Symbol(-token => 'T' );

my $rna = new Bio::Symbol::Alphabet( -symbols => [ $A, $U, $G, $T ] );
				     
ok($rna);
my @symbols = $rna->symbols;
ok(scalar @symbols, 4);

ok($rna->contains($A));
ok($rna->contains($T));
ok($rna->contains($U));
ok($rna->contains($G));


my $dna = new Bio::Symbol::DNAAlphabet();
ok($dna->isa('Bio::Symbol::AlphabetI'));
my $count = 0;

my @dnasymbols = sort qw( A B C D G H K M N R S T U V W X Y );
foreach my $s ( sort { $a->name cmp $b->name } $dna->symbols ) {
    ok($s->name, $dnasymbols[$count]);    
    ok($s->token, $dnasymbols[$count++]);    
}

my $prot = new Bio::Symbol::ProteinAlphabet();
ok($prot->isa('Bio::Symbol::AlphabetI'));

my @protsymbols = sort qw( * A B C D E F G H I K L M N P Q R S T V W X Y Z);
my %h = (
	 'Ala' => 'A',     'Asx' => 'B',
	 'Cys' => 'C',     'Asp' => 'D',
	 'Glu' => 'E',     'Phe' => 'F',
	 'Gly' => 'G',     'His' => 'H',
	 'Ile' => 'I',     'Lys' => 'K',
	 'Leu' => 'L',     'Met' => 'M',
	 'Asn' => 'N',     'Pro' => 'P',
	 'Gln' => 'Q',     'Arg' => 'R',
	 'Ser' => 'S',     'Thr' => 'T',
	 'Val' => 'V',     'Trp' => 'W',
	 'Tyr' => 'Y',     'Ter' => '*',
	 'Xaa' => 'X',     'Glx' => 'Z'    
	 );
my @protnms = sort { $h{$a} cmp $h{$b} } keys %h;
$count = 0;
foreach my $s ( sort { $a->token cmp $b->token } $prot->symbols ) {
    ok($s->name, $protnms[$count]);
    ok($s->token, $protsymbols[$count++]);    
}

