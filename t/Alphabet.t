# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 100);
	
	use_ok('Bio::Symbol::Alphabet');
	use_ok('Bio::Symbol::Symbol');
	use_ok('Bio::Symbol::DNAAlphabet');
	use_ok('Bio::Symbol::ProteinAlphabet');
}

my $A = Bio::Symbol::Symbol->new(-token => 'A' );
my $U = Bio::Symbol::Symbol->new(-token => 'U' );
my $G = Bio::Symbol::Symbol->new(-token => 'G' );
my $T = Bio::Symbol::Symbol->new(-token => 'T' );

my $rna = Bio::Symbol::Alphabet->new( -symbols => [ $A, $U, $G, $T ] );
				     
isa_ok($rna, 'Bio::Symbol::Alphabet');
my @symbols = $rna->symbols;
is(scalar @symbols, 4);

ok($rna->contains($A));
ok($rna->contains($T));
ok($rna->contains($U));
ok($rna->contains($G));

my $dna = Bio::Symbol::DNAAlphabet->new();
isa_ok($dna, 'Bio::Symbol::AlphabetI');
my $count = 0;

my @dnasymbols = sort qw( A B C D G H K M N R S T U V W X Y );
foreach my $s ( sort { $a->name cmp $b->name } $dna->symbols ) {
    is($s->name, $dnasymbols[$count]);    
    is($s->token, $dnasymbols[$count++]);    
}

my $prot = Bio::Symbol::ProteinAlphabet->new();
isa_ok($prot, 'Bio::Symbol::AlphabetI');

my @protsymbols = sort qw( * A B C D E F G H I J K L M N O P Q R S T U V W X Y Z);
my %h = (
	 'Ala' => 'A',     'Asx' => 'B',     'Cys' => 'C',     'Asp' => 'D',
	 'Glu' => 'E',     'Phe' => 'F',     'Gly' => 'G',     'His' => 'H',
	 'Ile' => 'I',     'Lys' => 'K',     'Leu' => 'L',     'Met' => 'M',
	 'Asn' => 'N',     'Pro' => 'P',     'Gln' => 'Q',     'Arg' => 'R',
	 'Ser' => 'S',     'Thr' => 'T',     'Val' => 'V',     'Trp' => 'W',
	 'Xaa' => 'X',     'Tyr' => 'Y',     'Glx' => 'Z',     'Ter' => '*',     
	 'Sec' => 'U',     'Pyl' => 'O',     'Xle' => 'J',     
	 );
my @protnms = sort { $h{$a} cmp $h{$b} } keys %h;
$count = 0;
foreach my $s ( sort { $a->token cmp $b->token } $prot->symbols ) {
    is($s->name, $protnms[$count]);
    is($s->token, $protsymbols[$count++]);    
}
