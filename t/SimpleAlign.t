

use strict;
BEGIN {     
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;

    plan tests => 6;
}
use Bio::SimpleAlign;
use Bio::Root::IO;
ok(1);

open(FH,Bio::Root::IO->catfile("t","data","test.mase")) || die "Could not open test.mase $!";
my $aln = Bio::SimpleAlign->new();
$aln->read_mase(\*FH);
close(FH);

ok( $aln );
open(OUT,">".Bio::Root::IO->catfile("t","data","out.aln_fasta")); 
$aln->write_fasta(\*OUT);
close(OUT);
ok(1);

$aln = Bio::SimpleAlign->new();
open(FH,Bio::Root::IO->catfile("t","data","test.pfam"));
$aln->read_Pfam(\*FH);
close(FH);

ok ( $aln );

open(OUT,">".Bio::Root::IO->catfile("t","data","out.pfam")); 
$aln->write_Pfam(\*OUT);
close(OUT);
ok(1);

$aln = Bio::SimpleAlign->new();
open(IN,Bio::Root::IO->catfile("t","data","out.pfam"));
$aln->read_Pfam(\*IN);
close(IN);

ok ( $aln );

unlink(Bio::Root::IO->catfile("t","data","out.pfam"), Bio::Root::IO->catfile("t","data","out.aln_fasta"));




