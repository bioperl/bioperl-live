


use strict;
BEGIN {     
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    plan tests => 8;
}

use Bio::SeqIO;
use Bio::Root::IO;

my $in = Bio::SeqIO->new(-file => "<".Bio::Root::IO->catfile("t","data", 
	"multifa.seq") , '-format' => 'Fasta');
ok $in;
my $c=0;
while ( my $seq = $in->next_seq() ) {
    ok($seq);
    $c++;
}
ok $c,6, " missing sequences in the file";
