
use Test;

BEGIN { plan tests => 2 };

use Bio::Tools::SiRNA;
use Bio::SeqIO;

# modules compile
print "ok 1\n";

my $input = Bio::SeqIO->new( -file 	=> 't/data/NM_002254.gb',
			     -format 	=> 'Genbank' );

my $seq = $input->next_seq;

my $sirna = Bio::Tools::SiRNA->new( -target 	=> $seq,
				    );

my @pairs = $sirna->design;

print "Got ",scalar(@pairs),"\n";

if( scalar(@pairs) == 65 ) {
    print "ok 2\n";
} else {
    print "not ok 2\n";
}
