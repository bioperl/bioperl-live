

use Test;

BEGIN { plan tests => 2 };

use Bio::Tools::SiRNA;
use Bio::SeqIO;
use Bio::Root::IO;

# modules compile
print "ok 1\n";

my $input = Bio::SeqIO->new( -file 	=> Bio::Root::IO->catfile("t","data","NM_002254.gb"),
			     -format 	=> 'Genbank' );

my $seq = $input->next_seq;

my $sirna = Bio::Tools::SiRNA->new( -target 	=> $seq,
				    -dust	=> 1,
				    );

my @pairs = $sirna->design;

print "Got ",scalar(@pairs),"\n";

if( scalar(@pairs) == 65 ) {
    print "ok 2\n";
} else {
    print "not ok 2\n";
}
