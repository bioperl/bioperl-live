use FindBin;
use File::Spec::Functions qw/catfile/;

use Test::More tests => 2;

use Bio::SearchIO;

my $searchio = Bio::SearchIO->new(
   -format => 'fasta',
   -file   => catfile($FindBin::Bin,'in.fasta'),
);

while ( my $result = $searchio->next_result() ) {
   while( my $hit = $result->next_hit ) {
      # process the Bio::Search::Hit::HitI object
      while( my $hsp = $hit->next_hsp ) {
         # process the Bio::Search::HSP::HSPI object
         is($hsp->num_identical, 30, "num_identical()");
         is($hsp->num_conserved, 30, "num_conserved()");
      }
   }
}

