use FindBin;
use File::Spec::Functions qw/catfile/;

use Test::More tests => 2;

use Bio::SearchIO;

my $searchio = Bio::SearchIO->new(
   -format => 'exonerate',
   -file   => catfile($FindBin::Bin,'exonerate.output.works'),
);
parse($searchio);


$searchio = Bio::SearchIO->new(
   -format => 'exonerate',
   -file   => catfile($FindBin::Bin,'exonerate.output.dontwork'),
);
parse($searchio);


sub parse {
   my ($searchio) = @_;
   while( my $r = $searchio->next_result ) {
      while(my $hit = $r->next_hit){
         while(my $hsp = $hit->next_hsp){
            # print $hsp->start. "\t". $hsp->end. "\n";
         }
      }
      ok($r->query_name, "query_name");
   }
}



