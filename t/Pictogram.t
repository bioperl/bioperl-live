# -*-Perl-*-
## Bioperl Test Harness Script for Modules


use strict;
BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    use vars qw($NTESTS);
    $NTESTS = 2;
    plan tests => $NTESTS;
}

END {
    for ( $Test::ntest..$NTESTS ) {
        skip("SVG module not found. Skipping. ",1);
    }
}
eval {
   require('SVG.pm');
};
if($@){
  warn("SVG not installed, skipping tests");
  exit;
}

use Bio::Graphics::Pictogram;
use Bio::SeqIO;
my $file =  Bio::Root::IO->catfile("t","data","pictogram.fa");
my $sio = Bio::SeqIO->new(-file=>$file,-format=>'fasta');
my @seq;
while(my $seq = $sio->next_seq){
  push @seq, $seq;
}
my $picto = Bio::Graphics::Pictogram->new(-width=>"800",-fontsize=>"80",-plot_bits=>1,-color=>{'A'=>'red',
                                                                                   'G'=>'blue',
                                                                                   'C'=>'green',
                                                                                   'T'=>'magenta'});
ok $picto->isa("Bio::Graphics::Pictogram");

my $svg = $picto->make_svg(\@seq);
ok $svg->xmlify;





