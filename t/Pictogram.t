# -*-Perl-*-
## Bioperl Test Harness Script for Modules


use strict;
use vars qw($NTESTS $SVG_AVAIL);

BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    $NTESTS = 3;
    plan tests => $NTESTS;
    eval {
	require Bio::Graphics::Pictogram;
    };
    $SVG_AVAIL = $@ ? 0 : 1;
}

END {
    for ( $Test::ntest..$NTESTS ) {
        skip("SVG module not found. Skipping. ",1);
    }
}

if(!$SVG_AVAIL){
  warn("SVG not installed, skipping tests");
  exit;
}

use Bio::SeqIO;
use Bio::Matrix::PSM::IO;

my $file =  Bio::Root::IO->catfile("t","data","pictogram.fa");
my $sio = Bio::SeqIO->new(-file=>$file,-format=>'fasta');
my @seq;
while(my $seq = $sio->next_seq){
  push @seq, $seq;
}
my $picto = Bio::Graphics::Pictogram->new(-width=>"800",
                                          -fontsize=>"80",
                                          -plot_bits=>1,
                                          -color=>{'A'=>'red',
                                                   'G'=>'blue',
                                                   'C'=>'green',
                                                   'T'=>'magenta'});
ok $picto->isa("Bio::Graphics::Pictogram");

my $svg = $picto->make_svg(\@seq);
ok $svg->xmlify;

my $psmIO =  new Bio::Matrix::PSM::IO(-format=>'meme', 
                                   -file=> Bio::Root::IO->catfile(qw(t data meme.dat)));
$picto = Bio::Graphics::Pictogram->new(-width=>"800",
                                          -normalize=>1,
                                          -fontsize=>"80",
                                          -plot_bits=>1,
                                          -color=>{'A'=>'red',
                                                   'G'=>'blue',
                                                   'C'=>'green',
                                                   'T'=>'magenta'});

my $psm = $psmIO->next_psm;
$svg = $picto->make_svg($psm);
ok $svg->xmlify;






