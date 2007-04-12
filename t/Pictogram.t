# -*-Perl-*-
# $Id$
## Bioperl Test Harness Script for Modules


use strict;
use vars qw($NTESTS $SVG_AVAIL);

BEGIN {
    eval { require Test::More; };
    if( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
    $NTESTS = 5;
    
    eval {
		 require Bio::Graphics::Pictogram;
		 require SVG;
	};
	if ($@) {
		plan skip_all => "SVG not installed, skipping tests";
	} else {
		plan tests => $NTESTS;
	}
	use_ok('Bio::SeqIO');
	use_ok('Bio::Matrix::PSM::IO');
}


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
isa_ok $picto,"Bio::Graphics::Pictogram";

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




