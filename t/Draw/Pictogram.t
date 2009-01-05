# -*-Perl-*- Test Harness script for Bioperl
# $Id: Pictogram.t 15112 2008-12-08 18:12:38Z sendu $

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 6,
			   -requires_module => 'SVG 2.26');
	
	use_ok('Bio::Draw::Pictogram');
	use_ok('Bio::SeqIO');
	use_ok('Bio::Matrix::PSM::IO');
}


my $file = test_input_file('pictogram.fa');
my $sio = Bio::SeqIO->new(-file=>$file,-format=>'fasta');
my @seq;
while(my $seq = $sio->next_seq){
  push @seq, $seq;
}
my $picto = Bio::Draw::Pictogram->new(-width=>"800",
									-fontsize=>"80",
									-plot_bits=>1,
									-color=>{'A'=>'red',
											 'G'=>'blue',
											 'C'=>'green',
											 'T'=>'magenta'});
isa_ok $picto,"Bio::Draw::Pictogram";

my $svg = $picto->make_svg(\@seq);
ok $svg->xmlify;

my $psmIO =  Bio::Matrix::PSM::IO->new(-format=>'meme', 
                                   -file=> test_input_file('meme.dat'));
$picto = Bio::Draw::Pictogram->new(-width=>"800",
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
