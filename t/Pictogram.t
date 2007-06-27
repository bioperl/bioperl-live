# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 5,
			   -requires_modules => [qw(Bio::Graphics::Pictogram SVG)]);
	
	use_ok('Bio::SeqIO');
	use_ok('Bio::Matrix::PSM::IO');
}


my $file = test_input_file('pictogram.fa');
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

my $psmIO =  Bio::Matrix::PSM::IO->new(-format=>'meme', 
                                   -file=> test_input_file('meme.dat'));
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
