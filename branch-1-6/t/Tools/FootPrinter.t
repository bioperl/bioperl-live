# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 27);
	
	use_ok('Bio::Tools::FootPrinter');
}

my $inputfilename= test_input_file('footprinter.out');
my $parser = Bio::Tools::FootPrinter->new(-file => $inputfilename);
my @sub;
my @species = qw(TETRAODON CHICKEN MOUSE HAMSTER HUMAN PIG);
while (my $feat = $parser->next_feature){
    is($feat->seq_id, shift @species);
    foreach my $sub ($feat->sub_SeqFeature){
      push @sub,$sub;
    }
}

is $sub[0]->seq_id, 'TETRAODON-motif1';
is $sub[0]->start,352;
is $sub[0]->end,362;
is $sub[0]->seq->seq,'tacaggatgca';
is $sub[1]->seq_id, 'TETRAODON-motif2';
is $sub[1]->start,363;
is $sub[1]->end,373;
is $sub[1]->seq->seq,'ccatatttgga';

is $sub[2]->seq_id, 'CHICKEN-motif1';
is $sub[2]->start,363;
is $sub[2]->end,373;
is $sub[2]->seq->seq,'cacaggatgta';

is $sub[3]->seq_id, 'CHICKEN-motif2';
is $sub[3]->start,376;
is $sub[3]->end,386;
is $sub[3]->seq->seq,'ccatataagga';

is $sub[4]->seq_id, 'MOUSE-motif1';
is $sub[4]->start,234;
is $sub[4]->end,243;
is $sub[4]->seq->seq,'cacaggatgt';
