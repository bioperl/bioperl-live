# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 4);

    use_ok('Bio::Tools::GuessSeqFormat');
}

my $file = 'example.vcf';
ok (my $guesser = Bio::Tools::GuessSeqFormat->new(-file => test_input_file($file)));
isa_ok ($guesser, 'Bio::Tools::GuessSeqFormat');
is ($guesser->guess, 'vcf');
