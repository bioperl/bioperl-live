# -*-Perl-*- Test Harness script for Bioperl

use strict;

BEGIN {
    use Bio::Root::Test;

    test_begin( -tests => 12 );

    use_ok('Bio::SeqIO::pir');
}

my $verbose = test_debug();

my $in = Bio::SeqIO->new(
    -file    => test_input_file('seqfile.pir'),
    -verbose => $verbose,
    -format  => 'pir'
);

ok( defined $in, 'new instance is defined ' );
isa_ok( $in, 'Bio::SeqIO' );

my $out = Bio::SeqIO->new(
    -format => 'pir',
    -fh     => \*STDOUT
);

while ( my $seq = $in->next_seq() ) {
    ok( $seq->length > 1, 'checked length' );
    $out->write_seq($seq) if $verbose > 0;
}

# Empty description line
$in = Bio::SeqIO->new(
    -file    => test_input_file('seqfile-no-desc.pir'),
    -verbose => $verbose,
    -format  => 'pir'
);
my $seq = $in->next_seq();
ok( $seq->seq =~ /^MGD/, 'Correct start' );
$seq = $in->next_seq();
ok( $seq->seq =~ /^GDV/, 'Correct start' );
$seq = $in->next_seq();
ok( $seq->seq =~ /^GDV/, 'Correct start' );
