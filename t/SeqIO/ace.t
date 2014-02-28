# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 7);

    use_ok 'Bio::SeqIO';
}

my $verbose = test_debug();

my $t_file = test_input_file('test.ace');
my $before;
{
    local $/ = undef;
    open my $BEFORE, '<', $t_file or die "Could not read file '$t_file': $!\n";
    $before = <$BEFORE>;
    close $BEFORE;
}

my $a_in = Bio::SeqIO->new( -FILE    => $t_file,
                            -verbose => $verbose,
                            -FORMAT  => 'ace' );
my @a_seq;
while (my $a = $a_in->next_seq) {
    push @a_seq, $a;
}

is @a_seq, 3, 'number of sequence objects';

my $esc_name = $a_seq[1]->display_id;
is $esc_name, 'Name; 4% strewn with \ various / escaped characters',
    "unescaping of characters, $esc_name";

is $a_seq[0]->alphabet, 'protein', 'alphabets detected';
is $a_seq[1]->alphabet, 'dna', 'alphabets detected';

my $o_file = test_output_file();
my $a_out = Bio::SeqIO->new( -FILE    => ">$o_file",
                             -verbose => $verbose,
                             -FORMAT  => 'ace' );
my $a_out_ok = 1;
for my $a (@a_seq) {
    $a_out->write_seq($a) or $a_out_ok = 0;
}
undef($a_out);  # Flush to disk
is $a_out_ok,1,'writing sequence';

my $after;
{
    local $/ = undef;
    open my $AFTER, '<', $o_file or die "Could not read file '$o_file': $!\n";
    $after = <$AFTER>;
    close $AFTER;
}

is( ($before and $after and ($before eq $after)), 1, 'test output');
