# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests           => 29,
               -requires_module => 'Data::Stag');

    use_ok('Bio::AlignIO');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# general filehandle tests
# not all parsers support output (noted as 0)
my %files = (
    # file                   format       I  O
    'testaln.phylip'     => ['phylip',     1, 1],
    'testaln.psi'         => ['psi',        1, 1],
    'testaln.arp'       => ['arp',        1, 0],
    'rfam_tests.stk'    => ['stockholm',  1, 1],
    'testaln.pfam'      => ['pfam',       1, 1],
    'testaln.msf'       => ['msf',        1, 1],
    'testaln.fasta'     => ['fasta',      1, 1],
    'testaln.selex'     => ['selex',      1, 1],
    'testaln.mase'      => ['mase',       1, 0],
    'testaln.prodom'    => ['prodom',     1, 0],
    'testaln.clustalw'  => ['clustalw',   1, 1],
    'testaln.metafasta' => ['metafasta',  1, 1],
    'testaln.nexus'     => ['nexus',      1, 1],
    'testaln.po'        => ['po',         1, 1],
    'testaln.xmfa'      => ['xmfa',       1, 1],
 );

# input file handles

$aln = Bio::AlignIO->new(
    -file  => test_input_file('longnames.aln'),
    -format=>'clustalw',
)->next_aln();
isa_ok($aln, 'Bio::AnnotatableI');

while (my ($file, $fdata) = each %files) {
    my ($format, $in, $out) = @{$fdata};
    if ($in) {
        my $fhin = Bio::AlignIO->newFh(
           '-file'  => test_input_file($file),
                           '-format' => $format);
        my $fhout = Bio::AlignIO->newFh(
           '-file' => ">".test_output_file(),
                        '-format' => 'clustalw');
        while ( $aln = <$fhin>) {
            cmp_ok($aln->num_sequences, '>=', 2, "input filehandle method test : $format");
            last;
        }
    }
}

# output file handles

while (my ($file, $fdata) = each %files) {
    my ($format, $in, $out) = @{$fdata};
    if ($out) {
        my $status = 0;
        my $fhin = Bio::AlignIO->newFh(
           '-file' => test_input_file('testaln.clustalw'),
                        '-format' => 'clustalw');
        my $fhout = Bio::AlignIO->newFh(
           '-file'  => '>'.test_output_file(),
                           '-format' => $format);
        while ( $aln = <$fhin> ) {
            $status = print $fhout $aln;
            last;
        }
        is $status, 1, "filehandle output test : $format";
    }
}
