# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error = 0;

use strict;
BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }

    use Test;
    plan tests => 5; 

}

# REQUIREMENTS:
# handle SP fuzzies (eg ?s)
# handle OX lines

if( $error == 1 ) {
    exit(0);
}
END {
   unlink(qw (swiss_unk.dat));
}
use Bio::SeqIO;
use Bio::Root::IO;
my $verbose = $ENV{'BIOPERLDEBUG'};
ok(1);

my $seqio =
  new Bio::SeqIO( -verbose => $verbose,
                  -format => 'swiss',
                  -file   => Bio::Root::IO->catfile('t','data', 
                                                    'swiss_unk.dat'));

ok($seqio);
my $seq = $seqio->next_seq;
my @gns = $seq->annotation->get_Annotations('gene_name');
$seqio =
  new Bio::SeqIO( -verbose => $verbose,
                  -format => 'swiss',
                  -file   => Bio::Root::IO->catfile('>swiss_unk.dat'));

$seqio->write_seq($seq);

# reads it in once again
$seqio =
  new Bio::SeqIO( -verbose => $verbose,
                  -format => 'swiss',
                  -file   => Bio::Root::IO->catfile('swiss_unk.dat'));

$seq = $seqio->next_seq;
ok($seq->species);
ok($seq->species->ncbi_taxid eq "6239");
my @gns2 = $seq->annotation->get_Annotations('gene_name');
# check gene name is preserved (was losing suffix in worm gene names)
ok($#gns2 == 0 && $gns[0]->value eq $gns2[0]->value);
