# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#
# modeled after the t/Allele.t test script

use strict;
use vars qw($DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'};
my $verbose = -1 unless $DEBUG;
BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan tests => 3;
}

# redirect STDERR to STDOUT
open (STDERR, ">&STDOUT");
use Bio::Assembly::IO;
use Bio::Assembly::Contig;
use Bio::Assembly::Singlet;
use Bio::Seq::SeqWithQuality;
use Bio::Seq::PrimaryQual;



my $aio = Bio::Assembly::IO->new(-file=>"<t/data/consed_project/edit_dir/test_project.fasta.screen.ace.1",
                                -forat=>'ace');

my $assembly = $aio->next_assembly();
my @contigs = $assembly->all_contigs();
my @singlets = $assembly->all_singlets();

print("Testing to see if the first contig is a Contig.\n");
my $first_contig = pop(@contigs);
ok(ref($first_contig) eq "Bio::Assembly::Contig");

print("Testing to see if the first singlet is a Singlet.\n");
my $first_singlet = pop(@singlets);
ok(ref($first_singlet) eq "Bio::Assembly::Singlet");

print("Testing to see if the Singlet ISA Contig.\n");
ok(UNIVERSAL::isa($first_singlet,'Bio::Assembly::Contig'));




