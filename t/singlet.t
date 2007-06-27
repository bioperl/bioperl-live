# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 4,
			   -requires_module => 'DB_File');
	
	use_ok('Bio::Assembly::IO');
}

my $aio = Bio::Assembly::IO->new(-file=>test_input_file(qw(consed_project edit_dir test_project.fasta.screen.ace.1)),
                                -format=>'ace');

my $assembly = $aio->next_assembly();
my @contigs = $assembly->all_contigs();
my @singlets = $assembly->all_singlets();

isa_ok $contigs[0], "Bio::Assembly::Contig",'Testing to see if the first contig is a Contig';

isa_ok $singlets[0], "Bio::Assembly::Singlet",'Testing to see if the first singlet is a Singlet';

isa_ok $singlets[0], 'Bio::Assembly::Contig', 'Testing to see if the Singlet ISA Contig';

# commented out entirely, no testing functionality - spiros denaxas
# this is what i really want to do:
#*** TODO?
# print("There were this many contigs: (".scalar(@contigs).")\n");
# print("There were this many singlets: (".scalar(@singlets).")\n");
# push @contigs,@singlets;
# print("This is a list of the ".scalar(@contigs)." contigs:\n");
# foreach my $contig (@contigs) {
#     # print &contig_dump($contig);
# }

###############################################

sub contig_dump {
     my ($contig) = @_;
     my $returner;
     #my $count = 1;
     my $prefix .= ("Contig: name(".$contig->id().") ");
     my @members = $contig->each_seq();
     if (!@members) { return $prefix." No Members\n"; }
     my $count = 1;
     foreach my $member (@members) {
          print("$prefix Member $count chromatfilename(".$member->{chromatfilename}.") phdfilenamename(".$member->{phdfilename}.") start(".$member->start().")\n");
          $count++;
     }
     return $returner;
}
