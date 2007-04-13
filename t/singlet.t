# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#
# modeled after the t/Allele.t test script

use strict;
use vars qw($DEBUG $NUMTESTS $HAVE_DB_FILE);
$DEBUG = $ENV{'BIOPERLDEBUG'};
my $verbose = -1 unless $DEBUG;
BEGIN {
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	# bioperl takes far too long to compile.
	
	unshift(@INC,'Bio');
	
	eval { require Test::More; };
	if( $@ ) {
		use lib 't/lib';
	}
	
	use Test::More;
	eval { require Bio::Assembly::Contig;
			 require DB_File;
			 $HAVE_DB_FILE = 1;
		 };
	if( $@ ) {
		warn "No DB_File installed which is needed for Bio::Assembly::Contig\n";
		$HAVE_DB_FILE=0;
	}
	plan tests => ($NUMTESTS = 7);
}

SKIP {
	skip 'No DB_File', $NUMTESTS unless $HAVE_DB_FILE ;
}

exit(0) unless $HAVE_DB_FILE;

# redirect STDERR to STDOUT
open (STDERR, ">&STDOUT");

use_ok('Bio::Assembly::IO');
use_ok('Bio::Assembly::Singlet');
use_ok('Bio::Seq::SeqWithQuality');
use_ok('Bio::Seq::PrimaryQual') ;
use Dumpvalue() ; 

my $dumper = new Dumpvalue();
$dumper->veryCompact(1);


my $aio = Bio::Assembly::IO->new(-file=>File::Spec->catfile(qw(t data consed_project edit_dir test_project.fasta.screen.ace.1)),
                                -format=>'ace');

my $assembly = $aio->next_assembly();
my @contigs = $assembly->all_contigs();
my @singlets = $assembly->all_singlets();

ok(ref($contigs[0]) eq "Bio::Assembly::Contig",'Testing to see if the first contig is a Contig');

ok(ref($singlets[0]) eq "Bio::Assembly::Singlet",'Testing to see if the first singlet is a Singlet');

ok(UNIVERSAL::isa($singlets[0],'Bio::Assembly::Contig') , 'Testing to see if the Singlet ISA Contig');

# commented out entirely, no testing functionality - spiros denaxas
# this is what i really want to do:
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
