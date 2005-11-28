# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

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
	plan tests => 10;
}
use Bio::Structure::Entry;
use Bio::Structure::IO;
use Bio::Root::IO;
ok(1);

# test reading PDB format

my $pdb_file = Bio::Root::IO->catfile("t","data","pdb1bpt.ent"); # BPTI
my $structin = Bio::Structure::IO->new(-file => $pdb_file, 
													-format => 'PDB');
ok(1);
my $struc = $structin->next_structure;
ok(1);
ok(ref($struc), "Bio::Structure::Entry");

# some basic checks, Structure objects are tested in Structure.t
my ($chain) = $struc->chain;
ok($struc->residue, 97);
my ($atom) = $struc->get_atom_by_serial(367);
ok($atom->id, "NZ");
ok($struc->parent($atom)->id, "LYS-46");
my ($ann) = $struc->annotation->get_Annotations("author");
ok($ann->as_text,
	"Value: D.HOUSSET,A.WLODAWER,F.TAO,J.FUCHS,C.WOODWARD              ");
my $pseq = $struc->seqres;
ok($pseq->subseq(1,20), "RPDFCLEPPYTGPCKARIIR");


# test writing PDB format

my $out_file = Bio::Root::IO->catfile("t","data","temp-pdb1bpt.ent");
my $structout = Bio::Structure::IO->new(-file => ">$out_file", 
                                        -format => 'PDB');
$structout->write_structure($struc);
ok(1);

END {
#	unlink $out_file if -e $out_file;
}
