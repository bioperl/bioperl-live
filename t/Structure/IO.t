# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 14);
	
	use_ok('Bio::Structure::IO');
}

#
# test reading PDB format - single model, single chain
#
my $pdb_file = test_input_file('1BPT.pdb');
ok my $structin = Bio::Structure::IO->new(-file => $pdb_file, 
										  -format => 'pdb');
ok my $struc = $structin->next_structure;
isa_ok $struc, "Bio::Structure::Entry";

# some basic checks, Structure objects are tested in Structure.t
my ($chain) = $struc->chain;
is($struc->residue, 97);
my ($atom) = $struc->get_atom_by_serial(367);
is($atom->id, "NZ");
is($struc->parent($atom)->id, "LYS-46");
my ($ann) = $struc->annotation->get_Annotations("author");
is($ann->as_text,
	"Value: D.HOUSSET,A.WLODAWER,F.TAO,J.FUCHS,C.WOODWARD              ");
($ann) = $struc->annotation->get_Annotations("header");
is($ann->as_text,
	"Value: PROTEINASE INHIBITOR (TRYPSIN)          11-DEC-91   1BPT");
my $pseq = $struc->seqres;
is($pseq->subseq(1,20), "RPDFCLEPPYTGPCKARIIR");

#
# test reading PDB format - single model, multiple chains
#
$pdb_file = test_input_file('1A3I.pdb');
$structin = Bio::Structure::IO->new(-file => $pdb_file, 
									-format => 'pdb');
$struc = $structin->next_structure;

my ($chaincount,$rescount,$atomcount);
for my $chain ($struc->get_chains) {
	$chaincount++;
   for my $res ($struc->get_residues($chain)) {
		$rescount++;
      for my $atom ($struc->get_atoms($res)) {
			$atomcount++;
		}
   }
}

is($chaincount, 4);  # 3 polypeptides and a group of hetero-atoms
is($rescount, 60);   # amino acid residues and solvent molecules
is($atomcount, 171); # ATOM and HETATM

#
# test reading PDB format - multiple models, single chain
#
#$pdb_file = test_input_file('1A11.pdb');
# TODO?

#
# test reading PDB format - chains with ATOMs plus HETATMs
#
#$pdb_file = test_input_file('8HVP.pdb');
# TODO?

#
# test writing PDB format
#
my $out_file = test_output_file();
my $structout = Bio::Structure::IO->new(-file => ">$out_file", 
                                        -format => 'PDB');
$structout->write_structure($struc);
ok -s $out_file;
