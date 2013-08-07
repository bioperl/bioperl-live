#
# BioPerl module for Bio::Tools::Protparam
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Richard Dobson, r.j.dobson at qmul dot ac dot uk
#
# Copyright Richard Dobson
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Protparam - submit to and parse output from protparam ;

=head1 SYNOPSIS

  use Bio::DB::GenBank;
  use Bio::Tools::Protparam;

  my $gb = new Bio::DB::GenBank(-retrievaltype => 'tempfile' , 
                                -format => 'Fasta');
  my @ids=qw(O14521 O43709 O43826);
  my $seqio = $gb->get_Stream_by_acc(\@ids );

  while( my $seq =  $seqio->next_seq ) {

	my $pp = Bio::Tools::Protparam->new(seq=>$seq->seq);

	print
	"ID : ", $seq->display_id,"\n",
	"Amino acid number : ",$pp->amino_acid_number(),"\n",
	"Number of negative amino acids : ",$pp->num_neg(),"\n",
	"Number of positive amino acids : ",$pp->num_pos(),"\n",
	"Molecular weight : ",$pp->molecular_weight(),"\n",
	"Theoretical pI : ",$pp->theoretical_pI(),"\n",
	"Total number of atoms : ", $pp->total_atoms(),"\n",
	"Number of carbon atoms : ",$pp->num_carbon(),"\n",
	"Number of hydrogen atoms : ",$pp->num_hydrogen(),"\n",
	"Number of nitrogen atoms : ",$pp->num_nitro(),"\n",
	"Number of oxygen atoms : ",$pp->num_oxygen(),"\n",
	"Number of sulphur atoms : ",$pp->num_sulphur(),"\n",
	"Half life : ", $pp->half_life(),"\n",
	"Instability Index : ", $pp->instability_index(),"\n",
	"Stability class : ", $pp->stability(),"\n",
	"Aliphatic_index : ",$pp->aliphatic_index(),"\n",
	"Gravy : ", $pp->gravy(),"\n",
	"Composition of A : ", $pp->AA_comp('A'),"\n",
	"Composition of R : ", $pp->AA_comp('R'),"\n",
	"Composition of N : ", $pp->AA_comp('N'),"\n",
	"Composition of D : ", $pp->AA_comp('D'),"\n",
	"Composition of C : ", $pp->AA_comp('C'),"\n",
	"Composition of Q : ", $pp->AA_comp('Q'),"\n",
	"Composition of E : ", $pp->AA_comp('E'),"\n",
	"Composition of G : ", $pp->AA_comp('G'),"\n",
	"Composition of H : ", $pp->AA_comp('H'),"\n",
	"Composition of I : ", $pp->AA_comp('I'),"\n",
	"Composition of L : ", $pp->AA_comp('L'),"\n",
	"Composition of K : ", $pp->AA_comp('K'),"\n",
	"Composition of M : ", $pp->AA_comp('M'),"\n",
	"Composition of F : ", $pp->AA_comp('F'),"\n",
	"Composition of P : ", $pp->AA_comp('P'),"\n",
	"Composition of S : ", $pp->AA_comp('S'),"\n",
	"Composition of T : ", $pp->AA_comp('T'),"\n",
	"Composition of W : ", $pp->AA_comp('W'),"\n",
	"Composition of Y : ", $pp->AA_comp('Y'),"\n",
	"Composition of V : ", $pp->AA_comp('V'),"\n",
	"Composition of B : ", $pp->AA_comp('B'),"\n",
	"Composition of Z : ", $pp->AA_comp('Z'),"\n",
	"Composition of X : ", $pp->AA_comp('X'),"\n";
}

=head1 DESCRIPTION

This module takes an amino acid sequence and submits it to the
Protparam program at www.expasy.org/cgi-bin/protparam.  Many
properties of the submitted sequence are returned.

=head1 AUTHOR

Richard Dobson, r.j.dobson at qmul dot ac dot uk

=cut

# Let the code begin...

package Bio::Tools::Protparam;

use strict;
use base qw(Bio::Root::Root);
use LWP 5.64;

=head2 new

  Title    : new
  Usage    : $pp = Protparam->new(seq=>$seq->seq);
  Function : Creates a new Protparam object
  Returns  : A Protparam object
  Args     : A sequence

=cut

sub new {
	my ($class,@args) = @_;
	@args=('-url'=>'http://web.expasy.org/cgi-bin/protparam/protparam','-form'=>'sequence',@args);
	my $self=$class->SUPER::new(@args);

	my ($url,$seq,$form)=$self->_rearrange([qw(URL SEQ FORM)],@args);

	my $browser = LWP::UserAgent->new;
	my $response;

	#send request to PROTPARAM @ Expasy
	$response = $browser->post($url,
                           [
                            $form => $seq
                           ],
                           'User-Agent' => 'Mozilla/4.76 [en] (Win2000; U)',
                          );

	#Check if successful
	$self->throw("$url error: ".$response->status_line) unless $response->is_success;
	$self->throw("Bad content type at $url ".$response->content_type) unless $response->content_type eq 'text/html';

	my $protParamOutput=$response->decoded_content;

	$self->{'output'}=$protParamOutput;

	return bless $self,$class;

}

=head2 num_neg

  Title    : num_neg
  Usage    : $pp->num_neg()
  Function : Retrieves the number of negative amino acids in a sequence
  Returns  : Returns the number of negative amino acids in a sequence
  Args     : none

=cut



sub num_neg{

	my $self=shift;

	($self->{'negAA'})=$self->{'output'}=~/<B>Total number of negatively charged residues.*?<\/B>\s*(\d*)/;

	return $self->{'negAA'};


}

=head2 num_pos

  Title    : num_pos
  Usage    : $pp->num_pos()
  Function : Retrieves the number of positive amino acids in a sequence
  Returns  : Returns the number of positive amino acids in a sequence
  Args     : none

=cut


sub num_pos{
	my $self=shift;
	($self->{'posAA'})=$self->{'output'}=~/<B>Total number of positively charged residues.*?<\/B>\s*(\d*)/;
	return $self->{'posAA'};
}

=head2 amino_acid_number

  Title    : amino_acid_number
  Usage    : $pp->amino_acid_number()
  Function : Retrieves the number of amino acids within a sequence
  Returns  : Returns the number of amino acids within a sequence
  Args     : none

=cut

sub amino_acid_number{
	my $self=shift;
	($self->{'numAA'})=$self->{'output'}=~/<B>Number of amino acids:<\/B> (\d+)/;
	return $self->{'numAA'};
}

=head2 total_atoms

  Title    : total_atoms
  Usage    : $pp->total_atoms()
  Function : Retrieves the total number of atoms within a sequence
  Returns  : Returns the total number of atoms within a sequence
  Args     : none

=cut


sub total_atoms{
	my $self=shift;
	$self->{'total_atoms'}=$self->{'output'}=~/<B>Total number of atoms:<\/B>\s*(\d*)/;
	return $self->{'total_atoms'};
}

=head2 molecular_weight

  Title    : molecular_weight
  Usage    : $pp->molecular_weight()
  Function : Retrieves the molecular weight of a sequence
  Returns  : Returns the molecular weight of a sequence
  Args     : none

=cut


sub molecular_weight{
	my $self=shift;
	($self->{'MolWt'})=$self->{'output'}=~/<B>Molecular weight:<\/B> (\d*\.{0,1}\d*)/;
	return $self->{'MolWt'};

}

=head2 theoretical_pI

  Title    : theoretical_pI
  Usage    : $pp->theoretical_pI()
  Function : Retrieve the theoretical pI for a sequence
  Returns  : Return the theoretical pI for a sequence
  Args     : none

=cut


sub theoretical_pI{
	my $self=shift;
	($self->{'TpI'})=$self->{'output'}=~/<B>Theoretical pI:<\/B> (-{0,1}\d*\.{0,1}\d*)/;
	return $self->{'TpI'};
}

=head2 num_carbon

  Title    : num_carbon
  Usage    : $pp->num_carbon()
  Function : Retrieves the number of carbon atoms in a sequence
  Returns  : Returns the number of carbon atoms in a sequence
  Args     : none

=cut


sub num_carbon{
	my $self=shift;
	($self->{'car'}) = $self->{'output'}=~/Carbon\s+C\s+(\d+)/;
	return $self->{'car'};
}

=head2 num_hydrogen

  Title    : num_hydrogen
  Usage    : $pp->num_hydrogen
  Function : Retrieves the number of hydrogen atoms in a sequence
  Returns  : Returns the number of hydrogen atoms in a sequence
  Args     : none

=cut


sub num_hydrogen{
	my $self=shift;
	($self->{'hyd'}) = $self->{'output'}=~/Hydrogen\s+H\s+(\d+)/;
	return $self->{'hyd'}
}

=head2 num_nitro

  Title    : num_nitro
  Usage    : $pp->num_nitro
  Function : Retrieves the number of nitrogen atoms in a sequence
  Returns  : Returns the number of nitrogen atoms in a sequence
  Args     : none

=cut


sub num_nitro{
	my $self=shift;
	($self->{'nitro'}) = $self->{'output'}=~/Nitrogen\s+N\s+(\d+)/;
	return $self->{'nitro'};
}

=head2 num_oxygen

  Title    : num_oxygen
  Usage    : $pp->num_oxygen()
  Function : Retrieves the number of oxygen atoms in a sequence
  Returns  : Returns the number of oxygen atoms in a sequence
  Args     : none

=cut


sub num_oxygen{
	my $self=shift;
	($self->{'oxy'}) = $self->{'output'}=~/Oxygen\s+O\s+(\d+)/;
	return $self->{'oxy'};
}

=head2 num_sulphur

  Title    : num_sulphur
  Usage    : $pp->num_sulphur()
  Function : Retrieves the number of sulphur atoms in a sequence
  Returns  : Returns the number of sulphur atoms in a sequence
  Args     : none

=cut


sub num_sulphur{
	my $self=shift;
	($self->{'sul'}) = $self->{'output'}=~/Sulfur\s+S\s+(\d+)/;
	return $self->{'sul'};
}

=head2 half_life

  Title    : half_life
  Usage    : $pp->half_life()
  Function : Retrieves the half life of a sequence
  Returns  : Returns the half life of a sequence
  Args     : none

=cut


sub half_life{
	my $self=shift;
	($self->{'half_life'}) = $self->{'output'}=~/The estimated half-life is.*?(-{0,1}\d*\.{0,1}\d*)\s*hours \(mammalian reticulocytes, in vitro\)/;
	return $self->{'half_life'};
}

=head2 instability_index

  Title    : instability_index
  Usage    : $pp->instability_index()
  Function : Retrieves the instability index of a sequence
  Returns  : Returns the instability index of a sequence
  Args     : none

=cut


sub instability_index{
	my $self=shift;
	($self->{'InstabilityIndex'})=$self->{'output'}=~/The instability index \(II\) is computed to be (-{0,1}\d*\.{0,1}\d*)/;
	return $self->{'InstabilityIndex'};
}

=head2 stability

  Title    : stability
  Usage    : $pp->stability()
  Function : Calculates whether the sequence is stable or unstable
  Returns  : 'stable' or 'unstable'
  Args     : none

=cut


sub stability{
	my $self=shift;
	($self->{'Stability'})=$self->{'output'}=~/This classifies the protein as\s(\w+)\./;
	return $self->{'Stability'};
}

=head2 aliphatic_index

  Title    : aliphatic_index
  Usage    : $pp->aliphatic_index()
  Function : Retrieves the aliphatic index of the sequence
  Returns  : Returns the aliphatic index of the sequence
  Args     : none

=cut


sub aliphatic_index{
	my $self=shift;
	($self->{'AliphaticIndex'})=$self->{'output'}=~/<B>Aliphatic index:<\/B>\s*(-{0,1}\d*\.{0,1}\d*)/;
	return $self->{'AliphaticIndex'};

}

=head2 gravy

  Title    : gravy
  Usage    : $pp->gravy()
  Function : Retrieves the grand average of hydropathicity (GRAVY) of a sequence
  Returns  : Returns the grand average of hydropathicity (GRAVY) of a sequence
  Args     : none

=cut


sub gravy{
	my $self=shift;
	($self->{'GRAVY'})=$self->{'output'}=~/<B>Grand average of hydropathicity \(GRAVY\):<\/B>\s*(-{0,1}\d*\.{0,1}\d*)/;
	return $self->{'GRAVY'};
}

=head2 AA_comp

  Title    : AA_comp
  Usage    : $pp->AA_comp('P')
  Function : Retrieves the percentage composition of a given amino acid for a sequence
  Returns  : Returns the percentage composition of a given amino acid for a sequence
  Args     : A single letter amino acid code eg A, R, G, P etc

=cut


sub AA_comp{
	my $self=shift;
	my $aa=shift;
	$aa=uc($aa);
	my $AA={qw(A Ala R Arg N Asn D Asp C Cys Q Gln E Glu G Gly H His I Ile L Leu K Lys M Met F Phe P Pro S Ser T Thr W Trp Y Tyr V Val B Asx Z Glx X Xaa)};
	($self->{$aa})= $self->{'output'}=~/$AA->{$aa} \($aa\)\s+\d+\s+(\d+\.\d+)%/;
	return $self->{$aa};
}


1;
