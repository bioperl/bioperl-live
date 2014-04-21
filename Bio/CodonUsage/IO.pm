#
# BioPerl module for Bio::CodonUsage::IO
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Richard Adams (richard.adams@ed.ac.uk)
#
# Copyright Richard Adams
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::CodonUsage::IO - for reading and writing codon usage tables to file

=head1 SYNOPSIS

  use Bio::CodonUsage::IO;

  ## read in a codon usage file
  my $io = Bio::CodonUsage::IO->new(-file => "in");
  my $cut = $io->next_data();

  ## write it out again
  my $out = Bio::CodonUsage::IO->new(-file => ">out");
  $out->write_data($cut);

=head1 DESCRIPTION

This class provides standard IO methods for reading and writing text files
of codon usage tables. These tables can initially be retrieved using
Bio::DB::CUTG. At present only this format is supported for read/write. 

Reading a CUTG will return a Bio::CodonUsage::Table object. 

=head1 SEE ALSO

L<Bio::Tools::CodonTable>, 
L<Bio::WebAgent>,
L<Bio::CodonUsage::Table>,
L<Bio::CodonUsage::IO>

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHORS

Richard Adams, Richard.Adams@ed.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin

package Bio::CodonUsage::IO;
use Bio::CodonUsage::Table;

use base qw(Bio::Root::IO);

=head2  new

 Title  : new
 Usage  : my $io = Bio::CodonUsage::IO->new(-file => "CUTfile");
 Purpose: To  read/write a Bio:CodonUsage::Table object  
 Returns: A  Bio:CodonUsage::IO object
 Args   : a file or file handle 

=cut

sub new  {
	my ($class , @args) = @_;
	my $self = $class->SUPER::new(@args);
}


=head2  next_data

 Title  : next_data
 Usage  : my $cut = $io->next_data();
 Purpose: To  obtain a Bio:CodonUsage::Table object 
 Returns: A  Bio:CodonUsage::Table object
 Args   : none

=cut

sub next_data {
	my $self = shift;
	my $cut = $self->_parse;
	return $cut;
}

=head2  write_data

 Title  : write_data
 Usage  : $io->write_data($cut);
 Purpose: To  write a CUT to file
 Returns: void
 Args   : a Bio::CodonUsage::Table object reference 

=cut


sub write_data {
	my ($self, $cut) = @_;
	if (!$cut || ! $cut->isa(Bio::CodonUsage::Table)) {
		$self->throw("must supply a Bio::CodonUsage::Table object for writing\n");
			}
	my $outstring = "Codon usage table\n\n";

	my $sp_string = $cut->species . "[" . $cut->_gb_db . "]  " .
					$cut->cds_count . "  CDS's\n\n";
	$outstring .= $sp_string;
	my $colhead = sprintf("%-9s%-9s%15s%12s%12s\n\n", "AmAcid",
							 "Codon", "Number", "/1000", "Fraction");
	$outstring .= $colhead;

	### now write bulk of codon data  ##
	my $ctable =  Bio::Tools::CodonTable->new;

	for my $f (qw(G A T C)) {
		for my $s (qw(G A T C)) {
			for my $t (qw(G A T C)) {
				$cod = $f . $s . $t;
				my $aa =$Bio::SeqUtils::THREECODE {$ctable->translate($cod)};
				my $codstr = sprintf("%-9s%-9s%15.2f%12.2f%12.2f\n",		

						$aa, $cod, my $tt = $cut->codon_count($cod)|| 0.00, 
						my $ll =$cut->{'_table'}{$aa}{$cod}{'per1000'}|| 0.00,
						my $ss = $cut->codon_rel_frequency($cod) || 0.00);
				$outstring .= $codstr;
			}
		$outstring .= "\n";
		}
	}
	$outstring .= "\n\n";

	## now append GC data
	$outstring .= "Coding GC ". $cut->get_coding_gc('all'). "%\n";
	$outstring .= "1st letter GC ". $cut->get_coding_gc(1). "%\n";
	$outstring .= "2nd letter GC ". $cut->get_coding_gc(2). "%\n";
	$outstring .= "3rd letter GC ". $cut->get_coding_gc(3). "%\n";
	$outstring .= "Genetic code " . $cut->genetic_code() ."\n\n\n";

$self->_print ($outstring);
$self->flush();

}

sub _parse {
	my $self = shift;
	my $cdtableobj = Bio::CodonUsage::Table->new();
	while (my $line = $self->_readline() ) {
		next if $line =~ /^$/ ;
		$line =~ s/End/Ter/;
		## now parse in species name, cds number

		if ($line =~ /^(.+?)\s*\[(\w+)\].+?(\d+)/) {
			$cdtableobj->species($1);
			$cdtableobj->{'_gb_db'} = $2;
			$cdtableobj->cds_count($3);
			}

		## now parse in bulk of codon usage table
		elsif ( $line =~ /^(\w\w\w)\s+(\w+)\s+(\d+\.\d+)
					\s+(\d+\.\d+)\s+(\d+\.\d+)/x){
			if (defined ($1)) {
				$cdtableobj->{'_table'}{$1}{$2} = {
						'abs_count'=>$3,
						 'per1000'=> $4,
						 'rel_freq'=> $5,
						};
				}
			}

		## now parse in gc data  ####
		if($line =~ /^Cod.+?(\d\d\.\d\d)/ ){
				$cdtableobj->{'_coding_gc'}{'all'} = $1;
					}
		elsif ($line =~ /^1st.+?(\d\d\.\d\d)/){ 
				$cdtableobj->{'_coding_gc'}{'1'} = $1;
			}
		elsif($line =~ /^2nd.+?(\d\d\.\d\d)/){ 
				$cdtableobj->{'_coding_gc'}{'2'} = $1;
				}
		elsif ($line =~ /^3rd.+?(\d\d\.\d\d)/){ 
				$cdtableobj->{'_coding_gc'}{'3'} = $1;
				}

		elsif	($line =~ /^Gen.+?(\d+)/){ 
				$cdtableobj->{'_genetic_code'} = $1;
			}
	}
		## check has been parsed ok ##
		if (scalar keys %{$cdtableobj->{'_table'}} != 21) {
			$cdtableobj->warn("probable parsing error - should be 21 entries for 20aa + stop codon");
		}
		return $cdtableobj;
		
}

1;

__END__

