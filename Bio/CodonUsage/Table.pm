# $Id$
#
# BioPerl module for Bio::CodonUsage::Table
#
# Cared for by Richard Adams (richard.adams@ed.ac.uk)
#
# Copyright Richard Adams
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::CodonUsage::Table - for access to the Codon usage Database
at http://www.kazusa.or.jp/codon.

=head1 SYNOPSIS

  use Bio::CodonUsage::Table;
  use Bio::DB::CUTG;

  ## get  a codon usage table from web database ##
  my $cdtable = Bio::DB::CUTG->new(-sp => 'Mus musculus'
                                         -gc => 1);

  ## or from local file

  my $io = Bio::CodonUsage::IO->new(-file=>"file");
  my $cdtable= $io->next_data(); 

  print "leu frequency is ", $cdtable->aa_frequency('LEU'), "\n";
  print "freqof ATG is ", $cdtable->codon_rel_frequency('ttc'), "\n";
  print "abs freq of ATG is ", $cdtable->codon_abs_frequency('ATG'), "\n";
  print "number of ATG codons is ", $cdtable->codon_count('ATG'), "\n";
  print "gc content at position 1 is ", $cdtable->get_coding_gc('1'), "\n";
  print "total CDSs for Mus musculus  is ", $cdtable->cds_count(), "\n";

=head1 DESCRIPTION


This class provides methods for accessing codon usage table data.

All of the methods at present are simple look-ups of the table or are
derived from simple calculations from the table. Future methods could
include measuring the codon usage of a sequence , for example, or
provide methods for examining codon usage in alignments.

=head1 SEE ALSO

L<Bio::Tools::CodonTable>, 
L<Bio::WebAgent>, 
L<Bio::CodonUsage::IO>, 
L<Bio::DB::CUTG>

=head1 FEEDBACK

=head2 Mailing Lists


User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                       - General discussion
  http://bio.perl.org/MailList.html           - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHORS

Richard Adams, Richard.Adams@ed.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...





package Bio::CodonUsage::Table;
use strict;
use vars qw(@ISA %STRICTAA @AA);
use Bio::DB::CUTG;
use Bio::Root::IO;
use Bio::SeqUtils;
use Bio::Tools::CodonTable;

@ISA = qw(Bio::Root::Root);

BEGIN{
 @AA = qw(A C D E F G H I K L M N P Q R S T V W Y);
 map {$STRICTAA{$_} = undef} @AA;
}

=head2 new

 Title   : new
 Usage   : not instantiated directly, made by Bio::DB::CUTG.pm 
           or by Bio::CodonUsage::IO.pm
 Returns : a reference to a new  Bio::CodonUsage::Table object
 Args    : none

=cut

sub new {
	my ($class, @args) = @_;
	my $self= $class->SUPER::new(@args);
	return $self;
}

=head2 all_aa_frequencies

 Title   : all_aa_frequencies
 Usage   : my $freq = $cdtable->all_aa_frequencies();
 Returns : a reference to a hash where each key is an amino acid
           and each value is its frequency in all proteins in that
           species.
 Args    : none

=cut

sub all_aa_frequencies {
	my $self = shift;
	my %aa_freqs =();
	for my $aa (keys %STRICTAA) {
		my $freq = $self->aa_frequency($aa);
		$aa_freqs{$aa} = $freq;
		}
	return \%aa_freqs;
}

=head2 codon_abs_frequency

 Title   : codon_abs_frequency
 Usage   : my $freq = $cdtable->codon_abs_frequency('CTG');
 Purpose : To return the frequency of that codon as a percentage
           of all codons in the organism. 
 Returns : a percentage frequency
 Args    : a non-ambiguous codon string

=cut

sub codon_abs_frequency {
	my ($self, $a) = @_;
	my $cod = uc $a;
	if ($self->_check_codon($cod))  {
		my $ctable =  Bio::Tools::CodonTable->new;
		$ctable->id($self->{'_gen_code'});
		my $aa =$Bio::SeqUtils::THREECODE {$ctable->translate($cod)};

		return $self->{'_table'}{$aa}{$cod}{'per1000'}/10 ;
		}
	else {return 0;}
}

=head2 codon_rel_frequency

 Title   : codon_rel_frequency
 Usage   : my $freq = $cdtable->codon_rel_frequency('CTG');
 Purpose : To return the frequency of that codon as a percentage
           of codons coding for the same amino acid. E.g., ATG and TGG
           would return 100 as those codons are unique.
 Returns : a percentage frequency
 Args    : a non-ambiguous codon string

=cut


sub codon_rel_frequency {
	my ($self, $a) = @_;
	my $cod = uc $a;
	if ($self->_check_codon($cod)) {
		my $ctable =  Bio::Tools::CodonTable->new;
		$ctable->id($self->{'_gen_code'});
		my $aa =$Bio::SeqUtils::THREECODE {$ctable->translate($cod)};
		return $self->{'_table'}{$aa}{$cod}{'rel_freq'};
	}
	else {
		return 0;
		}
}

=head2 codon_count

 Title   : codon_count
 Usage   : my $count = $cdtable->codon_count('CTG');
 Purpose : To obtain the absolute number of the codons in the
           organism. 
 Returns : an integer
 Args    : a non-ambiguous codon string

=cut

sub codon_count {
	my ($self, $a) = @_;
	my $cod = uc $a;
	if ($self->_check_codon($cod)) {
		my $ctable =  Bio::Tools::CodonTable->new;
		$ctable->id($self->{'_gen_code'});

		my $aa =$Bio::SeqUtils::THREECODE {$ctable->translate($cod)};
		return $self->{'_table'}{$aa}{$cod}{'abs_count'};
		}
	else { return 0;}
}

=head2 get_coding_gc

 Title   : get_coding_gc
 Usage   : my $count = $cdtable->get_coding_gc(1);
 Purpose : To return the percentage GC composition for the organism at
           codon positions 1,2 or 3, or an average for all coding sequence
          ('all').
 Returns : a number (%-age GC content) or 0 if these fields are undefined
 Args    : 1,2,3 or 'all'.

=cut

sub get_coding_gc {
	my $self  = shift;
	if (! @_) {
		$self->warn(" no parameters supplied must be  a codon position (1,2,3) or 'all'");
		return 0;
			}
	else{
		my $n = shift;
		##return request if valid ##
		if ( exists($self->{'_coding_gc'}{$n} ) ) {
			return $self->{'_coding_gc'}{$n};
			}
		##else return 'all' value if exists
		elsif (exists($self->{'_coding_gc'}{'all'} )) {
			$self->warn("coding gc doesn't have value for [$n], returning gc content for all CDSs");
			return $self->{'_coding_gc'}{'all'};
			}
		### else return 0, 
		else {
			$self->warn("coding gc values aren't defined, returning 0");
			return 0;
		}

	}#end of outer else
		
}

=head2 set_coding_gc

 Title   : set_coding_gc
 Usage   : my $count = $cdtable->set_coding_gc(-1=>55.78);
 Purpose : To set the percentage GC composition for the organism at
           codon positions 1,2 or 3, or an average for all coding sequence
           ('all').  
 Returns : void
 Args    : a hash where the key must be 1,2,3 or 'all' and the value the %age GC
           at that codon position..

=cut

sub set_coding_gc {
	my ($self, $key, $value) = @_;
	my @allowed = qw(1 2 3 all);
	$key =~ s/\-//;
	if (!grep {$key eq $_} @allowed ) {
		$self->warn ("invalid key! - must be one of [ ". (join " ", @allowed) . "]");
		return;
		}
	$self->{'_coding_gc'}{$key} = $value;
	

}

=head2 species

 Title     : species
 Usage     : my $sp = $cut->species();
 Purpose   : Get/setter for species name of codon table
 Returns   : Void or species name string
 Args      : None or species name string

=cut

sub species {
	my $self = shift;
	if (@_ ){
		$self->{'_species'} = shift;
		}
	return $self->{'_species'};
}

=head2 genetic_code

 Title     : genetic_code
 Usage     : my $sp = $cut->genetic_code();
 Purpose   : Get/setter for genetic_code name of codon table
 Returns   : Void or genetic_code id, 1 by default
 Args      : None or genetic_code id, 1 by default if invalid argument.

=cut

sub genetic_code {
	my $self = shift;
	if (@_ ){
		my $val = shift;
		if ($val < 0 || $val >16 || $val =~ /[^\d]/ 
				|| $val ==7 || $val ==8) {
			$self->warn ("invalid genetic code - must be 1-16 but not 7 or 8,setting to default [1]");
			$self->{'_genetic_code'} = 1;
			}
		else {
			$self->{'_genetic_code'} = shift;
			}
		}
	return $self->{'_genetic_code'} || 1;
}

=head2 cds_count

 Title   : cds_count
 Usage   : my $count = $cdtable->cds_count();
 Purpose : To retrieve the total number of CDSs used to generate the Codon Table
           for that organism. 
 Returns : an integer
 Args    : none (if retrieving the value) or an integer( if setting ). 

=cut

sub cds_count {
	my $self= shift;
	if (@_) {
		my $val = shift;
		if ($val < 0) {
			$self->warn("can't have negative count initializing to 1");
			$self->{'_cds_count'} = 0;
			}
		else{
			$self->{'_cds_count'} = $val;
		}
	}
	$self->warn("cds_count value is undefined, returning 0") 
		if !exists($self->{'_cds_count'});

	return $self->{'_cds_count'} || 0;
	}

=head2 aa_frequency

 Title   : aa_frequency
 Usage   : my $freq = $cdtable->aa_frequency('Leu');
 Purpose : To retrieve the frequency of an amino acid in the organism
 Returns : a percentage
 Args    : a 1 letter or 3 letter string representing the amino acid

=cut

	

sub aa_frequency {
	my ($self, $a) = @_;
	## process args ##

	## deal with cases ##
	 my $aa = lc $a;
	$aa =~ s/^(\w)/\U$1/;


	if (!exists($STRICTAA{$aa}) && !exists($Bio::SeqUtils::ONECODE{$aa}) ) {
		$self->warn("Invalid amino acid! must be a unique 1 letter or 3 letter identifier");
		return;
		}
	#translate to 3 letter code for Ctable #
	my $aa3 = $Bio::SeqUtils::THREECODE{$aa} || $aa;

	## return % of all amino acids in organism ## 
	my $freq = 0;
	map {$freq += $self->{'_table'}{$aa3}{$_}{'per1000'} } keys %{$self->{'_table'}{$aa3}};
	return $freq/10;
}

## internal sub that checks a codon is correct format
sub _check_codon {
	my ($self, $cod) = @_;
	if ($cod =~ /[^ATCG]/  || $cod !~ /\w\w\w/) {
		$self->warn(" impossible codon - must be 3 letters and just containing ATCG");
		return 0;
	}
	else {return 1;}
}
return 1;
