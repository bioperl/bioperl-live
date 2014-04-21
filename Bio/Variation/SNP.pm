# bioperl module for Bio::Variation::SNP
#
# Copyright Allen Day <allenday@ucla.edu>, Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

=head1 NAME

Bio::Variation::SNP - submitted SNP

=head1 SYNOPSIS

  $SNP = Bio::Variation::SNP->new ();

=head1 DESCRIPTION

Inherits from Bio::Variation::SeqDiff and Bio::Variation::Allele, with 
additional methods that are (db)SNP specific (ie, refSNP/subSNP IDs, batch
IDs, validation methods).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

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
the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Allen Day E<lt>allenday@ucla.eduE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Variation::SNP;

use strict;
use vars qw($AUTOLOAD);
use Bio::Root::Root;

use base qw(Bio::Variation::SeqDiff Bio::Variation::Allele);

=head2 get/set-able methods

 Usage   : $is = $snp->method()
 Function: for getting/setting attributes
 Returns : a value.  probably a scalar.
 Args    : if you're trying to set an attribute, pass in the new value.

 Methods:
 --------
 id
 type
 observed
 seq_5
 seq_3
 ncbi_build
 ncbi_chr_hits
 ncbi_ctg_hits
 ncbi_seq_loc
 ucsc_build 
 ucsc_chr_hits
 ucsc_ctg_hits
 heterozygous
 heterozygous_SE
 validated
 genotype
 handle
 batch_id
 method
 locus_id
 symbol
 mrna
 protein
 functional_class

=cut

#'
my %OK_AUTOLOAD = (
		id			=> '',
		type			=> '',
		observed		=> [],
		seq_5			=> '',
		seq_3			=> '',
		ncbi_build		=> '',
		ncbi_chr_hits		=> '',
		ncbi_ctg_hits		=> '',
		ncbi_seq_loc		=> '',
		ucsc_build		=> '',
		ucsc_chr_hits		=> '',
		ucsc_ctg_hits		=> '',
		heterozygous		=> '',
		heterozygous_SE		=> '',
		validated		=> '',
		genotype		=> '',
		handle			=> '',
		batch_id		=> '',
		method			=> '',
		locus_id		=> '',
		symbol			=> '',
		mrna			=> '',
		protein			=> '',
		functional_class	=> '',
		);

sub AUTOLOAD {
	my $self = shift;
	my $param = $AUTOLOAD;
	$param =~ s/.*:://;
	$self->throw(__PACKAGE__." doesn't implement $param") unless defined $OK_AUTOLOAD{$param};

	if( ref $OK_AUTOLOAD{$param} eq 'ARRAY' ) {
		push @{$self->{$param}}, shift if @_;
		return $self->{$param}->[scalar(@{$self->{$param}}) - 1];
	} else {
		$self->{$param} = shift if @_;
		return $self->{$param};
	}
}


#foreach my $slot (keys %RWSLOT){
#	no strict "refs"; #add class methods to package
#	*$slot = sub {
#		shift;
#		$RWSLOT{$slot} = shift if @_;
#		return $RWSLOT{$slot};
#	};
#}


=head2 is_subsnp

 Title   : is_subsnp
 Usage   : $is = $snp->is_subsnp()
 Function: returns 1 if $snp is a subSNP
 Returns : 1 or undef
 Args    : NONE

=cut

sub is_subsnp {
	return shift->{is_subsnp};
}

=head2 subsnp

 Title   : subsnp
 Usage   : $subsnp = $snp->subsnp()
 Function: returns the currently active subSNP of $snp
 Returns : Bio::Variation::SNP
 Args    : NONE

=cut

sub subsnp {
	my $self = shift;
	return $self->{subsnps}->[ scalar($self->each_subsnp) - 1 ];
}

=head2 add_subsnp

 Title   : add_subsnp
 Usage   : $subsnp = $snp->add_subsnp()
 Function: pushes the previous value returned by subsnp() onto a stack,
           accessible with each_subsnp().
           Sets return value of subsnp() to a new Bio::Variation::SNP
           object, and returns that object.
 Returns : Bio::Varitiation::SNP
 Args    : NONE

=cut

sub add_subsnp {
	my $self = shift;
	$self->throw("add_subsnp(): cannot add subSNP to subSNP, only to refSNP")
            if $self->is_subsnp;

	my $subsnp = Bio::Variation::SNP->new;
	push @{$self->{subsnps}}, $subsnp;
	$self->subsnp->{is_subsnp} = 1;
	return $self->subsnp;
}

=head2 each_subsnp

 Title   : each_subsnp
 Usage   : @subsnps = $snp->each_subsnp()
 Function: returns a list of the subSNPs of a refSNP
 Returns : list
 Args    : NONE

=cut

sub each_subsnp {
	my $self = shift;
	$self->throw("each_subsnp(): cannot be called on a subSNP")
            if $self->is_subsnp;
	return @{$self->{subsnps}};
}

1;
