# $Id$
# bioperl module for Bio::Variation::SNP
#
# Copyright Allen Day <allenday@ucla.edu>, Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

=head1 NAME

Bio::Variation::SNP - submitted SNP

=head1 SYNOPSIS

  $SNP = Bio::Variation::SNP->new ();
  # get an Allele object somehow
  $SNP->add_Variant($dnamut);

=head1 DESCRIPTION

Derived from Bio::Variation::SeqDiff, with additional methods that are
SNP specific (ie, refSNP/subSNP IDs, batch IDs, validation methods).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the 
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                         - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.  Bug reports can be submitted via
 email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR

Allen Day E<lt>allenday@ucla.eduE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Variation::SNP;
my $VERSION=1.0;

use strict;
use vars qw($VERSION @ISA $AUTOLOAD);
use Bio::Root::Root;
use Bio::Variation::SeqDiff;
use Bio::Variation::Allele;

@ISA = qw( Bio::Variation::SeqDiff );

my %OK_AUTOLOAD = map {$_=>1} qw(
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
		);

sub AUTOLOAD {
	my $self = shift;
	my $param = $AUTOLOAD;
	$param =~ s/.*:://;
	$self->throw(__PACKAGE__." doesn't implement $param") unless $OK_AUTOLOAD{$param};
	$self->{$param} = shift if @_;
	return $self->{$param};
}


#foreach my $slot (keys %RWSLOT){
#	no strict "refs"; #add class methods to package
#	*$slot = sub {
#		shift;
#		$RWSLOT{$slot} = shift if @_;
#		return $RWSLOT{$slot};
#	};
#}

sub is_subsnp {
	return shift->{is_subsnp};
}

sub subsnp {
	my $self = shift;
	return $self->{subsnps}->[ scalar($self->each_subsnp) - 1 ];
}

sub add_subsnp {
	my $self = shift;
	$self->throw("add_subsnp(): cannot add sunSNP to subSNP, only refSNP") if $self->is_subsnp;

	my $subsnp = Bio::Variation::SNP->new;
	push @{$self->{subsnps}}, $subsnp;
	$self->subsnp->{is_subsnp} = 1;
	return $self->subsnp;
}

sub each_subsnp {
	my $self = shift;
	$self->throw("each_subsnp(): cannot be called on a subSNP") if $self->is_subsnp;
	return @{$self->{subsnps}};
}

1;
