# $Id$
# bioperl module for Bio::Variation::SNP
#
# Copyright Allen Day <allenday@ucla.edu>, Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

=head1 NAME

Bio::Variation::SNP - submitted SNP

=head1 SYNOPSIS

  $SNP = Bio::Variation::SNP->new (
                                           -id => 1234,
                                           -gene_symbol => 'AR'
                                           -chromosome => 'X',
                                           );
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

Allen Day <allenday@ucla.edu>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Variation::SNP;
my $VERSION=1.0;

use strict;
use vars qw($VERSION @ISA);
use Bio::Root::Root;
use Bio::Variation::SeqDiff;

@ISA = qw( Bio::Variation::SeqDiff );

our %RWSLOT = 	(
		heterozygous => '',
		heterozygous_SE => '',
		snp_class => '',
		organism => '',
		snp_type => '',
		comment => '',
		proof => '',
		het => '',
		het_SE => '',
		genotype => '',
		date => '',
		handle => '',
		batch => '',
		orient => '',
		accession => '',
		);

foreach my $slot (keys %RWSLOT){
	no strict "refs"; #allow add methods to package
	*$slot = sub {
		shift;
		$RWSLOT{$slot} = shift if @_;
		return $RWSLOT{$slot};
	};
}

1;
