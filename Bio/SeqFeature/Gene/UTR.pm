# $Id$
#
# BioPerl module for Bio::SeqFeature::Gene::UTR
#
# Cared for by David Block <dblock@gene.pbi.nrc.ca>
#
# Copyright David Block
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Gene::UTR - A feature representing an untranslated region
          that is part of a transcription unit

=head1 SYNOPSIS

See documentation of methods

=head1 DESCRIPTION

A UTR is a Bio::SeqFeature::Gene::ExonI compliant object that is
non-coding, and can be either 5' or 3' in a transcript.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - David Block

Email dblock@gene.pbi.nrc.ca

=head1 CONTRIBUTORS

This is based on the Gene Structure scaffolding erected by Hilmar Lapp
(hlapp@gmx.net).

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Gene::UTR;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::SeqFeature::Gene::NC_Feature;

@ISA = qw(Bio::SeqFeature::Gene::NC_Feature);

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
}

1;
