
#
# BioPerl module for Bio::SeqFeatureProducerI
#
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeatureProducerI - Abstract interface of a Sequence Feature Producer

=head1 SYNOPSIS

Will provide interface for basic methods for producing sequence features.

=head1 DESCRIPTION

Interface for objects that create Sequence Features such as sequence analysis, 
BLAST report parsing, or gene prediction software.


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich with Himlar Lapp

Jason Stajich <jason@chg.mc.duke.edu>
Himlar Lapp <hilmar.lapp@pharma.novartis.com>

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeatureProducerI;
use vars qw(@ISA);
use Bio::Root::RootI;

use strict;

# Object preamble - inheriets from Bio::Root::RootI
use Carp;

@ISA = qw(Bio::Root::RootI);

# utility method
#
# Prints out a method like:
# Abstract method stop defined in interface Bio::SeqFeatureI not implemented by package You::BadFeature
sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller)[1];
  
  confess "Abstract method '$caller' defined in interface Bio::SeqFeatureProducerI not implemented by pacakge $package";
}

=head2 add_features

 Title   : add_features
 Usage   : $featprod->add_features($seq);
 Function: Adds features to the sequence based on
           already parsed sequence data
 Returns : none
 Args    : Bio::Seq object

=cut

sub add_features {
    my $self = shift;
   $self->_abstractDeath();
}

=head2 _parse_rpt

 Title   : _parse_rpt
 Usage   : $seqprod->_parse_rpt($filename);
 Function: Reads in rpt file
 Returns : none
 Args    : Bio::Seq object

=cut

sub _parse_rpt_file {
    my ($self,$rpt) = @_;
   $self->_abstractDeath();
}

1;
