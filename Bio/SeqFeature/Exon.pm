
#
# BioPerl module for Bio::SeqFeature::Exon
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Exon - Specialised Exon feature

=head1 SYNOPSIS

    foreach $exon ( $transcript->each_Exon() ) {
	print "Exon ", $exon->start, ":", $exon->end(), "\n";
    }

=head1 DESCRIPTION

Exon object for use with the transcript/gene set of information.
Exons are actually very simple objects. The ability to find out the
phase/utr'ness or protein coding abilities of an exon is encoded in
the relevant transcript and translation objects. Exon objects end up
being simply seqfeatures with hard coded primary keys.

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

=head1 AUTHOR - Ewan Birney

Email birney@sanger.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Exon;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::SeqFeature::Generic;

@ISA = qw(Bio::SeqFeature::Generic);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);
  
  $self->primary_tag('Exon');
  $self->source_tag("Bioperl");
# set stuff in self from @args
  return $make; # success - we hope!
}

1;







