#
# bioperl module for Bio::LiveSeq::Exon
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Joseph Insana <insana@ebi.ac.uk> <jinsana@gmx.net>
#
# Copyright Joseph Insana
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::LiveSeq::Exon - Range abstract class for LiveSeq

=head1 SYNOPSIS

  # documentation needed

=head1 DESCRIPTION

Class for EXON objects. They consist of a beginlabel, an endlabel (both
referring to a LiveSeq DNA object) and a strand.
The strand could be 1 (forward strand, default), -1 (reverse strand).

=head1 AUTHOR - Joseph A.L. Insana

Email:  Insana@ebi.ac.uk, jinsana@gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::LiveSeq::Exon;

use strict;
use base qw(Bio::LiveSeq::Range);

=head2 new

  Title   : new
  Usage   : $exon1 = Bio::LiveSeq::Exon-> new(-seq => $objref,
					      -start => $startlabel,
					      -end => $endlabel, -strand => 1);

  Function: generates a new Bio::LiveSeq::Exon
  Returns : reference to a new object of class Exon
  Errorcode -1
  Args    : two labels and an integer

=cut

=head2 get_Transcript

  Title   : get_Transcript
  Usage   : $transcript = $obj->get_Transcript()
  Function: retrieves the reference to the object of class Transcript (if any)
            attached to a LiveSeq object
  Returns : object reference
  Args    : none
  Note    : only Exons that compose a Transcript (i.e. those created out of
            a CDS Entry-Feature) will have an attached Transcript

=cut

sub get_Transcript {
  my $self=shift;
  return ($self->{'transcript'}); # this is set on all Exons a Transcript is made of when Transcript->new is called
}

# this checks if the attached Transcript has a Gene object attached
sub gene {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'gene'} = $value;
  }
  unless (exists $self->{'gene'}) {
    unless (exists $self->get_Transcript->{'gene'}) {
      return (0);
    } else {
      return ($self->get_Transcript->{'gene'});
    }
  } else {
    return $self->{'gene'};
  }
}

1;
