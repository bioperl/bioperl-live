# $Id$
#
# BioPerl module for Bio::AlternativeLocationHolderI
# Cared for by Lincoln Stein <lstein@cshl.org>
#
# Copyright Lincoln Stein
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::AlternativeLocationHolderI - Abstract interface of a Sequence that has multiple alternative locations

=head1 SYNOPSIS

    # get a Bio::AlternativeLocationHolderI somehow
    $holder->add_alternative_locations(Bio::Location::Simple->new(-start=>1,
                                                                -end=>10,
                                                                -seq_id=>'AF00128');
    @locations = $holder->alternative_locations();
    @locations = $holder->alternative_locations('AF00128');
    $holder->clear_alternative_locations();

=head1 DESCRIPTION

This Interface defines the methods for a Bio::AlternativeLocationHolderI,
an object which can have multiple alternative (secondary) locations in
addition to its primary one.  These alternative locations are to be
viewed as alternative coordinate systems for the *same* feature, such
as positions on alternative genomic assemblies.  The locations are
NOT:

  1) discontinuous locations on the same coordinate system, such as
     the disjunct ranges of a CDS.
nor

  2) the same feature type, such as an ALU, that is located in multiple
     places.

=head1 FEEDBACK

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org             - General discussion
  http://bio.perl.org/MailList.html - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Lincoln Stein

Email lstein@cshl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlternativeLocationHolderI;
use vars qw(@ISA);
use strict;
use Carp;

@ISA = qw(Bio::Root::RootI);

=head2 add_alternative_locations

 Title   : add_alternative_locations
 Usage   : $seqfeature->add_alternative_locations(@locationi)
 Function: adds new alternative location to the object
 Returns : void
 Args    : one or more LocationI-implementing object

 This adds one or more alternative locations to the feature.  These are
 to be viewed as alternative coordinate systems, such as
 assembly-to-assembly alignments, and not as alternative locations in
 the same coordinate space.

 Note: aliased in the interface to add_alternative_location(), for those
 who want to add a single location at a time.

=cut

sub add_alternative_locations {
    my $self = shift;
    my @locations = @_;
    $self->throw_not_implemented();
}
*add_alternative_location = \&add_alternative_locations;

=head2 alternative_locations

 Title   : alternative_locations
 Usage   : @locations = $seqfeature->alternative_locations([$seq_id])
 Function: returns alternative locations
 Returns : list of alternative locations
 Args    : optionally, a seq_id to filter on

=cut

sub alternative_locations {
    my $self = shift;
    my $seqid_filter = shift;
    $self->throw_not_implemented();
}

=head2 clear_alternative_locations

 Title   : clear_alternative_locations
 Usage   : $seqfeature->clear_alternative_locations([$seqid])
 Function: clears all alternative locations
 Returns : void
 Args    : optionally, a seq_id to clear locations on

=cut

sub clear_alternative_locations {
    my $self = shift;
    my $seqid_filter = shift;
    $self->throw_not_implemented();
}

1;

