# $Id $
#
# BioPerl module for Bio::PopGen::GenotypeI
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::PopGen::GenotypeI - A marker and alleles for a specific individual

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the interface here

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

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::PopGen::GenotypeI;
use vars qw(@ISA);
use strict;

use Bio::Root::RootI;

@ISA = qw( Bio::Root::RootI );

=head2 marker_name

 Title   : marker_name
 Usage   : my $name = $genotype->marker_name();
 Function: Get the marker name for a genotype result
 Returns : string
 Args    : none


=cut

sub marker_name{
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 individual_id

 Title   : individual_id
 Usage   : my $indid = $genotype->individual_id();
 Function: Gets the individual id associated with a genotype
           This is effectively a back reference since we will typically
           associate a genotype with an individual with an 
           individual HAS-A genotype relationship.
 Returns : unique id string for an individual
 Args    : none


=cut

sub individual_id{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 get_Alleles

 Title   : get_Alleles
 Usage   : my @alleles = $genotype->get_Alleles();
 Function: Get the alleles for a given marker and individual
 Returns : array of alleles (strings in many implementations)
 Args    : none


=cut

sub get_Alleles{
   my ($self) = @_;
   $self->throw_not_implemented();
}

1;
