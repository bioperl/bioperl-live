# $Id $
#
# BioPerl module for Bio::PopGen::IndividualI
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::PopGen::IndividualI - An individual who has Genotype or Sequence Results

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


package Bio::PopGen::IndividualI;
use vars qw(@ISA);
use strict;

use Bio::Root::RootI;

@ISA = qw( Bio::Root::RootI );


=head2 unique_id

 Title   : unique_id
 Usage   : my $id = $individual->unique_id
 Function: Unique Identifier
 Returns : string representing unique identifier
 Args    : string


=cut

sub unique_id{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 get_Genotypes

 Title   : get_Genotypes
 Usage   : my @genotypes = $ind->get_Genotypes(-marker => $markername);
 Function: Get the genotypes for an individual, based on a criteria
 Returns : Array of genotypes
 Args    : either none (return all genotypes) or 
           -marker => name of marker to return (exact match, case matters)


=cut

sub get_Genotypes{
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 has_Marker

 Title   : has_Marker
 Usage   : if( $ind->has_Marker($name) ) {}
 Function: Boolean test to see if an Individual has a genotype 
           for a specific marker
 Returns : Boolean (true or false)
 Args    : String representing a marker name


=cut

sub has_Marker{
   my ($self,$name) = @_;
   $self->throw_not_implemented();
}

=head2 get_Marker_Names

 Title   : get_Marker_Names
 Usage   : my @names = $individual->get_Marker_Names;
 Function: Returns the list of known marker names
 Returns : List of strings
 Args    : none


=cut

sub get_Marker_Names{
   my ($self) = @_;
   $self->throw_not_implemented();
}


1;
