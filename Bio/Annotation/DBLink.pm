
#
# BioPerl module for Bio::Annotation::Link
#
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Annotation::DBLink - Database Link

=head1 SYNOPSIS

    $annotation = $seq->annotation;
    foreach my $link ( $annotation->each_DBLink() ) {
	print "Linked to ",$link->primary_id()," in ",$link->database,"\n";
    }


=head1 DESCRIPTION

Provides an object which represents a link from one onbject to something
in another database without proscribing what is in the other database.
It gives back strings for the "primary_id", something like the accession
number or the main id for the database and "database" being a string
representation of the database.

=head1 CONTACT

Ewan Birney <birney@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Annotation::DBLink;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;


@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 database

 Title   : database
 Usage   : $self->database($newval)
 Function: 
 Example : 
 Returns : value of database
 Args    : newvalue (optional)


=cut

sub database{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'database'} = $value;
    }
    return $self->{'database'};

}

=head2 primary_id

 Title   : primary_id
 Usage   : $self->primary_id($newval)
 Function: 
 Example : 
 Returns : value of primary_id
 Args    : newvalue (optional)


=cut

sub primary_id{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'primary_id'} = $value;
    }
    return $self->{'primary_id'};

}

=head2 optional_id

 Title   : optional_id
 Usage   : $self->optional_id($newval)
 Function: 
 Example : 
 Returns : value of optional_id
 Args    : newvalue (optional)


=cut

sub optional_id{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'optional_id'} = $value;
    }
    return $self->{'optional_id'};

}

=head2 comment

 Title   : comment
 Usage   : $self->comment($newval)
 Function: 
 Example : 
 Returns : value of comment
 Args    : newvalue (optional)


=cut

sub comment{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'comment'} = $value;
    }
    return $self->{'comment'};

}

1;
