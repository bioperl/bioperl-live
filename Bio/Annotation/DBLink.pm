
#
# BioPerl module for Bio::Annotation::Link
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Annotation::DBLink - DESCRIPTION of Object

=head1 SYNOPSIS

$link1 = new Bio::Annotation::DBLink(-database => 'TSC',
				     -primary_id => 'TSC0000030'
				     );

#or 

$link2 = new Bio::Annotation::DBLink();
$link2->database('dbSNP');
$link2->primary_id('2367');

# $feat is Bio::Annotation object, Bio::SeqFeature::Generic inherits it
$feat->add_DBLink($link2);


=head1 DESCRIPTION

Provides an object which represents a link from one onbject to something
in another database without proscribing what is in the other database

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

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

  my ($database, $primary_id, $optional_id, $comment) =
      $self->_rearrange([qw(DATABASE
			    PRIMARY_ID
			    OPTIONAL_ID 
			    COMMENT
			    )], @args);
  
  my $make = $self->SUPER::_initialize(@args); 

  $database    && $self->database($database);
  $primary_id  && $self->primary_id($primary_id);
  $optional_id && $self->optional_id($optional_id);
  $comment     && $self->comment($comment);
  
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
