
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

use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);


sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  my ($database, $primary_id, $optional_id, $comment) =
      $self->_rearrange([qw(DATABASE
			    PRIMARY_ID
			    OPTIONAL_ID
			    COMMENT
			    )], @args);
  
  $database    && $self->database($database);
  $primary_id  && $self->primary_id($primary_id);
  $optional_id && $self->optional_id($optional_id);
  $comment     && $self->comment($comment);
  
  return $self;
}

=head2 database

 Title   : database
 Usage   : $self->database($newval)
 Function: set/get on the database string. Databases are just
           a string here which can then be interpretted elsewhere
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
 Function: set/get on the primary id (a string)
           The primary id is the main identifier used for this object in 
           the database. Good examples would be accession numbers. The id
           is meant to be the main, stable identifier for this object
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
 Function: get/set for the optional_id (a string)
           optional id is a slot for people to use as they wish. The main
           issue is that some databases do not have a clean single string
           identifier scheme. It is hoped that the primary_id can behave like
           a reasonably sane "single string identifier" of objects, and people
           can use/abuse optional ids to their heart's content to provide
           precise mappings. 
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
 Function: get/set of comments (comment object)
           Sets or gets comments of this dblink, which is sometimes relevant
 Example : 
 Returns : value of comment (Bio::Annotation::Comment)
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
