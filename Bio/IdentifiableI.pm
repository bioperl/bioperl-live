# $Id$

#
# This module is licensed under the same terms as Perl itself. You use,
# modify, and redistribute it under the terms of the Perl Artistic License.
#

=head1 NAME

Bio::IdentifiableI - interface for objects with identifiers

=head1 SYNOPSIS


    # to test this is an identifiable object

    $obj->isa("Bio::IdentifiableI") || 
      $obj->throw("$obj does not implement the Bio::IdentifiableI interface");

    # accessors

    $object_id = $obj->object_id();
    $namespace = $obj->namespace();
    $authority = $obj->authority();
    $version   = $obj->version();

    # utility function

    $lsid        = $obj->lsid_string();      # gives authority:namespace:object_id
    $ns_string   = $obj->namespace_string(); # gives namespace:object_id.version


=head1 DESCRIPTION

This interface describes methods expected on identifiable objects, ie
ones which have identifiers expected to make sense across a number of
instances and/or domains. This interface is modeled after pretty much
ubiquitous ideas for names in bioinformatics being 

 databasename:object_id.version

examples being

 swissprot:P012334.2

or

 GO:0007048

We also work well with LSID proposals which adds in the concept of an
authority, being the DNS name of the organisation assigning the namespace.
Helper functions are provided to make useful strings being


  lsid_string - string complying to the LSID standard
  namespace_string - string complying to the usual convention of 
     namespace:object_id.version

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                       - General discussion
  http://bio.perl.org/MailList.html           - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Ewan Birney

Email birney@sanger.ac.uk

=cut

package Bio::IdentifiableI;
use vars qw(@ISA );
use strict;
use Bio::Root::RootI;


@ISA = qw(Bio::Root::RootI);

=head1 Implementation Specific Functions

These functions are the ones that a specific implementation must
define.

=head2 object_id

 Title   : object_id
 Usage   : $string    = $obj->object_id()
 Function: a string which represents the stable primary identifier
           in this namespace of this object. For DNA sequences this
           is its accession_number, similarly for protein sequences

 Returns : A scalar
 Status  : Virtual

=cut

sub object_id {
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 version

 Title   : version
 Usage   : $version    = $obj->version()
 Function: a number which differentiates between versions of
           the same object. Higher numbers are considered to be
           later and more relevant, but a single object described
           the same identifier should represent the same concept

 Returns : A number
 Status  : Virtual

=cut

sub version {
   my ($self) = @_;
   $self->throw_not_implemented();
}


=head2 authority

 Title   : authority
 Usage   : $authority    = $obj->authority()
 Function: a string which represents the organisation which
           granted the namespace, written as the DNS name for  
           organisation (eg, wormbase.org)

 Returns : A scalar
 Status  : Virtual

=cut

sub authority {
   my ($self) = @_;
   $self->throw_not_implemented();
}


=head2 namespace

 Title   : namespace
 Usage   : $string    = $obj->namespace()
 Function: A string representing the name space this identifier
           is valid in, often the database name or the name
           describing the collection 

 Returns : A scalar
 Status  : Virtual

=cut

sub namespace {
   my ($self) = @_;
   $self->throw_not_implemented();
}



=head1 Implementation optional functions

These functions are helper functions that are provided by
the interface but can be overridden if so wished

=head2 lsid_string

 Title   : lsid_string
 Usage   : $string   = $obj->lsid_string()
 Function: a string which gives the LSID standard
           notation for the identifier of interest


 Returns : A scalar

=cut

sub lsid_string {
  my ($self) = @_;

  return $self->authority.":".$self->namespace.":".$self->object_id;
}



=head2 namespace_string

 Title   : namespace_string
 Usage   : $string   = $obj->namespace_string()
 Function: a string which gives the common notation of
           namespace:object_id.version

 Returns : A scalar

=cut

sub namespace_string {
  my ($self) = @_;

  return $self->namespace.":".$self->object_id .
      (defined($self->version()) ? ".".$self->version : '');
}

1;
