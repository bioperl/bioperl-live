# $Id$
#
# BioPerl module for Bio::Biblio::RefI
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::RefI - Abstract interface for a bibliographic reference
classes

=head1 SYNOPSIS

  # not instantiable

=head1 DESCRIPTION

Super class and interface class for bibliographic references. The
central class of the Bio::Biblio name space.

The class names and attributes in Martin Senger's java implementation
of Biblio objects are given in squate brackets in the documentation
for each method. The openBQS project pages are at
http://industry.ebi.ac.uk/openBQS/.

See Martin's the UML diagram at 
http://industry.ebi.ac.uk/openBQS/images/bibobjects_java.jpg
for an overview.


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

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Heikki Lehvaslaiho

Email heikki@ebi.ac.uk

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Biblio::RefI;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI

use Bio::Root::Root;


@ISA = qw(Bio::Root::Root );


=head2 id

 Title   : id
 Usage   : $obj->id();
 Function: 

           Sets and returns the name of the reference.
           [BibRef::identifer]

 Example : 
 Returns : string
 Args    : string

=cut


sub id {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_id'} = $value;
    }
    if ( ! exists $self->{'_id'} ) {
	return 0;
    } 
    return $self->{'_id'};
}


=head2 type

 Title   : type
 Usage   : $obj->type();
 Function: 

           Sets and returns the name of the reference.
           [BibRef::type]

 Example : 
 Returns : string
 Args    : string

=cut


sub type {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_type'} = $value;
    }
    if ( ! exists $self->{'_type'} ) {
	return 0;
    } 
    return $self->{'_type'};
}

=head2 title

 Title   : title
 Usage   : $obj->title();
 Function: 

           Sets and returns the name of the reference.
           [BibRef::title]

 Example : 
 Returns : string
 Args    : string

=cut


sub title {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_title'} = $value;
    }
    if ( ! exists $self->{'_title'} ) {
	return 0;
    } 
    return $self->{'_title'};
}


=head2 rights

 Title   : rights
 Usage   : $obj->rights();
 Function: 

           Sets and returns the name of the reference.
           [BibRef::rights]

 Example : 
 Returns : string
 Args    : string

=cut


sub rights {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_rights'} = $value;
    }
    if ( ! exists $self->{'_rights'} ) {
	return 0;
    } 
    return $self->{'_rights'};
}

=head2 language

 Title   : language
 Usage   : $obj->language();
 Function: 

           Sets and returns the name of the reference.
           [BibRef::language]

 Example : 
 Returns : string
 Args    : string

=cut


sub language {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_language'} = $value;
    }
    if ( ! exists $self->{'_language'} ) {
	return 0;
    } 
    return $self->{'_language'};
}


=head2 format

 Title   : format
 Usage   : $obj->format();
 Function: 

           Sets and returns the name of the reference.
           [BibRef::format]

 Example : 
 Returns : string
 Args    : string

=cut


sub format {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_format'} = $value;
    }
    if ( ! exists $self->{'_format'} ) {
	return 0;
    } 
    return $self->{'_format'};
}


=head2 date

 Title   : date
 Usage   : $obj->date();
 Function: 

           Sets and returns the name of the reference.
           [BibRef::date]

 Example : 
 Returns : string
 Args    : string

=cut


sub date {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_date'} = $value;
    }
    if ( ! exists $self->{'_date'} ) {
	return 0;
    } 
    return $self->{'_date'};
}


=head2 add_DBLink

 Title   : add_DBLink
 Usage   : $self->add_DBLink($allele)
 Function: 

	    Adds one Bio::Variation::DBLink into the list of alleles.
            Note that the method forces the convention that nucleotide
            sequence is in lower case and amino acds are in upper
            case.

 Example : 
 Returns : 1 when succeeds, 0 for failure.
 Args    : DBLink object

=cut


sub add_DBLink {
  my ($self,$value) = @_;
  if (defined $value) {
      if( ! $value->isa('Bio::Annotation::DBLink') ) {
	  my $com = ref $value;
	  $self->throw("Is not a DBLink object but a  [$com]");
      }
      push(@{$self->{'_dblinks'}},$value); 
      return 1;
  } else {
      return 0;
  }
}


=head2 each_DBLink

 Title   : alleles
 Usage   : $obj->each_DBLink();
 Function: 

	     Returns a list of Bio::Annotation::DBLink objects

 Example : 
 Returns : list of DBLinks
 Args    : none

=cut

sub each_DBLink{
   my ($self,@args) = @_;
   return @{$self->{'_dblinks'}};
}

=head2 Scope methods

=cut

=head2 spatial_location

 Title   : spatial_location
 Usage   : $obj->spatial_location();
 Function: 

           Sets and returns the name of the reference.
           [BiblioScope::spatialLocation]

 Example : 
 Returns : string
 Args    : string

=cut


sub spatial_location {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_spatial_location'} = $value;
    }
    if ( ! exists $self->{'_spatial_location'} ) {
	return 0;
    } 
    return $self->{'_spatial_location'};
}


=head2 temporal_period

 Title   : temporal_period
 Usage   : $obj->temporal_period();
 Function: 

           Sets and returns the name of the reference.
           [BiblioScope::temporalPeriod]

 Example : 
 Returns : string
 Args    : string

=cut


sub temporal_period {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_temporal_period'} = $value;
    }
    if ( ! exists $self->{'_temporal_period'} ) {
	return 0;
    } 
    return $self->{'_temporal_period'};
}

=head2 EntryStatus methods

=cut

=head2 last_modified

 Title   : last_modified
 Usage   : $obj->last_modified();
 Function: 

           Sets and returns the name of the reference.
           [BiblioentryStatus::lastModifiedDate]

 Example : 
 Returns : string
 Args    : string

=cut


sub last_modified {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_last_modified'} = $value;
    }
    if ( ! exists $self->{'_last_modified'} ) {
	return 0;
    } 
    return $self->{'_last_modified'};
}


=head2 repository_subset

 Title   : repository_subset
 Usage   : $obj->repository_subset();
 Function: 

           Sets and returns the name of the reference.
           [BiblioEntryStatus::repositorySubset]

 Example : 
 Returns : string
 Args    : string

=cut


sub repository_subset {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_repository_subset'} = $value;
    }
    if ( ! exists $self->{'_repository_subset'} ) {
	return 0;
    } 
    return $self->{'_repository_subset'};
}

=head2 Subject methods

[BiblioSubject needs to be added...]

=cut

=head2 Description methods

No need to implement a saparate abstract_language method [BiblioDescription::language].

=cut

=head2 abstract

 Title   : abstract
 Usage   : $obj->abstract();
 Function: 

           Sets and returns the name of the reference.
           [BiblioDescription::theAbstract]

 Example : 
 Returns : string
 Args    : string

=cut


sub abstract {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_abstract'} = $value;
    }
    if ( ! exists $self->{'_abstract'} ) {
	return 0;
    } 
    return $self->{'_abstract'};
}


=head2 abstract_type

 Title   : abstract_type
 Usage   : $obj->abstract_type();
 Function: 

           Sets and returns the name of the reference.
           [Biblio::abstractType]

 Example : 
 Returns : string
 Args    : string

=cut


sub abstract_type {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_abstract_type'} = $value;
    }
    if ( ! exists $self->{'_abstract_type'} ) {
	return 0;
    } 
    return $self->{'_abstract_type'};
}


=head2 toc

 Title   : toc
 Usage   : $obj->toc();
 Function: 

           Sets and returns the name of the reference.
           [Biblio::tableOfContents]

 Example : 
 Returns : string
 Args    : string

=cut


sub toc {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_toc'} = $value;
    }
    if ( ! exists $self->{'_toc'} ) {
	return 0;
    } 
    return $self->{'_toc'};
}

=head2 toc_type

 Title   : toc_type
 Usage   : $obj->toc_type();
 Function: 

           Sets and returns the name of the reference.
           [Biblio::tableOfContentsType]

 Example : 
 Returns : string
 Args    : string

=cut


sub toc_type {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_toc_type'} = $value;
    }
    if ( ! exists $self->{'_toc_type'} ) {
	return 0;
    } 
    return $self->{'_toc_type'};
}


=head2 publisher

 Title   : publisher
 Usage   : $publisher = $obj->publisher();
 Function: Returns or sets the reference to a L<Bio::Biblio::Provider> object.
           If there is no link, it will return undef
 Returns : an obj_ref or undef
 Args    : Bio::Biblio::Provider object

=cut

sub publisher {
  my ($self,$value) = @_;
  if (defined $value) {
      if( ! $value->isa('Bio::Biblio::ProviderI') ) {
	  $self->throw("Is not a Bio::Biblio::ProviderI object but a [$value]");
	  return (undef);
      }
      else {
	  $self->{'_publisher'} = $value;
      }
  }
  unless (exists $self->{'_publisher'}) {
      return (undef);
  } else {
      return $self->{'_publisher'};
  }
}


=head2 add_author

 Title   : add_author
 Usage   : $self->add_author($person)
 Function: 

           Adds one L<Bio::Biblio::ProviderI> object into the list of
           authors.

 Example : 
 Returns : 1 when succeeds, 0 for failure.
 Args    : Bio::Biblio::ProviderI object

=cut

sub add_author {
  my ($self,$value) = @_;
  if (defined $value) {
      if( ! $value->isa('Bio::Biblio::ProviderI') ) {
	  my $com = ref $value;
	  $self->throw("Is not a author object but a  [$com]");
	  return 0;
      } else {
	  push(@{$self->{'_authors'}},$value);
	  return 1;
      }
  } else {
      return 0;
  }
}


=head2 each_author

 Title   : each_author
 Usage   : $obj->each_author();
 Function: 

	     Returns a list of L<Bio::Biblio::ProviderI> objects

 Example : 
 Returns : list of authors
 Args    : none

=cut

sub each_author{
   my ($self,@args) = @_;
   return @{$self->{'_authors'}};
}



=head2 add_contributor

 Title   : add_contributor
 Usage   : $self->add_contributor($person)
 Function: 

           Adds one L<Bio::Biblio::ProviderI> object into the list of
           contributors.

 Example : 
 Returns : 1 when succeeds, 0 for failure.
 Args    : Bio::Biblio::ProviderI object

=cut

sub add_contributor {
  my ($self,$value) = @_;
  if (defined $value) {
      if( ! $value->isa('Bio::Biblio::ProviderI') ) {
	  my $com = ref $value;
	  $self->throw("Is not a contributor object but a  [$com]");
	  return 0;
      } else {
	  push(@{$self->{'_contributors'}},$value);
	  return 1;
      }
  } else {
      return 0;
  }
}


=head2 each_contributor

 Title   : each_contributor
 Usage   : $obj->each_contributor();
 Function: 

	     Returns a list of L<Bio::Biblio::ProviderI> objects

 Example : 
 Returns : list of contributors
 Args    : none

=cut

sub each_contributor{
   my ($self,@args) = @_;
   return @{$self->{'_contributors'}};
}



1;
