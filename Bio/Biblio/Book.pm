# $Id$
#
# BioPerl module for Bio::Biblio::Book
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::Book - A type of an author

=head1 SYNOPSIS

#

=head1 DESCRIPTION

#


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


package Bio::Biblio::Book;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Biblio::RefI;

@ISA = qw( Bio::Biblio::RefI Bio::Root::Root );

sub new {
    my($class,@args) = @_;
    my $self;
    $self = {};
    bless $self, $class;

    my ($id, $type, $title, $rights, $language, $format, $date, 
	$spatial_location, $temporal_period, $last_modified, 
	$repository_subset, $abstract, $abstract_type, $toc, 
	$toc_type, $publisher, $isbn, $volume, $edition, $series, $editor
	) =
	    $self->_rearrange([qw(ID
				  TYPE               
				  TITLE              
				  RIGHTS             
				  LANGUAGE           
				  FORMAT             
				  DATE               
				  SPATIAL_LOCATION   
				  TEMPORAL_PERIOD    
				  LAST_MODIFIED      
				  REPOSITORY_SUBSET  
				  ABSTRACT           
				  ABSTRACT_TYPE      
				  TOC 
				  TOC_TYPE           
				  PUBLISHER      
				  ISBN
				  VOLUME
				  EDITION  
				  SERIES
				  )],
			      @args);

    $id && $self->id($id);
    $type && $self->type($type);
    $title && $self->title($title);
    $rights && $self->rights($rights);
    $language && $self->language($language);
    $format && $self->format($format);
    $date && $self->date($date);
    $spatial_location && $self->spatial_location($spatial_location);
    $temporal_period && $self->temporal_period($temporal_period);
    $last_modified && $self->last_modified($last_modified);
    $repository_subset && $self->repository_subset($repository_subset);
    $abstract && $self->abstract($abstract);
    $abstract_type && $self->abstract_type($abstract_type);
    $toc && $self->toc($toc);
    $toc_type && $self->toc_type($toc_type);
    $publisher && $self->publisher($publisher);

    $isbn && $self->isbn($isbn);
    $volume && $self->volume($volume);
    $edition && $self->edition($edition);
    $series && $self->series($series);
    $editor && $self->editor($editor);

    return $self; # success - we hope!

}


=head2 isbn

 Title   : isbn
 Usage   : $obj->isbn();
 Function: 

           Sets and returns the Web ISBN
           [BiblioBook::docNumber]

 Example : 
 Returns : string
 Args    : string

=cut


sub isbn {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_isbn'} = $value;
    }
    if ( ! exists $self->{'_isbn'} ) {
	return 0;
    } 
    return $self->{'_isbn'};
}


=head2 volume

 Title   : volume
 Usage   : $obj->volume();
 Function: 

           Sets and returns the estimated volume of the resource in characters.
           [BiblioBook::docOffice]

 Example : 
 Returns : integer
 Args    : integer

=cut


sub volume {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_volume'} = $value;
    }
    if ( ! exists $self->{'_volume'} ) {
	return 0;
    } 
    return $self->{'_volume'};
}

=head2 edition

 Title   : edition
 Usage   : $obj->edition();
 Function: 

           Sets and returns the name of the reference.
           [BibRefBook::docType]

 Example : 
 Returns : string
 Args    : string

=cut


sub edition {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_edition'} = $value;
    }
    if ( ! exists $self->{'_edition'} ) {
	return 0;
    } 
    return $self->{'_edition'};
}


=head2 series

 Title   : series
 Usage   : $obj->series();
 Function: 

           Sets and returns the name of the reference.
           [BibRefPatent::docType]

 Example : 
 Returns : string
 Args    : string

=cut


sub series {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_series'} = $value;
    }
    if ( ! exists $self->{'_series'} ) {
	return 0;
    } 
    return $self->{'_series'};
}


=head2 editor

 Title   : editor
 Usage   : $editor = $obj->editor();
 Function: Returns or sets the reference to a Bio::Biblio::Provider object.
           If there is no link, it will return undef
 Returns : an obj_ref or undef
 Args    : Bio::Biblio::Provider object

See L<Bio::Biblio::Provider> for more information.

=cut

sub editor {
  my ($self,$value) = @_;
  if (defined $value) {
      if( ! $value->isa('Bio::Biblio::ProviderI') ) {
	  $self->throw("Is not a Bio::Biblio::ProviderI object but a [$value]");
	  return (undef);
      }
      else {
	  $self->{'_editor'} = $value;
      }
  }
  unless (exists $self->{'_editor'}) {
      return (undef);
  } else {
      return $self->{'_editor'};
  }
}



1;
