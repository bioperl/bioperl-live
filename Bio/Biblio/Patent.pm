# $Id$
#
# BioPerl module for Bio::Biblio::Patent
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::Patent - A type of an author

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


package Bio::Biblio::Patent;
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
	$toc_type, $publisher, $doc_number, $doc_office, $doc_type, 
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
				  DOC_NUMBER
				  DOC_OFFICE
				  DOC_TYPE    
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

    $doc_number && $self->doc_number($doc_number);
    $doc_office && $self->doc_office($doc_office);
    $doc_type && $self->doc_type($doc_type);

    return $self; # success - we hope!

}


=head2 doc_number

 Title   : doc_number
 Usage   : $obj->doc_number();
 Function: 

           Sets and returns the Web DOC_NUMBER
           [BiblioPatent::docNumber]

 Example : 
 Returns : string
 Args    : string

=cut


sub doc_number {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_doc_number'} = $value;
    }
    if ( ! exists $self->{'_doc_number'} ) {
	return 0;
    } 
    return $self->{'_doc_number'};
}


=head2 doc_office

 Title   : doc_office
 Usage   : $obj->doc_office();
 Function: 

           Sets and returns the estimated doc_office of the resource in characters.
           [BiblioPatent::docOffice]

 Example : 
 Returns : integer
 Args    : integer

=cut


sub doc_office {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_doc_office'} = $value;
    }
    if ( ! exists $self->{'_doc_office'} ) {
	return 0;
    } 
    return $self->{'_doc_office'};
}

=head2 doc_type

 Title   : doc_type
 Usage   : $obj->doc_type();
 Function: 

           Sets and returns the name of the reference.
           [BibRefPatent::docType]

 Example : 
 Returns : string
 Args    : string

=cut


sub doc_type {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_doc_type'} = $value;
    }
    if ( ! exists $self->{'_doc_type'} ) {
	return 0;
    } 
    return $self->{'_doc_type'};
}


=head2 add_applicant

 Title   : add_applicant
 Usage   : $self->add_applicant($person)
 Function: 

           Adds one Bio::Biblio::ProviderI object into the list of
           applicants.

 Example : 
 Returns : 1 when succeeds, 0 for failure.
 Args    : Bio::Biblio::ProviderI object

See L<Bio::Biblio::Provider> for more information.

=cut

sub add_applicant {
  my ($self,$value) = @_;
  if (defined $value) {
      if( ! $value->isa('Bio::Biblio::ProviderI') ) {
	  my $com = ref $value;
	  $self->throw("Is not a applicant object but a  [$com]");
	  return 0;
      } else {
	  push(@{$self->{'_applicants'}},$value);
	  return 1;
      }
  } else {
      return 0;
  }
}


=head2 each_applicant

 Title   : each_applicant
 Usage   : $obj->each_applicant();
 Function: 

	     Returns a list of L<Bio::Biblio::ProviderI> objects

 Example : 
 Returns : list of applicants
 Args    : none

=cut

sub each_applicant{
   my ($self,@args) = @_;
   return @{$self->{'_applicants'}};
}



1;
