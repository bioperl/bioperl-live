# $Id$
#
# BioPerl module for Bio::Biblio::WebResource
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::WebResource - A type of an author

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


package Bio::Biblio::WebResource;
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
	$toc_type, $publisher, $url, $size, $cost
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
				  URL
				  SIZE
				  COST    
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

    $url && $self->url($url);
    $size && $self->size($size);
    $cost && $self->cost($cost);

    return $self; # success - we hope!

}


=head2 url

 Title   : url
 Usage   : $obj->url();
 Function: 

           Sets and returns the Web URL
           [BiblioWebResource::url]

 Example : 
 Returns : string
 Args    : string

=cut


sub url {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_url'} = $value;
    }
    if ( ! exists $self->{'_url'} ) {
	return 0;
    } 
    return $self->{'_url'};
}


=head2 size

 Title   : size
 Usage   : $obj->size();
 Function: 

           Sets and returns the estimated size of the resource in characters.
           [BiblioWebResource::estimatedSize]

 Example : 
 Returns : integer
 Args    : integer

=cut


sub size {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_size'} = $value;
    }
    if ( ! exists $self->{'_size'} ) {
	return 0;
    } 
    return $self->{'_size'};
}

=head2 cost

 Title   : cost
 Usage   : $obj->cost();
 Function: 

           Sets and returns the name of the reference.
           [BibRefWebResource::cost]

 Example : 
 Returns : string
 Args    : string

=cut


sub cost {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_cost'} = $value;
    }
    if ( ! exists $self->{'_cost'} ) {
	return 0;
    } 
    return $self->{'_cost'};
}



1;
