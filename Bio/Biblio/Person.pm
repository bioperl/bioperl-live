# $Id$
#
# BioPerl module for Bio::Biblio::Person
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::Person - A type of an author

=head1 SYNOPSIS



=head1 DESCRIPTION

Person as an author.

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


package Bio::Biblio::Person;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Biblio::ProviderI;

@ISA = qw( Bio::Biblio::ProviderI Bio::Root::Root );

sub new {
    my($class,@args) = @_;
    my $self;
    $self = {};
    bless $self, $class;

    my ($name, $surname, $first_name, $mid_initials, $email, $address, $affiliation) =
	$self->_rearrange([qw(NAME
			      SURNAME
			      FIRST_NAME
			      MID_INITIALS
			      EMAIL
			      ADDRESS
			      AFFILIATION
				  )],
			      @args);

    $name && $self->name($name);
    $surname && $self->surname($surname);
    $first_name && $self->first_name($first_name);
    $mid_initials && $self->mid_initials($mid_initials);
    $email && $self->email($email);
    $address && $self->address($address);
    $affiliation && $self->affiliation($affiliation);
    
    return $self; # success - we hope!

}

=head2 name

 Title   : name
 Usage   : $obj->name();
 Function: 

           Sets and returns the surname of the author.
           [BiblioPerson::surname]

 Example : 
 Returns : string
 Args    : string

=cut


=head2 surname

 Title   : surname
 Usage   : $obj->surname();
 Function: 

           Sets and returns the surname of the author.
           Calls $self->name().
           [BiblioPerson::surname]

 Example : 
 Returns : string
 Args    : string

=cut


sub surname {
    my ($self,$value) = @_;
    $self->name($value);
}


=head2 first_name

 Title   : first_name
 Usage   : $obj->first_name();
 Function: 

           Sets and returns the first name of the provider.
           [BiblioPerson::firstName]

 Example : 
 Returns : string
 Args    : string

=cut


sub first_name {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_first_name'} = $value;
    }
    if ( ! exists $self->{'_first_name'} ) {
	return 0;
    } 
    return $self->{'_first_name'};
}


=head2 mid_initials

 Title   : mid_initials
 Usage   : $obj->mid_initials();
 Function: 

           Sets and returns the mid initials of the provider.
           [BiblioPerson::midInitials]

 Example : 
 Returns : string
 Args    : string

=cut


sub mid_initials {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_mid_initials'} = $value;
    }
    if ( ! exists $self->{'_mid_initials'} ) {
	return 0;
    } 
    return $self->{'_mid_initials'};
}

=head2 email

 Title   : email
 Usage   : $obj->email();
 Function: 

           Sets and returns the email address of the provider.
           [BiblioPerson::email]

 Example : 
 Returns : string
 Args    : string

=cut


sub email {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_email'} = $value;
    }
    if ( ! exists $self->{'_email'} ) {
	return 0;
    } 
    return $self->{'_email'};
}


=head2 address

 Title   : address
 Usage   : $obj->address();
 Function: 

           Sets and returns the mid initials of the provider.
           [BiblioPerson::postalAddress]

 Example : 
 Returns : string
 Args    : string

=cut


sub address {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_address'} = $value;
    }
    if ( ! exists $self->{'_address'} ) {
	return 0;
    } 
    return $self->{'_address'};
}


=head2 affiliation

 Title   : affiliation
 Usage   : $obj->affiliation();
 Function: 

           Sets and returns the affiliation of the provider.
           [BiblioPerson::affiliation]

 Example : 
 Returns : string
 Args    : string

=cut


sub affiliation {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_affiliation'} = $value;
    }
    if ( ! exists $self->{'_affiliation'} ) {
	return 0;
    } 
    return $self->{'_affiliation'};
}

1;
