# $Id$
#
# BioPerl module for Bio::Biblio::ProviderI
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::ProviderI - Abstract interface for a author classes

=head1 SYNOPSIS

# not instantiable

=head1 DESCRIPTION

Super class and interface class for bibliographic reference providers
- "authors".

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


package Bio::Biblio::ProviderI;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;


=head2 name

 Title   : name
 Usage   : $obj->name();
 Function: 

           Sets and returns the name of the provider.
           [BiblioProvider::name]

 Example : 
 Returns : string
 Args    : string

=cut


sub name {
    my ($self,$value) = @_;
    if ( defined $value) {
	$self->{'_name'} = $value;
    }
    if ( ! exists $self->{'_name'} ) {
	return 0;
    } 
    return $self->{'_name'};
}

1;
