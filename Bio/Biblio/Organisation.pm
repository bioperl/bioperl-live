# $Id$
#
# BioPerl module for Bio::Biblio::Organisation
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::Organisation - A type of an author

=head1 SYNOPSIS



=head1 DESCRIPTION

Organisation as an author.

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


package Bio::Biblio::Organisation;
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

    my ($name) =
	$self->_rearrange([qw(NAME
				  )],
			      @args);

    $name && $self->name($name);
    
    return $self; # success - we hope!

}
