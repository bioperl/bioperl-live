
# BioPerl module for Bio::Tools::Primer::Assessor::Base
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Primer::Assessor::Base - base class for common assessor things

=head1 SYNOPSIS

    use Bio::Tools::Primer::Assessor::Base

    $base->weight(10);



=head1 DESCRIPTION

Base class for assessors, probably only defining the weight function

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
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


use vars qw(@ISA);

use Bio::Root::Root;

package Bio::Tools::Primer::Assessor::Base;

@ISA = qw(Bio::Root::Root);

sub new {
    my ( $caller, @args) = @_;   
    my ($self) = $caller->SUPER::new(@args); 

    my ($weight) = $self->_rearrange([qw(WEIGHT)],@args);

    if( !defined $weight ) {
	$weight = 10;
    }

    $self->weight($weight);

    # done - we hope
    return $self;
}

sub weight {
    my $self   = shift;
    my $weight = shift;

    if( defined $weight ) {
	$self->{'weight'} = $weight;
    }

    return $self->{'weight'};
}

1;
