
# BioPerl module for Bio::Tools::Primer::AssessorI
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Primer::AssessorI - interface for assessing primer pairs

=head1 SYNOPSIS

    use Bio::Tools::Primer::AssessorI;

    if( $obj->isa('Bio::Tools::Primer::AssessorI') ) {
	my $score = $obj->assess($primer_pair);
    }


=head1 DESCRIPTION

The Primer Assessor interface provides a interface for scoring
functions of primer pairs to comply to. It is mainly used by
Bio::Tools::Primer::Design module

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

package Bio::Tools::Primer::AssessorI.pm;

use Bio::Root::RootI;

use vars qw(@ISA);

@ISA = qw(Bio::Root::RootI);

sub assess {
    my ($self) = shift;
    $self->throw_not_implemented();
}

1;
