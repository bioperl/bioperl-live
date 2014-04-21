#
# BioPerl module for Bio::Tools::Primer::Assessor::Base
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Ewan Birney

Email birney-at-ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...




package Bio::Tools::Primer::Assessor::Base;

use base qw(Bio::Root::Root);

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
