#-----------------------------------------------------------------
# $Id$
#
# BioPerl module for Bio::Factory::BlastResultFactory
#
# Cared for by Steve Chervitz <sac@bioperl.org>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

# POD documentation - main docs before the code

=head1 NAME

Bio::Factory::BlastResultFactory - Factory for Bio::Search::Result::BlastResult objects

=head1 SYNOPSIS

    use Bio::Factory::BlastResultFactory;

    my $result_fact = Bio::Factory::BlastResultFactory->new();

    my $result = $result_fact->create_result( %parameters );

See documentation for create_result() for information about C<%parameters>.

=head1 DESCRIPTION

This module encapsulates code for creating Bio::Search::Result::BlastResult
and Bio::Search::HSP::BlastHSP objects from traditional BLAST report
data (i.e., non-XML formatted).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                - General discussion
  http://bio.perl.org/MailList.html    - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR 

Steve Chervitz <sac@bioperl.org>

See L<the FEEDBACK section | FEEDBACK> for where to send bug reports and comments.

=head1 COPYRIGHT

Copyright (c) 2001 Steve Chervitz. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

#'

package Bio::Factory::BlastResultFactory;

use strict;
use Bio::Root::Root;
use Bio::Factory::ResultFactoryI;
use Bio::Search::Result::BlastResult;

use vars qw(@ISA);

@ISA = qw(Bio::Root::Root Bio::Factory::ResultFactoryI); 

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    return $self;
}

=head2 create_result

 Title   : create_result
 Usage   : $result = $factory->create_result( %params );
 Function: Creates a new Bio::Search::Result::BlastResult object.
 Returns : A single Bio::Search::Result::BlastResult object
 Args    : none

=cut

sub create_result {
    my ($self, @args) = @_;

    my $result = Bio::Search::Result::BlastResult->new( @args );

    return $result;
}



1;
