#
# BioPerl module for Bio::Search::Result::INFERNALResult.pm
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Paul Cantalupo
#
# Copyright Paul Cantalupo
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Result::INFERNALResult - A Result object for INFERNAL results

=head1 SYNOPSIS

    # typically one gets Results from a SearchIO stream
    use Bio::SearchIO;
    my $io = Bio::SearchIO->new(-format => 'infernal',
                                -file   => 't/data/cmsearch_output.txt');
    while( my $result = $io->next_result ) {
        while( my $hit = $result->next_hit ) {
          print join(" ", $result->query_name, $result->algorithm, $result->num_hits), "\n";
        }
    }

=head1 DESCRIPTION

This object is a specialization of L<Bio::Search::Result::GenericResult>. There
is one extra method called L<cm_name>.

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

=head1 AUTHOR - Paul Cantalupo

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Search::Result::INFERNALResult;
use strict;
use warnings;

use base qw(Bio::Search::Result::GenericResult);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::Result::INFERNALResult->new();
 Function: Builds a new Bio::Search::Result::INFERNALResult object
 Returns : Bio::Search::Result::INFERNALResult
 Args    : -cm_name    => string, name of covariance model (CM) file.
           plus Bio::Search::Result::GenericResult parameters

=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($cm) = $self->_rearrange([qw(CM_NAME)], @args);
  if (defined $cm) { $self->cm_name($cm) }

  return $self;
}

=head2 cm_name

 Title   : cm_name
 Usage   : $obj->cm_name($newvalue)
 Function: Get/Set value of the covariance model file name (cm_name)
 Returns : value of cm_name
 Args    : newvalue (optional)

=cut

sub cm_name {
  my ($self, $value) = @_;
  if (defined $value) {
    $self->{'_cm_name'} = $value;
  }
  return $self->{'_cm_name'};
}

1;
