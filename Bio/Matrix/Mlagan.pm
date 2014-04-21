#
# BioPerl module for Bio::Matrix::Mlagan
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Sendu Bala
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Matrix::Mlagan - A generic matrix with mlagan fields

=head1 SYNOPSIS

  # See L<Bio::Matrix::Generic> for most methods.
  # These are relevant for mlagan IO:
  $matrix->gap_open(-400);
  $matrix->gap_continue(-25);

=head1 DESCRIPTION

This is based on Bio::Matrix::Generic, differing by storing gap_open and
gap_continue data members to allow mlagan IO (see Bio::Matrix::IO::mlagan).
(Those values are 'outside' the matrix.)

It also limits the structure to a 6x6 matrix with row & column names 'A', 'C',
'G', 'T', '.' and 'N'.

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

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Matrix::Mlagan;
use strict;

use base qw(Bio::Matrix::Generic);


=head2 new

 Title   : new
 Usage   : my $obj = Bio::Matrix::Generic->new();
 Function: Builds a new Bio::Matrix::Generic object 
 Returns : an instance of Bio::Matrix::Generic
 Args    : -values            => arrayref of arrayrefs of data initialization
           -matrix_id         => id of the matrix
           -matrix_name       => name of the matrix
           -matrix_init_value => default value to initialize empty cells
           -gap_open          => gap open penalty (int)
           -gap_continue      => gap continue penalty (int)

           NB: -rownames and -colnames should not be given here, since they are
           always being set to 'A', 'C', 'G', 'T', '.' and 'N'.

=cut

sub new {
    my($class, @args) = @_;
    my %args = (@args, -rownames => [qw(A C G T . N)],
                       -colnames => [qw(A C G T . N)]);
    my $self = $class->SUPER::new(%args);
    
    $self->_set_from_args(\@args, -methods => [qw(gap_open gap_continue)]);
    
    return $self;
}

=head2 gap_open

 Title   : gap_open
 Usage   : $obj->gap_open(-400);
 Function: Get/set the gap open amount.
 Returns : int
 Args    : none to get, OR int to set

=cut

sub gap_open {
    my $self = shift;
    if (@_) { $self->{gap_open} = shift }
    return $self->{gap_open} || return;
}

=head2 gap_continue

 Title   : gap_continue
 Usage   : $obj->gap_continue(-25);
 Function: Get/set the gap continue amount.
 Returns : int
 Args    : none to get, OR int to set

=cut

sub gap_continue {
    my $self = shift;
    if (@_) { $self->{gap_continue} = shift }
    return $self->{gap_continue} || return;
}

=head2 add_row

 Title   : add_row
 Usage   : Do not use
 Function: This generic method is not suitable for mlagan, where the number of
           rows is fixed.
 Returns : Warning
 Args    : none

=cut

sub add_row {
    shift->warn("Mlagan matricies are fixed at 6x6");
}

=head2 remove_row

 Title   : remove_row
 Usage   : Do not use
 Function: This generic method is not suitable for mlagan, where the number of
           rows is fixed.
 Returns : Warning
 Args    : none

=cut

sub remove_row {
    shift->warn("Mlagan matricies are fixed at 6x6");
}

=head2 add_column

 Title   : add_column
 Usage   : Do not use
 Function: This generic method is not suitable for mlagan, where the number of
           columns is fixed.
 Returns : Warning
 Args    : none

=cut

sub add_column {
    shift->warn("Mlagan matricies are fixed at 6x6");
}

=head2 remove_column

 Title   : remove_column
 Usage   : Do not use
 Function: This generic method is not suitable for mlagan, where the number of
           columns is fixed.
 Returns : Warning
 Args    : none

=cut

sub remove_column {
    shift->warn("Mlagan matricies are fixed at 6x6");
}

1;
