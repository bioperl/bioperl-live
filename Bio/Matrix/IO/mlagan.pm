#
# BioPerl module for Bio::Matrix::IO::mlagan
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

Bio::Matrix::IO::mlagan - A parser for the mlagan substitution matrix

=head1 SYNOPSIS

  use Bio::Matrix::IO;
  my $parser = Bio::Matrix::IO->new(-format => 'mlagan',
                                   -file   => 'nucmatrix.txt');
  my $matrix = $parser->next_matrix;
  my $gap_open = $parser->gap_open;
  my $gap_continue = $parser->gap_continue;

=head1 DESCRIPTION

Use to read in and write out substitution matrix files suitable for use by
mlagan.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Matrix::IO::mlagan;
use strict;

use Bio::Matrix::Mlagan;
use base qw(Bio::Matrix::IO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Matrix::IO::mlagan->new();
 Function: Builds a new Bio::Matrix::IO::mlagan object 
 Returns : an instance of Bio::Matrix::IO::mlagan
 Args    :

=cut

=head2 next_matrix

 Title   : next_matrix
 Usage   : my $matrix = $obj->next_matrix();
 Function: parses a matrix file
 Returns : L<Bio::Matrix::Mlagan>
 Args    : none

=cut

sub next_matrix {
    my $self = shift;
    
    my (@matrix, $gap_open, $gap_cont);
    while (defined ($_ = $self->_readline)) {
        if (/^[ACGTN\.]/) {
            my (undef, @values) = split;
            push(@matrix, \@values);
        }
        elsif (/^[-\d]/) {
            ($gap_open, $gap_cont) = split;
            last;
        }
    }
    
    @matrix == 6 || $self->throw("Something wrong with file, was it the correct format?");
    
    my $matrix = Bio::Matrix::Mlagan->new(-values => \@matrix,
                                          -gap_open => $gap_open,
                                          -gap_continue => $gap_cont);
    
    return $matrix;
}

=head2 write_matrix

 Title   : write_matrix
 Usage   : $obj->write_matrix($matrix)
 Function: Write out a matrix in mlagan format
 Returns : n/a
 Args    : L<Bio::Matrix::Generic>

=cut

sub write_matrix {
    my ($self, $matrix) = @_;
    $matrix || $self->throw("Matrix required as input");
    my $gap_open = $matrix->gap_open;
    my $gap_continue = $matrix->gap_continue;
    
    unless (defined $gap_open && defined $gap_continue) {
        $self->throw("gap_open() and gap_continue() in the supplied matrix object must both be set");
    }
    
    $self->_print("     A    C    G    T    .    N\n");
    
    foreach my $char (qw(A C G T . N)) {
        my @row = $matrix->get_row($char);
        my $row = $char;
        foreach my $val (@row) {
            $row .= " " x (5 - length($val)) . $val;
        }
        
        $self->_print($row."\n");
    }
    
    $self->_print("\n$gap_open $gap_continue");
    
    return;
}

1;
