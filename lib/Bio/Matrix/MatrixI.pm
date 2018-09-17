# $Id $
#
# BioPerl module for Bio::Matrix::MatrixI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Matrix::MatrixI - An interface for describing a Matrix 

=head1 SYNOPSIS

  # Get a Matrix object

=head1 DESCRIPTION

This is an interface describing how one should be able to interact
with a matrix.  One can have a lot of information I suppose and this
outline won't really work for PWM or PSSMs.  We will have to derive a
particular interface for those.

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
email or the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Matrix::MatrixI;
use strict;


use base qw(Bio::Root::RootI);


=head2 matrix_id

 Title   : matrix_id
 Usage   : my $id = $matrix->matrix_id
 Function: Get the matrix ID
 Returns : string value
 Args    : 


=cut

sub matrix_id{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 matrix_name

 Title   : matrix_name
 Usage   : my $name = $matrix->matrix_name();
 Function: Get the matrix name
 Returns : string value
 Args    :


=cut

sub matrix_name{
   my ($self) = @_;

   $self->throw_not_implemented();
}

=head2 get_entry

 Title   : get_entry
 Usage   : my $entry = $matrix->get_entry($rowname,$columname)
 Function: Get the entry for a given row,column pair
 Returns : scalar
 Args    : $row name
           $column name 


=cut

sub get_entry{
   my ($self) = @_;

    $self->throw_not_implemented();
}


=head2 get_column

 Title   : get_column
 Usage   : my @row = $matrix->get_column('ALPHA');
 Function: Get a particular column
 Returns : Array (in array context) or arrayref (in scalar context)
           of values
 Args    : name of the column


=cut

sub get_column{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 get_row

 Title   : get_row
 Usage   : my @row = $matrix->get_row('ALPHA');
 Function: Get a particular row
 Returns : Array (in array context) or arrayref (in scalar context)
           of values
 Args    : name of the row

=cut

sub get_row{
   my ($self) = @_;
   $self->throw_not_implemented();
}


=head2 get_diagonal

 Title   : get_diagonal
 Usage   : my @diagonal = $matrix->get_diagonal; 
 Function: Get the diagonal of the matrix
 Returns : Array (in array context) or arrayref (in scalar context)
 Args    : none


=cut

sub get_diagonal{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 column_num_for_name

 Title   : column_num_for_name
 Usage   : my $num = $matrix->column_num_for_name($name)
 Function: Gets the column number for a particular column name
 Returns : integer
 Args    : string


=cut

sub column_num_for_name{
   my ($self) = @_;

    $self->throw_not_implemented();
}

=head2 row_num_for_name

 Title   : row_num_for_name
 Usage   : my $num = $matrix->row_num_for_name($name)
 Function: Gets the row number for a particular row name
 Returns : integer
 Args    : string


=cut

sub row_num_for_name{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 num_rows

 Title   : num_rows
 Usage   : my $rowcount = $matrix->num_rows;
 Function: Get the number of rows
 Returns : integer
 Args    : none


=cut

sub num_rows{
   my ($self) = @_;
   $self->throw_not_implemented();
}


=head2 num_columns

 Title   : num_columns
 Usage   : my $colcount = $matrix->num_columns
 Function: Get the number of columns
 Returns : integer
 Args    : none


=cut

sub num_columns{
   my ($self) = @_;
   $self->throw_not_implemented();
}


# inverse?
=head2 reverse

 Title   : reverse
 Usage   : my $matrix = $matrix->reverse
 Function: Get the reverse of a matrix
 Returns : 
 Args    :


=cut

sub reverse{
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 row_names

 Title   : row_names
 Usage   : my @rows = $matrix->row_names
 Function: The names of all the rows
 Returns : array in array context, arrayref in scalar context
 Args    : none


=cut

sub row_names{
   my ($self) = @_;
   $self->throw_not_implemented();
}


=head2 column_names

 Title   : column_names
 Usage   : my @columns = $matrix->column_names
 Function: The names of all the columns
 Returns : array in array context, arrayref in scalar context
 Args    : none


=cut

sub column_names{
   my ($self) = @_;
   $self->throw_not_implemented();
}

1;
