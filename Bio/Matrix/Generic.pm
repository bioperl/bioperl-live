#
# BioPerl module for Bio::Matrix::Generic
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Matrix::Generic - A generic matrix implementation

=head1 SYNOPSIS

  # A matrix has columns and rows 
  my $matrix = Bio::Matrix::Generic->new;
  $matrix->add_column(1,$column1);
  $matrix->add_column(2,$column2);

  my $element = $matrix->entry_by_num(1,2);
  $matrix->entry_by_num(1,2,$newval);

  my $entry = $matrix->entry('human', 'mouse');

  $matrix->entry('human','mouse', $newval);


=head1 DESCRIPTION

This is a general purpose matrix object for dealing with row+column
data which is typical when enumerating all the pairwise combinations
and desiring to get slices of the data.

Data can be accessed by column and row names or indexes.  Matrix
indexes start at 0.

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

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Matrix::Generic;
use strict;


use base qw(Bio::Root::Root Bio::Matrix::MatrixI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Matrix::Generic->new();
 Function: Builds a new Bio::Matrix::Generic object 
 Returns : an instance of Bio::Matrix::Generic
 Args    : -values     => arrayref of arrayrefs of data initialization 
           -rownames   => arrayref of row names
           -colnames   => arrayref of col names
           -matrix_id  => id of the matrix
           -matrix_name=> name of the matrix
           -matrix_init_value => default value to initialize empty cells

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($values, $rownames, $colnames,
      $id,$name,$init_val) = 
	  $self->_rearrange([qw(VALUES ROWNAMES COLNAMES 
			        MATRIX_ID MATRIX_NAME 
                                MATRIX_INIT_VALUE)],@args);
  $self->matrix_id($id) if  defined $id;
  $self->matrix_name($name) if defined $name;
  if( defined $rownames && defined $colnames ) {
      if( ref($rownames) !~ /ARRAY/i ) {
	  $self->throw("need an arrayref for the -rownames option");
      }
      # insure we copy the values
      $self->{'_rownames'} = [ @$rownames ];
      my $count = 0;
      %{$self->{'_rownamesmap'}} = map { $_ => $count++ } @$rownames; 

      if( ref($colnames) !~ /ARRAY/i ) {
	  $self->throw("need an arrayref for the -colnames option");
      }
      # insure we copy the values
      $self->{'_colnames'} = [ @$colnames ];
      $count = 0;
      %{$self->{'_colnamesmap'}} = map { $_ => $count++ } @$colnames; 

      $self->{'_values'} = [];
      if( defined $values ) {
	  if( ref($values) !~ /ARRAY/i ) {
	      $self->throw("Need an arrayref of arrayrefs (matrix) for -values option");
	  }	  
	  for my $v ( @$values ) {
	      if( ref($v) !~ /ARRAY/i ) {
		  $self->throw("Need and array of arrayrefs (matrix) for -values option");
	      }
	      push @{$self->{'_values'}}, [@$v];
	  }
      } else {
	  my @fill = ($init_val) x scalar @$colnames; # undef init_val will be default
	  for ( @$rownames ) {
	      push @{$self->{'_values'}}, [@fill];
	  }
      }
  } elsif( ! defined $rownames && ! defined $colnames && ! defined $values ) {
      $self->{'_values'}   = [];
      $self->{'_rownames'} = [];
      $self->{'_colnames'} = [];
  } else { 
      $self->throw("Must have either provided no values/colnames/rownames or provided all three");
  }

  return $self;
}


=head2 matrix_id

 Title   : matrix_id
 Usage   : my $id = $matrix->matrix_id
 Function: Get/Set the matrix ID
 Returns : scalar value
 Args    : [optional] new id value to store


=cut

sub matrix_id{
   my $self = shift;
   return $self->{'_matid'} = shift if @_;
   return $self->{'_matid'};

   
}

=head2 matrix_name

 Title   : matrix_name
 Usage   : my $name = $matrix->matrix_name();
 Function: Get/Set the matrix name
 Returns : scalar value
 Args    : [optional] new matrix name value


=cut

sub matrix_name{
   my $self = shift;
   return $self->{'_matname'} = shift if @_;
   return $self->{'_matname'};
}


=head2 entry

 Title   : entry
 Usage   : my $entry = $matrix->entry($row,$col)
 Function: Get the value for a specific cell as specified
           by the row and column names
 Returns : scalar value or undef if row or col does not
           exist
 Args    : $rowname - name of the row
           $colname - column name

=cut

sub entry{
   my ($self,$row,$column,$newvalue) = @_;
   if( ! defined $row || ! defined $column ) {
       $self->throw("Need at least 2 ids");
   }

   my ($rownum) = $self->row_num_for_name($row);
   my ($colnum) = $self->column_num_for_name($column);
   return $self->entry_by_num($rownum,$colnum,$newvalue);
}

=head2 get_entry

 Title   : get_entry
 Usage   : my $entry = $matrix->get_entry($rowname,$columname)
 Function: Get the entry for a given row,column pair
 Returns : scalar
 Args    : $row name
           $column name 


=cut

sub get_entry{ $_[0]->entry($_[1],$_[2]) }

=head2 entry_by_num

 Title   : entry_by_num
 Usage   : my $entry = $matrix->entry_by_num($rownum,$colnum)
 Function: Get an entry by row and column numbers instead of by name
           (rows and columns start at 0)
 Returns : scalar value or undef if row or column name does not
           exist
 Args    : $row - row number
           $col - column number
           [optional] $newvalue to store at this cell

=cut

sub entry_by_num {
    my ($self,$row,$col,$newvalue) = @_;
    if( ! defined $row || ! defined $col || 
	$row !~ /^\d+$/ ||
	$col !~ /^\d+$/ ) {
	$self->warn("expected to get 2 number for entry_by_num");
	return;
    }
    
    if( defined $newvalue ) {
       return $self->_values->[$row][$col] = $newvalue;
   } else { 
       return $self->_values->[$row][$col];
   }
}

sub get_element { 
    my $self = shift;
    $self->entry(@_);
}


=head2 column

 Title   : column
 Usage   : my @col = $matrix->column('ALPHA');
           OR
           $matrix->column('ALPHA', \@col);
 Function: Get/Set a particular column
 Returns : Array (in array context) or arrayref (in scalar context)
           of values.  
           For setting will warn if the new column is of a different
           length from the rest of the columns.
 Args    : name of the column
           [optional] new column to store here 

=cut

sub column{
    my ($self,$column,$newcol) = @_;

    if( ! defined $column ) {
	$self->warn("Need at least a column id");
	return;
    }
    my $colnum  = $self->column_num_for_name($column);
    if( ! defined $colnum ) { 
	$self->warn("could not find column number for $column");
	return;
    }
    return $self->column_by_num($colnum,$newcol);
}


=head2 get_column

 Title   : get_column
 Usage   : my @row = $matrix->get_column('ALPHA');
 Function: Get a particular column
 Returns : Array (in array context) or arrayref (in scalar context)
           of values
 Args    : name of the column


=cut

sub get_column { $_[0]->column($_[1]) }


=head2 column_by_num

 Title   : column_by_num
 Usage   : my @col = $matrix->column_by_num(1);
           OR
           $matrix->column_by_num(1,\@newcol);
 Function: Get/Set a column by its number instead of name
           (cols/rows start at 0)
 Returns : Array (in array context) or arrayref (in scalar context)
           of values
 Args    : name of the column
           [optional] new value to store for a particular column

=cut

sub column_by_num{
    my ($self,$colnum,$newcol) = @_;
    if( ! defined $colnum ) {
	$self->warn("need at least a column number");
	return;
    }
    my $rowcount = $self->num_rows;
    my $colcount = $self->num_columns;
    my $ret;
    
    if( defined $newcol ) {
	if( ref($newcol) !~ /ARRAY/i) {
	    $self->warn("expected a valid arrayref for resetting a column");
	    return;
	}
	if( scalar @$newcol != $rowcount ) {
	    $self->warn("new column is not the correct length ($rowcount) - call add or remove row to shrink or grow the number of rows first");
	    return;
	}
	for(my $i=0; $i < $rowcount; $i++) {
	    $self->entry_by_num($i,$colnum,$newcol->[$i]);
	}
	$ret = $newcol;
    } else { 
	$ret = [];
	for(my $i=0; $i < $rowcount; $i++) {
	    push @$ret,$self->entry_by_num($i,$colnum);
	}
    }
    if( wantarray ) { return @$ret } 
    return $ret;

}

=head2 row

 Title   : row
 Usage   : my @row = $matrix->row($rowname);
             OR
           $matrix->row($rowname,\@rowvalues);
 Function: Get/Set the row of the matrix
 Returns : Array (in array context) or arrayref (in scalar context)
 Args    : rowname
           [optional] new value of row to store


=cut

sub row {
    my ($self,$row,$newrow) = @_;
    if( ! defined $row) {
	$self->warn("Need at least a row id");
	return;
    }
    my $rownum = $self->row_num_for_name($row);
    return $self->row_by_num($rownum,$newrow);
}


=head2 get_row

 Title   : get_row
 Usage   : my @row = $matrix->get_row('ALPHA');
 Function: Get a particular row
 Returns : Array (in array context) or arrayref (in scalar context)
           of values
 Args    : name of the row

=cut

sub get_row { $_[0]->row($_[1]) }

=head2 row_by_num

 Title   : row_by_num
 Usage   : my @row = $matrix->row_by_num($rownum);
             OR
           $matrix->row($rownum,\@rowvalues);
 Function: Get/Set the row of the matrix
 Returns : Array (in array context) or arrayref (in scalar context)
 Args    : rowname
           [optional] new value of row to store

=cut

sub row_by_num{
   my ($self,$rownum,$newrow) = @_;
   if( ! defined $rownum ) {
       $self->warn("need at least a row number");
       return;
   }
    my $colcount = $self->num_columns;
    my $ret;
    if( defined $newrow ) {
	if( ref($newrow) !~ /ARRAY/i) {
	    $self->warn("expected a valid arrayref for resetting a row");
	    return;
	}
	if( scalar @$newrow != $colcount ) {
	    $self->warn("new row is not the correct length ($colcount) - call add or remove column to shrink or grow the number of columns first");
	    return;
	}
	for(my $i=0; $i < $colcount; $i++) {
	    $self->entry_by_num($rownum,$i, $newrow->[$i]);
	}
	$ret = $newrow;
    } else { 
	$ret = [];
	for(my $i=0; $i < $colcount; $i++) {
	    # we're doing this to explicitly 
	    # copy the entire row
	    push @$ret, $self->entry_by_num($rownum,$i);
	}
    }
    if( wantarray ) { return @$ret } 
    return $ret;


}


=head2 diagonal

 Title   : diagonal
 Usage   : my @diagonal = $matrix->get_diagonal()
 Function: Get the diagonal of a matrix
 Returns : Array (in array context) or arrayref (in scalar context)
           of values which lie along the diagonal
 Args    : none


=cut

sub get_diagonal{
   my ($self) = @_;
   my @diag;
   my $rowcount = $self->num_rows;
   my $colcount = $self->num_columns;
   for(my $i = 0; $i < $rowcount; $i++ ) {
       push @diag, $self->entry_by_num($i,$i);
   }
   return @diag;
}


=head2 add_row

 Title   : add_row
 Usage   : $matrix->add_row($index,\@newrow);
 Function: Adds a row at particular location in the matrix.
           If $index < the rowcount will shift all the rows down
           by the number of new rows.
           To add a single empty row, simply call
           $matrix->add_row($index,undef);
 Returns : the updated number of total rows in the matrix
 Args    : index to store
           name of the row (header)
           newrow to add, if this is undef will add a single
                     row with all values set to undef 

=cut

sub add_row{
   my ($self,$index,$name,$newrow) = @_;
   if( !defined $index || 
       $index !~ /^\d+$/ ) {
       $self->warn("expected a valid row index in add_row");
       return;
   } elsif( ! defined $name) {
       $self->warn("Need a row name or heading");
       return;
   } elsif( defined $self->row_num_for_name($name) ) {
       $self->warn("Need a unqiue name for the column heading, $name is already used");
       return;
   }
   my $colcount = $self->num_columns;
   my $rowcount = $self->num_rows;

   if( $index >  $rowcount ) { 
       $self->warn("cannot add a row beyond 1+last row at the end ($rowcount) not $index - adding at $rowcount instead");
       $index = $rowcount;
   }

   if( ! defined $newrow ) {
       $newrow = [];
       $newrow->[$colcount] = undef;
   } elsif( ref($newrow) !~ /ARRAY/i ) {
       $self->throw("Expected either undef or a valid arrayref for add_row");
   }
   # add this row to the matrix by carving out space for it with 
   # splice
   splice(@{$self->{'_values'}}, $index,0,[]);
   for( my $i = 0; $i < $colcount; $i++ ) {
       $self->entry_by_num($index,$i,$newrow->[$i]);
   }
   splice(@{$self->{'_rownames'}}, $index,0,$name);
   # Sadly we have to remap these each time (except for the case
   # when we're adding a new column to the end, but I don't think
   # the speedup for that case warrants the extra code at this time.
   my $ct = 0;
   %{$self->{'_rownamesmap'}} = map { $_ => $ct++} @{$self->{'_rownames'}};
   return $self->num_rows;
}

=head2 remove_row

 Title   : remove_row
 Usage   : $matrix->remove_row($colnum)
 Function: remove a row from the matrix shifting all the rows
           up by one
 Returns : Updated number of rows in the matrix
 Args    : row index


=cut

sub remove_row{
   my ($self,$rowindex) = @_;
   my $rowcount = $self->num_rows;
   
   if( $rowindex > $rowcount ) {
       $self->warn("rowindex $rowindex is greater than number of rows $rowcount, cannot process");
       return 0;
   } else { 
       splice(@{$self->_values},$rowindex,1);
       delete $self->{'_rownamesmap'}->{$self->{'_rownames'}->[$rowindex]};
       splice(@{$self->{'_rownames'}},$rowindex,1);
   }
   my $ct = 0;
   %{$self->{'_rownamesmap'}} = map { $_ => $ct++} @{$self->{'_rownames'}};
   return $self->num_rows;
}

=head2 add_column

 Title   : add_column
 Usage   : $matrix->add_column($index,$colname,\@newcol);
 Function: Adds a column at particular location in the matrix.
           If $index < the colcount will shift all the columns right
           by the number of new columns.
           To add a single empty column, simply call
           $matrix->add_column($index,undef);
 Returns : the updated number of total columns in the matrix
 Args    : index to store
           name of the column (header)
           newcolumn to add, if this is undef will add a single
                 column with all values set to undef 


=cut


sub add_column{
   my ($self,$index,$name,$newcol) = @_;
   if( !defined $index ||
       $index !~ /^\d+$/ ) {
       $self->warn("expected a valid col index in add_column");
       return;
   } elsif( ! defined $name) {
       $self->warn("Need a column name or heading");
       return;
   } elsif( defined $self->column_num_for_name($name) ) {
       $self->warn("Need a unqiue name for the column heading, $name is already used");
       return;
   }
   my $colcount = $self->num_columns;
   my $rowcount = $self->num_rows;
   if( $index > $colcount ) { 
       $self->warn("cannot add a column beyond 1+last column at the end ($colcount) not $index - adding at $colcount instead");
       $index = $colcount;
   }

   if( ! defined $newcol ) {
       $newcol = [];
       $newcol->[$rowcount] = undef; # make the array '$rowcount' long
   } elsif( ref($newcol) !~ /ARRAY/i ) {
       $self->throw("Expected either undef or a valid arrayref for add_row");
   }
   for( my $i = 0; $i < $rowcount; $i++ ) {
       # add this column to each row
       splice(@{$self->_values->[$i]},$index,0,[]);
       $self->entry_by_num($i,$index,$newcol->[$i]);
   }
   splice(@{$self->{'_colnames'}}, $index,0,$name);
   # Sadly we have to remap these each time (except for the case
   # when we're adding a new column to the end, but I don't think
   # the speedup for that case warrants the extra code at this time.
   my $ct = 0;
   %{$self->{'_colnamesmap'}} = map {$_ => $ct++} @{$self->{'_colnames'}};
   return $self->num_columns;
}

=head2 remove_column

 Title   : remove_column
 Usage   : $matrix->remove_column($colnum)
 Function: remove a column from the matrix shifting all the columns
           to the left by one
 Returns : Updated number of columns in the matrix
 Args    : column index

=cut

sub remove_column{
   my ($self,$colindex) = @_;

   my $colcount = $self->num_columns;
   my $rowcount = $self->num_rows;
   if( $colindex > $colcount ) {
		$self->warn("colindex $colindex is greater than number of columns ($colcount), cannot process");
		return 0;
   } else { 
		for(my $i = 0; $i < $rowcount; $i++ ) {
			splice(@{$self->_values->[$i]},$colindex,1);
		}
		delete $self->{'_colnamesmap'}->{$self->{'_colnames'}->[$colindex]};
		splice(@{$self->{'_colnames'}},$colindex,1);
   }
   my $ct = 0;
   %{$self->{'_colnamesmap'}} = map {$_ => $ct++} @{$self->{'_colnames'}};
   return $self->num_columns;
}

=head2 column_num_for_name

 Title   : column_num_for_name
 Usage   : my $num = $matrix->column_num_for_name($name)
 Function: Gets the column number for a particular column name
 Returns : integer
 Args    : string


=cut

sub column_num_for_name{
   my ($self,$name) = @_;
   
   return $self->{'_colnamesmap'}->{$name};
}

=head2 row_num_for_name

 Title   : row_num_for_name
 Usage   : my $num = $matrix->row_num_for_name
 Function: Gets the row number for a particular row name
 Returns : integer
 Args    : string


=cut

sub row_num_for_name{
   my ($self,$name) = @_;
   return $self->{'_rownamesmap'}->{$name}
}


=head2 column_header

 Title   : column_header
 Usage   : my $name = $matrix->column_header(0)
 Function: Gets the column header for a particular column number
 Returns : string
 Args    : integer


=cut

sub column_header{
   my ($self,$num) = @_;
   return $self->{'_colnames'}->[$num];
}


=head2 row_header

 Title   : row_header
 Usage   : my $name = $matrix->row_header(0)
 Function: Gets the row header for a particular row number
 Returns : string
 Args    : integer


=cut

sub row_header{
   my ($self,$num) = @_;
   return $self->{'_rownames'}->[$num];
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
   return scalar @{$self->_values};
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
   return scalar @{$self->_values->[0] || []};
}


=head2 row_names

 Title   : row_names
 Usage   : my @rows = $matrix->row_names
 Function: The names of all the rows
 Returns : array in array context, arrayref in scalar context
 Args    : none


=cut

sub row_names{
   if( wantarray ) { 
       return @{shift->{'_rownames'}};
   } else { 
       return shift->{'_rownames'};
   }
}


=head2 column_names

 Title   : column_names
 Usage   : my @columns = $matrix->column_names
 Function: The names of all the columns
 Returns : array in array context, arrayref in scalar context
 Args    : none


=cut

sub column_names{
   if( wantarray ) { 
       return @{shift->{'_colnames'}};
   } else { 
       return shift->{'_colnames'};
   }
}

=head2 private methods

Private methods for a Generic Matrix

=head2 _values

 Title   : _values
 Usage   : $matrix->_values();
 Function: get/set for array ref of the matrix containing
           distance values 
 Returns : an array reference 
 Args    : an array reference


=cut

sub _values{
   my ($self,$val) = @_;
   if( $val ){
       $self->{'_values'} = $val;
   }
   return $self->{'_values'};
}

1;
