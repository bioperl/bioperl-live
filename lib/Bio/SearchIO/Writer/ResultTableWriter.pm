
=head1 NAME

Bio::SearchIO::Writer::ResultTableWriter - Outputs tab-delimited data for each Bio::Search::Result::ResultI object.

=head1 SYNOPSIS

=head2 Example 1: Using the default columns

    use Bio::SearchIO;
    use Bio::SearchIO::Writer::ResultTableWriter;

    my $in = Bio::SearchIO->new();

    my $writer = Bio::SearchIO::Writer::ResultTableWriter->new();

    my $out = Bio::SearchIO->new( -writer => $writer );

    while ( my $result = $in->next_result() ) {
        $out->write_result($result, ($in->report_count - 1 ? 0 : 1) );
    }

=head2 Example 2: Specifying a subset of columns 

    use Bio::SearchIO;
    use Bio::SearchIO::Writer::ResultTableWriter;

    my $in = Bio::SearchIO->new();

    my $writer = Bio::SearchIO::Writer::ResultTableWriter->new( 
                                  -columns => [qw(
                                                  query_name
                                                  query_length
                                                  num_hits
                                                  )]  );

    my $out = Bio::SearchIO->new( -writer => $writer,
				  -file   => ">result.out" );

    while ( my $result = $in->next_result() ) {
        $out->write_result($result, ($in->report_count - 1 ? 0 : 1) );
    }

=head2 Custom Labels

You can also specify different column labels if you don't want to use
the defaults.  Do this by specifying a C<-labels> hash reference
parameter when creating the ResultTableWriter object.  The keys of the
hash should be the column number (left-most column = 1) for the label(s)
you want to specify. Here's an example:

    my $writer = Bio::SearchIO::Writer::ResultTableWriter->new( 
                               -columns => [qw( query_name 
                                                query_length
                                                query_description 
                                                num_hits)],
                               -labels  => { 1 => 'QUERY_GI',
  	                                     2 => 'QUERY_LENGTH' } );


=head1 DESCRIPTION

Bio::SearchIO::Writer::ResultTableWriter outputs data in tab-delimited
format for each search result, one row per search result. This is a very
coarse-grain level of information since it only includes data
stored in the Bio::Search::Result::ResultI object itself and does not
include any information about hits or HSPs.

You most likely will never use this object but instead will use one of
its subclasses: Bio::SearchIO::Writer::HitTableWriter or
Bio::SearchIO::Writer::HSPTableWriter.

=head2 Available Columns

Here are the columns that can be specified in the C<-columns>
parameter when creating a ResultTableWriter object.  If a C<-columns> parameter
is not specified, this list, in this order, will be used as the default.

    query_name
    query_length
    query_description

For more details about these columns, see the documentation for the
corresponding method in L<Bio::Search::Result::ResultI|Bio::Search::Result::ResultI>.

=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other
Bioperl modules.  Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues           

=head1 AUTHOR 

Steve Chervitz E<lt>sac@bioperl.orgE<gt>

See L<the FEEDBACK section | FEEDBACK> for where to send bug reports
and comments.

=head1 COPYRIGHT

Copyright (c) 2001 Steve Chervitz. All Rights Reserved.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 SEE ALSO

L<Bio::SearchIO::Writer::HitTableWriter>,
L<Bio::SearchIO::Writer::HSPTableWriter>

=head1 METHODS

=cut


package Bio::SearchIO::Writer::ResultTableWriter;

use strict;

use base qw(Bio::Root::Root Bio::SearchIO::SearchWriterI);

# Array fields: column, object, method[/argument], printf format, column label
# Methods are defined in Bio::Search::Result::ResultI.
# Tech note: If a bogus method is supplied, it will result in all values to be zero.
#            Don't know why this is.
my %column_map = (
                  'query_name'        => ['1', 'result', 'query_name', 's', 'QUERY' ],
                  'query_length'      => ['2', 'result', 'query_length', 'd', 'LEN_Q'],
                  'query_description' => ['3', 'result', 'query_description', 's', 'DESC_Q'],
                  'num_hits'          => ['4', 'result', 'num_hits', 'd', 'NUM_HITS'],
                 );

sub column_map { return %column_map }

sub new {
    my ($class, @args) = @_; 
    my $self = $class->SUPER::new(@args);

    my( $col_spec, $label_spec,
	$filters ) = $self->_rearrange( [qw(COLUMNS 
					    LABELS
					    FILTERS)], @args);
    
    $self->_set_cols( $col_spec );
    $self->_set_labels( $label_spec ) if $label_spec;
    $self->_set_printf_fmt();
    $self->_set_row_data_func();
    $self->_set_column_labels();
    
    if( defined $filters ) {
	if( !ref($filters) =~ /HASH/i ) { 
	    $self->warn("Did not provide a hashref for the FILTERS option, ignoring.");
	} else { 
	    while( my ($type,$code) = each %{$filters} ) {
		$self->filter($type,$code);
	    }
	}
    }


    return $self;
}


# Purpose : Stores the column spec internally. Also performs QC on the 
#           user-supplied column specification.
#
sub _set_cols {
    my ($self, $col_spec_ref) = @_;
    return if defined $self->{'_cols'};  # only set columns once

    my %map = $self->column_map;

    if( not defined $col_spec_ref) {
        print STDERR "\nUsing default column map.\n";
	$col_spec_ref = [ map { $_ } sort { $map{$a}->[0] <=> $map{$b}->[0] } keys %map ];
    }

    if( ref($col_spec_ref) eq 'ARRAY') {
        # printf "%d columns to process\n", scalar(@$col_spec_ref);
        my @col_spec = @{$col_spec_ref};
        while( my $item = shift @col_spec ) {
            $item = lc($item);
            if( not defined ($map{$item}) ) {
                $self->throw(-class =>'Bio::Root::BadParameter',
                             -text => "Unknown column name: $item"
                            );
            }
            push @{$self->{'_cols'}}, $item;
            #print "pushing on to col $col_num, $inner: $item\n";
        }
    }
    else {
        $self->throw(-class =>'Bio::Root::BadParameter',
                     -text => "Can't set columns: not a ARRAY ref",
                     -value => $col_spec_ref
                    );
    }
}

sub _set_printf_fmt {
    my ($self) = @_;

    my @cols = $self->columns();
    my %map = $self->column_map;

    my $printf_fmt = '';

    foreach my $col ( @cols ) {
	$printf_fmt .= "\%$map{$col}->[3]\t";
    }

    $printf_fmt =~ s/\\t$//;

    $self->{'_printf_fmt'} = $printf_fmt;
}

sub printf_fmt { shift->{'_printf_fmt'} }

# Sets the data to be used for the labels.
sub _set_labels {
    my ($self, $label_spec) = @_;
    if( ref($label_spec) eq 'HASH') {
        foreach my $col ( sort { $a <=> $b } keys %$label_spec ) {
#            print "LABEL: $col $label_spec->{$col}\n";
            $self->{'_custom_labels'}->{$col} = $label_spec->{$col};
        }
    }
    else {
        $self->throw(-class =>'Bio::Root::BadParameter',
                     -text => "Can't set labels: not a HASH ref: $label_spec"
                    );
    }
}

sub _set_column_labels {
    my $self = shift;

    my @cols = $self->columns;
    my %map = $self->column_map;
    my $printf_fmt = '';
    my (@data, $label, @underbars);

    my $i = 0;
    foreach my $col( @cols ) {
	$i++;
        $printf_fmt .= "\%s\t";

        if(defined $self->{'_custom_labels'}->{$i}) {
	    $label = $self->{'_custom_labels'}->{$i};
        }
	else {
	    $label = $map{$col}->[4];
	}
	push @data, $label;
        push @underbars, '-' x length($label);

    }
    $printf_fmt =~ s/\\t$//;

    my $str = sprintf "$printf_fmt\n", @data;

    $str =~ s/\t\n/\n/;
    $str .= sprintf "$printf_fmt\n", @underbars;

    $str =~ s/\t\n/\n/gs;
    $self->{'_column_labels'} = $str;
}

# Purpose : Generate a function that will call the appropriate
# methods on the result, hit, and hsp objects to retrieve the column data 
# specified in the column spec.
#
# We should only have to go through the column spec once
# for a given ResultTableWriter. To permit this, we'll generate code 
# for a method that returns an array of the data for a row of output
# given a result, hit, and hsp object as arguments.
#
sub _set_row_data_func {
    my $self = shift;

    # Now we need to generate a string that can be eval'd to get the data.
    my @cols = $self->columns();
    my %map = $self->column_map;
    my @data;
    while( my $col = shift @cols ) {
	my $object = $map{$col}->[1];
	my $method = $map{$col}->[2];
        my $arg = '';
        if( $method =~ m!(\w+)/(\w+)! ) {
            $method = $1;
            $arg = "\"$2\"";
        }
        push @data, "\$$object->$method($arg)";
    }
    my $code = join( ",", @data);

    if( $self->verbose > 0 ) {
## Begin Debugging	
	$self->debug( "Data to print:\n");
	foreach( 0..$#data) { $self->debug( " [". ($_+ 1) . "] $data[$_]\n");}
	$self->debug( "CODE:\n$code\n");
	$self->debug("Printf format: ". $self->printf_fmt. "\n");
## End Debugging
    }

    my $func = sub {
        my ($result, $hit, $hsp) = @_;
        my @r = eval $code;
        # This should reduce the occurrence of those opaque "all zeros" bugs.
	if( $@ ) { $self->throw("Trouble in ResultTableWriter::_set_row_data_func() eval: $@\n\n"); 
               }
	return @r;
    };
    $self->{'_row_data_func'} = $func;
}

sub row_data_func { shift->{'_row_data_func'} }


=head2 to_string()

Note: this method is not intended for direct use. The
SearchIO::write_result() method calls it automatically if the writer
is hooked up to a SearchIO object as illustrated in L<the SYNOPSIS section | SYNOPSIS>.

 Title     : to_string()
           :
 Usage     : print $writer->to_string( $result_obj, [$include_labels] );
           :
 Argument  : $result_obj = A Bio::Search::Result::ResultI object
           : $include_labels = boolean, if true column labels are included (default: false)
           :
 Returns   : String containing tab-delimited set of data for each hit 
           : in a ResultI object. Some data is summed across multiple HSPs.
           :
 Throws    : n/a

=cut

#----------------
sub to_string {
#----------------
    my ($self, $result, $include_labels) = @_;

    my $str = $include_labels ? $self->column_labels() : '';
    my $resultfilter = $self->filter('RESULT');
    if( ! defined $resultfilter ||
        &{$resultfilter}($result) ) {	
	my @row_data  = &{$self->{'_row_data_func'}}( $result );
	$str .= sprintf "$self->{'_printf_fmt'}\n", @row_data;
	$str =~ s/\t\n/\n/gs;
    }
    return $str;
}



sub columns {
    my $self = shift;
    my @cols;
    if( ref $self->{'_cols'} ) {
        @cols = @{$self->{'_cols'}};
    }
    else {
        my %map = $self->column_map;
        @cols = sort { $map{$a}->[0] <=> $map{$b}->[0] } keys %map;
   }
    return @cols;
}


=head2 column_labels

 Usage     : print $result_obj->column_labels();
 Purpose   : Get column labels for to_string().
 Returns   : String containing column labels. Tab-delimited.
 Argument  : n/a
 Throws    : n/a

=cut

sub column_labels { shift->{'_column_labels'} }

=head2 end_report

 Title   : end_report
 Usage   : $self->end_report()
 Function: The method to call when ending a report, this is
           mostly for cleanup for formats which require you to 
           have something at the end of the document.  Nothing for
           a text message.
 Returns : string
 Args    : none

=cut

sub end_report {
    return '';
}

=head2 filter

 Title   : filter
 Usage   : $writer->filter('hsp', \&hsp_filter);
 Function: Filter out either at HSP,Hit,or Result level
 Returns : none
 Args    : string => data type,
           CODE reference


=cut


# Is this really needed?
#=head2 signif_format
#
# Usage     : $writer->signif_format( [FMT] );
# Purpose   : Allows retrieval of the P/Expect exponent values only
#           : or as a two-element list (mantissa, exponent).
# Usage     : $writer->signif_format('exp');
#           : $writer->signif_format('parts');
# Returns   : String or '' if not set.
# Argument  : String, FMT = 'exp' (return the exponent only)
#           :             = 'parts'(return exponent + mantissa in 2-elem list)
#           :              = undefined (return the raw value)
# Comments  : P/Expect values are still stored internally as the full,
#           : scientific notation value.
#
#=cut
#
##-------------
#sub signif_format {
##-------------
#    my $self = shift;
#    if(@_) { $self->{'_signif_format'} = shift; }
#    return $self->{'_signif_format'};
#}

1;
