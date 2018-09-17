
=head1 NAME

Bio::SearchIO::Writer::HSPTableWriter - Tab-delimited data for Bio::Search::HSP::HSPI objects

=head1 SYNOPSIS

=head2 Example 1: Using the default columns

    use Bio::SearchIO;
    use Bio::SearchIO::Writer::HSPTableWriter;

    my $in = Bio::SearchIO->new();

    my $writer = Bio::SearchIO::Writer::HSPTableWriter->new();

    my $out = Bio::SearchIO->new( -writer => $writer );

    while ( my $result = $in->next_result() ) {
        $out->write_result($result, ($in->report_count - 1 ? 0 : 1) );
    }

=head2 Example 2: Specifying a subset of columns 

    use Bio::SearchIO;
    use Bio::SearchIO::Writer::HSPTableWriter;

    my $in = Bio::SearchIO->new();

    my $writer = Bio::SearchIO::Writer::HSPTableWriter->new( 
                                  -columns => [qw(
                                                  query_name
                                                  query_length
                                                  hit_name
                                                  hit_length
                                                  rank
                                                  frac_identical_query
                                                  expect
                                                  )]  );

    my $out = Bio::SearchIO->new( -writer => $writer,
				  -file   => ">searchio.out" );

    while ( my $result = $in->next_result() ) {
        $out->write_result($result, ($in->report_count - 1 ? 0 : 1) );
    }

=head2 Custom Labels

You can also specify different column labels if you don't want to use
the defaults.  Do this by specifying a C<-labels> hash reference
parameter when creating the HSPTableWriter object.  The keys of the
hash should be the column number (left-most column = 1) for the label(s)
you want to specify. Here's an example:

    my $writer = Bio::SearchIO::Writer::HSPTableWriter->new( 
                               -columns => [qw( query_name 
                                                query_length
                                                hit_name
                                                hit_length  )],
                               -labels  => { 1 => 'QUERY_GI',
  	                                     3 => 'HIT_IDENTIFIER' } );


=head1 DESCRIPTION

Bio::SearchIO::Writer::HSPTableWriter generates output at the finest
level of granularity for data within a search result. Data for each HSP
within each hit in a search result is output in tab-delimited format,
one row per HSP.

=head2 Available Columns

Here are the columns that can be specified in the C<-columns>
parameter when creating a HSPTableWriter object.  If a C<-columns> parameter
is not specified, this list, in this order, will be used as the default.

    query_name             # Sequence identifier of the query.
    query_length           # Full length of the query sequence
    hit_name               # Sequence identifier of the hit
    hit_length             # Full length of the hit sequence
    round                  # Round number for hit (PSI-BLAST)
    rank
    expect                 # Expect value for the alignment
    score                  # Score for the alignment (e.g., BLAST score)
    bits                   # Bit score for the alignment
    frac_identical_query   # fraction of identical substitutions in query
    frac_identical_hit     # fraction of identical substitutions in hit
    frac_conserved_query   # fraction of conserved substitutions in query
    frac_conserved_hit     # fraction of conserved substitutions in hit
    length_aln_query       # Length of the aligned portion of the query sequence
    length_aln_hit         # Length of the aligned portion of the hit sequence
    gaps_query             # Number of gap characters in the aligned query sequence
    gaps_hit               # Number of gap characters in the aligned hit sequence
    gaps_total             # Number of gap characters in the aligned query and hit sequences
    start_query            # Starting coordinate of the aligned portion of the query sequence
    end_query              # Ending coordinate of the aligned portion of the query sequence
    start_hit              # Starting coordinate of the aligned portion of the hit sequence
    end_hit                # Ending coordinate of the aligned portion of the hit sequence
    strand_query           # Strand of the aligned query sequence
    strand_hit             # Strand of the aligned hit sequence
    frame                  # Reading frame of the aligned query sequence 
    hit_description        # Full description of the hit sequence
    query_description      # Full description of the query sequence
    frac_identical_total   # fraction of total identical substitutions
    frac_conserved_total   # fraction of total conserved substitutions

For more details about these columns, see the documentation for the
corresponding method in Bio::Search::HSP::HSPI.

=head1 TODO

Figure out the best way to incorporate algorithm-specific score columns.
The best route is probably to have algorith-specific subclasses 
(e.g., BlastHSPTableWriter, FastaHSPTableWriter).

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

    Bio::SearchIO::Writer::HitTableWriter
    Bio::SearchIO::Writer::ResultTableWriter

=head1 METHODS

=cut

package Bio::SearchIO::Writer::HSPTableWriter;

use strict;

use base qw(Bio::SearchIO::Writer::ResultTableWriter);


# Array fields: column, object, method[/argument], printf format, column label
# Methods for result object are defined in Bio::Search::Result::ResultI.
# Methods for hit object are defined in Bio::Search::Hit::HitI.
# Methods for hsp object are defined in Bio::Search::HSP::HSPI.
# Tech note: If a bogus method is supplied, it will result in all values to be zero.
#            Don't know why this is.
# TODO (maybe): Allow specification of signif_format (i.e., separate mantissa/exponent)
my %column_map = (
                  'query_name'            => ['1', 'result', 'query_name', 's', 'QUERY' ],
                  'query_length'          => ['2', 'result', 'query_length', 'd', 'LEN_Q'],
                  'hit_name'              => ['3', 'hit', 'name', 's', 'HIT'],
                  'hit_length'            => ['4', 'hit', 'hit_length', 'd', 'LEN_H'],
                  'round'                 => ['5', 'hit', 'iteration', 'd', 'ROUND', 'hit'],
                  'rank'                  => ['6', 'hsp', 'rank', 'd', 'RANK'],
                  'expect'                => ['7', 'hsp', 'expect', '.1e', 'EXPCT'],
                  'score'                 => ['8', 'hsp', 'score', 'd', 'SCORE'],
                  'bits'                  => ['9', 'hsp', 'bits', 'd', 'BITS'],
                  'frac_identical_query'  => ['10', 'hsp', 'frac_identical/query', '.2f', 'FR_IDQ'],
                  'frac_identical_hit'    => ['11', 'hsp', 'frac_identical/hit', '.2f', 'FR_IDH'],
                  'frac_conserved_query'  => ['12', 'hsp', 'frac_conserved/query', '.2f', 'FR_CNQ'],
                  'frac_conserved_hit'    => ['13', 'hsp', 'frac_conserved/hit', '.2f', 'FR_CNH'],
                  'length_aln_query'      => ['14', 'hsp', 'length/query', 'd', 'LN_ALQ'],
                  'length_aln_hit'        => ['15', 'hsp', 'length/hit', 'd', 'LN_ALH'],
                  'gaps_query'            => ['16', 'hsp', 'gaps/query', 'd', 'GAPS_Q'],
                  'gaps_hit'              => ['17', 'hsp', 'gaps/hit', 'd', 'GAPS_H'],
                  'gaps_total'            => ['18', 'hsp', 'gaps/total', 'd', 'GAPS_QH'],
                  'start_query'           => ['19', 'hsp', 'start/query', 'd', 'START_Q'],
                  'end_query'             => ['20', 'hsp', 'end/query', 'd', 'END_Q'],
                  'start_hit'             => ['21', 'hsp', 'start/hit', 'd', 'START_H'],
                  'end_hit'               => ['22', 'hsp', 'end/hit', 'd', 'END_H'],
                  'strand_query'          => ['23', 'hsp', 'strand/query', 'd', 'STRND_Q'],
                  'strand_hit'            => ['24', 'hsp', 'strand/hit', 'd', 'STRND_H'],
                  'frame_hit'             => ['25', 'hsp', 'frame/hit', 's', 'FRAME_H'],
                  'frame_query'           => ['26', 'hsp', 'frame/query', 's', 'FRAME_Q'],
                  'hit_description'       => ['27', 'hit', 'hit_description', 's', 'DESC_H'],
                  'query_description'     => ['28', 'result', 'query_description', 's', 'DESC_Q'],
                  'frac_identical_total'  => ['29', 'hsp', 'frac_identical/total', '.2f', 'FR_IDT'],
                  'frac_conserved_total'  => ['30', 'hsp', 'frac_conserved/total', '.2f', 'FR_CNT'],
                 );

sub column_map { return %column_map }

=head2 to_string()

Note: this method is not intended for direct use. 
The SearchIO::write_result() method calls it automatically 
if the writer is hooked up to a SearchIO object as illustrated in
L<the SYNOPSIS section | SYNOPSIS>.

 Title     : to_string()
           :
 Usage     : print $writer->to_string( $result_obj, [$include_labels] );
           :
 Argument  : $result_obj = A Bio::Search::Result::ResultI object
           : $include_labels = boolean, if true column labels are included (default: false)
           :
 Returns   : String containing tab-delimited set of data for each HSP
           : in each Hit of the supplied ResultI object. 
           :
 Throws    : n/a

=cut

sub to_string {
    my ($self, $result, $include_labels) = @_;
    
    my $str = $include_labels ? $self->column_labels() : '';
    my ($resultfilter,$hitfilter,
	$hspfilter) = ( $self->filter('RESULT'),
			$self->filter('HIT'),
			$self->filter('HSP'));
    if( ! defined $resultfilter || &{$resultfilter}($result) ) {
	my $func_ref = $self->row_data_func;
	my $printf_fmt = $self->printf_fmt;
	$result->can('rewind') && 
	    $result->rewind(); # insure we're at the beginning
	while( my $hit = $result->next_hit) {
	    next if( defined $hitfilter && ! &{$hitfilter}($hit) );
	    $hit->can('rewind') && $hit->rewind;# insure we're at the beginning
	    while(my $hsp = $hit->next_hsp) {
            next if ( defined $hspfilter && ! &{$hspfilter}($hsp));
            my @row_data  = &{$func_ref}($result, $hit, $hsp);
            $str .= sprintf("$printf_fmt\n", map {$_ || ($printf_fmt eq 's' ? '' : 0)} @row_data);
	    }
	}
    }
    $str =~ s/\t\n/\n/gs;
    return $str;
}

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


1;
