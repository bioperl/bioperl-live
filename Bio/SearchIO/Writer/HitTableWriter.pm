
=head1 NAME

Bio::SearchIO::Writer::HitTableWriter - Tab-delimited data for Bio::Search::Hit::HitI objects

=head1 SYNOPSIS

=head2 Example 1: Using the default columns

    use Bio::SearchIO;
    use Bio::SearchIO::Writer::HitTableWriter;

    my $in = Bio::SearchIO->new();

    my $writer = Bio::SearchIO::Writer::HitTableWriter->new();

    my $out = Bio::SearchIO->new( -writer => $writer );

    while ( my $result = $in->next_result() ) {
        $out->write_result($result, ($in->report_count - 1 ? 0 : 1) );
    }

=head2 Example 2: Specifying a subset of columns 

    use Bio::SearchIO;
    use Bio::SearchIO::Writer::HitTableWriter;

    my $in = Bio::SearchIO->new();

    my $writer = Bio::SearchIO::Writer::HitTableWriter->new( 
                                  -columns => [qw(
                                                  query_name
                                                  query_length
                                                  hit_name
                                                  hit_length
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
parameter when creating the HitTableWriter object. The keys of the
hash should be the column number (left-most column = 1) for the label(s)
you want to specify. Here's an example:

    my $writer = Bio::SearchIO::Writer::HitTableWriter->new( 
                               -columns => [qw( query_name 
                                                query_length
                                                hit_name
                                                hit_length  )],
                               -labels  => { 1 => 'QUERY_GI',
  	                                     3 => 'HIT_IDENTIFIER' } );


=head1 DESCRIPTION

Bio::SearchIO::Writer::HitTableWriter outputs summary data 
for each Hit within a search result. Output is in tab-delimited format,
one row per Hit. 

The reason why this is considered summary data is that if a hit
contains multiple HSPs, the HSPs will be tiled and 
the data represents a summary across all HSPs.
See below for which columns are affected.
See the docs in L<Bio::Search::Hit::BlastHit|Bio::Search::Hit::BlastHit>
 for more details on HSP tiling.

=head2 Available Columns

Here are the columns that can be specified in the C<-columns>
parameter when creating a HitTableWriter object.  If a C<-columns> parameter
is not specified, this list, in this order, will be used as the default.

    query_name             # Sequence identifier of the query.
    query_length           # Full length of the query sequence
    hit_name               # Sequence identifier of the hit
    hit_length             # Full length of the hit sequence
    round                  # Round number for hit (PSI-BLAST)
    expect                 # Expect value for the alignment
    score                  # Score for the alignment (e.g., BLAST score)
    bits                   # Bit score for the alignment
    num_hsps               # Number of HSPs (not the "N" value)
    frac_identical_query*  # fraction of identical substitutions in query
    frac_identical_hit*    # fraction of identical substitutions in hit
    frac_conserved_query*  # fraction of conserved substitutions in query
    frac_conserved_hit*    # fraction of conserved substitutions in hit
    frac_aligned_query*    # fraction of the query sequence that is aligned
    frac_aligned_hit*      # fraction of the hit sequence that is aligned
    length_aln_query*      # Length of the aligned portion of the query sequence
    length_aln_hit*        # Length of the aligned portion of the hit sequence
    gaps_query*            # Number of gap characters in the aligned query sequence
    gaps_hit*              # Number of gap characters in the aligned hit sequence
    gaps_total*            # Number of gap characters in the aligned query and hit sequences
    start_query*           # Starting coordinate of the aligned portion of the query sequence
    end_query*             # Ending coordinate of the aligned portion of the query sequence
    start_hit*             # Starting coordinate of the aligned portion of the hit sequence
    end_hit*               # Ending coordinate of the aligned portion of the hit sequence
    strand_query           # Strand of the aligned query sequence
    strand_hit             # Strand of the aligned hit sequence
    frame                  # Frame of the alignment (0,1,2)
    ambiguous_aln          # Ambiguous alignment indicator ('qs', 'q', 's')
    hit_description        # Full description of the hit sequence
    query_description      # Full description of the query sequence
    rank                   # The rank order of the hit
    num_hits               # Number of hits for the query finding this hit

Items marked with a C<*> report data summed across all HSPs
after tiling them to avoid counting data from overlapping regions
multiple times.

For more details about these columns, see the documentation for the
corresponding method in Bio::Search::Result::BlastHit.

=head1 TODO

Figure out the best way to incorporate algorithm-specific score columns.
The best route is probably to have algorithm-specific subclasses 
(e.g., BlastHitTableWriter, FastaHitTableWriter).

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

Copyright (c) 2001, 2002 Steve Chervitz. All Rights Reserved.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 SEE ALSO

L<Bio::SearchIO::Writer::HitTableWriter>, 
L<Bio::SearchIO::Writer::ResultTableWriter>

=head1 METHODS

=cut

package Bio::SearchIO::Writer::HitTableWriter;

use strict;

use base qw(Bio::SearchIO::Writer::ResultTableWriter);


# Array fields: column, object, method[/argument], printf format,
# column label Methods for result object are defined in
# Bio::Search::Result::ResultI.  Methods for hit object are defined in
# Bio::Search::Hit::HitI.  Tech note: If a bogus method is supplied,
# it will result in all values to be zero.  Don't know why this is.

# TODO (maybe): Allow specification of separate mantissa/exponent for
# significance data.

my %column_map = (
                  'query_name'            => ['1', 'result', 'query_name', 's', 'QUERY' ],
                  'query_length'          => ['2', 'result', 'query_length', 'd', 'LEN_Q'],
                  'hit_name'              => ['3', 'hit', 'name', 's', 'HIT'],
                  'hit_length'            => ['4', 'hit', 'length', 'd', 'LEN_H'],
                  'round'                 => ['5', 'hit', 'iteration', 'd', 'ROUND'],
                  'expect'                => ['6', 'hit', 'significance', '.1e', 'EXPCT'],
                  'score'                 => ['7', 'hit', 'raw_score', 'd', 'SCORE'],
                  'bits'                  => ['8', 'hit', 'bits', 'd', 'BITS'],
                  'num_hsps'              => ['9', 'hit', 'num_hsps', 'd', 'HSPS'],
                  'frac_identical_query'  => ['10', 'hit', 'frac_identical/query', '.2f', 'FR_IDQ'],
                  'frac_identical_hit'    => ['11', 'hit', 'frac_identical/hit', '.2f', 'FR_IDH'],
                  'frac_conserved_query'  => ['12', 'hit', 'frac_conserved/query', '.2f', 'FR_CNQ'],
                  'frac_conserved_hit'    => ['13', 'hit', 'frac_conserved/hit', '.2f', 'FR_CNH'],
                  'frac_aligned_query'    => ['14', 'hit', 'frac_aligned_query', '.2f', 'FR_ALQ'],
                  'frac_aligned_hit'      => ['15', 'hit', 'frac_aligned_hit', '.2f', 'FR_ALH'],
                  'length_aln_query'      => ['16', 'hit', 'length_aln/query', 'd', 'LN_ALQ'],
                  'length_aln_hit'        => ['17', 'hit', 'length_aln/hit', 'd', 'LN_ALH'],
                  'gaps_query'            => ['18', 'hit', 'gaps/query', 'd', 'GAPS_Q'],
                  'gaps_hit'              => ['19', 'hit', 'gaps/hit', 'd', 'GAPS_H'],
                  'gaps_total'            => ['20', 'hit', 'gaps/total', 'd', 'GAPS_QH'],
                  'start_query'           => ['21', 'hit', 'start/query', 'd', 'START_Q'],
                  'end_query'             => ['22', 'hit', 'end/query', 'd', 'END_Q'],
                  'start_hit'             => ['23', 'hit', 'start/hit', 'd', 'START_H'],
                  'end_hit'               => ['24', 'hit', 'end/hit', 'd', 'END_H'],
                  'strand_query'          => ['25', 'hit', 'strand/query', 's', 'STRND_Q'],
                  'strand_hit'            => ['26', 'hit', 'strand/hit', 's', 'STRND_H'],
                  'frame'                 => ['27', 'hit', 'frame', 'd', 'FRAME'],
                  'ambiguous_aln'         => ['28', 'hit', 'ambiguous_aln', 's', 'AMBIG'],
                  'hit_description'       => ['29', 'hit', 'description', 's', 'DESC_H'],
                  'query_description'     => ['30', 'result', 'query_description', 's', 'DESC_Q'],
                  'rank'                  => ['31', 'hit', 'rank', 's', 'RANK'],
                  'num_hits'              => ['32', 'result', 'num_hits', 's', 'NUM_HITS'],
                 );

sub column_map { return %column_map }


=head2 to_string()

Note: this method is not intended for direct use. The
SearchIO::write_result() method calls it automatically if the writer
is hooked up to a SearchIO object as illustrated in 
L<the SYNOPSIS section | SYNOPSIS>.

 Title     : to_string()
           :
 Usage     : print $writer->to_string( $result_obj, [$include_labels] );
           :
 Argument  : $result_obj = A Bio::Search::Result::BlastResult object
           : $include_labels = boolean, if true column labels are included (default: false)
           :
 Returns   : String containing tab-delimited set of data for each hit 
           : in a BlastResult object. Some data is summed across multiple HSPs.
           :
 Throws    : n/a

=cut

#----------------
sub to_string {
#----------------
    my ($self, $result, $include_labels) = @_;

    my $str = $include_labels ? $self->column_labels() : '';
    my $func_ref = $self->row_data_func;
    my $printf_fmt = $self->printf_fmt;
    
    my ($resultfilter,$hitfilter) = ( $self->filter('RESULT'),
				      $self->filter('HIT') );
    if( ! defined $resultfilter ||
        &{$resultfilter}($result) ) {
	$result->can('rewind') && 
	    $result->rewind(); # insure we're at the beginning
	foreach my $hit($result->hits) {	    
	    next if( defined $hitfilter && ! &{$hitfilter}($hit));
	    my @row_data  = map { defined $_ ? $_ : 0 } &{$func_ref}($result, $hit);
	    $str .= sprintf "$printf_fmt\n", @row_data;
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
