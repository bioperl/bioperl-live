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
parameter when creating the HitTableWriter object.  The keys of the
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

    query_name                 # Sequence identifier of the query.
    query_length               # Full length of the query sequence.
    hit_name                   # Sequence identifier of the hit
    hit_length                 # Full length of the hit sequence
    round                      # Round number for hit (PSI-BLAST).
    expect                     # Expect value for the alignment.
    score                      # Score for the alignment.
    bits
    num_hsps
    frac_identical_query*
    frac_identical_hit*
    frac_conserved_query*
    frac_conserved_hit*
    frac_aligned_query*
    frac_aligned_hit*
    length_aln_query*
    length_aln_hit*
    gaps_query*
    gaps_hit*
    gaps_total*
    start_query*
    end_query*
    start_hit*
    end_hit*
    strand_query
    strand_hit
    ambiguous_aln
    hit_description
    query_description

Items marked with a C<*> report data summed across all HSPs
after tiling them to avoid counting data from overlapping regions
multiple times.

For more details about these columns, see the documentation for the
corresponding method in Bio::Search::Result::BlastHit.

=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other
Bioperl modules.  Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

    bioperl-l@bioperl.org              - General discussion
    http://bio.perl.org/MailList.html  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via email
or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           

=head1 AUTHOR 

Steve Chervitz E<lt>sac@bioperl.orgE<gt>

See L<the FEEDBACK section | FEEDBACK> for where to send bug reports
and comments.

=head1 SEE ALSO

L<Bio::SearchIO::Writer::HitTableWriter>, 
L<Bio::SearchIO::Writer::ResultTableWriter>

=head1 METHODS

=cut

package Bio::SearchIO::Writer::HitTableWriter;

use strict;
use Bio::SearchIO::Writer::ResultTableWriter;

use vars qw( @ISA );
@ISA = qw( Bio::SearchIO::Writer::ResultTableWriter );


# Array fields: column, object, method[/argument], printf format, column label
# Methods for result object are defined in Bio::Search::Result::ResultI.
# Methods for hit object are defined in Bio::Search::Hit::HitI.
# Tech note: If a bogus method is supplied, it will result in all values to be zero.
#            Don't know why this is.
# TODO (maybe): Allow specification of separate mantissa/exponent for significance data.
my %column_map = (
                  'query_name'            => ['1', 'result', 'query_name', 's', 'QUERY' ],
                  'query_length'          => ['2', 'result', 'query_length', 'd', 'LEN_Q'],
                  'hit_name'              => ['3', 'hit', 'name', 's', 'HIT'],
                  'hit_length'            => ['4', 'hit', 'length', 'd', 'LEN_H'],
                  'round'                 => ['5', 'hit', 'psiblast_round', 'd', 'ROUND'],
                  'expect'                => ['6', 'hit', 'expect', '.1e', 'EXPCT'],
                  'score'                 => ['7', 'hit', 'score', 'd', 'SCORE'],
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
                  'ambiguous_aln'         => ['27', 'hit', 'ambiguous_aln', 's', 'AMBIG'],
                  'hit_description'       => ['28', 'hit', 'description', 's', 'DESC_H'],
                  'query_description'     => ['29', 'result', 'query_description', 's', 'DESC_Q'],
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

    foreach my $hit($result->hits) {
        my @row_data  = &{$func_ref}($result, $hit);
        $str .= sprintf "$printf_fmt\n", @row_data;
    }
    $str =~ s/\t\n/\n/gs;
    return $str;
}


1;
__END__

  TODO: Integrate the column descriptions into the POD

           : Left-to-Right order of fields is customizable (see new()).
           : Here is the default order:
           :
           : 1 QUERY_NAME           # Sequence identifier of the query.
           : 2 QUERY_LENGTH         # Full length of the query sequence.
           : 3 HIT_NAME             # Sequence identifier of the hit
           : 4 HIT_LENGTH           # Full length of the hit sequence.
           : 5 ROUND                # Round number for hit (PSI-BLAST).
           : 6 EXPECT               # Expect value for the alignment.
           : 7 SCORE                # Result score for the alignment.
           : 8 BITS                 # Bit score for the alignment.
           : 9 NUM_HSPS             # Number of HSPs (not the "N" value).
           : 10 FRAC_IDENTICAL_QUERY*      # fraction of identical substitutions in query.
           : 11 FRAC_IDENTICAL_HIT*      # fraction of identical substitutions in query.
           : 12 FRAC_CONSERVED_QUERY*     # fraction of conserved ("positive") substitutions in query.
           : 13 FRAC_CONSERVED_HIT*     # fraction of conserved ("positive") substitutions in hit.
           : 14 FRAC_ALN_QUERY*     # fraction of the query sequence that is aligned.
           : 15 FRAC_ALN_HIT*       # fraction of the hit sequence that is aligned.
           : 16 LENGTH_ALN_QUERY*   # Length of the aligned portion of the query sequence.
           : 17 LENGTH_ALN_HIT*     # Length of the aligned portion of the hit sequence.
           : 18 GAPS_QUERY*         # Number of gaps in the aligned query sequence.
           : 19 GAPS_HIT*           # Number of gaps in the aligned hit sequence.
           : 20 GAPS_TOTAL*         # Total number of gaps in the aligned query and hit sequences.
           : 21 START_QUERY*        # Starting coordinate of the query sequence.
           : 22 END_QUERY*          # Ending coordinate of the query sequence.
           : 23 START_HIT*          # Starting coordinate of the hit sequence.
           : 24 END_HIT*            # Ending coordinate of the hit sequence.
           : 25 AMBIGUOUS_ALN       # Ambiguous alignment indicator ('qs', 'q', 's').
           : 26 DESC_HIT            # Full description of the hit sequence.
           : 27 DESC_QUERY          # Full description of the query sequence.
           :
           : * Items marked with a "*" report data summed across all HSPs
           :   after tiling them to avoid counting data from overlapping regions
           :   multiple times.
