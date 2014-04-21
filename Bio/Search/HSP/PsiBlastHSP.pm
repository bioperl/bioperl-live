#-----------------------------------------------------------------
#
# BioPerl module Bio::Search::HSP::PsiBlastHSP
#
# (This module was originally called Bio::Tools::Blast::HSP)
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Steve Chervitz <sac@bioperl.org>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

## POD Documentation:

=head1 NAME

Bio::Search::HSP::PsiBlastHSP - Bioperl BLAST High-Scoring Pair object

=head1 SYNOPSIS

See L<Bio::Search::Hit::BlastHit>.

=head1 DESCRIPTION

A Bio::Search::HSP::PsiBlastHSP object provides an interface to data
obtained in a single alignment section of a Blast report (known as a
"High-scoring Segment Pair"). This is essentially a pairwise
alignment with score information.

PsiBlastHSP objects are accessed via L<Bio::Search::Hit::BlastHit>
objects after parsing a BLAST report using the L<Bio::SearchIO>
system.

The construction of PsiBlastHSP objects is performed by
Bio::Factory::BlastHitFactory in a process that is
orchestrated by the Blast parser (L<Bio::SearchIO::psiblast>).
The resulting PsiBlastHSPs are then accessed via
L<Bio::Search::Hit::BlastHit>). Therefore, you do not need to
use L<Bio::Search::HSP::PsiBlastHSP>) directly. If you need to construct
PsiBlastHSPs directly, see the new() function for details.

For L<Bio::SearchIO> BLAST parsing usage examples, see the
C<examples/search-blast> directory of the Bioperl distribution.


=head2 Start and End coordinates

Sequence endpoints are swapped so that start is always less than
end. This affects For TBLASTN/X hits on the minus strand. Strand
information can be recovered using the strand() method. This
normalization step is standard Bioperl practice. It also facilitates
use of range information by methods such as match().

=over 1

=item * Supports BLAST versions 1.x and 2.x, gapped and ungapped.

=back

Bio::Search::HSP::PsiBlastHSP.pm has the ability to extract a list of all
residue indices for identical and conservative matches along both
query and sbjct sequences. Since this degree of detail is not always
needed, this behavior does not occur during construction of the PsiBlastHSP
object.  These data will automatically be collected as necessary as
the PsiBlastHSP.pm object is used.

=head1 DEPENDENCIES

Bio::Search::HSP::PsiBlastHSP.pm is a concrete class that inherits from
L<Bio::SeqFeature::SimilarityPair> and L<Bio::Search::HSP::HSPI>.
L<Bio::Seq> and L<Bio::SimpleAlign> are employed for creating
sequence and alignment objects, respectively.

=head2 Relationship to L<Bio::SimpleAlign> and L<Bio::Seq>

PsiBlastHSP.pm can provide the query or sbjct sequence as a L<Bio::Seq>
object via the L<seq()|seq> method. The PsiBlastHSP.pm object can also create a
two-sequence L<Bio::SimpleAlign> alignment object using the the query
and sbjct sequences via the L<get_aln()|get_aln> method. Creation of alignment
objects is not automatic when constructing the PsiBlastHSP.pm object since
this level of functionality is not always required and would generate
a lot of extra overhead when crunching many reports.


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

Steve Chervitz E<lt>sac-at-bioperl.orgE<gt>

See L<the FEEDBACK section | FEEDBACK> for where to send bug reports and comments.

=head1 ACKNOWLEDGEMENTS

This software was originally developed in the Department of Genetics
at Stanford University. I would also like to acknowledge my
colleagues at Affymetrix for useful feedback.

=head1 SEE ALSO

 Bio::Search::Hit::BlastHit.pm          - Blast hit object.
 Bio::Search::Result::BlastResult.pm    - Blast Result object.
 Bio::Seq.pm                            - Biosequence object

=head2 Links:

 http://bio.perl.org/                       - Bioperl Project Homepage

=head1 COPYRIGHT

Copyright (c) 1996-2001 Steve Chervitz. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut


# END of main POD documentation.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Search::HSP::PsiBlastHSP;

use strict;
use Bio::SeqFeature::Similarity;

use vars qw($GAP_SYMBOL %STRAND_SYMBOL);

use overload
    '""' => \&to_string;

use base qw(Bio::SeqFeature::SimilarityPair Bio::Search::HSP::HSPI);

$GAP_SYMBOL    = '-';          # Need a more general way to handle gap symbols.
%STRAND_SYMBOL = ('Plus' => 1, 'Minus' => -1 );


=head2 new

 Usage     : $hsp = Bio::Search::HSP::PsiBlastHSP->new( %named_params );
           : Bio::Search::HSP::PsiBlastHSP.pm objects are constructed
           : automatically by Bio::SearchIO::BlastHitFactory.pm,
           : so there is no need for direct instantiation.
 Purpose   : Constructs a new PsiBlastHSP object and Initializes key variables
           : for the HSP.
 Returns   : A Bio::Search::HSP::PsiBlastHSP object
 Argument  : Named parameters:
           : Parameter keys are case-insensitive.
           :      -RAW_DATA  => array ref containing raw BLAST report data for
           :                    for a single HSP. This includes all lines
           :                    of the HSP alignment from a traditional BLAST
                                or PSI-BLAST (non-XML) report,
           :      -RANK         => integer (1..n).
           :      -PROGRAM      => string ('TBLASTN', 'BLASTP', etc.).
           :      -QUERY_NAME   => string, id of query sequence
           :      -HIT_NAME     => string, id of hit sequence
           :
 Comments  : Having the raw data allows this object to do lazy parsing of
           : the raw HSP data (i.e., not parsed until needed).
           :
           : Note that there is a fair amount of basic parsing that is
           : currently performed in this module that would be more appropriate
           : to do within a separate factory object.
           : This parsing code will likely be relocated and more initialization
           : parameters will be added to new().
           :
See Also   : L<Bio::SeqFeature::SimilarityPair::new()>, L<Bio::SeqFeature::Similarity::new()>

=cut

#----------------
sub new {
#----------------
    my ($class, @args ) = @_;

    my $self = $class->SUPER::new( @args );
    # Initialize placeholders
    $self->{'_queryGaps'} = $self->{'_sbjctGaps'} = 0;
    my ($raw_data, $qname, $hname, $qlen, $hlen);

    ($self->{'_prog'}, $self->{'_rank'}, $raw_data,
     $qname, $hname) =
      $self->_rearrange([qw( PROGRAM
                             RANK
                             RAW_DATA
                             QUERY_NAME
                             HIT_NAME
                           )], @args );

    # _set_data() does a fair amount of parsing.
    # This will likely change (see comment above.)
    $self->_set_data( @{$raw_data} );
    # Store the aligned query as sequence feature
    my ($qb, $hb) = ($self->start());
    my ($qe, $he) = ($self->end());
    my ($qs, $hs) = ($self->strand());
    my ($qf,$hf) = ($self->query->frame(),
                    $self->hit->frame);

    $self->query( Bio::SeqFeature::Similarity->new (-start   =>$qb,
                                                    -end     =>$qe,
                                                    -strand  =>$qs,
                                                    -bits    =>$self->bits,
                                                    -score   =>$self->score,
                                                    -frame   =>$qf,
                                                    -seq_id  => $qname,
                                                    -source  =>$self->{'_prog'} ));

    $self->hit( Bio::SeqFeature::Similarity->new (-start   =>$hb,
                                                  -end     =>$he,
                                                  -strand  =>$hs,
                                                  -bits    =>$self->bits,
                                                  -score   =>$self->score,
                                                  -frame   =>$hf,
                                                  -seq_id  => $hname,
                                                  -source  =>$self->{'_prog'} ));

    # set lengths
    $self->query->seqlength($qlen); # query
    $self->hit->seqlength($hlen); # subject

    $self->query->frac_identical($self->frac_identical('query'));
    $self->hit->frac_identical($self->frac_identical('hit'));
    return $self;
}

#sub DESTROY {
#    my $self = shift;
#    #print STDERR "--->DESTROYING $self\n";
#}


# Title   : _id_str;
# Purpose : Intended for internal use only to provide a string for use
#           within exception messages to help users figure out which
#           query/hit caused the problem.
# Returns : Short string with name of query and hit seq
sub _id_str {
    my $self = shift;
    if( not defined $self->{'_id_str'}) {
        my $qname = $self->query->seqname;
        my $hname = $self->hit->seqname;
        $self->{'_id_str'} = "QUERY=\"$qname\" HIT=\"$hname\"";
    }
    return $self->{'_id_str'};
}

#=================================================
# Begin Bio::Search::HSP::HSPI implementation
#=================================================

=head2 algorithm

 Title   : algorithm
 Usage   : $alg = $hsp->algorithm();
 Function: Gets the algorithm specification that was used to obtain the hsp
           For BLAST, the algorithm denotes what type of sequence was aligned
           against what (BLASTN: dna-dna, BLASTP prt-prt, BLASTX translated
           dna-prt, TBLASTN prt-translated dna, TBLASTX translated
           dna-translated dna).
 Returns : a scalar string
 Args    : none

=cut

#----------------
sub algorithm {
#----------------
    my ($self,@args) = @_;
    return $self->{'_prog'};
}




=head2 signif()

 Usage     : $hsp_obj->signif()
 Purpose   : Get the P-value or Expect value for the HSP.
 Returns   : Float (0.001 or 1.3e-43)
           : Returns P-value if it is defined, otherwise, Expect value.
 Argument  : n/a
 Throws    : n/a
 Comments  : Provided for consistency with BlastHit::signif()
           : Support for returning the significance data in different
           : formats (e.g., exponent only), is not provided for HSP objects.
           : This is only available for the BlastHit or Blast object.

See Also   : L</p>, L</expect>, L<Bio::Search::Hit::BlastHit::signif()|Bio::Search::Hit::BlastHit>

=cut

#-----------
sub signif {
#-----------
    my $self = shift;
    my $val ||= defined($self->{'_p'}) ? $self->{'_p'} : $self->{'_expect'};
    $val;
}



=head2 evalue

 Usage     : $hsp_obj->evalue()
 Purpose   : Get the Expect value for the HSP.
 Returns   : Float (0.001 or 1.3e-43)
 Argument  : n/a
 Throws    : n/a
 Comments  : Support for returning the expectation data in different
           : formats (e.g., exponent only), is not provided for HSP objects.
           : This is only available for the BlastHit or Blast object.

See Also   : L</p>

=cut

#----------
sub evalue { shift->{'_expect'} }
#----------


=head2 p

 Usage     : $hsp_obj->p()
 Purpose   : Get the P-value for the HSP.
 Returns   : Float (0.001 or 1.3e-43) or undef if not defined.
 Argument  : n/a
 Throws    : n/a
 Comments  : P-value is not defined with NCBI Blast2 reports.
           : Support for returning the expectation data in different
           : formats (e.g., exponent only) is not provided for HSP objects.
           : This is only available for the BlastHit or Blast object.

See Also   : L</expect>

=cut

#-----
sub p { my $self = shift; $self->{'_p'}; }
#-----

# alias
sub pvalue { shift->p(@_); }

=head2 length

 Usage     : $hsp->length( [seq_type] )
 Purpose   : Get the length of the aligned portion of the query or sbjct.
 Example   : $hsp->length('query')
 Returns   : integer
 Argument  : seq_type: 'query' | 'hit' or 'sbjct' | 'total'  (default = 'total')
             ('sbjct' is synonymous with 'hit')
 Throws    : n/a
 Comments  : 'total' length is the full length of the alignment
           : as reported in the denominators in the alignment section:
           : "Identical = 34/120 Positives = 67/120".

See Also   : L</gaps>

=cut

#-----------
sub length {
#-----------
## Developer note: when using the built-in length function within
##                 this module, call it as CORE::length().
    my( $self, $seqType ) = @_;
    $seqType  ||= 'total';
    $seqType = 'sbjct' if $seqType eq 'hit';

    $seqType ne 'total' and $self->_set_seq_data() unless $self->{'_set_seq_data'};

    ## Sensitive to member name format.
    $seqType = "_\L$seqType\E";
    $self->{$seqType.'Length'};
}



=head2 gaps

 Usage     : $hsp->gaps( [seq_type] )
 Purpose   : Get the number of gap characters in the query, sbjct, or total alignment.
           : Also can return query gap chars and sbjct gap chars as a two-element list
           : when in array context.
 Example   : $total_gaps      = $hsp->gaps();
           : ($qgaps, $sgaps) = $hsp->gaps();
           : $qgaps           = $hsp->gaps('query');
 Returns   : scalar context: integer
           : array context without args: (int, int) = ('queryGaps', 'sbjctGaps')
 Argument  : seq_type: 'query' or 'hit' or 'sbjct' or 'total'
           :  ('sbjct' is synonymous with 'hit')
           : (default = 'total', scalar context)
           : Array context can be "induced" by providing an argument of 'list' or 'array'.
 Throws    : n/a

See Also   : L</length>, L</matches>

=cut

#---------
sub gaps {
#---------
    my( $self, $seqType ) = @_;

    $self->_set_seq_data() unless $self->{'_set_seq_data'};

    $seqType  ||= (wantarray ? 'list' : 'total');
    $seqType = 'sbjct' if $seqType eq 'hit';

    if($seqType =~ /list|array/i) {
        return (($self->{'_queryGaps'} || 0), ($self->{'_sbjctGaps'} || 0));
    }

    if($seqType eq 'total') {
        return ($self->{'_queryGaps'} + $self->{'_sbjctGaps'}) || 0;
    } else {
        ## Sensitive to member name format.
        $seqType = "_\L$seqType\E";
        return $self->{$seqType.'Gaps'} || 0;
    }
}


=head2 frac_identical

 Usage     : $hsp_object->frac_identical( [seq_type] );
 Purpose   : Get the fraction of identical positions within the given HSP.
 Example   : $frac_iden = $hsp_object->frac_identical('query');
 Returns   : Float (2-decimal precision, e.g., 0.75).
 Argument  : seq_type: 'query' or 'hit' or 'sbjct' or 'total'
           :  ('sbjct' is synonymous with 'hit')
           : default = 'total' (but see comments below).
 Throws    : n/a
 Comments  : Different versions of Blast report different values for the total
           : length of the alignment. This is the number reported in the
           : denominators in the stats section:
           : "Identical = 34/120 Positives = 67/120".
           : NCBI-BLAST uses the total length of the alignment (with gaps)
           : WU-BLAST uses the length of the query sequence (without gaps).
           : Therefore, when called without an argument or an argument of 'total',
           : this method will report different values depending on the
           : version of BLAST used.
           :
           : To get the fraction identical among only the aligned residues,
           : ignoring the gaps, call this method with an argument of 'query'
           : or 'sbjct' ('sbjct' is synonymous with 'hit').

See Also   : L</frac_conserved>, L</num_identical>, L</matches>

=cut

#-------------------
sub frac_identical {
#-------------------
# The value is calculated as opposed to storing it from the parsed results.
# This saves storage and also permits flexibility in determining for which
# sequence (query or sbjct) the figure is to be calculated.

    my( $self, $seqType ) = @_;
    $seqType ||= 'total';
    $seqType = 'sbjct' if $seqType eq 'hit';

    if($seqType ne 'total') {
      $self->_set_seq_data() unless $self->{'_set_seq_data'};
    }
    ## Sensitive to member name format.
    $seqType = "_\L$seqType\E";

    sprintf( "%.2f", $self->{'_numIdentical'}/$self->{$seqType.'Length'});
}


=head2 frac_conserved

 Usage     : $hsp_object->frac_conserved( [seq_type] );
 Purpose   : Get the fraction of conserved positions within the given HSP.
           : (Note: 'conservative' positions are called 'positives' in the
           : Blast report.)
 Example   : $frac_cons = $hsp_object->frac_conserved('query');
 Returns   : Float (2-decimal precision, e.g., 0.75).
 Argument  : seq_type: 'query' or 'hit' or 'sbjct' or 'total'
           :  ('sbjct' is synonymous with 'hit')
           : default = 'total' (but see comments below).
 Throws    : n/a
 Comments  : Different versions of Blast report different values for the total
           : length of the alignment. This is the number reported in the
           : denominators in the stats section:
           : "Identical = 34/120 Positives = 67/120".
           : NCBI-BLAST uses the total length of the alignment (with gaps)
           : WU-BLAST uses the length of the query sequence (without gaps).
           : Therefore, when called without an argument or an argument of 'total',
           : this method will report different values depending on the
           : version of BLAST used.
           :
           : To get the fraction conserved among only the aligned residues,
           : ignoring the gaps, call this method with an argument of 'query'
           : or 'sbjct'.

See Also   : L</frac_conserved>, L</num_conserved>, L</matches>

=cut

#--------------------
sub frac_conserved {
#--------------------
# The value is calculated as opposed to storing it from the parsed results.
# This saves storage and also permits flexibility in determining for which
# sequence (query or sbjct) the figure is to be calculated.

    my( $self, $seqType ) = @_;
    $seqType ||= 'total';
    $seqType = 'sbjct' if $seqType eq 'hit';

    if($seqType ne 'total') {
      $self->_set_seq_data() unless $self->{'_set_seq_data'};
    }

    ## Sensitive to member name format.
    $seqType = "_\L$seqType\E";

    sprintf( "%.2f", $self->{'_numConserved'}/$self->{$seqType.'Length'});
}

=head2 query_string

 Title   : query_string
 Usage   : my $qseq = $hsp->query_string;
 Function: Retrieves the query sequence of this HSP as a string
 Returns : string
 Args    : none


=cut

#----------------
sub query_string{ shift->seq_str('query'); }
#----------------

=head2 hit_string

 Title   : hit_string
 Usage   : my $hseq = $hsp->hit_string;
 Function: Retrieves the hit sequence of this HSP as a string
 Returns : string
 Args    : none


=cut

#----------------
sub hit_string{ shift->seq_str('hit'); }
#----------------


=head2 homology_string

 Title   : homology_string
 Usage   : my $homo_string = $hsp->homology_string;
 Function: Retrieves the homology sequence for this HSP as a string.
         : The homology sequence is the string of symbols in between the
         : query and hit sequences in the alignment indicating the degree
         : of conservation (e.g., identical, similar, not similar).
 Returns : string
 Args    : none

=cut

#----------------
sub homology_string{ shift->seq_str('match'); }
#----------------

#=================================================
# End Bio::Search::HSP::HSPI implementation
#=================================================

# Older method delegating to method defined in HSPI.

=head2 expect

See L<Bio::Search::HSP::HSPI::expect()|Bio::Search::HSP::HSPI>

=cut

#----------
sub expect { shift->evalue( @_ ); }
#----------


=head2 rank

 Usage     : $hsp->rank( [string] );
 Purpose   : Get the rank of the HSP within a given Blast hit.
 Example   : $rank = $hsp->rank;
 Returns   : Integer (1..n) corresponding to the order in which the HSP
             appears in the BLAST report.

=cut

#'

#----------
sub rank { shift->{'_rank'} }
#----------

# For backward compatibility
#----------
sub name { shift->rank }
#----------

=head2 to_string

 Title   : to_string
 Usage   : print $hsp->to_string;
 Function: Returns a string representation for the Blast HSP.
           Primarily intended for debugging purposes.
 Example : see usage
 Returns : A string of the form:
           [PsiBlastHSP] <rank>
           e.g.:
           [BlastHit] 1
 Args    : None

=cut

#----------
sub to_string {
#----------
    my $self = shift;
    return "[PsiBlastHSP] " . $self->rank();
}


=head2 _set_data

 Usage     : called automatically during object construction.
 Purpose   : Parses the raw HSP section from a flat BLAST report and
             sets the query sequence, sbjct sequence, and the "match" data
           : which consists of the symbols between the query and sbjct lines
           : in the alignment.
 Argument  : Array (all lines for a single, complete HSP, from a raw,
             flat (i.e., non-XML) BLAST report)
 Throws    : Propagates any exceptions from the methods called ("See Also")

See Also   : L</_set_seq>, L</_set_score_stats>, L</_set_match_stats>

=cut

#--------------
sub _set_data {
#--------------
    my $self = shift;
    my @data = @_;
    my @queryList  = ();  # 'Query' = SEQUENCE USED TO QUERY THE DATABASE.
    my @sbjctList  = ();  # 'Sbjct' = HOMOLOGOUS SEQUENCE FOUND IN THE DATABASE.
    my @matchList  = ();
    my $matchLine  = 0;   # Alternating boolean: when true, load 'match' data.
    my @linedat = ();

    #print STDERR "PsiBlastHSP: set_data()\n";

    my($line, $aln_row_len, $length_diff);
    $length_diff = 0;

    # Collecting data for all lines in the alignment
    # and then storing the collections for possible processing later.
    #
    # Note that "match" lines may not be properly padded with spaces.
    # This loop now properly handles such cases:
    # Query: 1141 PSLVELTIRDCPRLEVGPMIRSLPKFPMLKKLDLAVANIIEEDLDVIGSLEELVIXXXXX 1200
    #             PSLVELTIRDCPRLEVGPMIRSLPKFPMLKKLDLAVANIIEEDLDVIGSLEELVI
    # Sbjct: 1141 PSLVELTIRDCPRLEVGPMIRSLPKFPMLKKLDLAVANIIEEDLDVIGSLEELVILSLKL 1200

    foreach $line( @data ) {
        next if $line =~ /^\s*$/;

        if( $line =~ /^ ?Score/ ) {
            $self->_set_score_stats( $line );
        } elsif( $line =~ /^ ?(Identities|Positives|Strand)/ ) {
            $self->_set_match_stats( $line );
        } elsif( $line =~ /^ ?Frame = ([\d+-]+)/ ) {
          # Version 2.0.8 has Frame information on a separate line.
          # Storing frame according to SeqFeature::Generic::frame()
          # which does not contain strand info (use strand()).
          my $frame = abs($1) - 1;
          $self->frame( $frame );
        } elsif( $line =~ /^(Query:?[\s\d]+)([^\s\d]+)/ ) {
            push @queryList, $line;
            $self->{'_match_indent'} = CORE::length $1;
            $aln_row_len = (CORE::length $1) + (CORE::length $2);
            $matchLine = 1;
        } elsif( $matchLine ) {
            # Pad the match line with spaces if necessary.
            $length_diff = $aln_row_len - CORE::length $line;
            $length_diff and $line .= ' 'x $length_diff;
            push @matchList, $line;
            $matchLine = 0;
        } elsif( $line =~ /^Sbjct/ ) {
            push @sbjctList, $line;
        }
    }
    # Storing the query and sbjct lists in case they are needed later.
    # We could make this conditional to save memory.
    $self->{'_queryList'} = \@queryList;
    $self->{'_sbjctList'} = \@sbjctList;

    # Storing the match list in case it is needed later.
    $self->{'_matchList'} = \@matchList;

    if(not defined ($self->{'_numIdentical'})) {
        my $id_str = $self->_id_str;
        $self->throw( -text  => "Can't parse match statistics. Possibly a new or unrecognized Blast format. ($id_str)");
    }

    if(!scalar @queryList or !scalar @sbjctList) {
        my $id_str = $self->_id_str;
        $self->throw( "Can't find query or sbjct alignment lines. Possibly unrecognized Blast format. ($id_str)");
    }
}


=head2 _set_score_stats

 Usage     : called automatically by _set_data()
 Purpose   : Sets various score statistics obtained from the HSP listing.
 Argument  : String with any of the following formats:
           : blast2:  Score = 30.1 bits (66), Expect = 9.2
           : blast2:  Score = 158.2 bits (544), Expect(2) = e-110
           : blast1:  Score = 410 (144.3 bits), Expect = 1.7e-40, P = 1.7e-40
           : blast1:  Score = 55 (19.4 bits), Expect = 5.3, Sum P(3) = 0.99
 Throws    : Exception if the stats cannot be parsed, probably due to a change
           : in the Blast report format.

See Also   : L</_set_data>

=cut

#--------------------
sub _set_score_stats {
#--------------------
    my ($self, $data) = @_;

    my ($expect, $p);

    if($data =~ /Score = +([\d.e+-]+) bits \(([\d.e+-]+)\), +Expect = +([\d.e+-]+)/) {
        # blast2 format n = 1
        $self->bits($1);
        $self->score($2);
        $expect            = $3;
    } elsif($data =~ /Score = +([\d.e+-]+) bits \(([\d.e+-]+)\), +Expect\((\d+)\) = +([\d.e+-]+)/) {
        # blast2 format n > 1
        $self->bits($1);
        $self->score($2);
        $self->{'_n'}      = $3;
        $expect            = $4;

    } elsif($data =~ /Score = +([\d.e+-]+) \(([\d.e+-]+) bits\), +Expect = +([\d.e+-]+), P = +([\d.e-]+)/) {
        # blast1 format, n = 1
        $self->score($1);
        $self->bits($2);
        $expect            = $3;
        $p                 = $4;

    } elsif($data =~ /Score = +([\d.e+-]+) \(([\d.e+-]+) bits\), +Expect = +([\d.e+-]+), +Sum P\((\d+)\) = +([\d.e-]+)/) {
        # blast1 format, n > 1
        $self->score($1);
        $self->bits($2);
        $expect            = $3;
        $self->{'_n'}      = $4;
        $p                 = $5;

    } else {
        my $id_str = $self->_id_str;
        $self->throw(-class => 'Bio::Root::Exception',
                     -text => "Can't parse score statistics: unrecognized format. ($id_str)",
                     -value => $data);
    }
    $expect = "1$expect" if $expect =~ /^e/i;
    $p      = "1$p"      if defined $p and $p=~ /^e/i;

    $self->{'_expect'} = $expect;
    $self->{'_p'}      = $p || undef;
    $self->significance( $p || $expect );
}


=head2 _set_match_stats

 Usage     : Private method; called automatically by _set_data()
 Purpose   : Sets various matching statistics obtained from the HSP listing.
 Argument  : blast2: Identities = 23/74 (31%), Positives = 29/74 (39%), Gaps = 17/74 (22%)
           : blast2: Identities = 57/98 (58%), Positives = 74/98 (75%)
           : blast1: Identities = 87/204 (42%), Positives = 126/204 (61%)
           : blast1: Identities = 87/204 (42%), Positives = 126/204 (61%), Frame = -3
           : WU-blast: Identities = 310/553 (56%), Positives = 310/553 (56%), Strand = Minus / Plus
 Throws    : Exception if the stats cannot be parsed, probably due to a change
           : in the Blast report format.
 Comments  : The "Gaps = " data in the HSP header has a different meaning depending
           : on the type of Blast: for BLASTP, this number is the total number of
           : gaps in query+sbjct; for TBLASTN, it is the number of gaps in the
           : query sequence only. Thus, it is safer to collect the data
           : separately by examining the actual sequence strings as is done
           : in _set_seq().

See Also   : L</_set_data>, L</_set_seq>

=cut

#--------------------
sub _set_match_stats {
#--------------------
    my ($self, $data) = @_;

    if($data =~ m!Identities = (\d+)/(\d+)!) {
      # blast1 or 2 format
      $self->{'_numIdentical'} = $1;
      $self->{'_totalLength'}  = $2;
    }

    if($data =~ m!Positives = (\d+)/(\d+)!) {
      # blast1 or 2 format
      $self->{'_numConserved'} = $1;
      $self->{'_totalLength'}  = $2;
    }

    if($data =~ m!Frame = ([\d+-]+)!) {
      $self->frame($1);
    }

    # Strand data is not always present in this line.
    # _set_seq() will also set strand information.
    if($data =~ m!Strand = (\w+) / (\w+)!) {
        $self->{'_queryStrand'} = $1;
        $self->{'_sbjctStrand'} = $2;
    }

#    if($data =~ m!Gaps = (\d+)/(\d+)!) {
#         $self->{'_totalGaps'} = $1;
#    } else {
#         $self->{'_totalGaps'} = 0;
#    }
}



=head2 _set_seq_data

 Usage     : called automatically when sequence data is requested.
 Purpose   : Sets the HSP sequence data for both query and sbjct sequences.
           : Includes: start, stop, length, gaps, and raw sequence.
 Argument  : n/a
 Throws    : Propagates any exception thrown by _set_match_seq()
 Comments  : Uses raw data stored by _set_data() during object construction.
           : These data are not always needed, so it is conditionally
           : executed only upon demand by methods such as gaps(), _set_residues(),
           : etc. _set_seq() does the dirty work.

See Also   : L</_set_seq>

=cut

#-----------------
sub _set_seq_data {
#-----------------
    my $self = shift;

    $self->_set_seq('query', @{$self->{'_queryList'}});
    $self->_set_seq('sbjct', @{$self->{'_sbjctList'}});

    # Liberate some memory.
    @{$self->{'_queryList'}} = @{$self->{'_sbjctList'}} = ();
    undef $self->{'_queryList'};
    undef $self->{'_sbjctList'};

    $self->{'_set_seq_data'} = 1;
}



=head2 _set_seq

 Usage     : called automatically by _set_seq_data()
           : $hsp_obj->($seq_type, @data);
 Purpose   : Sets sequence information for both the query and sbjct sequences.
           : Directly counts the number of gaps in each sequence (if gapped Blast).
 Argument  : $seq_type = 'query' or 'sbjct'
           : @data = all seq lines with the form:
           : Query: 61  SPHNVKDRKEQNGSINNAISPTATANTSGSQQINIDSALRDRSSNVAAQPSLSDASSGSN 120
 Throws    : Exception if data strings cannot be parsed, probably due to a change
           : in the Blast report format.
 Comments  : Uses first argument to determine which data members to set
           : making this method sensitive data member name changes.
           : Behavior is dependent on the type of BLAST analysis (TBLASTN, BLASTP, etc).
 Warning   : Sequence endpoints are normalized so that start < end. This affects HSPs
           : for TBLASTN/X hits on the minus strand. Normalization facilitates use
           : of range information by methods such as match().

See Also   : L</_set_seq_data>, L</matches>, L</range>, L</start>, L</end>

=cut

#-------------
sub _set_seq {
#-------------
    my $self      = shift;
    my $seqType   = shift;
    my @data      = @_;
    my @ranges    = ();
    my @sequence  = ();
    my $numGaps   = 0;

    foreach( @data ) {
        if( m/(\d+) *([^\d\s]+) *(\d+)/) {
            push @ranges, ( $1, $3 ) ;
            push @sequence, $2;
        #print STDERR "_set_seq found sequence \"$2\"\n";
        } else {
            $self->warn("Bad sequence data: $_");
        }
    }

    if( !(scalar(@sequence) and scalar(@ranges))) {
        my $id_str = $self->_id_str;
        $self->throw("Can't set sequence: missing data. Possibly unrecognized Blast format. ($id_str)");
   }

    # Sensitive to member name changes.
    $seqType = "_\L$seqType\E";
    $self->{$seqType.'Start'} = $ranges[0];
    $self->{$seqType.'Stop'}  = $ranges[ $#ranges ];
    $self->{$seqType.'Seq'}   = \@sequence;

    $self->{$seqType.'Length'} = abs($ranges[ $#ranges ] - $ranges[0]) + 1;

    # Adjust lengths for BLASTX, TBLASTN, TBLASTX sequences
    # Converting nucl coords to amino acid coords.

    my $prog = $self->algorithm;
    if($prog eq 'TBLASTN' and $seqType eq '_sbjct') {
        $self->{$seqType.'Length'} /= 3;
    } elsif($prog eq 'BLASTX' and $seqType eq '_query') {
        $self->{$seqType.'Length'} /= 3;
    } elsif($prog eq 'TBLASTX') {
        $self->{$seqType.'Length'} /= 3;
    }

    if( $prog ne 'BLASTP' ) {
        $self->{$seqType.'Strand'} = 'Plus' if $prog =~ /BLASTN/;
        $self->{$seqType.'Strand'} = 'Plus' if ($prog =~ /BLASTX/ and $seqType eq '_query');
        # Normalize sequence endpoints so that start < end.
        # Reverse complement or 'minus strand' HSPs get flipped here.
        if($self->{$seqType.'Start'} > $self->{$seqType.'Stop'}) {
            ($self->{$seqType.'Start'}, $self->{$seqType.'Stop'}) =
                ($self->{$seqType.'Stop'}, $self->{$seqType.'Start'});
            $self->{$seqType.'Strand'} = 'Minus';
        }
    }

    ## Count number of gaps in each seq. Only need to do this for gapped Blasts.
#    if($self->{'_gapped'}) {
        my $seqstr = join('', @sequence);
        $seqstr =~ s/\s//g;
        my $num_gaps = CORE::length($seqstr) - $self->{$seqType.'Length'};
        $self->{$seqType.'Gaps'} = $num_gaps if $num_gaps > 0;
#    }
}


=head2 _set_residues

 Usage     : called automatically when residue data is requested.
 Purpose   : Sets the residue numbers representing the identical and
           : conserved positions. These data are obtained by analyzing the
           : symbols between query and sbjct lines of the alignments.
 Argument  : n/a
 Throws    : Propagates any exception thrown by _set_seq_data() and _set_match_seq().
 Comments  : These data are not always needed, so it is conditionally
           : executed only upon demand by methods such as seq_inds().
           : Behavior is dependent on the type of BLAST analysis (TBLASTN, BLASTP, etc).

See Also   : L</_set_seq_data>, L</_set_match_seq>, L</seq_inds>

=cut

#------------------
sub _set_residues {
#------------------
    my $self      = shift;
    my @sequence  = ();

    $self->_set_seq_data() unless $self->{'_set_seq_data'};

    # Using hashes to avoid saving duplicate residue numbers.
    my %identicalList_query = ();
    my %identicalList_sbjct = ();
    my %conservedList_query = ();
    my %conservedList_sbjct = ();

    my $aref = $self->_set_match_seq() if not ref $self->{'_matchSeq'};
    $aref  ||= $self->{'_matchSeq'};
    my $seqString = join('', @$aref );

    my $qseq = join('',@{$self->{'_querySeq'}});
    my $sseq = join('',@{$self->{'_sbjctSeq'}});
    my $resCount_query = $self->{'_queryStop'} || 0;
    my $resCount_sbjct = $self->{'_sbjctStop'} || 0;

    my $prog = $self->algorithm;
    if($prog !~ /^BLASTP|^BLASTN/) {
        if($prog eq 'TBLASTN') {
            $resCount_sbjct /= 3;
        } elsif($prog eq 'BLASTX') {
            $resCount_query /= 3;
        } elsif($prog eq 'TBLASTX') {
            $resCount_query /= 3;
            $resCount_sbjct /= 3;
        }
    }

    my ($mchar, $schar, $qchar);
    while( $mchar = chop($seqString) ) {
        ($qchar, $schar) = (chop($qseq), chop($sseq));
        if( $mchar eq '+' ) {
            $conservedList_query{ $resCount_query } = 1;
            $conservedList_sbjct{ $resCount_sbjct } = 1;
        } elsif( $mchar ne ' ' ) {
            $identicalList_query{ $resCount_query } = 1;
            $identicalList_sbjct{ $resCount_sbjct } = 1;
        }
        $resCount_query-- if $qchar ne $GAP_SYMBOL;
        $resCount_sbjct-- if $schar ne $GAP_SYMBOL;
    }
    $self->{'_identicalRes_query'} = \%identicalList_query;
    $self->{'_conservedRes_query'} = \%conservedList_query;
    $self->{'_identicalRes_sbjct'} = \%identicalList_sbjct;
    $self->{'_conservedRes_sbjct'} = \%conservedList_sbjct;

}




=head2 _set_match_seq

 Usage     : $hsp_obj->_set_match_seq()
 Purpose   : Set the 'match' sequence for the current HSP (symbols in between
           : the query and sbjct lines.)
 Returns   : Array reference holding the match sequences lines.
 Argument  : n/a
 Throws    : Exception if the _matchList field is not set.
 Comments  : The match information is not always necessary. This method
           : allows it to be conditionally prepared.
           : Called by _set_residues>() and seq_str().

See Also   : L</_set_residues>, L</seq_str>

=cut

#-------------------
sub _set_match_seq {
#-------------------
    my $self = shift;

    if( ! ref($self->{'_matchList'}) ) {
        my $id_str = $self->_id_str;
        $self->throw("Can't set HSP match sequence: No data ($id_str)");
    }

    my @data = @{$self->{'_matchList'}};

    my(@sequence);
    foreach( @data ) {
        chomp($_);
        ## Remove leading spaces; (note: aln may begin with a space
        ## which is why we can't use s/^ +//).
        s/^ {$self->{'_match_indent'}}//;
        push @sequence, $_;
    }
    # Liberate some memory.
    @{$self->{'_matchList'}} = undef;
    $self->{'_matchList'} = undef;

    $self->{'_matchSeq'} = \@sequence;

    return $self->{'_matchSeq'};
}


=head2 n

 Usage     : $hsp_obj->n()
 Purpose   : Get the N value (num HSPs on which P/Expect is based).
           : This value is not defined with NCBI Blast2 with gapping.
 Returns   : Integer or null string if not defined.
 Argument  : n/a
 Throws    : n/a
 Comments  : The 'N' value is listed in parenthesis with P/Expect value:
           : e.g., P(3) = 1.2e-30  ---> (N = 3).
           : Not defined in NCBI Blast2 with gaps.
           : This typically is equal to the number of HSPs but not always.
           : To obtain the number of HSPs, use Bio::Search::Hit::BlastHit::num_hsps().

See Also   : L<Bio::SeqFeature::SimilarityPair::score()|Bio::SeqFeature::SimilarityPair>

=cut

#-----
sub n { my $self = shift; $self->{'_n'} || ''; }
#-----


=head2 matches

 Usage     : $hsp->matches([seq_type], [start], [stop]);
 Purpose   : Get the total number of identical and conservative matches
           : in the query or sbjct sequence for the given HSP. Optionally can
           : report data within a defined interval along the seq.
           : (Note: 'conservative' matches are called 'positives' in the
           : Blast report.)
 Example   : ($id,$cons) = $hsp_object->matches('hit');
           : ($id,$cons) = $hsp_object->matches('query',300,400);
 Returns   : 2-element array of integers
 Argument  : (1) seq_type = 'query' or 'hit' or 'sbjct' (default = query)
           :  ('sbjct' is synonymous with 'hit')
           : (2) start = Starting coordinate (optional)
           : (3) stop  = Ending coordinate (optional)
 Throws    : Exception if the supplied coordinates are out of range.
 Comments  : Relies on seq_str('match') to get the string of alignment symbols
           : between the query and sbjct lines which are used for determining
           : the number of identical and conservative matches.

See Also   : L</length>, L</gaps>, L</seq_str>, L<Bio::Search::Hit::BlastHit::_adjust_contigs()|Bio::Search::Hit::BlastHit>

=cut

#-----------
sub matches {
#-----------
    my( $self, %param ) = @_;
    my(@data);
    my($seqType, $beg, $end) = ($param{-SEQ}, $param{-START}, $param{-STOP});
    $seqType ||= 'query';
    $seqType = 'sbjct' if $seqType eq 'hit';

    my($start,$stop);

    if(!defined $beg && !defined $end) {
        ## Get data for the whole alignment.
        push @data, ($self->{'_numIdentical'}, $self->{'_numConserved'});
    } else {
        ## Get the substring representing the desired sub-section of aln.
        $beg ||= 0;
        $end ||= 0;
        ($start,$stop) = $self->range($seqType);
        if($beg == 0) { $beg = $start; $end = $beg+$end; }
        elsif($end == 0) { $end = $stop; $beg = $end-$beg; }

        if($end >= $stop) { $end = $stop; } ##ML changed from if (end >stop)
        else { $end += 1;}   ##ML moved from commented position below, makes
                             ##more sense here
#        if($end > $stop) { $end = $stop; }
        if($beg < $start) { $beg = $start; }
#        else { $end += 1;}

#        my $seq = substr($self->seq_str('match'), $beg-$start, ($end-$beg));

        ## ML: START fix for substr out of range error ------------------
        my $seq = "";
        my $prog = $self->algorithm;
        if (($prog eq 'TBLASTN') and ($seqType eq 'sbjct'))
        {
            $seq = substr($self->seq_str('match'),
                          int(($beg-$start)/3), int(($end-$beg+1)/3));

        } elsif (($prog eq 'BLASTX') and ($seqType eq 'query'))
        {
            $seq = substr($self->seq_str('match'),
                          int(($beg-$start)/3), int(($end-$beg+1)/3));
        } else {
            $seq = substr($self->seq_str('match'),
                          $beg-$start, ($end-$beg));
        }
        ## ML: End of fix for  substr out of range error -----------------


        ## ML: debugging code
        ## This is where we get our exception.  Try printing out the values going
        ## into this:
        ##
#         print STDERR
#             qq(*------------MY EXCEPTION --------------------\nSeq: ") ,
#             $self->seq_str("$seqType"), qq("\n),$self->rank,",(  index:";
#         print STDERR  $beg-$start, ", len: ", $end-$beg," ), (HSPRealLen:",
#             CORE::length $self->seq_str("$seqType");
#         print STDERR ", HSPCalcLen: ", $stop - $start +1 ," ),
#             ( beg: $beg, end: $end ), ( start: $start, stop: stop )\n";
         ## ML: END DEBUGGING CODE----------

        if(!CORE::length $seq) {
            my $id_str = $self->_id_str;
            $self->throw("Undefined $seqType sub-sequence ($beg,$end). Valid range = $start - $stop ($id_str)");
        }
        ## Get data for a substring.
#        printf "Collecting HSP subsection data: beg,end = %d,%d; start,stop = %d,%d\n%s<---\n", $beg, $end, $start, $stop, $seq;
#        printf "Original match seq:\n%s\n",$self->seq_str('match');
        $seq =~ s/ //g;  # remove space (no info).
        my $len_cons = CORE::length $seq;
        $seq =~ s/\+//g;  # remove '+' characters (conservative substitutions)
        my $len_id = CORE::length $seq;
        push @data, ($len_id, $len_cons);
#        printf "  HSP = %s\n  id = %d; cons = %d\n", $self->rank, $len_id, $len_cons; <STDIN>;
    }
    @data;
}


=head2 num_identical

 Usage     : $hsp_object->num_identical();
 Purpose   : Get the number of identical positions within the given HSP.
 Example   : $num_iden = $hsp_object->num_identical();
 Returns   : integer
 Argument  : n/a
 Throws    : n/a

See Also   : L</num_conserved>, L</frac_identical>

=cut

#-------------------
sub num_identical {
#-------------------
    my( $self) = shift;

    $self->{'_numIdentical'};
}


=head2 num_conserved

 Usage     : $hsp_object->num_conserved();
 Purpose   : Get the number of conserved positions within the given HSP.
 Example   : $num_iden = $hsp_object->num_conserved();
 Returns   : integer
 Argument  : n/a
 Throws    : n/a

See Also   : L</num_identical>, L</frac_conserved>

=cut

#-------------------
sub num_conserved {
#-------------------
    my( $self) = shift;

    $self->{'_numConserved'};
}



=head2 range

 Usage     : $hsp->range( [seq_type] );
 Purpose   : Gets the (start, end) coordinates for the query or sbjct sequence
           : in the HSP alignment.
 Example   : ($query_beg, $query_end) = $hsp->range('query');
           : ($hit_beg, $hit_end) = $hsp->range('hit');
 Returns   : Two-element array of integers
 Argument  : seq_type = string, 'query' or 'hit' or 'sbjct'  (default = 'query')
           :  ('sbjct' is synonymous with 'hit')
 Throws    : n/a

See Also   : L</start>, L</end>

=cut

#----------
sub range {
#----------
    my ($self, $seqType) = @_;

    $self->_set_seq_data() unless $self->{'_set_seq_data'};

    $seqType ||= 'query';
    $seqType = 'sbjct' if $seqType eq 'hit';

    ## Sensitive to member name changes.
    $seqType = "_\L$seqType\E";

    return ($self->{$seqType.'Start'},$self->{$seqType.'Stop'});
}

=head2 start

 Usage     : $hsp->start( [seq_type] );
 Purpose   : Gets the start coordinate for the query, sbjct, or both sequences
           : in the HSP alignment.
           : NOTE: Start will always be less than end.
           : To determine strand, use $hsp->strand()
 Example   : $query_beg = $hsp->start('query');
           : $hit_beg = $hsp->start('hit');
           : ($query_beg, $hit_beg) = $hsp->start();
 Returns   : scalar context: integer
           : array context without args: list of two integers
 Argument  : In scalar context: seq_type = 'query' or 'hit' or 'sbjct' (default= 'query')
           :  ('sbjct' is synonymous with 'hit')
           : Array context can be "induced" by providing an argument of 'list' or 'array'.
 Throws    : n/a

See Also   : L</end>, L</range>

=cut

#----------
sub start {
#----------
    my ($self, $seqType) = @_;

    $seqType ||= (wantarray ? 'list' : 'query');
    $seqType = 'sbjct' if $seqType eq 'hit';

    $self->_set_seq_data() unless $self->{'_set_seq_data'};

    if($seqType =~ /list|array/i) {
        return ($self->{'_queryStart'}, $self->{'_sbjctStart'});
    } else {
        ## Sensitive to member name changes.
        $seqType = "_\L$seqType\E";
        return $self->{$seqType.'Start'};
    }
}

=head2 end

 Usage     : $hsp->end( [seq_type] );
 Purpose   : Gets the end coordinate for the query, sbjct, or both sequences
           : in the HSP alignment.
           : NOTE: Start will always be less than end.
           : To determine strand, use $hsp->strand()
 Example   : $query_end = $hsp->end('query');
           : $hit_end = $hsp->end('hit');
           : ($query_end, $hit_end) = $hsp->end();
 Returns   : scalar context: integer
           : array context without args: list of two integers
 Argument  : In scalar context: seq_type = 'query' or 'hit' or 'sbjct' (default= 'query')
           :  ('sbjct' is synonymous with 'hit')
           : Array context can be "induced" by providing an argument of 'list' or 'array'.
 Throws    : n/a

See Also   : L</start>, L</range>, L</strand>

=cut

#----------
sub end {
#----------
    my ($self, $seqType) = @_;

    $seqType ||= (wantarray ? 'list' : 'query');
    $seqType = 'sbjct' if $seqType eq 'hit';

    $self->_set_seq_data() unless $self->{'_set_seq_data'};

    if($seqType =~ /list|array/i) {
        return ($self->{'_queryStop'}, $self->{'_sbjctStop'});
    } else {
        ## Sensitive to member name changes.
        $seqType = "_\L$seqType\E";
        return $self->{$seqType.'Stop'};
    }
}



=head2 strand

 Usage     : $hsp_object->strand( [seq_type] )
 Purpose   : Get the strand of the query or sbjct sequence.
 Example   : print $hsp->strand('query');
           : ($query_strand, $hit_strand) = $hsp->strand();
 Returns   : -1, 0, or 1
           : -1 = Minus strand, +1 = Plus strand
           : Returns 0 if strand is not defined, which occurs
           : for BLASTP reports, and the query of TBLASTN
           : as well as the hit if BLASTX reports.
           : In scalar context without arguments, returns queryStrand value.
           : In array context without arguments, returns a two-element list
           :    of strings (queryStrand, sbjctStrand).
           : Array context can be "induced" by providing an argument of 'list' or 'array'.
 Argument  : seq_type: 'query' or 'hit' or 'sbjct' or undef
           :  ('sbjct' is synonymous with 'hit')
 Throws    : n/a

See Also   : L</_set_seq>, L</_set_match_stats>

=cut

#-----------
sub strand {
#-----------
    my( $self, $seqType ) = @_;

    $seqType  ||= (wantarray ? 'list' : 'query');
    $seqType = 'sbjct' if $seqType eq 'hit';

    ## Sensitive to member name format.
    $seqType = "_\L$seqType\E";

    # $seqType could be '_list'.
    $self->{'_queryStrand'} or $self->_set_seq_data() unless $self->{'_set_seq_data'};

    my $prog = $self->algorithm;

    if($seqType  =~ /list|array/i) {
        my ($qstr, $hstr);
        if( $prog eq 'BLASTP') {
            $qstr = 0;
            $hstr = 0;
        }
        elsif( $prog eq 'TBLASTN') {
            $qstr = 0;
            $hstr = $STRAND_SYMBOL{$self->{'_sbjctStrand'}};
        }
        elsif( $prog eq 'BLASTX') {
            $qstr = $STRAND_SYMBOL{$self->{'_queryStrand'}};
            $hstr = 0;
        }
        else {
            $qstr = $STRAND_SYMBOL{$self->{'_queryStrand'}} if defined $self->{'_queryStrand'};
            $hstr = $STRAND_SYMBOL{$self->{'_sbjctStrand'}} if defined $self->{'_sbjctStrand'};
        }
        $qstr ||= 0;
        $hstr ||= 0;
        return ($qstr, $hstr);
    }
    local $^W = 0;
    $STRAND_SYMBOL{$self->{$seqType.'Strand'}} || 0;
}


=head2 seq

 Usage     : $hsp->seq( [seq_type] );
 Purpose   : Get the query or sbjct sequence as a Bio::Seq.pm object.
 Example   : $seqObj = $hsp->seq('query');
 Returns   : Object reference for a Bio::Seq.pm object.
 Argument  : seq_type = 'query' or 'hit' or 'sbjct' (default = 'query').
           :  ('sbjct' is synonymous with 'hit')
 Throws    : Propagates any exception that occurs during construction
           : of the Bio::Seq.pm object.
 Comments  : The sequence is returned in an array of strings corresponding
           : to the strings in the original format of the Blast alignment.
           : (i.e., same spacing).

See Also   : L</seq_str>, L</seq_inds>, L<Bio::Seq>

=cut

#-------
sub seq {
#-------
    my($self,$seqType) = @_;
    $seqType ||= 'query';
    $seqType = 'sbjct' if $seqType eq 'hit';
    my $str = $self->seq_str($seqType);

    require Bio::Seq;

    Bio::Seq->new(-ID   => $self->to_string,
                  -SEQ  => $str,
                  -DESC => "$seqType sequence",
                  );
}

=head2 seq_str

 Usage     : $hsp->seq_str( seq_type );
 Purpose   : Get the full query, sbjct, or 'match' sequence as a string.
           : The 'match' sequence is the string of symbols in between the
           : query and sbjct sequences.
 Example   : $str = $hsp->seq_str('query');
 Returns   : String
 Argument  : seq_Type = 'query' or 'hit' or 'sbjct' or 'match'
           :  ('sbjct' is synonymous with 'hit')
 Throws    : Exception if the argument does not match an accepted seq_type.
 Comments  : Calls _set_seq_data() to set the 'match' sequence if it has
           : not been set already.

See Also   : L</seq>, L</seq_inds>, L</_set_match_seq>

=cut

#------------
sub seq_str {
#------------
    my($self,$seqType) = @_;

    $seqType ||= 'query';
    $seqType = 'sbjct' if $seqType eq 'hit';
    ## Sensitive to member name changes.
    $seqType = "_\L$seqType\E";

    $self->_set_seq_data() unless $self->{'_set_seq_data'};

    if($seqType =~ /sbjct|query/) {
        my $seq = join('',@{$self->{$seqType.'Seq'}});
        $seq =~ s/\s+//g;
        return $seq;

    } elsif( $seqType =~ /match/i) {
        # Only need to call _set_match_seq() if the match seq is requested.
        my $aref = $self->_set_match_seq() unless ref $self->{'_matchSeq'};
        $aref =  $self->{'_matchSeq'};

        return join('',@$aref);

    } else {
        my $id_str = $self->_id_str;
        $self->throw(-class => 'Bio::Root::BadParameter',
                     -text => "Invalid or undefined sequence type: $seqType ($id_str)\n" .
                               "Valid types: query, sbjct, match",
                     -value => $seqType);
    }
}

=head2 seq_inds

 Usage     : $hsp->seq_inds( seq_type, class, collapse );
 Purpose   : Get a list of residue positions (indices) for all identical
           : or conserved residues in the query or sbjct sequence.
 Example   : @s_ind = $hsp->seq_inds('query', 'identical');
           : @h_ind = $hsp->seq_inds('hit', 'conserved');
           : @h_ind = $hsp->seq_inds('hit', 'conserved', 1);
 Returns   : List of integers
           : May include ranges if collapse is true.
 Argument  : seq_type  = 'query' or 'hit' or 'sbjct'  (default = query)
           :  ('sbjct' is synonymous with 'hit')
           : class     = 'identical' or 'conserved' (default = identical)
           :              (can be shortened to 'id' or 'cons')
           :              (actually, anything not 'id' will evaluate to 'conserved').
           : collapse  = boolean, if true, consecutive positions are merged
           :             using a range notation, e.g., "1 2 3 4 5 7 9 10 11"
           :             collapses to "1-5 7 9-11". This is useful for
           :             consolidating long lists. Default = no collapse.
 Throws    : n/a.
 Comments  : Calls _set_residues() to set the 'match' sequence if it has
           : not been set already.

See Also   : L</seq>, L</_set_residues>, L<Bio::Search::BlastUtils::collapse_nums()|Bio::Search::BlastUtils>, L<Bio::Search::Hit::BlastHit::seq_inds()|Bio::Search::Hit::BlastHit>

=cut

#---------------
sub seq_inds {
#---------------
    my ($self, $seqType, $class, $collapse) = @_;

    $seqType  ||= 'query';
    $class ||= 'identical';
    $collapse ||= 0;
    $seqType = 'sbjct' if $seqType eq 'hit';

    $self->_set_residues() unless defined $self->{'_identicalRes_query'};

    $seqType  = ($seqType !~ /^q/i ? 'sbjct' : 'query');
    $class = ($class !~ /^id/i ? 'conserved' : 'identical');

    ## Sensitive to member name changes.
    $seqType  = "_\L$seqType\E";
    $class = "_\L$class\E";

    my @ary = sort { $a <=> $b } keys %{ $self->{"${class}Res$seqType"}};

    require Bio::Search::BlastUtils if $collapse;

    return $collapse ? &Bio::Search::BlastUtils::collapse_nums(@ary) : @ary;
}


=head2 get_aln

 Usage     : $hsp->get_aln()
 Purpose   : Get a Bio::SimpleAlign object constructed from the query + sbjct
           : sequences of the present HSP object.
 Example   : $aln_obj = $hsp->get_aln();
 Returns   : Object reference for a Bio::SimpleAlign.pm object.
 Argument  : n/a.
 Throws    : Propagates any exception ocurring during the construction of
           : the Bio::SimpleAlign object.
 Comments  : Requires Bio::SimpleAlign.
           : The Bio::SimpleAlign object is constructed from the query + sbjct
           : sequence objects obtained by calling seq().
           : Gap residues are included (see $GAP_SYMBOL).

See Also   : L</seq>, L<Bio::SimpleAlign>

=cut

#------------
sub get_aln {
#------------
    my $self = shift;

    require Bio::SimpleAlign;
    require Bio::LocatableSeq;
    my $qseq = $self->seq('query');
    my $sseq = $self->seq('sbjct');

    my $type = $self->algorithm =~ /P$|^T/ ? 'amino' : 'dna';
    my $aln = Bio::SimpleAlign->new();
    $aln->add_seq(Bio::LocatableSeq->new(-seq => $qseq->seq(),
                                        -id  => 'query_'.$qseq->display_id(),
                                        -start => 1,
                                        -end   => CORE::length($qseq)));

    $aln->add_seq(Bio::LocatableSeq->new(-seq => $sseq->seq(),
                                        -id  => 'hit_'.$sseq->display_id(),
                                        -start => 1,
                                        -end   => CORE::length($sseq)));

    return $aln;
}


1;
__END__


=head1 FOR DEVELOPERS ONLY

=head2 Data Members

Information about the various data members of this module is provided for those
wishing to modify or understand the code. Two things to bear in mind:

=over 4

=item 1 Do NOT rely on these in any code outside of this module.

All data members are prefixed with an underscore to signify that they are private.
Always use accessor methods. If the accessor doesn't exist or is inadequate,
create or modify an accessor (and let me know, too!).

=item 2 This documentation may be incomplete and out of date.

It is easy for these data member descriptions to become obsolete as
this module is still evolving. Always double check this info and search
for members not described here.

=back

An instance of Bio::Search::HSP::PsiBlastHSP.pm is a blessed reference to a hash containing
all or some of the following fields:

 FIELD           VALUE
 --------------------------------------------------------------
 (member names are mostly self-explanatory)

 _score              :
 _bits               :
 _p                  :
 _n                  : Integer. The 'N' value listed in parenthesis with P/Expect value:
                     : e.g., P(3) = 1.2e-30  ---> (N = 3).
                     : Not defined in NCBI Blast2 with gaps.
                     : To obtain the number of HSPs, use Bio::Search::Hit::BlastHit::num_hsps().
 _expect             :
 _queryLength        :
 _queryGaps          :
 _queryStart         :
 _queryStop          :
 _querySeq           :
 _sbjctLength        :
 _sbjctGaps          :
 _sbjctStart         :
 _sbjctStop          :
 _sbjctSeq           :
 _matchSeq           : String. Contains the symbols between the query and sbjct lines
                       which indicate identical (letter) and conserved ('+') matches
                       or a mismatch (' ').
 _numIdentical       :
 _numConserved       :
 _identicalRes_query :
 _identicalRes_sbjct :
 _conservedRes_query :
 _conservedRes_sbjct :
 _match_indent       : The number of leading space characters on each line containing
                       the match symbols. _match_indent is 13 in this example:
                         Query:   285 QNSAPWGLARISHRERLNLGSFNKYLYDDDAG
                                      Q +APWGLARIS       G+ + Y YD+ AG
                         ^^^^^^^^^^^^^


=cut

1;

