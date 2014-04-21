#
# BioPerl module for Bio::Search::HSP::HMMERHSP
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::HSP::HMMERHSP - A HSP object for HMMER results

=head1 SYNOPSIS

    use Bio::Search::HSP::HMMERHSP;
    # use it just like a Bio::Search::HSP::GenericHSP object

=head1 DESCRIPTION

This object is a specialization of L<Bio::Search::HSP::GenericHSP>.

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

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Search::HSP::HMMERHSP;
use strict;

use base qw(Bio::Search::HSP::GenericHSP);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::HSP::HMMERHSP->new();
 Function: Builds a new Bio::Search::HSP::HMMERHSP object 
 Returns : Bio::Search::HSP::HMMERHSP
 Args    :

Plus Bio::Search::HSP::GenericHSP methods

    -algorithm    => algorithm used (BLASTP, TBLASTX, FASTX, etc)
    -evalue       => evalue
    -pvalue       => pvalue
    -bits         => bit value for HSP
    -score        => score value for HSP (typically z-score but depends on
					           analysis)
    -hsp_length   => Length of the HSP (including gaps)
    -identical    => # of residues that that matched identically
    -conserved    => # of residues that matched conservatively 
                     (only protein comparisions - 
                     conserved == identical in nucleotide comparisons)
    -hsp_gaps     => # of gaps in the HSP
    -query_gaps   => # of gaps in the query in the alignment
    -hit_gaps     => # of gaps in the subject in the alignment    
    -query_name   => HSP Query sequence name (if available)
    -query_start  => HSP Query start (in original query sequence coords)
    -query_end    => HSP Query end (in original query sequence coords)
    -hit_name     => HSP Hit sequence name (if available)
    -hit_start    => HSP Hit start (in original hit sequence coords)
    -hit_end      => HSP Hit end (in original hit sequence coords)
    -hit_length   => total length of the hit sequence
    -query_length => total length of the query sequence
    -query_seq    => query sequence portion of the HSP
    -hit_seq      => hit sequence portion of the HSP
    -homology_seq => homology sequence for the HSP
    -hit_frame    => hit frame (only if hit is translated protein)
    -query_frame  => query frame (only if query is translated protein)

=cut

=head2 Bio::Search::HSP::HSPI methods

Implementation of Bio::Search::HSP::HSPI methods follow

=head2 algorithm

 Title   : algorithm
 Usage   : my $r_type = $hsp->algorithm
 Function: Obtain the name of the algorithm used to obtain the HSP
 Returns : string (e.g., BLASTP)
 Args    : [optional] scalar string to set value

=cut

=head2 pvalue

 Title   : pvalue
 Usage   : my $pvalue = $hsp->pvalue();
 Function: Returns the P-value for this HSP or undef 
 Returns : float or exponential (2e-10)
           P-value is not defined with NCBI Blast2 reports.
 Args    : [optional] numeric to set value

=cut

=head2 evalue

 Title   : evalue
 Usage   : my $evalue = $hsp->evalue();
 Function: Returns the e-value for this HSP
 Returns : float or exponential (2e-10)
 Args    : [optional] numeric to set value

=cut

=head2 frac_identical

 Title   : frac_identical
 Usage   : my $frac_id = $hsp->frac_identical( ['query'|'hit'|'total'] );
 Function: Returns the fraction of identitical positions for this HSP 
 Returns : Float in range 0.0 -> 1.0
 Args    : arg 1:  'query' = num identical / length of query seq (without gaps)
                   'hit'   = num identical / length of hit seq (without gaps)
                   'total' = num identical / length of alignment (with gaps)
                   default = 'total' 
           arg 2: [optional] frac identical value to set for the type requested

=cut

=head2 frac_conserved

 Title    : frac_conserved
 Usage    : my $frac_cons = $hsp->frac_conserved( ['query'|'hit'|'total'] );
 Function : Returns the fraction of conserved positions for this HSP.
            This is the fraction of symbols in the alignment with a 
            positive score.
 Returns : Float in range 0.0 -> 1.0
 Args    : arg 1: 'query' = num conserved / length of query seq (without gaps)
                  'hit'   = num conserved / length of hit seq (without gaps)
                  'total' = num conserved / length of alignment (with gaps)
                  default = 'total' 
           arg 2: [optional] frac conserved value to set for the type requested

=cut

=head2 gaps

 Title    : gaps
 Usage    : my $gaps = $hsp->gaps( ['query'|'hit'|'total'] );
 Function : Get the number of gaps in the query, hit, or total alignment.
 Returns  : Integer, number of gaps or 0 if none
 Args     : arg 1: 'query' = num gaps in query seq
                   'hit'   = num gaps in hit seq
                   'total' = num gaps in whole alignment 
                   default = 'total' 
            arg 2: [optional] integer gap value to set for the type requested

=cut

=head2 query_string

 Title   : query_string
 Usage   : my $qseq = $hsp->query_string;
 Function: Retrieves the query sequence of this HSP as a string
 Returns : string
 Args    : [optional] string to set for query sequence


=cut

=head2 hit_string

 Title   : hit_string
 Usage   : my $hseq = $hsp->hit_string;
 Function: Retrieves the hit sequence of this HSP as a string
 Returns : string
 Args    : [optional] string to set for hit sequence


=cut


=head2 homology_string

 Title   : homology_string
 Usage   : my $homo_string = $hsp->homology_string;
 Function: Retrieves the homology sequence for this HSP as a string.
         : The homology sequence is the string of symbols in between the 
         : query and hit sequences in the alignment indicating the degree
         : of conservation (e.g., identical, similar, not similar).
 Returns : string
 Args    : [optional] string to set for homology sequence

=cut

=head2 length

 Title    : length
 Usage    : my $len = $hsp->length( ['query'|'hit'|'total'] );
 Function : Returns the length of the query or hit in the alignment 
            (without gaps) 
            or the aggregate length of the HSP (including gaps;
            this may be greater than either hit or query )
 Returns  : integer
 Args     : arg 1: 'query' = length of query seq (without gaps)
                   'hit'   = length of hit seq (without gaps)
                   'total' = length of alignment (with gaps)
                   default = 'total' 
            arg 2: [optional] integer length value to set for specific type

=cut

=head2 percent_identity

 Title   : percent_identity
 Usage   : my $percentid = $hsp->percent_identity()
 Function: Returns the calculated percent identity for an HSP
 Returns : floating point between 0 and 100 
 Args    : none


=cut


=head2 frame

 Title   : frame
 Usage   : my ($qframe, $hframe) = $hsp->frame('list',$queryframe,$subjectframe)
 Function: Set the Frame for both query and subject and insure that
           they agree.
           This overrides the frame() method implementation in
           FeaturePair.
 Returns : array of query and subject frame if return type wants an array
           or query frame if defined or subject frame if not defined
 Args    : 'hit' or 'subject' or 'sbjct' to retrieve the frame of the subject (default)
           'query' to retrieve the query frame 
           'list' or 'array' to retrieve both query and hit frames together
 Note    : Frames are stored in the GFF way (0-2) not 1-3
           as they are in BLAST (negative frames are deduced by checking
                                 the strand of the query or hit)

=cut

=head2 get_aln

 Title   : get_aln
 Usage   : my $aln = $hsp->gel_aln
 Function: Returns a Bio::SimpleAlign representing the HSP alignment
 Returns : Bio::SimpleAlign
 Args    : none

=cut

sub get_aln {
    my ($self) = shift;
    $self->warn("Inappropriate to build a Bio::SimpleAlign from a HMMER HSP object");
    return;
}

=head2 num_conserved

 Title   : num_conserved
 Usage   : $obj->num_conserved($newval)
 Function: returns the number of conserved residues in the alignment
 Returns : inetger
 Args    : integer (optional)


=cut

=head2 num_identical

 Title   : num_identical
 Usage   : $obj->num_identical($newval)
 Function: returns the number of identical residues in the alignment
 Returns : integer
 Args    : integer (optional)


=cut

=head2 seq_inds

 Title   : seq_inds
 Purpose   : Get a list of residue positions (indices) for all identical 
           : or conserved residues in the query or sbjct sequence.
 Example   : @s_ind = $hsp->seq_inds('query', 'identical');
           : @h_ind = $hsp->seq_inds('hit', 'conserved');
           : @h_ind = $hsp->seq_inds('hit', 'conserved', 1);
 Returns   : List of integers 
           : May include ranges if collapse is true.
 Argument  : seq_type  = 'query' or 'hit' or 'sbjct'  (default = query)
           :  ('sbjct' is synonymous with 'hit') 
           : class     = 'identical' or 'conserved' or 'nomatch' or 'gap'
           :              (default = identical)
           :              (can be shortened to 'id' or 'cons')
           :              
           : collapse  = boolean, if true, consecutive positions are merged
           :             using a range notation, e.g., "1 2 3 4 5 7 9 10 11" 
           :             collapses to "1-5 7 9-11". This is useful for 
           :             consolidating long lists. Default = no collapse.
 Throws    : n/a.
 Comments  : 

See Also   : L<Bio::Search::BlastUtils::collapse_nums()|Bio::Search::BlastUtils>, 
L<Bio::Search::Hit::HitI::seq_inds()|Bio::Search::Hit::HitI>

=cut

=head2 Inherited from Bio::SeqFeature::SimilarityPair

These methods come from Bio::SeqFeature::SimilarityPair

=head2 query

 Title   : query
 Usage   : my $query = $hsp->query
 Function: Returns a SeqFeature representing the query in the HSP
 Returns : Bio::SeqFeature::Similarity
 Args    : [optional] new value to set


=head2 hit

 Title   : hit
 Usage   : my $hit = $hsp->hit
 Function: Returns a SeqFeature representing the hit in the HSP
 Returns : Bio::SeqFeature::Similarity
 Args    : [optional] new value to set


=head2 significance

 Title   : significance
 Usage   : $evalue = $obj->significance();
           $obj->significance($evalue);
 Function: Get/Set the significance value
 Returns : numeric
 Args    : [optional] new value to set


=head2 score

 Title   : score
 Usage   : my $score = $hsp->score();
 Function: Returns the score for this HSP or undef 
 Returns : numeric           
 Args    : [optional] numeric to set value

=cut 

=head2 bits

 Title   : bits
 Usage   : my $bits = $hsp->bits();
 Function: Returns the bit value for this HSP or undef 
 Returns : numeric
 Args    : none

=cut


1;
