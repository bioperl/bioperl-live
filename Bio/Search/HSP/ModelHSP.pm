#
# BioPerl module for Bio::Search::HSP::ModelHSP
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chris Fields <cjfields at bioperl dot org>
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::HSP::ModelHSP - A HSP object for model-based searches

=head1 SYNOPSIS

    use Bio::Search::HSP::ModelHSP;
    # us it just like a Bio::Search::HSP::ModelHSP object

=head1 DESCRIPTION

This object is a specialization of L<Bio::Search::HSP::ModelHSP> and is used
for searches which involve a query model, such as a Hidden Markov Model (HMM),
covariance model (CM), descriptor, or anything else besides a sequence. Note
that results from any HSPI class methods which rely on the query being a
sequence are unreliable and have thus been overridden with warnings indicating
they have not been implemented at this time.

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

=head1 AUTHOR - Chris Fields

Email cjfields at bioperl dot org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Search::HSP::ModelHSP;
use strict;
use Bio::Seq::Meta;

use base qw(Bio::Search::HSP::GenericHSP);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::HSP::ModelHSP->new();
 Function: Builds a new Bio::Search::HSP::ModelHSP object 
 Returns : Bio::Search::HSP::ModelHSP
 Args    :

Plus Bio::Seach::HSP::ModelHSP methods

           -algorithm => algorithm used (Infernal, RNAMotif, ERPIN, etc)
           -evalue    => evalue
           -pvalue    => pvalue
           -bits      => bit value for HSP
           -score     => score value for HSP (typically z-score but depends on
					      analysis)
           -hsp_length=> Length of the HSP (including gaps)
           -identical => # of residues that that matched identically
           -conserved => # of residues that matched conservatively 
                           (only protein comparisions; 
			    conserved == identical in nucleotide comparisons)
           -hsp_gaps   => # of gaps in the HSP
           -query_gaps => # of gaps in the query in the alignment
           -hit_gaps   => # of gaps in the subject in the alignment    
           -query_name  => HSP Query sequence name (if available)
           -query_start => HSP Query start (in original query sequence coords)
           -query_end   => HSP Query end (in original query sequence coords)
           -hit_name    => HSP Hit sequence name (if available)
           -hit_start   => HSP Hit start (in original hit sequence coords)
           -hit_end     => HSP Hit end (in original hit sequence coords)
           -hit_length  => total length of the hit sequence
           -query_length=> total length of the query sequence
           -query_seq   => query sequence portion of the HSP
           -hit_seq     => hit sequence portion of the HSP
           -homology_seq=> homology sequence for the HSP
           -hit_frame   => hit frame (only if hit is translated protein)
           -query_frame => query frame (only if query is translated protein)
           -meta        => optional meta data (sec structure, markup, etc)
           -custom_score=> custom score data

=cut

=head2 meta

 Title   : meta
 Usage   : my $meta = $hsp->meta();
 Function: Returns meta data for this HSP or undef
 Returns : string of meta data or undef
 Args    : [optional] string to set value
 Note    : At some point very soon this will likely be a Bio::AnnotationI.
           Don't get used to a simple string!

=cut

sub meta {
    my ($self,$value) = @_;
    my $previous = $self->{'META'};
    if( defined $value  ) {
        $self->{'META'} = $value;
    }
    return $previous;
}

=head2 custom_score

 Title   : custom_score
 Usage   : my $data = $hsp->custom_score();
 Function: Returns custom_score data for this HSP, or undef
 Returns : custom_score data or undef
 Args    : [optional] custom_score
 Note    : This is a Get/Set used to deal with odd score-like data generated
           from RNAMotif (and other programs) where the score section
           can be customized to include non-standard data, including sequence
           data, user-based scores, and other values.

=cut

sub custom_score {
    my ($self,$value) = @_;
    my $previous = $self->{'CUSTOMSCORE'};
    if( defined $value  ) {
        $self->{'CUSTOMSCORE'} = $value;
    }
    return $previous;
}

=head2 Bio::Search::HSP::HSPI methods

Implementation of Bio::Search::HSP::HSPI methods follow

=head2 algorithm

 Title   : algorithm
 Usage   : my $r_type = $hsp->algorithm
 Function: Obtain the name of the algorithm used to obtain the HSP
 Returns : string (e.g., BLASTP)
 Args    : [optional] scalar string to set value

=cut

=head2 strand

 Title   : strand
 Usage   : $hsp->strand('hit')
 Function: Retrieves the strand for the HSP component requested
 Returns : +1 or -1 (0 if unknown)
 Args    : 'hit' or 'subject' or 'sbjct' to retrieve the strand of the subject.
           There is no strand available for 'query', as the query is a model
           and not a true sequence.

=cut

# overrides HSPI::seq()

=head2 seq

 Usage     : $hsp->seq( [seq_type] );
 Purpose   : Get the query or sbjct sequence as a Bio::Seq.pm object.
 Example   : $seqObj = $hsp->seq('sbjct');
 Returns   : Object reference for a Bio::Seq.pm object.
 Argument  : seq_type = 'query' or 'hit' or 'sbjct' (default = 'sbjct').
           :  ('sbjct' is synonymous with 'hit') 
           : default is 'sbjct'
           : Note: if there is no sequence available (eg for a model-based
           : search), this returns a LocatableSeq object w/o a sequence
 Throws    : Propagates any exception that occurs during construction
           : of the Bio::Seq.pm object.
 Comments  : The sequence is returned in an array of strings corresponding
           : to the strings in the original format of the Blast alignment.
           : (i.e., same spacing).

See Also   : L<seq_str()|seq_str>, L<Bio::Seq>

=cut

#-------
sub seq {
#-------
    my($self,$seqType) = @_; 
    $seqType ||= 'sbjct';
    $seqType = 'sbjct' if $seqType eq 'hit';
    my $str = $self->seq_str($seqType);
    if( $seqType =~ /^(m|ho)/i ) {
        $self->throw("cannot call seq on the homology match string, it isn't really a sequence, use get_aln to convert the HSP to a Bio::AlignIO and generate a consensus from that.");
    }
    require Bio::LocatableSeq;
    my $id = $seqType =~ /^q/i ? $self->query->seq_id : $self->hit->seq_id;
    $str =~ s{\*\[\s*(\d+)\s*\]\*}{'N' x $1}ge;
    $str =~ s{\s+}{}g;
    my $seq = Bio::LocatableSeq->new (-ID    => $id,
                           -START => $self->start($seqType),
                           -END   => $self->end($seqType),
                           -STRAND=> $self->strand($seqType),
                           -DESC  => "$seqType sequence ",
                           );
    $seq->seq($str) if $str;
    $seq;
}

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
    my ($self) = @_;
    require Bio::LocatableSeq;
    require Bio::SimpleAlign;
    my $aln = Bio::SimpleAlign->new;
    my %hsp = (hit =>  $self->hit_string,
               midline => $self->homology_string,
               query => $self->query_string,
               meta  => $self->meta);
    
    # this takes care of infernal issues
    if ($hsp{meta} && $hsp{meta} =~ m{~+}) {
        $self->_postprocess_hsp(\%hsp);
    }
    
    if (!$hsp{query}) {
        $self->warn("Missing query string, can't build alignment");
        return;
    }
    
    my $seqonly = $hsp{query};
    $seqonly =~ s/[\-\s]//g;
    my ($q_nm,$s_nm) = ($self->query->seq_id(),
                        $self->hit->seq_id());
    unless( defined $q_nm && CORE::length ($q_nm) ) {
        $q_nm = 'query';
    }
    unless( defined $s_nm && CORE::length ($s_nm) ) {
        $s_nm = 'hit';
    }
    my $query = Bio::LocatableSeq->new('-seq'   => $hsp{query},
                                      '-id'    => $q_nm,
                                      '-start' => $self->query->start,
                                      '-end'   => $self->query->end,
                                      );
    $seqonly = $hsp{hit};
    $seqonly =~ s/[\-\s]//g;
    my $hit =  Bio::LocatableSeq->new('-seq'    => $hsp{hit},
                                      '-id'    => $s_nm,
                                      '-start' => $self->hit->start,
                                      '-end'   => $self->hit->end,
                                      );
    $aln->add_seq($query);
    $aln->add_seq($hit);
    if ($hsp{meta}) {
        my $meta_obj = Bio::Seq::Meta->new();
        $meta_obj->named_meta('ss_cons', $hsp{meta});
        $aln->consensus_meta($meta_obj);
    }
    return $aln;
}

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

=head1 ModelHSP methods overridden in ModelHSP

The following methods have been overridden due to their current reliance on
sequence-based queries. They may be implemented in future versions of this class.

=head2 seq_inds

=cut

sub seq_inds {
    my $self = shift;
    $self->warn('$hsp->seq_inds not implemented for Model-based searches');
    return;    
}

=head2 frac_identical

=cut

sub frac_identical {
    my $self = shift;
    $self->warn('$hsp->frac_identical not implemented for Model-based searches');
    return;
}

=head2 frac_conserved

=cut

sub frac_conserved {
    my $self = shift;
    $self->warn('$hsp->frac_conserved not implemented for Model-based searches');
    return;
}

=head2 matches

=cut

sub matches {
    my $self = shift;
    $self->warn('$hsp->matches not implemented for Model-based searches');
    return;
}

=head2 num_conserved

=cut

sub num_conserved {
    my $self = shift;
    $self->warn('$hsp->num_conserved not implemented for Model-based searches');
    return;
}

=head2 num_identical

=cut

sub num_identical {
    my $self = shift;
    $self->warn('$hsp->num_identical not implemented for Model-based searches');
    return;
}

=head2 cigar_string

=cut


sub cigar_string {
    my $self = shift;
    $self->warn('$hsp->cigar_string not implemented for Model-based searches');
    return;
}

=head2 generate_cigar_string

=cut

sub generate_cigar_string {
    my $self = shift;
    $self->warn('$hsp->generate_cigar_string not implemented for Model-based searches');
    return;    
}

=head2 percent_identity

=cut

sub percent_identity {
    my $self = shift;
    $self->warn('$hsp->percent_identity not implemented for Model-based searches');
    return;
}

############## PRIVATE ##############

# the following method postprocesses HSP data in cases where the sequences
# aren't complete (which can trigger a validation error)

{
	my $SEQ_REGEX = qr/\*\[\s*(\d+)\s*\]\*/;
    my $META_REGEX = qr/(~+)/;

sub _postprocess_hsp {
	my ($self, $hsp) = @_;
	$self->throw('Must pass a hash ref for HSP processing') unless ref($hsp) eq 'HASH';
	my @ins;
	for my $type (qw(query hit meta)) {
        $hsp->{$type} =~ s{\s+$}{};
		my $str = $hsp->{$type};
		my $regex = $type eq 'meta' ? $META_REGEX : $SEQ_REGEX;
		my $ind = 0;		
		while ($str =~ m{$regex}g) {
			$ins[$ind]->{$type} = {pos => pos($str) - length($1), str => $1};
            $ind++;
		}
	}
	for my $chunk (reverse @ins) {
        my ($max, $min) = ($chunk->{hit}->{str} >= $chunk->{query}->{str}) ?
            ('hit', 'query') : ('query', 'hit');
        my %rep;
        $rep{$max} = 'N' x $chunk->{$max}->{str};
        $rep{$min} = 'N' x $chunk->{$min}->{str}.
            ('-'x($chunk->{$max}->{str}-$chunk->{$min}->{str}));
        $rep{'meta'} = '~' x $chunk->{$max}->{str};
        $rep{'midline'} = ' ' x $chunk->{$max}->{str};
        for my $t (qw(hit query meta midline)) {
            substr($hsp->{$t}, $chunk->{meta}->{pos}, length($chunk->{meta}->{str}) , $rep{$t});
        }
	}
}

}

1;
