# $Id$
#
# BioPerl module for Bio::Search::HSP::ModelHSP
#
# Cared for by Chris Fields <cjfields at uiuc dot edu>
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::HSP::ModelHSP - A HSP object for model-based searches

=head1 SYNOPSIS

    use Bio::Search::HSP::ModelHSP;
    # us it just like a Bio::Search::HSP::GenericHSP object

=head1 DESCRIPTION

This object is a specialization of L<Bio::Search::HSP::GenericHSP> and is used
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

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Chris Fields

Email cjfields at uiuc dot edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::Search::HSP::ModelHSP;
use strict;

use base qw(Bio::Search::HSP::GenericHSP);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::HSP::ModelHSP->new();
 Function: Builds a new Bio::Search::HSP::ModelHSP object 
 Returns : Bio::Search::HSP::ModelHSP
 Args    :

Plus Bio::Seach::HSP::GenericHSP methods

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

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  
  my ($meta, $cs) = $self->_rearrange([qw(META
                                       CUSTOM_SCORE)], @args);

  defined $meta  && $self->meta($meta);
  defined $cs    && $self->custom_score($cs);

  return $self;
}

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

sub strand {
    my $self = shift;
    my $val = shift;
    $val = 'hit' unless defined $val;
    return $self->SUPER::strand($val);
}

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
 Usage   : $hsp->frame($queryframe,$subjectframe)
 Function: Set the Frame for both query and subject and insure that
           they agree.
           This overrides the frame() method implementation in
           FeaturePair.
 Returns : array of query and subjects if return type wants an array
           or query frame if defined or subject frame
 Args    : none
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
    my $hs = $self->hit_string();
    my $qs = $self->query_string();
    if (!$qs) {
        $self->warn("Missing query string, can't build alignment");
        return;
    }
    my $seqonly = $qs;
    $seqonly =~ s/[\-\s]//g;
    my ($q_nm,$s_nm) = ($self->query->seq_id(),
                        $self->hit->seq_id());
    unless( defined $q_nm && CORE::length ($q_nm) ) {
        $q_nm = 'query';
    }
    unless( defined $s_nm && CORE::length ($s_nm) ) {
        $s_nm = 'hit';
    }
    my $query = Bio::LocatableSeq->new('-seq'   => $qs,
                                      '-id'    => $q_nm,
                                      '-start' => $self->query->start,
                                      '-end'   => $self->query->end,
                                      );
    $seqonly = $hs;
    $seqonly =~ s/[\-\s]//g;
    my $hit =  Bio::LocatableSeq->new('-seq'    => $hs,
                                      '-id'    => $s_nm,
                                      '-start' => $self->hit->start,
                                      '-end'   => $self->hit->end,
                                      );
    $aln->add_seq($query);
    $aln->add_seq($hit);
    return $aln;
}

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

See Also   : L<Bio::Search::BlastUtils::collapse_nums()|Bio::Search::BlastUtils>, L<Bio::Search::Hit::HitI::seq_inds()|Bio::Search::Hit::HitI>

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

=head1 GenericHSP methods overridden in ModelHSP

The following methods have been overridden due to their current reliance on
sequence-based queries. They may be implemented in future versions of this class.

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

# the following subs override several GenericHSP private methods
# to allow for Model-based queries

# needed before seqfeatures can be made

sub _pre_seq_feature {
    my $self = shift;
    my $algo = $self->{ALGORITHM};
    my ($queryfactor, $hitfactor) = (0,1);
    # no exceptions made for algorithms yet
    $self->{_query_factor} = $queryfactor;
    $self->{_hit_factor} = $hitfactor;
}

# make query seq feature
sub _query_seq_feature {
    my $self = shift;
    $self->{_making_qff} = 1;
    my $qs = $self->{QUERY_START};
    my $qe = $self->{QUERY_END};
    unless (defined $self->{_query_factor}) {
        $self->_pre_seq_feature;
    }
    my $queryfactor = $self->{_query_factor};
    unless( defined $qe && defined $qs ) { $self->throw("Did not specify a Query End or Query Begin"); }
    my $sim1 = $self->{_sim1} || Bio::SeqFeature::Similarity->new(-verbose => $self->verbose);
    $sim1->start($qs);
    $sim1->end($qe);
    $sim1->significance($self->{EVALUE});
    $sim1->bits($self->{BITS});
    $sim1->score($self->{SCORE});
    # models do not have strandedness
    $sim1->strand(0);
    $sim1->seq_id($self->{QUERY_NAME});
    $sim1->seqlength($self->{QUERY_LENGTH});
    $sim1->source_tag($self->{ALGORITHM});
    $sim1->seqdesc($self->{QUERY_DESC});
    $sim1->add_tag_value('meta', $self->{META}) if $self->meta;
    
    $self->Bio::Search::HSP::HSPI::feature1($sim1);

    $self->{QUERY_FRAME} = 0; # no frame for a model

    $self->{_created_qff} = 1;
    $self->{_making_qff} = 0;
    $self->_pre_frame;
}

# make subject seq feature
sub _subject_seq_feature {
    my $self = shift;
    $self->{_making_sff} = 1;
    my $hs = $self->{HIT_START};
    my $he = $self->{HIT_END};
    unless (defined $self->{_hit_factor}) {
        $self->_pre_seq_feature;
    }
    my $hitfactor = $self->{_hit_factor};

    unless( defined $he && defined $hs ) { $self->throw("Did not specify a Hit End or Hit Begin"); }

    my $strand;
    if ($he > $hs) { # normal subject
        if ($hitfactor) {
            $strand = 1;
        }
        else {
            $strand = undef;
        }
    }
    else {
        if ($hitfactor) {
            $strand = -1;
        }
        else {
            $strand = undef;
        }
        ($hs,$he) = ( $he,$hs); # reverse subject: start bigger than end
    }

    my $sim2 = $self->{_sim2} || Bio::SeqFeature::Similarity->new(-verbose => $self->verbose);
    $sim2->start($hs);
    $sim2->end($he);
    $sim2->significance($self->{EVALUE});
    $sim2->bits($self->{BITS});
    $sim2->score($self->{SCORE});
    $sim2->strand($strand);
    $sim2->seq_id($self->{HIT_NAME});
    $sim2->seqlength($self->{HIT_LENGTH});
    $sim2->source_tag($self->{ALGORITHM});
    $sim2->seqdesc($self->{HIT_DESC});
    $sim2->add_tag_value('meta', $self->{META}) if $self->meta;
    $self->Bio::Search::HSP::HSPI::feature2($sim2);

    my $hframe = $self->{HIT_FRAME};
    if (defined $strand && ! defined $hframe && $hitfactor) {
        $hframe = ( $hs % 3 ) * $strand;
    }
    elsif (! defined $strand) {
        $hframe = 0;
    }
    $self->{HIT_FRAME} = $hframe;

    $self->{_created_sff} = 1;
    $self->{_making_sff} = 0;
    $self->_pre_frame;
}

1;
