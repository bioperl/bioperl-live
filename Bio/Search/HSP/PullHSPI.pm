#
# BioPerl module for Bio::Search::HSP::PullHSPI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::HSP::PullHSPI - Bio::Search::HSP::HSPI interface for pull parsers.

=head1 SYNOPSIS

  # This is an interface and cannot be instantiated

  # generally we use Bio::SearchIO to build these objects
  use Bio::SearchIO;
  my $in = Bio::SearchIO->new(-format => 'hmmer_pull',
                              -file   => 'result.hmmer');

  while (my $result = $in->next_result) {
      while (my $hit = $result->next_hit) {
          while (my $hsp = $hit->next_hsp) {
              $r_type = $hsp->algorithm;
              $pvalue = $hsp->p();
              $evalue = $hsp->evalue();
              $frac_id = $hsp->frac_identical( ['query'|'hit'|'total'] );
              $frac_cons = $hsp->frac_conserved( ['query'|'hit'|'total'] );
              $gaps = $hsp->gaps( ['query'|'hit'|'total'] );
              $qseq = $hsp->query_string;
              $hseq = $hsp->hit_string;
              $homo_string = $hsp->homology_string;
              $len = $hsp->length( ['query'|'hit'|'total'] );
              $len = $hsp->length( ['query'|'hit'|'total'] );
              $rank = $hsp->rank;
          }
      }
  }


=head1 DESCRIPTION

PullHSP is for fast implementations that only do parsing work on the hsp
data when you actually request information by calling one of the HSPI
methods.

Many methods of HSPI are implemented in a way suitable for inheriting classes
that use Bio::PullParserI. It only really makes sense for PullHSP modules to be
created by (and have as a -parent) PullHit modules.

In addition to the usual -chunk and -parent, -hsp_data is all you should supply
when making a PullHSP object. This will store that data and make it accessible
via _raw_hsp_data, which you can access in your subclass. It would be best to
simply provide the data as the input -chunk instead, if the raw data is large
enough.

=head1 SEE ALSO

This module inherits methods from these other modules:

L<Bio::SeqFeatureI>,
L<Bio::SeqFeature::FeaturePair>
L<Bio::SeqFeature::SimilarityPair>

Please refer to these modules for documentation of the 
many additional inherited methods.

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

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 COPYRIGHT

Copyright (c) 2006 Sendu Bala. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Search::HSP::PullHSPI;


use strict;

use base qw(Bio::Search::HSP::HSPI Bio::PullParserI);

=head2 _setup

 Title   : _setup
 Usage   : $self->_setup(@args)
 Function: Implementers should call this to setup common fields and deal with
           common arguments to new().
 Returns : n/a
 Args    : @args received in new().

=cut

sub _setup {
    my ($self, @args) = @_;
	
	# fields most subclasses probably will want
	$self->_fields( { ( hsp_length => undef,
                        identical => undef,
                        percent_identity => undef,
                        conserved => undef,
                        hsp_gaps => undef,
                        query_gaps => undef,
                        hit_gaps => undef,
						evalue => undef,
						pvalue => undef,
						score => undef,
						query_start => undef,
						query_end => undef,
						query_string => undef,
						hit_start => undef,
						hit_end => undef,
						hit_string => undef,
						homology_string => undef,
						rank => undef,
                        seq_inds => undef,
                        hit_identical_inds => undef,
                        hit_conserved_inds => undef,
                        hit_nomatch_inds => undef,
                        hit_gap_inds => undef,
                        query_identical_inds => undef,
                        query_conserved_inds => undef,
                        query_nomatch_inds => undef,
                        query_gap_inds => undef ) } );
	
	my ($parent, $chunk, $hsp_data) = $self->_rearrange([qw(PARENT
														    CHUNK
															HSP_DATA)], @args);
	
    $self->throw("Need -parent or -chunk to be defined") unless defined $parent || $chunk;
    
	$self->parent($parent) if $parent;
    
    if ($chunk) {
        my ($io, $start, $end) = (undef, 0, undef);
        if (ref($chunk) eq 'ARRAY') {
            ($io, $start, $end) = @{$chunk};
        }
        else {
            $io = $chunk;
        }
        $self->chunk($io, -start => $start, -end => $end);
    }
    
	$self->_raw_hsp_data($hsp_data) if $hsp_data;
	
    return $self;
}

sub _raw_hsp_data {
	my $self = shift;
	if (@_) {
		$self->{_raw_hsp_data} = shift;
	}
	return $self->{_raw_hsp_data};
}

#
# Some of these methods are written explitely to avoid HSPI throwing not
# implemented or the wrong ancestor class being used to answer the method;
# if it didn't do that then PullParserI AUTOLOAD would have cought them.
#

=head2 algorithm

 Title   : algorithm
 Usage   : my $r_type = $hsp->algorithm
 Function: Obtain the name of the algorithm used to obtain the HSP
 Returns : string (e.g., BLASTP)
 Args    : none

=cut

sub algorithm {
	return shift->get_field('algorithm');
}

=head2 pvalue

 Title   : pvalue
 Usage   : my $pvalue = $hsp->pvalue();
 Function: Returns the P-value for this HSP or undef 
 Returns : float or exponential (2e-10)
 Args    : none

=cut

sub pvalue {
	return shift->get_field('pvalue');
}

=head2 evalue

 Title   : evalue
 Usage   : my $evalue = $hsp->evalue();
 Function: Returns the e-value for this HSP
 Returns : float or exponential (2e-10)
 Args    : none

=cut

sub evalue {
	return shift->get_field('evalue');
}

*expect = \&evalue;

=head2 frac_identical

 Title   : frac_identical
 Usage   : my $frac_id = $hsp->frac_identical( ['query'|'hit'|'total'] );
 Function: Returns the fraction of identitical positions for this HSP 
 Returns : Float in range 0.0 -> 1.0
 Args    : 'query' = num identical / length of query seq (without gaps)
           'hit'   = num identical / length of hit seq (without gaps)
           'total' = num identical / length of alignment (with gaps)
           default = 'total' 

=cut

sub frac_identical {
	my ($self, $type) = @_;
	
	$type = lc $type if defined $type;
	$type = 'hit' if (defined $type && $type =~ /subject|sbjct/);
	$type = 'total' if (! defined $type || $type eq 'hsp' || $type !~ /query|hit|subject|sbjct|total/);
	
	my $ratio = $self->num_identical($type) / $self->length($type);
    return sprintf( "%.4f", $ratio);
}

=head2 frac_conserved

 Title    : frac_conserved
 Usage    : my $frac_cons = $hsp->frac_conserved( ['query'|'hit'|'total'] );
 Function : Returns the fraction of conserved positions for this HSP.
            This is the fraction of symbols in the alignment with a 
            positive score.
 Returns : Float in range 0.0 -> 1.0
 Args    : 'query' = num conserved / length of query seq (without gaps)
           'hit'   = num conserved / length of hit seq (without gaps)
           'total' = num conserved / length of alignment (with gaps)
           default = 'total' 

=cut

sub frac_conserved {
    my ($self, $type) = @_;
	
	$type = lc $type if defined $type;
	$type = 'hit' if (defined $type && $type =~ /subject|sbjct/);
	$type = 'total' if (! defined $type || $type eq 'hsp' || $type !~ /query|hit|subject|sbjct|total/);
	
	my $ratio = $self->num_conserved($type) / $self->length($type);
    return sprintf( "%.4f", $ratio);
}

=head2 num_identical

 Title   : num_identical
 Usage   : $obj->num_identical($newval)
 Function: returns the number of identical residues in the alignment
 Returns : integer
 Args    : integer (optional)

=cut

sub num_identical {
    my $self = shift;
	return scalar($self->seq_inds('hit', 'identical'));
}

=head2 num_conserved

 Title   : num_conserved
 Usage   : $obj->num_conserved($newval)
 Function: returns the number of conserved residues in the alignment
 Returns : inetger
 Args    : integer (optional)

=cut

sub num_conserved {
    my $self = shift;
	return scalar($self->seq_inds('hit', 'conserved-not-identical'));
}

=head2 gaps

 Title    : gaps
 Usage    : my $gaps = $hsp->gaps( ['query'|'hit'|'total'] );
 Function : Get the number of gap characters in the query, hit, or total alignment.
 Returns  : Integer, number of gap characters or 0 if none
 Args     : 'query', 'hit' or 'total'; default = 'total' 

=cut

sub gaps {
    my ($self, $type) = @_;
    $type = lc $type if defined $type;
    $type = 'total' if (! defined $type || $type eq 'hsp' || $type !~ /query|hit|subject|sbjct|total/); 
    $type = 'hit' if $type =~ /sbjct|subject/;
	
	if ($type eq 'total') {
		return scalar($self->seq_inds('hit', 'gap')) + scalar($self->seq_inds('query', 'gap'));
	}
	return scalar($self->seq_inds($type, 'gap'));
}

=head2 query_string

 Title   : query_string
 Usage   : my $qseq = $hsp->query_string;
 Function: Retrieves the query sequence of this HSP as a string
 Returns : string
 Args    : none

=cut

sub query_string {
	return shift->get_field('query_string');
}

=head2 hit_string

 Title   : hit_string
 Usage   : my $hseq = $hsp->hit_string;
 Function: Retrieves the hit sequence of this HSP as a string
 Returns : string
 Args    : none

=cut

sub hit_string {
	return shift->get_field('hit_string');
}

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

sub homology_string {
	return shift->get_field('homology_string');
}

=head2 length

 Title    : length
 Usage    : my $len = $hsp->length( ['query'|'hit'|'total'] );
 Function : Returns the length of the query or hit in the alignment (without gaps) 
            or the aggregate length of the HSP (including gaps;
            this may be greater than either hit or query )
 Returns  : integer
 Args     : 'query' = length of query seq (without gaps)
            'hit'   = length of hit seq (without gaps)
            'total' = length of alignment (with gaps)
            default = 'total' 
 Args    : none

=cut

sub length {
    my ($self, $type) = @_;
    $type = 'total' unless defined $type;
    $type = lc $type;

    if ($type =~ /^q/i) {
        return $self->query->length;
    }
	elsif ($type =~ /^(hit|subject|sbjct)/) {
        return $self->hit->length;
    }
	else { 
        return $self->hit->length + $self->gaps('hit');
	}
}

=head2 hsp_length

 Title   : hsp_length
 Usage   : my $len = $hsp->hsp_length()
 Function: shortcut  length('hsp')
 Returns : floating point between 0 and 100 
 Args    : none

=cut

sub hsp_length {
	return shift->length('total');
}

=head2 percent_identity

 Title   : percent_identity
 Usage   : my $percentid = $hsp->percent_identity()
 Function: Returns the calculated percent identity for an HSP
 Returns : floating point between 0 and 100 
 Args    : none

=cut

sub percent_identity{
	my ($self) = @_;
	return $self->frac_identical('hsp') * 100;   
}

=head2 get_aln

 Title   : get_aln
 Usage   : my $aln = $hsp->get_aln
 Function: Returns a Bio::SimpleAlign representing the HSP alignment
 Returns : Bio::SimpleAlign
 Args    : none

=cut

sub get_aln {
	my $self = shift;
	
    require Bio::LocatableSeq;
    require Bio::SimpleAlign;
    my $aln = Bio::SimpleAlign->new();
    my $hs = $self->seq('hit');
    my $qs = $self->seq('query');
	if ($hs && $qs) {
		$aln->add_seq($hs);
		$aln->add_seq($qs);
		return $aln;
	}
	return;
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
              ('sbjct' is synonymous with 'hit') 
             class     = 'identical' or 'conserved' or 'nomatch' or 'gap'
                          (default = identical)
                          (can be shortened to 'id' or 'cons')
		                  Note that 'conserved' includes identical unless you
		                  use 'conserved-not-identical'

             collapse  = boolean, if true, consecutive positions are merged
                         using a range notation, e.g., "1 2 3 4 5 7 9 10 11" 
                         collapses to "1-5 7 9-11". This is useful for 
                         consolidating long lists. Default = no collapse.
 Throws    : n/a.
 Comments  : 

See Also   : L<Bio::Search::BlastUtils::collapse_nums()|Bio::Search::BlastUtils>, L<Bio::Search::Hit::HitI::seq_inds()|Bio::Search::Hit::HitI>

=cut

sub seq_inds {
    my ($self, $seqType, $class, $collapse) = @_;
    
    $seqType ||= 'query';
    $class ||= 'identical';
    $collapse ||= 0;
    $seqType = lc($seqType);
    $class = lc($class);
    $seqType = 'hit' if $seqType eq 'sbjct';
    
    my $t = substr($seqType,0,1);
    if ($t eq 'q') {
        $seqType = 'query';
    }
    elsif ($t eq 's' || $t eq 'h') {
        $seqType = 'hit';
    }
    else { 
        $self->warn("unknown seqtype $seqType using 'query'");
        $seqType = 'query';
    }
    
    $t = substr($class,0,1);
    if ($t eq 'c') {
        if ($class eq 'conserved-not-identical') {
            $class = 'conserved';
        }
        else { 
            $class = 'conservedall';
        }
    }
    elsif ($t eq 'i') {
        $class = 'identical';
    }
    elsif ($t eq 'n') {
        $class = 'nomatch';
    }
    elsif ($t eq 'g') {
        $class = 'gap';
    }
    else { 
        $self->warn("unknown sequence class $class using 'identical'");
        $class = 'identical';
    }
    
    $seqType .= '_';
    $class .= '_inds';
    
    my @ary;
    if ($class eq 'conservedall_inds') {
		my %tmp = map { $_, 1 } @{$self->get_field($seqType.'conserved_inds')},
								@{$self->get_field($seqType.'identical_inds')};
		@ary = sort {$a <=> $b} keys %tmp;
    }
    else { 
        @ary = @{$self->get_field($seqType.$class)};
    }
    
    return $collapse ? &Bio::Search::SearchUtils::collapse_nums(@ary) : @ary;
}

=head2 Inherited from L<Bio::SeqFeature::SimilarityPair>

These methods come from L<Bio::SeqFeature::SimilarityPair>

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
 Function: Get/Set the significance value (see Bio::SeqFeature::SimilarityPair)
 Returns : significance value (scientific notation string)
 Args    : significance value (sci notation string)

=cut

sub significance {
	return shift->get_field('evalue');
}

=head2 score

 Title   : score
 Usage   : my $score = $hsp->score();
 Function: Returns the score for this HSP or undef 
 Returns : numeric           
 Args    : [optional] numeric to set value

=cut

sub score {
	return shift->get_field('score');
}

=head2 bits

 Title   : bits
 Usage   : my $bits = $hsp->bits();
 Function: Returns the bit value for this HSP or undef 
 Returns : numeric
 Args    : none

=cut

sub bits {
	return shift->get_field('bits');
}

# override 

=head2 strand

 Title   : strand
 Usage   : $hsp->strand('query')
 Function: Retrieves the strand for the HSP component requested
 Returns : +1 or -1 (0 if unknown)
 Args    : 'hit' or 'subject' or 'sbjct' to retrieve the strand of the subject
           'query' to retrieve the query strand (default)
           'list' or 'array' to retreive both query and hit together

=cut

sub strand {
    my $self = shift;
    my $val = shift;
    $val = 'query' unless defined $val;
    $val =~ s/^\s+//;

    if ($val =~ /^q/i) {
        return $self->query->strand(@_);
    }
    elsif ($val =~ /^hi|^s/i) {
        return $self->hit->strand(@_);
    }
    elsif ($val =~ /^list|array/i) {
        return ($self->query->strand(@_), $self->hit->strand(@_) );
    }
    else { 
        $self->warn("unrecognized component '$val' requested\n");
    }
    return 0;
}

=head2 start

 Title   : start
 Usage   : $hsp->start('query')
 Function: Retrieves the start for the HSP component requested
 Returns : integer, or list of two integers (query start and subject start) in
           list context
 Args    : 'hit' or 'subject' or 'sbjct' to retrieve the start of the subject
           'query' to retrieve the query start (default)

=cut

sub start {
    my $self = shift;
    my $val = shift;
    $val = (wantarray ? 'list' : 'query') unless defined $val;
    $val =~ s/^\s+//;

    if ($val =~ /^q/i) { 
        return $self->query->start(@_);
    }
    elsif ($val =~ /^(hi|s)/i) {
        return $self->hit->start(@_);
    }
    elsif ($val =~ /^list|array/i) {
        return ($self->query->start(@_), $self->hit->start(@_) );
    }
    else { 
        $self->warn("unrecognized component '$val' requested\n");
    }
    return 0;
}

=head2 end

 Title   : end
 Usage   : $hsp->end('query')
 Function: Retrieves the end for the HSP component requested
 Returns : integer, or list of two integers (query end and subject end) in
           list context
 Args    : 'hit' or 'subject' or 'sbjct' to retrieve the end of the subject
           'query' to retrieve the query end (default)

=cut

sub end {
    my $self = shift;
    my $val = shift;
    $val = (wantarray ? 'list' : 'query') unless defined $val;
    $val =~ s/^\s+//;

    if ($val =~ /^q/i) { 
        return $self->query->end(@_);
    }
    elsif ($val =~ /^(hi|s)/i) {
        return $self->hit->end(@_);
    }
    elsif ($val =~ /^list|array/i) {
        return ($self->query->end(@_), $self->hit->end(@_) );
    }
    else {
        $self->warn("unrecognized end component '$val' requested\n");
    }
    return 0;
}

=head2 seq

 Usage     : $hsp->seq( [seq_type] );
 Purpose   : Get the query or sbjct sequence as a Bio::Seq.pm object.
 Example   : $seqObj = $hsp->seq('query');
 Returns   : Object reference for a Bio::LocatableSeq object.
 Argument  : seq_type = 'query' or 'hit' or 'sbjct' (default = 'query').
           : ('sbjct' is synonymous with 'hit') 
           : default is 'query'

=cut

sub seq {
    my ($self, $seqType) = @_; 
    $seqType ||= 'query';
    $seqType = 'hit' if $seqType eq 'sbjct';
    if ($seqType =~ /^(m|ho)/i ) {
        $self->throw("cannot call seq on the homology match string, it isn't really a sequence, use get_aln to convert the HSP to a Bio::AlignIO and generate a consensus from that.");
    }
	
    my $str = $self->seq_str($seqType) || return;
    require Bio::LocatableSeq;
    my $id = ($seqType =~ /^q/i) ? $self->query->seq_id : $self->hit->seq_id;
    return Bio::LocatableSeq->new(  -ID        => $id,
                                    -SEQ       => $str,
                                    -START     => $self->start($seqType),
                                    -END       => $self->end($seqType),
                                    -STRAND    => $self->strand($seqType),
                                    -FORCE_NSE => $id ? 0 : 1,
                                    -DESC      => "$seqType sequence " );

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
           : default is 'query'
 Throws    : Exception if the argument does not match an accepted seq_type.
 Comments  : 

See Also   : L<seq()|seq>, L<seq_inds()|seq_inds>, B<_set_match_seq()>

=cut

sub seq_str {  
    my $self = shift;
    my $type = shift || 'query';

    if ($type =~ /^q/i) {
        return $self->query_string(@_);
    }
    elsif ($type =~ /^(s)|(hi)/i) {
        return $self->hit_string(@_);
    }
    elsif ($type =~ /^(ho)|(ma)/i) {
        return $self->homology_string(@_);
    }
    else { 
        $self->warn("unknown sequence type $type");
    }
    return '';
}

=head2 rank

 Usage     : $hsp->rank( [string] );
 Purpose   : Get the rank of the HSP within a given Blast hit.
 Example   : $rank = $hsp->rank;
 Returns   : Integer (1..n) corresponding to the order in which the HSP
             appears in the BLAST report.

=cut

sub rank {
    return shift->get_field('rank');
}

=head2 matches

 Usage     : $hsp->matches(-seq   => 'hit'|'query', 
                           -start => $start, 
                           -stop  => $stop);
 Purpose   : Get the total number of identical and conservative matches 
           : in the query or sbjct sequence for the given HSP. Optionally can
           : report data within a defined interval along the seq.
 Example   : ($id,$cons) = $hsp_object->matches(-seq   => 'hit');
           : ($id,$cons) = $hsp_object->matches(-seq   => 'query',
                                                -start => 300,
                                                -stop  => 400);
 Returns   : 2-element array of integers 
 Argument  : (1) seq_type = 'query' or 'hit' or 'sbjct' (default = query)
           :  ('sbjct' is synonymous with 'hit') 
           : (2) start = Starting coordinate (optional)
           : (3) stop  = Ending coordinate (optional)

=cut

sub matches {
    my ($self, @args) = @_;
    my($seqType, $beg, $end) = $self->_rearrange([qw(SEQ START STOP)], @args);
    $seqType ||= 'query';
    $seqType = 'hit' if $seqType eq 'sbjct';
	
    my @data;
    if ((!defined $beg && !defined $end) || ! $self->seq_str('match')) {
        push @data, ($self->num_identical, $self->num_conserved);
    }
	else {
        $beg ||= 0;
        $end ||= 0;
        my ($start, $stop) = $self->range($seqType);
		
        if ($beg == 0) {
			$beg = $start;
			$end = $beg+$end;
		}
        elsif ($end == 0) {
			$end = $stop;
			$beg = $end-$beg;
		}
		
        if ($end >= $stop) {
			$end = $stop;
		}
        else {
			$end += 1;
		}
        if ($beg < $start) {
			$beg = $start;
		}
        
        my $seq = substr($self->seq_str('homology'), $beg-$start, ($end-$beg));
        
        if (!CORE::length $seq) {
            $self->throw("Undefined sub-sequence ($beg,$end). Valid range = $start - $stop");
        }
        ## Get data for a substring.
        $seq =~ s/ //g;  # remove space (no info).
        my $len_cons = CORE::length $seq;
        $seq =~ s/\+//g;  # remove '+' characters (conservative substitutions)
        my $len_id = CORE::length $seq;
        push @data, ($len_id, $len_cons);
    }
	
    return @data;
}

=head2 n

 Usage     : $hsp_obj->n()
 Purpose   : Get the N value (num HSPs on which P/Expect is based).
 Returns   : Integer or null string if not defined.
 Argument  : n/a
 Throws    : n/a
 Comments  : The 'N' value is listed in parenthesis with P/Expect value:
           : e.g., P(3) = 1.2e-30  ---> (N = 3).
           : Not defined in NCBI Blast2 with gaps.
           : This typically is equal to the number of HSPs but not always.

=cut

sub n {
    return shift->get_field('num_hsps');
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
 Comments  : This is a convenience method for constructions such as
             ($hsp->query->start, $hsp->query->end)

=cut

sub range {
    my ($self, $seqType) = @_;
    $seqType ||= 'query';
	
    my ($start, $end);
    if ($seqType eq 'query') {
        $start = $self->query->start;
        $end = $self->query->end;
    }
    else {
        $start = $self->hit->start;
        $end = $self->hit->end;
    }
    return ($start, $end);
}

#*** would want cigar stuff from GenericHSP - move to HSPI?

1;

