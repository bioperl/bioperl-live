#-----------------------------------------------------------------
# $Id$
#
# BioPerl module for Bio::Search::HSP::HSPI
#
# Cared for by Steve Chervitz <sac@bioperl.org>
# and Jason Stajich <jason@bioperl.org>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::HSP::HSPI - Interface for a High Scoring Pair in a similarity search result

=head1 SYNOPSIS

    $r_type = $hsp->algorithm

    $pvalue = $hsp->pvalue();

    $evalue = $hsp->evalue();

    $frac_id = $hsp->frac_identical( ['query'|'hit'|'total'] );

    $frac_cons = $hsp->frac_conserved( ['query'|'hit'|'total'] );

    $gaps = $hsp->gaps( ['query'|'hit'|'total'] );

    $qseq = $hsp->query_string;

    $hseq = $hsp->hit_string;

    $homology_string = $hsp->homology_string;

    $len = $hsp->length( ['query'|'hit'|'total'] );

    $rank = $hsp->rank;

=head1 DESCRIPTION

Bio::Search::HSP::HSPI objects cannot be instantiated since this
module defines a pure interface.

Given an object that implements the Bio::Search::HSP::HSPI  interface,
you can do the following things with it:

=head1 SEE ALSO

This interface inherits methods from these other modules:

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

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Steve Chervitz, Jason Stajich

Email sac@bioperl.org
Email jason@bioperl.org

=head1 COPYRIGHT

Copyright (c) 2001 Steve Chervitz, Jason Stajich. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::HSP::HSPI;
use vars qw(@ISA);

use Bio::Root::RootI;
use Bio::SeqFeature::SimilarityPair;

use strict;
use Carp;

@ISA = qw(Bio::SeqFeature::SimilarityPair Bio::Root::RootI);


=head2 algorithm

 Title   : algorithm
 Usage   : my $r_type = $hsp->algorithm
 Function: Obtain the name of the algorithm used to obtain the HSP
 Returns : string (e.g., BLASTP)
 Args    : none

=cut

sub algorithm{
   my ($self,@args) = @_;
   $self->throw_not_implemented;
}

=head2 pvalue

 Title   : pvalue
 Usage   : my $pvalue = $hsp->pvalue();
 Function: Returns the P-value for this HSP or undef 
 Returns : float or exponential (2e-10)
           P-value is not defined with NCBI Blast2 reports.
 Args    : none

=cut

sub pvalue {
   my ($self) = @_;
   $self->throw_not_implemented;
}

=head2 evalue

 Title   : evalue
 Usage   : my $evalue = $hsp->evalue();
 Function: Returns the e-value for this HSP
 Returns : float or exponential (2e-10)
 Args    : none

=cut

sub evalue {
   my ($self) = @_;
   $self->throw_not_implemented;
}

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
   $self->throw_not_implemented;
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
    $self->throw_not_implemented;
}

=head2 num_identical

 Title   : num_identical
 Usage   : $obj->num_identical($newval)
 Function: returns the number of identical residues in the alignment
 Returns : integer
 Args    : integer (optional)


=cut

sub num_identical{
    shift->throw_not_implemented;
}

=head2 num_conserved

 Title   : num_conserved
 Usage   : $obj->num_conserved($newval)
 Function: returns the number of conserved residues in the alignment
 Returns : inetger
 Args    : integer (optional)


=cut

sub num_conserved{
    shift->throw_not_implemented();
}

=head2 gaps

 Title    : gaps
 Usage    : my $gaps = $hsp->gaps( ['query'|'hit'|'total'] );
 Function : Get the number of gaps in the query, hit, or total alignment.
 Returns  : Integer, number of gaps or 0 if none
 Args     : 'query' = num conserved / length of query seq (without gaps)
            'hit'   = num conserved / length of hit seq (without gaps)
            'total' = num conserved / length of alignment (with gaps)
            default = 'total' 

=cut

sub gaps        {
    my ($self, $type) = @_;
    $self->throw_not_implemented;
}

=head2 query_string

 Title   : query_string
 Usage   : my $qseq = $hsp->query_string;
 Function: Retrieves the query sequence of this HSP as a string
 Returns : string
 Args    : none


=cut

sub query_string{
   my ($self) = @_;
   $self->throw_not_implemented;
}

=head2 hit_string

 Title   : hit_string
 Usage   : my $hseq = $hsp->hit_string;
 Function: Retrieves the hit sequence of this HSP as a string
 Returns : string
 Args    : none


=cut

sub hit_string{
   my ($self) = @_;
   $self->throw_not_implemented;
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

sub homology_string{
   my ($self) = @_;
   $self->throw_not_implemented;
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

sub length{
    shift->throw_not_implemented();
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
 Usage   : my $aln = $hsp->gel_aln
 Function: Returns a Bio::SimpleAlign representing the HSP alignment
 Returns : Bio::SimpleAlign
 Args    : none

=cut

sub get_aln {
   my ($self) = @_;
   $self->throw_not_implemented;
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

             collapse  = boolean, if true, consecutive positions are merged
                         using a range notation, e.g., "1 2 3 4 5 7 9 10 11" 
                         collapses to "1-5 7 9-11". This is useful for 
                         consolidating long lists. Default = no collapse.
 Throws    : n/a.
 Comments  : 

See Also   : L<Bio::Search::BlastUtils::collapse_nums()|Bio::Search::BlastUtils>, L<Bio::Search::Hit::HitI::seq_inds()|Bio::Search::Hit::HitI>

=cut

sub seq_inds {
    shift->throw_not_implemented();
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
 Function: Get/Set the significance value (see Bio::SeqFeature::SimilarityPair)
 Returns : significance value (scientific notation string)
 Args    : significance value (sci notation string)


=head2 score

 Title   : score
 Usage   : my $score = $hsp->score();
 Function: Returns the score for this HSP or undef 
 Returns : numeric           
 Args    : [optional] numeric to set value

=head2 bits

 Title   : bits
 Usage   : my $bits = $hsp->bits();
 Function: Returns the bit value for this HSP or undef 
 Returns : numeric
 Args    : none

=cut

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

    if( $val =~ /^q/i ) { 
	return $self->query->strand(shift);
    } elsif( $val =~ /^(hi|s)/i ) {
	return $self->hit->strand(shift);
    } elsif ( $val =~ m/^(list|array)/) {
	return ($self->query->strand(shift), $self->hit->strand(shift));
    } else { 
	$self->warn("unrecognized component $val requested\n");
    }
    return 0;
}

=head2 start

 Title   : start
 Usage   : $hsp->start('query')
 Function: Retrieves the start for the HSP component requested
 Returns : integer
 Args    : 'hit' or 'subject' or 'sbjct' to retrieve the start of the subject
           'query' to retrieve the query start (default)

=cut

sub start {
    my $self = shift;
    my $val = shift;
    $val = 'query' unless defined $val;
    $val =~ s/^\s+//;

    if( $val =~ /^q/i ) { 
	return $self->query->start(shift);
    } elsif( $val =~ /^(hi|s)/i ) {
	return $self->hit->start(shift);
    } else { 
	$self->warn("unrecognized component $val requested\n");
    }
    return 0;
}

=head2 end

 Title   : end
 Usage   : $hsp->end('query')
 Function: Retrieves the end for the HSP component requested
 Returns : integer
 Args    : 'hit' or 'subject' or 'sbjct' to retrieve the end of the subject
           'query' to retrieve the query end (default)

=cut

sub end {
    my $self = shift;
    my $val = shift;
    $val = 'query' unless defined $val;
    $val =~ s/^\s+//;

    if( $val =~ /^q/i ) { 
	return $self->query->end(shift);
    } elsif( $val =~ /^(hi|s)/i ) {
	return $self->hit->end(shift);
    } else { 
	$self->warn("unrecognized component $val requested\n");
    }
    return 0;
}

sub seq_str {  
    my ($self,$type) = @_;
    if( $type =~ /^q/i ) { return $self->query_string(shift) }
    elsif( $type =~ /^s/i || $type =~ /^hi/i ) { return $self->hit_string(shift)}
    elsif ( $type =~ /^ho/i ) { return $self->hit_string(shift) }
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

sub rank { shift->throw_not_implemented }

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

See Also   : L<length()|length>, L<gaps()|gaps>, L<seq_str()|seq_str>, L<Bio::Search::Hit::BlastHit::_adjust_contigs()|Bio::Search::Hit::BlastHit>

=cut

#-----------
sub matches {
#-----------
    my( $self, %param ) = @_;
    my(@data);
    my($seqType, $beg, $end) = ($param{-SEQ}, $param{-START}, $param{-STOP});
    $seqType ||= 'query';
   $seqType = 'sbjct' if $seqType eq 'hit';

    if(!defined $beg && !defined $end) {
	## Get data for the whole alignment.
	push @data, ($self->num_identical, $self->num_conserved);
    } else {
	## Get the substring representing the desired sub-section of aln.
	$beg ||= 0;
	$end ||= 0;
	my($start,$stop) = $self->range($seqType);
	if($beg == 0) { $beg = $start; $end = $beg+$end; }
	elsif($end == 0) { $end = $stop; $beg = $end-$beg; }

	if($end >= $stop) { $end = $stop; } ##ML changed from if (end >stop)
	else { $end += 1;}   ##ML moved from commented position below, makes
                             ##more sense here
#	if($end > $stop) { $end = $stop; }
	if($beg < $start) { $beg = $start; }
#	else { $end += 1;}

#	my $seq = substr($self->seq_str('match'), $beg-$start, ($end-$beg));

	## ML: START fix for substr out of range error ------------------
	my $seq = "";
	if (($self->algorithm eq 'TBLASTN') and ($seqType eq 'sbjct'))
	{
	    $seq = substr($self->seq_str('match'),
			  int(($beg-$start)/3), int(($end-$beg+1)/3));

	} elsif (($self->algorithm eq 'BLASTX') and ($seqType eq 'query'))
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
#	 print STDERR 
#	     qq(*------------MY EXCEPTION --------------------\nSeq: ") , 
#	     $self->seq_str("$seqType"), qq("\n),$self->rank,",(  index:";
#	 print STDERR  $beg-$start, ", len: ", $end-$beg," ), (HSPRealLen:", 
#	     CORE::length $self->seq_str("$seqType");
#	 print STDERR ", HSPCalcLen: ", $stop - $start +1 ," ), 
#	     ( beg: $beg, end: $end ), ( start: $start, stop: stop )\n";
	 ## ML: END DEBUGGING CODE----------

	if(!CORE::length $seq) {
	    $self->throw("Undefined sub-sequence ($beg,$end). Valid range = $start - $stop");
	}
	## Get data for a substring.
#	printf "Collecting HSP subsection data: beg,end = %d,%d; start,stop = %d,%d\n%s<---\n", $beg, $end, $start, $stop, $seq;
#	printf "Original match seq:\n%s\n",$self->seq_str('match');
	$seq =~ s/ //g;  # remove space (no info).
	my $len_cons = CORE::length $seq;
	$seq =~ s/\+//g;  # remove '+' characters (conservative substitutions)
	my $len_id = CORE::length $seq;
	push @data, ($len_id, $len_cons);
#	printf "  HSP = %s\n  id = %d; cons = %d\n", $self->rank, $len_id, $len_cons; <STDIN>;
    }
    @data;
}
1;
