###############################################################################
# Bio::Tools::BPlite::HSP
###############################################################################
# HSP = High Scoring Pair (to all non-experts as I am)
#
# The original BPlite.pm module has been written by Ian Korf !
# see http://sapiens.wustl.edu/~ikorf
#
# You may distribute this module under the same terms as perl itself

package Bio::Tools::BPlite::HSP;

use vars qw(@ISA);
use strict;

# to disable overloading comment this out:
#use overload '""' => '_overload';

# Object preamble - inheriets from Bio::SeqFeature::SimilarityPair

use Bio::SeqFeature::SimilarityPair;
use Bio::SeqFeature::Similarity;

@ISA = qw(Bio::SeqFeature::SimilarityPair);

sub new {
    my ($class, @args) = @_;

    # workaround to make sure frame is not set before strand is
    # interpreted from query/subject info 
    # this workaround removes the key from the hash
    # so the superclass does not try and work with it
    # we'll take care of setting it in this module later on

    my %newargs = @args;
    foreach ( keys %newargs ) {
	if( /frame$/i ) {
	    delete $newargs{$_};
	} 
    }
    # done with workaround

    my $self = $class->SUPER::new(%newargs);
    
    my ($score,$bits,$match,$hsplength,$positive,$gaps,$p,$qb,$qe,$sb,$se,$qs,
	$ss,$hs,$qname,$sname,$qlength,$slength,$qframe,$sframe,$blasttype) = 
	    $self->_rearrange([qw(SCORE
				  BITS
				  MATCH
				  HSPLENGTH
				  POSITIVE
				  GAPS				  
				  P
				  QUERYBEGIN
				  QUERYEND
				  SBJCTBEGIN
				  SBJCTEND
				  QUERYSEQ
				  SBJCTSEQ
				  HOMOLOGYSEQ
				  QUERYNAME
				  SBJCTNAME
				  QUERYLENGTH
				  SBJCTLENGTH
				  QUERYFRAME
				  SBJCTFRAME
				  BLASTTYPE
				  )],@args);
    
	# Determine strand meanings
	my ($queryfactor, $sbjctfactor);
	if ($blasttype eq 'BLASTN' || 
	    $blasttype eq 'BLASTX' || 
	    $blasttype eq 'TBLASTX')  { $queryfactor = 1; }
	else                          { $queryfactor = 0; }
	if ($blasttype eq 'TBLASTN' || 
	    $blasttype eq 'TBLASTX')  { $sbjctfactor = 1; }
	else                          { $sbjctfactor = 0; }

	# Set BLAST type
	$self->{'BLAST_TYPE'} = $blasttype;
	
    # Store the aligned query as sequence feature
    my $strand;
    if ($qe > $qb) {		# normal query: start < end
		if ($queryfactor) { $strand = 1; } else { $strand = undef; }
		$self->query( Bio::SeqFeature::Similarity->new
		      (-start=>$qb, -end=>$qe, -strand=>$strand, 
		       -source=>"BLAST" ) ) }
    else {			# reverse query (i dont know if this is possible, but feel free to correct)
		if ($queryfactor) { $strand = -1; } else { $strand = undef; }
		$self->query( Bio::SeqFeature::Similarity->new
		      (-start=>$qe, -end=>$qb, -strand=>$strand,
		       -source=>"BLAST" ) ) }

    # store the aligned subject as sequence feature
    if ($se > $sb) {		# normal subject
		if ($sbjctfactor) { $strand = 1; } else { $strand = undef; }
		$self->subject( Bio::SeqFeature::Similarity->new
			(-start=>$sb, -end=>$se, -strand=>$strand,
			 -source=>"BLAST" ) ) }
    else {			# reverse subject: start bigger than end
		if ($sbjctfactor) { $strand = -1; } else { $strand = undef; }
		$self->subject( Bio::SeqFeature::Similarity->new
			(-start=>$se, -end=>$sb, -strand=>$strand,
			 -source=>"BLAST" ) ) }
    
    # name the sequences
    $self->query->seqname($qname); # query
    $self->subject->seqname($sname); # subject

    # set lengths
    $self->query->seqlength($qlength); # query
    $self->subject->seqlength($slength); # subject

    # set object vars
    $self->score($score);
    $self->bits($bits);
    $self->significance($p);
    $self->query->frac_identical($match);
    $self->subject->frac_identical($match);
    $self->{'HSPLENGTH'} = $hsplength;
    $self->{'PERCENT'} = int((1000 * $match)/$hsplength)/10;
    $self->{'POSITIVE'} = $positive;
    $self->{'GAPS'} = $gaps;
    $self->{'QS'} = $qs;
    $self->{'SS'} = $ss;
    $self->{'HS'} = $hs;
    
    $self->frame($qframe, $sframe);
    return $self;		# success - we hope!
}

# to disable overloading comment this out:
sub _overload {
	my $self = shift;
	return $self->start."..".$self->end." ".$self->bits;
}

=head2 P

 Title    : P
 Usage    : $hsp->P();
 Function : returns the P (significance) value for a HSP 
 Example  : 
 Returns  : (double) significance value
 Args     :

=cut

sub P {
	my ($self, @args) = @_;
	my $float = $self->significance(@args);
	my $match = '([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?'; # Perl Cookbook 2.1
	if ($float =~ /^$match$/) {
	    # Is a C float
	    return $float;
	} elsif ("1$float" =~ /^$match$/) {
	    # Almost C float, Jitterbug 974
	    return "1$float";
	} else {
		$self->warn("[HSP::P()] '$float' is not a known number format. Returning zero (0) instead.");
		return 0;
	}
}

=head2 percent

 Title    : percent
 Usage    : $hsp->percent();
 Function : returns the percent matching 
 Returns  : (double) percent matching
 Args     : none

=cut

sub percent         {shift->{'PERCENT'}}


=head2 match

 Title    : match
 Usage    : $hsp->match();
 Function : returns the match
 Example  : 
 Returns  : (double) frac_identical 
 Args     :

=cut

sub match           {shift->query->frac_identical(@_)}

=head2 hsplength

 Title    : hsplength
 Usage    : $hsp->hsplength();
 Function : returns the HSP length (including gaps)
 Returns  : (integer) HSP length
 Args     : none

=cut

sub hsplength              {shift->{'HSPLENGTH'}}

=head2 positive

 Title    : positive
 Usage    : $hsp->positive();
 Function : returns the number of positive hits 
 Returns  : (int) number of positive residue hits 
 Args     : none

=cut

sub positive        {shift->{'POSITIVE'}}

=head2 gaps

 Title    : gaps
 Usage    : $hsp->gaps();
 Function : returns the number of gaps or 0 if none
 Returns  : (int) number of gaps or 0 if none
 Args     : none

=cut

sub gaps        {shift->{'GAPS'}}

=head2 querySeq

 Title    : querySeq
 Usage    : $hsp->querySeq();
 Function : returns the query sequence
 Returns  : (string) the Query Sequence 
 Args     : none

=cut

sub querySeq        {shift->{'QS'}}

=head2 sbjctSeq

 Title    : sbjctSeq
 Usage    : $hsp->sbjctSeq();
 Function : returns the Sbjct sequence 
 Returns  : (string) the Sbjct Sequence 
 Args     : none

=cut

sub sbjctSeq        {shift->{'SS'}}

=head2 homologySeq

 Title    : homologySeq
 Usage    : $hsp->homologySeq();
 Function : returns the homologous sequence 
 Returns  : (string) homologous sequence 
 Args     : none

=cut

sub homologySeq     {shift->{'HS'}}

=head2 qs

 Title    : qs
 Usage    : $hsp->qs();
 Function : returns the Query Sequence (same as querySeq)
 Returns  : (string) query Sequence 
 Args     : none

=cut

sub qs              {shift->{'QS'}}

=head2 ss

 Title    : ss
 Usage    : $hsp->ss();
 Function : returns the subject sequence ( same as sbjctSeq) 
 Returns  : (string) Sbjct Sequence
 Args     : none

=cut

sub ss              {shift->{'SS'}}

=head2 hs

 Title    : hs
 Usage    : $hsp->hs();
 Function : returns the Homologous Sequence (same as homologySeq ) 
 Returns  : (string) Homologous Sequence
 Args     : none

=cut

sub hs              {shift->{'HS'}}


sub frame {
    my ($self, $qframe, $sframe) = @_;
    if( defined $qframe ) {
	  if( $qframe == 0 ) {
	    $qframe = undef;
	  } elsif( $qframe !~ /^([+-])?([1-3])/ ) {	    
	    $self->warn("Specifying an invalid query frame ($qframe)");
	    $qframe = undef;
	  } else { 
	    if( ($1 eq '-' && $self->query->strand >= 0) || ($1 eq '+' && $self->query->strand <= 0) ) {
			$self->warn("Query frame ($qframe) did not match strand of query (". $self->query->strand() . ")");
	    }

	    # Set frame to GFF [0-2]
	    $qframe = $2 - 1;
	  }
	  $self->{'QFRAME'} = $qframe;
    }
    if( defined $sframe ) {
	  if( $sframe == 0 ) {
	    $sframe = undef;
	  } elsif( $sframe !~ /^([+-])?([1-3])/ ) {	    
	    $self->warn("Specifying an invalid subject frame ($sframe)");
	    $sframe = undef;
	  } else { 
	    if( ($1 eq '-' && $self->subject->strand >= 0) || ($1 eq '+' && $self->subject->strand <= 0) ) {
			$self->warn("Subject frame ($sframe) did not match strand of subject (". $self->subject->strand() . ")");
	    }

	    # Set frame to GFF [0-2]
	    $sframe = $2 - 1;
	  }
      $self->{'SFRAME'} = $sframe;
    }
    (defined $qframe && $self->SUPER::frame($qframe) && ($self->{'FRAME'} = $qframe)) || (defined $sframe && $self->SUPER::frame($sframe) && ($self->{'FRAME'} = $sframe));
    if    (wantarray() && 
           $self->{'BLAST_TYPE'} eq 'TBLASTX') { return ($self->{'QFRAME'}, $self->{'SFRAME'}); } 
    elsif (wantarray())                        { (defined $self->{'QFRAME'} && return ($self->{'QFRAME'}, undef)) || (defined $self->{'SFRAME'} && return (undef, $self->{'SFRAME'})); }
    else                                       { (defined $self->{'QFRAME'} && return $self->{'QFRAME'}) || (defined $self->{'SFRAME'} && return $self->{'SFRAME'}); }
}
1;
