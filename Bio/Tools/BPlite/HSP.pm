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
    
    my ($score,$bits,$match,$positive,$p,$qb,$qe,$sb,$se,$qs,
	$ss,$hs,$qname,$sname,$qlength,$slength, $frame) = 
	    $self->_rearrange([qw(SCORE
				  BITS
				  MATCH
				  POSITIVE
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
				  FRAME
				  )],@args);
    
    # Store the aligned query as sequence feature
    if ($qe > $qb) {		# normal query: start < end
	$self->query( Bio::SeqFeature::Similarity->new
		      (-start=>$qb, -end=>$qe, -strand=>1, 
		       -source=>"BLAST" ) ) }
    else {			# reverse query (i dont know if this is possible, but feel free to correct)
	$self->query( Bio::SeqFeature::Similarity->new
		      (-start=>$qe, -end=>$qb, -strand=>-1,
		       -source=>"BLAST" ) ) }

    # store the aligned subject as sequence feature
    if ($se > $sb) {		# normal subject
	$self->subject( Bio::SeqFeature::Similarity->new
			(-start=>$sb, -end=>$se, -strand=>1,
			 -source=>"BLAST" ) ) }
    else {			# reverse subject: start bigger than end
	$self->subject( Bio::SeqFeature::Similarity->new
			(-start=>$se, -end=>$sb, -strand=>-1, 
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
    $self->{'PERCENT'} = int((1000 * $match)/
			     $self->query->length)/10;
    $self->{'POSITIVE'} = $positive;
    $self->{'QS'} = $qs;
    $self->{'SS'} = $ss;
    $self->{'HS'} = $hs;
    
    defined $frame && $self->frame($frame);
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

sub P               {shift->significance(@_)}

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

=head2 positive

 Title    : positive
 Usage    : $hsp->positive();
 Function : returns the number of positive hits 
 Returns  : (int) number of positive residue hits 
 Args     : none

=cut

sub positive        {shift->{'POSITIVE'}}

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
    my ($self, $frame) = @_;
    if( defined $frame ) {
	if( $frame == 0 ) {
	    $frame = undef;
	} elsif( $frame !~ /^([+-])?([1-3])/ ) {	    
	    $self->warn("Specifying an invalid frame ($frame)");
	    $frame = undef;
	} else { 
	    # JB 949 - Creates too many warnings for blastx report.
	    #          Future enhancement to BPLite::_parseHeader needed to set report type
	    #          so that subject strand is used with tblastn and query strand with blastx
	    # if( ($1 eq '-' && $self->subject->strand >= 0) ||
		# ($1 eq '+' && $self->subject->strand <= 0) ) {
		# $self->warn("Frame ($frame) did not match strand of query match (".
		# 	    $self->subject->strand().")");
	    # }
	    
	    # Set frame to GFF [0-2]
	    $frame = $2 - 1;
	}
    }
    return $self->SUPER::frame($frame);
}
1;
