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
use overload '""' => '_overload';

# Object preamble - inheriets from Bio::SeqFeature::SimilarityPair

use Bio::SeqFeature::SimilarityPair;
use Bio::SeqFeature::Similarity;

@ISA = qw(Bio::SeqFeature::SimilarityPair);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  my $make = $self->SUPER::_initialize;

  my ($score,$bits,$match,$positive,$p,$qb,$qe,$sb,$se,$qs,
      $ss,$hs,$qname,$sname,$qlength,$slength) = 
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
			)],@args);

  # Store the aligned query as sequence feature
  if ($qe > $qb) { # normal query: start < end
    $self->query( Bio::SeqFeature::Similarity->new
      (-start=>$qb, -end=>$qe, -strand=>1, -source=>"BLAST" ) ) }
  else { # reverse query (i dont know if this is possible, but feel free to correct)
    $self->query( Bio::SeqFeature::Similarity->new
      (-start=>$qe, -end=>$qb, -strand=>-1, -source=>"BLAST" ) ) }
  
  # store the aligned subject as sequence feature
  if ($se > $sb) { # normal subject
    $self->subject( Bio::SeqFeature::Similarity->new
      (-start=>$sb, -end=>$se, -strand=>1, -source=>"BLAST" ) ) }
  else { # reverse subject: start bigger than end
    $self->subject( Bio::SeqFeature::Similarity->new
      (-start=>$se, -end=>$sb, -strand=>-1, -source=>"BLAST" ) ) }

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
  $self->{PERCENT} = int((1000 * $match)/
			 $self->query->length)/10;
  $self->{POSITIVE} = $positive;
  $self->{QS} = $qs;
  $self->{SS} = $ss;
  $self->{HS} = $hs;

  return $make; # success - we hope!
}

# to disable overloading comment this out:
sub _overload {
	my $self = shift;
	return $self->start."..".$self->end." ".$self->bits;
}

sub P               {shift->significance(@_)}
sub percent         {shift->{PERCENT}}
sub match           {shift->query->frac_identical(@_)}
sub positive        {shift->{POSITIVE}}
sub querySeq        {shift->{QS}}
sub sbjctSeq        {shift->{SS}}
sub homologySeq     {shift->{HS}}
sub qs              {shift->{QS}}
sub ss              {shift->{SS}}
sub hs              {shift->{HS}}
