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

# Object preamble - inheriets from Bio::SeqFeature::FeaturePair

use Bio::SeqFeature::FeaturePair;
use Bio::SeqFeature::Generic; # we want to use Generic Features

@ISA = qw(Bio::SeqFeature::FeaturePair);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  my $make = $self->SUPER::_initialize;

  my ($score,$bits,$match,$positive,$length,
      $p,$qb,$qe,$sb,$se,$qs,$ss,$hs,$qname,$sname) = 
      $self->_rearrange([qw(SCORE
			    BITS
			    MATCH
			    POSITIVE
			    LENGTH
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
			    )],@args);
  # set object vars
  $self->{SCORE} = $score;
  $self->{BITS} = $bits;
  $self->{PERCENT} = int((1000 * $match)/$length)/10;
  $self->{MATCH} = $match;
  $self->{POSITIVE} = $positive;
  $self->{LENGTH} = $length;
  $self->{P} = $p;
  $self->{QS} = $qs;
  $self->{SS} = $ss;
  $self->{HS} = $hs;

  # Store the aligned query as sequence feature
  if ($qe > $qb) { # normal query: start < end
    $self->feature1( Bio::SeqFeature::Generic->new
      (-start=>$qb, -end=>$qe, -strand=>1, -primary=>"similarity", source=>"BLAST" ) ) }
  else { # reverse query (i dont know if this is possible, but feel free to correct)
    $self->feature1( Bio::SeqFeature::Generic->new
      (-start=>$qe, -end=>$qb, -strand=>-1, -primary=>"similarity", source=>"BLAST" ) ) }
  
  # store the aligned subject as sequence feature
  if ($se > $sb) { # normal subject
    $self->feature2( Bio::SeqFeature::Generic->new
      (-start=>$sb, -end=>$se, -strand=>1, -primary=>"similarity", source=>"BLAST" ) ) }
  else { # reverse subject: start bigger than end
    $self->feature2( Bio::SeqFeature::Generic->new
      (-start=>$se, -end=>$sb, -strand=>-1, -primary=>"similarity", source=>"BLAST" ) ) }

  # name the sequences
  $self->feature1->seqname($qname); # query
  $self->feature2->seqname($sname); # sbjct

  return $make; # success - we hope!
}

# to disable overloading comment this out:
sub _overload {
	my $self = shift;
	return $self->start."..".$self->end." ".$self->bits;
}

# aliases for features 1 & 2
sub query {shift->feature1(@_)}
sub subject {shift->feature2(@_)}
sub sbjct {shift->feature2(@_)}

sub score           {shift->{SCORE}}
sub hscore          {shift->{SCORE}} # for compatibility with FeaturePair
sub bits            {shift->{BITS}}
sub percent         {shift->{PERCENT}}
sub match           {shift->{MATCH}}
sub positive        {shift->{POSITIVE}}
sub length          {shift->{LENGTH}}
sub P               {shift->{P}}
sub querySeq        {shift->{QS}}
sub sbjctSeq        {shift->{SS}}
sub homologySeq     {shift->{HS}}
sub qs              {shift->{QS}}
sub ss              {shift->{SS}}
sub hs              {shift->{HS}}
