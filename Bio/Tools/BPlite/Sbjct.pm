###############################################################################
# Bio::Tools::BPlite::Sbjct
###############################################################################
#
# The original BPlite.pm module has been written by Ian Korf !
# see http://sapiens.wustl.edu/~ikorf
#
# You may distribute this module under the same terms as perl itself

package Bio::Tools::BPlite::Sbjct;

use strict;

use Bio::Root::Object; # root object to inherit from
use Bio::Tools::BPlite::HSP; # we want to use HSP

use overload '""' => 'name';
use vars qw(@ISA);

@ISA = qw(Bio::Root::Object);

sub _initialize {
  my ($self, @param) = @_; 
  my $make = $self->SUPER::_initialize;
  
  ($self->{NAME},$self->{FH},$self->{LASTLINE},$self->{PARENT}) = @param;
  $self->{HSP_ALL_PARSED} = 0;
    
  return $make; # success - we hope!
}

sub name {shift->{NAME}}

sub nextFeaturePair {shift->nextHSP}; # just another name

sub nextHSP {
  my ($self) = @_;
  return 0 if $self->{HSP_ALL_PARSED};
  
  ############################
  # get and parse scorelines #
  ############################
  my $scoreline = $self->{LASTLINE};
  my $FH = $self->{FH};
  my $nextline = <$FH>;
  return undef if not defined $nextline;
  $scoreline .= $nextline;
  my ($score, $bits);
  if ($scoreline =~ /\d bits\)/) {
    ($score, $bits) = $scoreline =~
      /Score = (\d+) \((\S+) bits\)/; # WU-BLAST
  }
  else {
    ($bits, $score) = $scoreline =~
      /Score =\s+(\S+) bits \((\d+)/; # NCBI-BLAST
  }
  
  my ($match, $length) = $scoreline =~ /Identities = (\d+)\/(\d+)/;
  my ($positive) = $scoreline =~ /Positives = (\d+)/;
  $positive = $match if not defined $positive;
  my ($p)        = $scoreline =~ /[Sum ]*P[\(\d+\)]* = (\S+)/;
  if (not defined $p) {($p) = $scoreline =~ /Expect = (\S+)/}
  
  $self->throw("Unable to parse '$scoreline'") if not defined $score;
  
  #######################
  # get alignment lines #
  #######################
  my @hspline;
  while(<$FH>) {
    if ($_ =~ /^WARNING:|^NOTE:/) {
      while(<$FH>) {last if $_ !~ /\S/}
    }
    elsif ($_ !~ /\S/)            {next}
    elsif ($_ =~ /Strand HSP/)    {next} # WU-BLAST non-data
    elsif ($_ =~ /^\s*Strand/)    {next} # NCBI-BLAST non-data
    elsif ($_ =~ /^\s*Score/)     {$self->{LASTLINE} = $_; last}
    elsif ($_ =~ /^>|^Parameters|^\s+Database:|^CPU\stime/)   {
      $self->{LASTLINE} = $_;
      $self->{PARENT}->{LASTLINE} = $_;
      $self->{HSP_ALL_PARSED} = 1;
      last;
    }
    else {
      push @hspline, $_;           #      store the query line
      my $l1 = <$FH>; push @hspline, $l1; # grab/store the alignment line
      my $l2 = <$FH>; push @hspline, $l2; # grab/store the sbjct line
    }
  }
  
  #########################
  # parse alignment lines #
  #########################
  my ($ql, $sl, $as) = ("", "", "");
  my ($qb, $qe, $sb, $se) = (0,0,0,0);
  my (@QL, @SL, @AS); # for better memory management
  
  for(my $i=0;$i<@hspline;$i+=3) {
    #warn $hspline[$i], $hspline[$i+2];
    $hspline[$i]   =~ /^Query:\s+(\d+)\s+(\S+)\s+(\d+)/;
    $ql = $2; $qb = $1 unless $qb; $qe = $3;
    
    my $offset = index($hspline[$i], $ql);
    $as = substr($hspline[$i+1], $offset, CORE::length($ql));
    
    $hspline[$i+2] =~ /^Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)/;
    $sl = $2; $sb = $1 unless $sb; $se = $3;
    
    push @QL, $ql; push @SL, $sl; push @AS, $as;
  }
  
  ##################
  # the HSP object #
  ##################
  # score 
  # bits
  # match
  # positive
  # percent
  # length of alignment ? 
  # p: p value
  # qb: query begin
  # qe: query end
  # sb: subject begin
  # se: subject end
  # ql: qery alignment
  # sl: subject alignment
  # as: alignment sequence
  $ql = join("", @QL);
  $sl = join("", @SL);
  $as = join("", @AS);
  my $hsp = new Bio::Tools::BPlite::HSP(
    -score=>$score,-bits=>$bits,-match=>$match,
    -positive=>$positive,-length=>$length,-p=>$p,
    -queryBegin=>$qb,-queryEnd=>$qe,-sbjctBegin=>$sb,
    -sbjctEnd=>$se,-queryLine=>$ql,-sbjctLine=>$sl,
    -homologyLine=>$as);
  return $hsp;
}
