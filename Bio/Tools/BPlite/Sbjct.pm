# $Id$
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

use Bio::Root::RootI;        # root object to inherit from
use Bio::Tools::BPlite::HSP; # we want to use HSP
#use overload '""' => 'name';
use vars qw(@ISA);

@ISA = qw(Bio::Root::RootI);

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    ($self->{'NAME'},$self->{'LENGTH'},$self->{'FH'},
     $self->{'LASTLINE'},$self->{'PARENT'}) =
	 $self->_rearrange([qw(NAME
			       LENGTH
			       FH
			       LASTLINE
			       PARENT
			       )],@args);
    
    $self->{'HSP_ALL_PARSED'} = 0;
    
  return $self;
}

=head2 name

 Title    : name
 Usage    : $name = $obj->name();
 Function : returns the name of the Sbjct 
 Example  : 
 Returns  : name of the Sbjct 
 Args     :

=cut

sub name {shift->{'NAME'}}

=head2 nextFeaturePair

 Title    : nextFeaturePair
 Usage    : $name = $obj->nextFeaturePair();
 Function : same as the nextHSP function 
 Example  : 
 Returns  : next FeaturePair 
 Args     :

=cut

sub nextFeaturePair {shift->nextHSP}; # just another name

=head2 nextHSP

 Title    : nextHSP
 Usage    : $hsp = $obj->nextHSP();
 Function : returns the next available High Scoring Pair
 Example  : 
 Returns  : Bio::Tools::HSP  or null if finished
 Args     :

=cut

sub nextHSP {
  my ($self) = @_;
  return undef if $self->{'HSP_ALL_PARSED'};
  
  ############################
  # get and parse scorelines #
  ############################

  my $scoreline = $self->{'LASTLINE'};
  my $FH = $self->{'FH'};
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
  
  my ($match, $hsplength) = ($scoreline =~ /Identities = (\d+)\/(\d+)/);
  my ($positive) = ($scoreline =~ /Positives = (\d+)/);
  my ($gaps) = ($scoreline =~ /Gaps = (\d+)/);  
  my $frame = '0';
  $positive = $match if not defined $positive;
  $gaps = '0' if not defined $gaps;
  my ($p)        = $scoreline =~ /[Sum ]*P[\(\d+\)]* = (\S+)/;
  if (not defined $p) {($p) = $scoreline =~ /Expect =\s+(\S+)/}
  
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
    elsif ($_ =~ /^\s*Score/)     {$self->{'LASTLINE'} = $_; last}
    elsif ($_ =~ /^>|^Parameters|^\s+Database:|^CPU\stime|^\s*Lambda/)   {    #ps 5/28/01
#    elsif ($_ =~ /^>|^Parameters|^\s+Database:|^CPU\stime/)   {
      $self->{'LASTLINE'} = $_;
      $self->{'PARENT'}->{'LASTLINE'} = $_;
      $self->{'HSP_ALL_PARSED'} = 1;
      last;
    }
    elsif( $_ =~ /^\s*Frame\s*=\s*([-\+]\d+)/ ) {
	$frame = $1;
    }

    else {
      push @hspline, $_;           #      store the query line
      $nextline = <$FH> ;
# Skip "pattern" line when parsing PHIBLAST reports, otherwise store the alignment line
      my $l1 = ($nextline =~ /^\s*pattern/) ? <$FH> : $nextline;
#     my $l1 =  push @hspline, $l1; # grab/store the alignment line
      push @hspline, $l1; # store the alignment line
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
    # warn $hspline[$i], $hspline[$i+2];
    $hspline[$i]   =~ /^Query:\s+(\d+)\s*([\D\S]+)\s+(\d+)/;
    $ql = $2; $qb = $1 unless $qb; $qe = $3;
    
    my $offset = index($hspline[$i], $ql);
    $as = substr($hspline[$i+1], $offset, CORE::length($ql));
    
    $hspline[$i+2] =~ /^Sbjct:\s+(\d+)\s*([\D\S]+)\s+(\d+)/;
    $sl = $2; $sb = $1 unless $sb; $se = $3;
    
    push @QL, $ql; push @SL, $sl; push @AS, $as;
  }
  
  ##################
  # the HSP object #
  ##################
  $ql = join("", @QL);
  $sl = join("", @SL);
  $as = join("", @AS);
# Query name and length are not in the report for a bl2seq report so {'PARENT'}->query and
# {'PARENT'}->qlength will not be available.
  my ($qname, $qlength) = ('unknown','unknown');
  if ($self->{'PARENT'}->can('query')) {
	$qname = $self->{'PARENT'}->query;
	$qlength = $self->{'PARENT'}->qlength;
  }	
  my $hsp = new Bio::Tools::BPlite::HSP
      ('-score'      => $score, 
       '-bits'       => $bits, 
       '-match'      => $match,
       '-positive'   => $positive, 
       '-gaps'       => $gaps,
       '-hsplength'  => $hsplength,
       '-p'          => $p,
       '-queryBegin' => $qb, 
       '-queryEnd'   => $qe, 
       '-sbjctBegin' => $sb,
       '-sbjctEnd'   => $se, 
       '-querySeq'   => $ql, 
       '-sbjctSeq'   => $sl,
       '-homologySeq'=> $as, 
       '-queryName'  => $qname,
#					'-queryName'=>$self->{'PARENT'}->query,
       '-sbjctName'  => $self->{'NAME'},
       '-queryLength'=> $qlength,
#					'-queryLength'=>$self->{'PARENT'}->qlength,
       '-sbjctLength'=> $self->{'LENGTH'},
       '-frame'      => $frame);
  return $hsp;
}

1;
