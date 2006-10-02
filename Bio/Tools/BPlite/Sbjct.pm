# $Id$
###############################################################################
# Bio::Tools::BPlite::Sbjct
###############################################################################
#
# The original BPlite.pm module has been written by Ian Korf !
# see http://sapiens.wustl.edu/~ikorf
#
# You may distribute this module under the same terms as perl itself


#
# BioPerl module for Bio::Tools::BPlite::Sbjct
#
# Cared for by Peter Schattner <schattner@alum.mit.edu>
#
# Copyright Peter Schattner
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::BPlite::Sbjct - A Blast Subject (database search Hit)

=head1 SYNOPSIS

  use Bio::Tools::BPlite;
  my $report = new Bio::Tools::BPlite(-fh=>\*STDIN);
  while(my $sbjct = $report->nextSbjct) {
      $sbjct->name;    # access to the hit name
      "$sbjct";        # overloaded to return name
      $sbjct->nextHSP; # gets the next HSP from the sbjct
      while (my $hsp = $sbjct->nextHSP) {
 	 # canonical form is again a while loop
      }
  }

=head1 DESCRIPTION

See L<Bio::Tools::BPlite> for a more detailed information about the
BPlite BLAST parsing objects.

The original BPlite.pm module has been written by Ian Korf!
See http://sapiens.wustl.edu/~ikorf

The Sbjct object encapsulates a Hit in a Blast database
search.  The Subjects are the "Hits" for a particular query.  A
Subject may be made up of multiple High Scoring Pairs (HSP) which are
accessed through the nextHSP method.

If you are searching for the P-value or percent identity that is
specific to each HSP and you will need to use the nextHSP method to
get access to that data.

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

=head1 AUTHOR - Peter Schattner

Email: schattner@alum.mit.edu

=head1 CONTRIBUTORS

Jason Stajich, jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tools::BPlite::Sbjct;

use strict;

use Bio::Tools::BPlite::HSP; # we want to use HSP
#use overload '""' => 'name';

use base qw(Bio::Root::Root);

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    ($self->{'NAME'},$self->{'LENGTH'},
     $self->{'PARENT'}) =
	 $self->_rearrange([qw(NAME
			       LENGTH
			       PARENT
			       )],@args);
    $self->report_type($self->{'PARENT'}->{'BLAST_TYPE'} || 'UNKNOWN');
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

=head2 report_type

 Title    : report_type
 Usage    : $type = $sbjct->report_type()
 Function : Returns the type of report from which this hit was obtained.
            This usually pertains only to BLAST and friends reports, for which
            the report type denotes what type of sequence was aligned against
            what (BLASTN: dna-dna, BLASTP prt-prt, BLASTX translated dna-prt, 
            TBLASTN prt-translated dna, TBLASTX translated dna-translated dna).
 Example  : 
 Returns  : A string (BLASTN, BLASTP, BLASTX, TBLASTN, TBLASTX, UNKNOWN)
 Args     : a string on set (you should know what you are doing)

=cut

sub report_type {
    my ($self, $rpt) = @_;
    if($rpt) {
	$self->{'_report_type'} = $rpt;
    }
    return $self->{'_report_type'};
}

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
  return  if $self->{'HSP_ALL_PARSED'};
  
  ############################
  # get and parse scorelines #
  ############################
  my ($qframe, $sframe);
  my $scoreline = $self->_readline();
  my $nextline = $self->_readline();
  return if not defined $nextline;
  $scoreline .= $nextline;
  my ($score, $bits);
  if ($scoreline =~ /\d bits\)/) {
      ($score, $bits) = ( $scoreline =~
			  /Score = (\d+) \((\S+) bits\)/); # WU-BLAST
  } else {
      ($bits, $score) = ( $scoreline =~
			  /Score =\s+(\S+) bits \((\d+)/); # NCBI-BLAST
  }
  unless( defined $bits && defined $score ) { 
      $self->warn("Weird scoreline ($scoreline) bailing\n");
      return;
  }
  my ($match, $hsplength) = ($scoreline =~ /Identities = (\d+)\/(\d+)/);
  my ($positive) = ($scoreline =~ /Positives = (\d+)/);
  my ($gaps) = ($scoreline =~ /Gaps = (\d+)/);
  if($self->report_type() eq 'TBLASTX') {
      ($qframe, $sframe) = $scoreline =~ /Frame =\s+([+-]\d)\s+\/\s+([+-]\d)/;
  } elsif ($self->report_type() eq 'TBLASTN')  {
      ($sframe) = $scoreline =~ /Frame =\s+([+-]\d)/;
  } else {
      ($qframe) = $scoreline =~ /Frame =\s+([+-]\d)/;
  }
  $positive = $match if not defined $positive;
  $gaps = '0' if not defined $gaps;
  my ($p)        = ($scoreline =~ /[Sum ]*P[\(\d+\)]* = (\S+)/);
  unless (defined $p) {(undef, $p) = $scoreline =~ /Expect(\(\d+\))? =\s+(\S+)/}
  my ($exp) = ($scoreline =~ /Expect(?:\(\d+\))? =\s+([^\s,]+)/);
  $exp = -1 unless( defined $exp );

  $self->throw("Unable to parse '$scoreline'") unless defined $score;
  
  #######################
  # get alignment lines #
  #######################
  my (@hspline);
  local $_;
  while( defined($_ = $self->_readline()) ) {
      if (/^WARNING:|^NOTE:/) {
	  while(defined($_ = $self->_readline())) {last if $_ !~ /\S/}
      }
      elsif ( ! /\S/o)         {next}
      elsif (/Strand HSP/o)    {next} # WU-BLAST non-data
      elsif (/^\s*Strand/o)    {next} # NCBI-BLAST non-data
      elsif (/^\s*Score/o)     {$self->_pushback($_); last}

      elsif (/^>|^Histogram|^Searching|^Parameters|^\s+Database:|^CPU\stime|^\s*Lambda|^\s+Subset/o)   
      {    
	  #ps 5/28/01	
	  # elsif ($_ =~ /^>|^Parameters|^\s+Database:|^CPU\stime/)   {
	  $self->_pushback($_);

	  $self->{'HSP_ALL_PARSED'} = 1;
	  last;
      } elsif( /^BLAST/ ) {
	  $self->_pushback($_);
	  $self->{'HSP_ALL_PARSED'} = 1;
	  last;
      } elsif( $_ =~ /^\s*Frame/ ) {
	  if ($self->report_type() eq 'TBLASTX') {
	      ($qframe, $sframe) = $_ =~ /Frame = ([\+-]\d)\s+\/\s+([\+-]\d)/;
	  } elsif ($self->report_type() eq 'TBLASTN') {
	      ($sframe) = $_ =~ /Frame = ([\+-]\d)/;
	  } else {
	      ($qframe) = $_ =~ /Frame = ([\+-]\d)/;
	  }
      }
      else {
	  push @hspline, $_;	#      store the query line
	  $nextline = $self->_readline();
	  # Skip "pattern" line when parsing PHIBLAST reports, otherwise store the alignment line
	  my $l1 = ($nextline =~ /^\s*pattern/) ? $self->_readline() : $nextline;
	  push @hspline, $l1;	# store the alignment line
	  my $l2 = $self->_readline(); push @hspline, $l2; # grab/store the sbjct line
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
    $hspline[$i]   =~ /^(?:Query|Trans):\s+(\d+)\s*([\D\S]+)\s+(\d+)/o;
    $ql = $2; $qb = $1 unless $qb; $qe = $3;
    
    my $offset = index($hspline[$i], $ql);
    $as = substr($hspline[$i+1], $offset, CORE::length($ql));
    
    $hspline[$i+2] =~ /^Sbjct:\s+(\d+)\s*([\D\S]+)\s+(\d+)/o;
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
	$qname   = $self->{'PARENT'}->query;
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
       '-exp'        => $exp,
       '-queryBegin' => $qb, 
       '-queryEnd'   => $qe, 
       '-sbjctBegin' => $sb,
       '-sbjctEnd'   => $se, 
       '-querySeq'   => $ql, 
       '-sbjctSeq'   => $sl,
       '-homologySeq'=> $as, 
       '-queryName'  => $qname,
#			'-queryName'=>$self->{'PARENT'}->query,
       '-sbjctName'  => $self->{'NAME'},
       '-queryLength'=> $qlength,
#		       	'-queryLength'=>$self->{'PARENT'}->qlength,
       '-sbjctLength'=> $self->{'LENGTH'},
       '-queryFrame' => $qframe,
       '-sbjctFrame' => $sframe,
       '-blastType'  => $self->report_type());
  return $hsp;
}

=head2 _readline

 Title   : _readline
 Usage   : $obj->_readline
 Function: Reads a line of input.

           Note that this method implicitely uses the value of $/ that is
           in effect when called.

           Note also that the current implementation does not handle pushed
           back input correctly unless the pushed back input ends with the
           value of $/.
 Example :
 Returns : 

=cut

sub _readline{
   my ($self) = @_;
   return $self->{'PARENT'}->_readline();
}

=head2 _pushback

 Title   : _pushback
 Usage   : $obj->_pushback($newvalue)
 Function: puts a line previously read with _readline back into a buffer
 Example :
 Returns :
 Args    : newvalue

=cut

sub _pushback {
   my ($self, $arg) = @_;   
   return $self->{'PARENT'}->_pushback($arg);    
}

1;
