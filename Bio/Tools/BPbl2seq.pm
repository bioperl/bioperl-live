# $Id$
#
# Bioperl module Bio::Tools::BPbl2seq
#	based closely on the Bio::Tools::BPlite modules
#	Ian Korf (ikorf@sapiens.wustl.edu, http://sapiens.wustl.edu/~ikorf), 
#	Lorenz Pollak (lorenz@ist.org, bioperl port)
#
#
# Copyright Peter Schattner
#
# You may distribute this module under the same terms as perl itself
# _history
# October 20, 2000
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::BPbl2seq - Lightweight BLAST parser for pair-wise sequence 
alignment using the BLAST algorithm.

=head1 SYNOPSIS

use Bio::Tools::BPbl2seq;
my $report = Bio::Tools::BPbl2seq->new(-file => 't/bl2seq.out');
 $report->query;
 $report->score;
 $report->bits;
 $report->percent;
 $report->P;
 $report->match;
 $report->positive;
 $report->length;
 $report->querySeq;
 $report->sbjctSeq;
 $report->homologySeq;
 $report->subject->start;
 $report->subject->end;
 $report->subject->seqname;

=head1 DESCRIPTION

BPbl2seq is a package for parsing BLAST bl2seq reports. BLAST bl2seq is a
program for comparing and aligning two sequences using BLAST.  Although
the report format is similar to that of a conventional BLAST, there are a
few differences so that the standard bioperl BLAST parsers Blast.pm and
BPlite are unable to read bl2seq reports directly.

From the user's perspective, the main difference between bl2seq and
other blast reports is that the bl2seq report does not print out the
name of the first of the two aligned sequences.  (The second sequence
name is given in the report as the name of the "hit").  Consequently,
BPbl2seq has no way of identifying the name of the initial sequence
unless it is passed to constructor as a second argument as in:

	my $report = Bio::Tools::BPbl2seq->new(\*FH, "ALEU_HORVU");

If the name of the first sequence is not passed to BPbl2seq.pm in this
manner, the name of the first sequence will be left as "unknown".
(Note that to preserve a common interface with the other BLAST
programs the two sequences being compared are referred to in bl2seq as
"query" and "subject" although this is perhaps a bit misleading when
simply comparing 2 sequences as opposed to querying a database.)


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.  Bug reports can be submitted via
 email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Peter Schattner

Email: schattner@alum.mit.edu

Based on work of:
Ian Korf (ikorf@sapiens.wustl.edu, http://sapiens.wustl.edu/~ikorf), 
Lorenz Pollak (lorenz@ist.org, bioperl port)

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT

BPlite.pm is copyright (C) 1999 by Ian Korf. 

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut
#'
package Bio::Tools::BPbl2seq;

use vars qw(@ISA);
use strict;

# to enable overloading uncomment the following line:
#use overload '""' => '_overload';

# Object preamble - inherits from Bio::SeqFeature::SimilarityPair

use Bio::SeqFeature::SimilarityPair;
use Bio::SeqFeature::Similarity;

@ISA = qw(Bio::SeqFeature::SimilarityPair);

=head2 new

 Title   : new
 Function: Create a new Bio::Tools::BPbl2seq object
 Returns : Bio::Tools::BPbl2seq
 Args    : -file     input file (alternative to -fh)
           -fh       input stream (alternative to -file)
           -queryname    name of query sequence
=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($file, $fh, $query) = $self->_rearrange([qw(FILE FH QUERYNAME)], @args);
    $query = 'unknown' if( ! defined $query );

    if( defined $file && defined $fh ) {
	$self->throw("Cannot define both a file and fh for input");
    }
    if( defined $file ) {
	$fh = Symbol::gensym();
	open ($fh,$file) || $self->throw("Could not open file [$file] $!");
    } elsif( defined $fh ) {
	if (ref $fh !~ /GLOB/)
	{ $self->throw("Expecting a GLOB reference, not $fh!"); }
    }
    $self->fh($fh);
    my ($score,$bits,$match,$positive,$p,$qb,$qe,$sb,$se,$qs,
	$ss,$hs,$qname,$sname,$qlength,$slength) = $self->_parsebl2seq($query);
    unless ( $positive ) {
	$self->throw("No match found or other problem parsing bl2seq report");
    }

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
    $self->{'PERCENT'} = int((1000 * $match)/ $self->query->length)/10;
    $self->{'POSITIVE'} = $positive;
    $self->{'QS'} = $qs;
    $self->{'SS'} = $ss;
    $self->{'HS'} = $hs;
    return $self;
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
 Example  : 
 Returns  : (double) percent matching
 Args     :

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
 Example  : 
 Returns  : (int) number of positive residue hits 
 Args     :

=cut

sub positive        {shift->{'POSITIVE'}}

=head2 querySeq

 Title    : querySeq
 Usage    : $hsp->querySeq();
 Function : returns the query sequence
 Example  : 
 Returns  : (string) the Query Sequence 
 Args     :

=cut

sub querySeq        {shift->{'QS'}}

=head2 sbjctSeq

 Title    : sbjctSeq
 Usage    : $hsp->sbjctSeq();
 Function : returns the Sbjct sequence 
 Example  : 
 Returns  : (string) the Sbjct Sequence 
 Args     :

=cut

sub sbjctSeq        {shift->{'SS'}}

=head2 homologySeq

 Title    : homologySeq
 Usage    : $hsp->homologySeq();
 Function : returns the homologous sequence 
 Example  : 
 Returns  : (string) homologous sequence 
 Args     :

=cut

sub homologySeq     {shift->{'HS'}}

=head2 qs

 Title    : qs
 Usage    : $hsp->qs();
 Function : returns the Query Sequence (same as querySeq)
 Example  : 
 Returns  : (string) query Sequence 
 Args     :

=cut

sub qs              {shift->{'QS'}}

=head2 ss

 Title    : ss
 Usage    : $hsp->ss();
 Function : returns the subject sequence ( same as sbjctSeq) 
 Example  : 
 Returns  : (string) Sbjct Sequence
 Args     :

=cut

sub ss              {shift->{'SS'}}

=head2 hs

 Title    : hs
 Usage    : $hsp->hs();
 Function : returns the Homologous Sequence (same as homologySeq ) 
 Example  : 
 Returns  : (string) Homologous Sequence
 Args     :

=cut

sub hs              {shift->{'HS'}}

sub _parsebl2seq {
  my ($self,$query) = @_;
  my $def = "";

  ############################
  # get seq2 (the "hit") name & lrngth  
  ############################
  my $FH = $self->fh;
  while(<$FH>) {
    if    ($_ !~ /\w/)            {next}
    elsif ($_ =~ /^\s*Length/) {$def .= $_; last}
    else                          {$def .= $_}
  }
  $def =~ s/\s+/ /g;
  $def =~ s/\s+$//g;
  $def =~ s/Length = ([\d,]+)$//g;
  my $hlength = $1;
  return 0 unless $def =~ />/;
  $def =~ s/(.*)>//;

  ############################
  # get and parse scorelines #
  ############################
  my $nextline;

 BLANKS: while (defined($nextline = <$FH>)) {
     last BLANKS if ($nextline =~ /\w/) ;
 }
  return undef if not defined $nextline;
  my $scoreline = $nextline;
  my ($score, $bits);
  if ($scoreline =~ /\d bits\)/) {
    ($score, $bits) = $scoreline =~
      /Score = (\d+) \((\S+) bits\)/; # WU-BLAST
  }
  else {
    ($bits, $score) = $scoreline =~
      /Score =\s+(\S+) bits \((\d+)/; # NCBI-BLAST
  }
  my ($p)        = $scoreline =~ /[Sum ]*P[\(\d+\)]* = (\S+)/;
  if (not defined $p) {($p) = $scoreline =~ /Expect = (\S+)/}

# Read second score line
  my $scoreline2 = <$FH> ;
  my ($match, $length) = $scoreline2 =~ /Identities = (\d+)\/(\d+)/;
  my ($positive) = $scoreline2 =~ /Positives = (\d+)/;
  $positive = $match if not defined $positive;

  $self->throw("Unable to parse $scoreline") if not defined $score;
  $self->throw("Unable to parse $scoreline2") if not defined $match;

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
    elsif ($_ =~ /^>|^Parameters|^\s+Database:|^CPU\stime:|^\s*Lambda/)   {
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
  # return the parsed data
  ##################
  $ql = join("", @QL);
  $sl = join("", @SL);
  $as = join("", @AS);
  my ($queryid, $querylength);
  
  if  (!defined $query || 
       $query eq 'unknown') { $queryid  = 'unknown'; $querylength = 0;}
  elsif( $query->can('id') && $query->can('length') ) {
	$queryid  =  $query->id;
	$querylength = $query->length;
  }
  my @returnarray = ($score,  $bits,   $match, $positive,  $p, $qb,  $qe,
			  $sb, $se,  $ql,  $sl, $as,  $queryid, $def,
			 $querylength, $hlength);
  return @returnarray;
}

=head2 fh

 Title    : fh
 Usage    : do not use
 Function : Get/Set method for filehandle
 Example  : my $fh = $self->fh
 Returns  : filehandle reference
 Args     : filehandle reference [optional]

=cut
sub fh {
    my ($self, $value) = @_;    
    if( defined $value && ref($value) =~ /GLOB/i ) {
	$self->{'FH'} = $value;
    } 
    return $self->{'FH'};
}
1;
