###############################################################################
# Bio::Tools::BPlite::HSP
###############################################################################
# HSP = High Scoring Pair (to all non-experts as I am)
#
# The original BPlite.pm module has been written by Ian Korf !
# see http://sapiens.wustl.edu/~ikorf
#
# You may distribute this module under the same terms as perl itself


#
# BioPerl module for Bio::Tools::BPlite::HSP
#
# Cared for by Peter Schattner <schattner@alum.mit.edu>
#
# Copyright Peter Schattner
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::BPlite::HSP - Blast report High Scoring Pair (HSP)

=head1 SYNOPSIS

 use Bio::Tools::BPlite;
 my $report = new Bio::Tools::BPlite(-fh=>\*STDIN);
 {
    while(my $sbjct = $report->nextSbjct) {
	while (my $hsp = $sbjct->nextHSP) {
	    $hsp->score;
	    $hsp->bits;
	    $hsp->percent;
	    $hsp->P;
	    $hsp->match;
	    $hsp->positive;
	    $hsp->length;
	    $hsp->querySeq;
	    $hsp->sbjctSeq;
	    $hsp->homologySeq;
	    $hsp->query->start;
	    $hsp->query->end;
	    $hsp->hit->start;
	    $hsp->hit->end;
	    $hsp->hit->seq_id;
	    $hsp->hit->overlaps($exon);
	}
    }

    # the following line takes you to the next report in the stream/file
    # it will return 0 if that report is empty,
    # but that is valid for an empty blast report.
    # Returns -1 for EOF.

    last if ($report->_parseHeader == -1));

 redo
 }

=head1 DESCRIPTION

This object handles the High Scoring Pair data for a Blast report.
This is where the percent identity, query and hit sequence length,
P value, etc are stored and where most of the necessary information is located when building logic around parsing a Blast report.

See L<Bio::Tools::BPlite> for more detailed information on the entire
BPlite Blast parsing system.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org            - General discussion
http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Peter Schattner

Email: schattner@alum.mit.edu

=head1 CONTRIBUTORS

Jason Stajich, jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

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
    # interpreted from query/hit info 
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
    
    my ($score,$bits,$match,$hsplength,$positive,$gaps,$p,$exp,$qb,$qe,$sb,
	$se,$qs,$ss,$hs,$qname,$sname,$qlength,$slength,$qframe,$sframe,
	$blasttype) = 
	    $self->_rearrange([qw(SCORE
				  BITS
				  MATCH
				  HSPLENGTH
				  POSITIVE
				  GAPS				  
				  P
				  EXP
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
    
    $blasttype = 'UNKNOWN' unless $blasttype;
    $self->report_type($blasttype);
    # Determine strand meanings
    my ($queryfactor, $sbjctfactor) = (1,0); # default
    if ($blasttype eq 'BLASTP' || $blasttype eq 'TBLASTN' ) {
	$queryfactor = 0;
    }
    if ($blasttype eq 'TBLASTN' || $blasttype eq 'TBLASTX' || 
	$blasttype eq 'BLASTN' )  {
	$sbjctfactor = 1;
    }
    
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

    # store the aligned hit as sequence feature
    if ($se > $sb) {		# normal hit
	if ($sbjctfactor) { $strand = 1; } else { $strand = undef; }
	$self->hit( Bio::SeqFeature::Similarity->new
			(-start=>$sb, -end=>$se, -strand=>$strand,
			 -source=>"BLAST" ) ) }
    else { # reverse hit: start bigger than end
	if ($sbjctfactor) { $strand = -1; } else { $strand = undef; }
	$self->hit( Bio::SeqFeature::Similarity->new
			(-start=>$se, -end=>$sb, -strand=>$strand,
			 -source=>"BLAST" ) ) }
    
    # name the sequences
    $self->query->seq_id($qname); # query name
    $self->hit->seq_id($sname);   # hit name

    # set lengths
    $self->query->seqlength($qlength); # query length
    $self->hit->seqlength($slength);   # hit length

    # set object vars
    $self->score($score);
    $self->bits($bits);

    $self->significance($p);
    $self->{'EXP'} = $exp;
    
    $self->query->frac_identical($match);
    $self->hit->frac_identical($match);
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

=head2 EXP

 Title   : EXP
 Usage   : my $exp = $hsp->EXP;
 Function: returns the EXP value for the HSP
 Returns : string value
 Args    : none
 Note    : Patch provided by Sami Ashour for BTK parsing


=cut

sub EXP{
    return $_[0]->{'EXP'};
}


=head2 P

 Title    : P
 Usage    : $hsp->P();
 Function : returns the P (significance) value for a HSP 
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
 Function : returns the number of positive matches (symbols in the alignment
            with a positive score)
 Returns  : (int) number of positive matches in the alignment
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
	    if( ($1 eq '-' && $self->query->strand >= 0) || 
		($1 eq '+' && $self->query->strand <= 0) ) {
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
	    $self->warn("Specifying an invalid hit frame ($sframe)");
	    $sframe = undef;
	  } else { 
	      if( ($1 eq '-' && $self->hit->strand >= 0) || 
		  ($1 eq '+' && $self->hit->strand <= 0) ) 
	      {
		  $self->warn("Hit frame ($sframe) did not match strand of hit (". $self->hit->strand() . ")");
	      }
	      
	      # Set frame to GFF [0-2]
	      $sframe = $2 - 1;
	  }
	  $self->{'SFRAME'} = $sframe;
      }

    (defined $qframe && $self->SUPER::frame($qframe) && 
     ($self->{'FRAME'} = $qframe)) || 
    (defined $sframe && $self->SUPER::frame($sframe) && 
     ($self->{'FRAME'} = $sframe));

    if (wantarray() && 
	$self->{'BLAST_TYPE'} eq 'TBLASTX') 
    { 
	return ($self->{'QFRAME'}, $self->{'SFRAME'}); 
    } elsif (wantarray())  { 
	(defined $self->{'QFRAME'} && 
	 return ($self->{'QFRAME'}, undef)) || 
	     (defined $self->{'SFRAME'} && 
	      return (undef, $self->{'SFRAME'})); 
    } else { 
	(defined $self->{'QFRAME'} && 
	 return $self->{'QFRAME'}) || 
	(defined $self->{'SFRAME'} && 
	 return $self->{'SFRAME'}); 
    }
}

1;
