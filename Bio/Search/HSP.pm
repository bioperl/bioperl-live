# $Id$
#
# BioPerl module for Bio::Search::HSP
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::HSP - Implementation of Bio::Search::HSPI to handle High Scoring Pairs in memory

=head1 SYNOPSIS

    # typically one gets a HSP From a SearchIO report or Subject Object
    use Bio::Search::HSP;
    my $hsp = new Bio::Search::HSP
       (
	'-report_type'   => 'BLASTP',
	'-score'         => '101',
	'-bits'          => '15.3',
	'-match'         => 19,
	'-hsp_length'    => 212,
	'-positive'      => 30,
	'-gaps'          => 10,
	'-evalue'        => 0.2,
	'-query_begin'   => 20,
	'-query_end'     => 252,
	'-subject_begin' => 37,
	'-subject_end'   => 190,
	'-query_seq'     => 'MVTYW',
	'-subject_seq'   => 'MVXX-A',
	'-homology_seq'  => 'MV    ',
	'-query_length'  => '232',
	'-subject_length'=> '153',
	'-query_name'    => 'seqa',
	'-subject_name'  => 'seqb',
	'-query_frame'   => '0',
	'-subject_frame' => '0',
	);

    print (join(',',( $hsp->report_type,
		      $hsp->score,
		      $hsp->bits,

		      $hsp->P,
		      $hsp->evalue,
		      $hsp->percent_identity,
		      $hsp->gaps,
		      $hsp->positive,

		      $hsp->hsp_length,
		      $hsp->query->start,
		      $hsp->query->end,
		      $hsp->query->strand,
		      $hsp->query->strand,
		      $hsp->query->frame,
		      $hsp->query->length,
		      $hsp->subject->start,
		      $hsp->subject->end,
		      $hsp->subject->strand,
		      $hsp->subject->frame,
		      $hsp->subject->length,
		      $hsp->query_seq,
		      $hsp->subject_seq,
		      $hsp->homology_seq,
		      ))), "\n";

=head1 DESCRIPTION

This object provides an in-memory copy of a High Scoring Pair
implementation and appropriately initializes values like
percent_identity from the values in identity and hsp_length.

An HSP is-a FeaturePair so it has methods query() and subject() which
provides access to the paired features that make up a HSP().

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
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::HSP;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::Search::HSPI;
use Bio::SeqFeature::Similarity;

@ISA = qw(Bio::Root::Root Bio::Search::HSPI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Search::HSP();
 Function: Builds a new Bio::Search::HSP object
 Returns : Bio::Search::HSP
 Args    :
           -query_seq   => query sequence
           -subject_seq => subject sequence
           -homology_seq=> homology sequence
           -hsp_length  => hsp length
           -score       => HSP Z score
           -bits        => HSP bit score
           -P           => P value
           -match       => # of bases that matched exactly
           -positive    => # of positive matches
           -gaps        => # of gaps
           -query_begin   => start of HSP in query coordinates
           -query_end     => end of HSP in query coordinates
           -subject_begin => start of HSP in subject coordinates
           -subject_end   => end of HSP in subject coordinates
           -query_length  => total length of query sequence
           -subject_length=> total length of subject sequence
           -query_frame   => query frame (if any)
           -subject_frame => subject frame (if any)
           -report_type   => Type of report ([t]blast[pnx] or hmmer, etc)

=cut

sub new {
    my($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
    my ($score,$bits,$match,$hsplength,$positive,$gaps,$e,$qb,$qe,$sb,$se,$qs,
	$ss,$hs,$qname,$sname,$qlength,$slength,$qframe,$sframe,$reporttype) =
	    $self->_rearrange([qw(SCORE
				  BITS
				  MATCH
				  HSP_LENGTH
				  POSITIVE
				  GAPS
				  EVALUE
				  QUERY_BEGIN
				  QUERY_END
				  SUBJECT_BEGIN
				  SUBJECT_END
				  QUERY_SEQ
				  SUBJECT_SEQ
				  HOMOLOGY_SEQ
				  QUERY_NAME
				  SUBJECT_NAME
				  QUERY_LENGTH
				  SUBJECT_LENGTH
				  QUERY_FRAME
				  SUBJECT_FRAME
				  REPORT_TYPE
				  )],@args);

    $reporttype = 'SEARCHREPORT' unless $reporttype;
    $self->{'_reporttype'} = $reporttype;
    # Determine strand meanings
    my ($queryfactor, $sbjctfactor) = (1,0); # default
    
    if ($reporttype eq 'BLASTP' || $reporttype eq 'TBLASTN'
	|| $reporttype eq 'FASTX' || $reporttype eq 'FASTY') {
	$queryfactor = 0;
    } elsif ($reporttype eq 'TBLASTN' || $reporttype eq 'TBLASTX' ||
	$reporttype eq 'BLASTN' || $reporttype eq 'TFASTX' ||
	$reporttype eq 'FASTA' )  {
	$sbjctfactor = 1;
    } elsif( $reporttype eq 'RPSBLAST' ) {
	$queryfactor = $sbjctfactor = 0;
	$qframe = $sframe = 0;
    }

    # Store the aligned query as sequence feature
    my $strand;
    if( ! $qe || ! $qb ) { $self->throw("Did not specify a Query End or Query Begin"); }
    if ($qe > $qb) {		# normal query: start < end
	if ($queryfactor) { $strand = 1; } else { $strand = undef; }
	$self->query( Bio::SeqFeature::Similarity->new
		      ('-start' => $qb,
		       '-end'   => $qe,
		       '-strand'=> $strand,
		       '-source'=> $reporttype ) ) }
    else {			# reverse query (i dont know if this is possible, but feel free to correct)
	if ($queryfactor) { $strand = -1; } else { $strand = undef; }
	$self->query( Bio::SeqFeature::Similarity->new
		      ('-start' => $qe,
		       '-end'   => $qb,
		       '-strand'=> $strand,
		       '-source'=> $reporttype ) );
    }
    $qframe = 0 unless defined $strand;
    # store the aligned subject as sequence feature
    if ($se > $sb) {		# normal subject
	if ($sbjctfactor) { $strand = 1; } else { $strand = undef; }
	$self->subject( Bio::SeqFeature::Similarity->new
			('-start' => $sb,
			 '-end'   => $se,
			 '-strand'=> $strand,
			 '-source'=> $reporttype ) ) }
    else { # reverse subject: start bigger than end
	if ($sbjctfactor) { $strand = -1; } else { $strand = undef; }
	$self->subject( Bio::SeqFeature::Similarity->new
			('-start' => $se,
			 '-end'   => $sb,
			 '-strand'=> $strand,
			 '-source'=> $reporttype ) );
    }
    $sframe = 0 unless defined $strand;
    # name the sequences
    $self->query->seqname($qname); # query
    $self->subject->seqname($sname); # subject

    # set lengths
    $self->query->seqlength($qlength); # query
    $self->subject->seqlength($slength); # subject

    # set object vars
    $self->score($score);
    $self->bits($bits);
    $self->significance($e);
    $self->query->frac_identical($match);
    $self->subject->frac_identical($match);
    $self->{'_hsp_length'} = $hsplength;
    if( $hsplength ) {
	$self->{'_percent_id'} = int((1000 * $match)/$hsplength)/10;
    } else {
	$self->{'_percent_id'} = 0;
    }
    $self->{'_positive'} = $positive;
    $self->{'_gaps'} = defined $gaps ? $gaps : 0;
    $self->{'_query_seq'} = $qs;
    $self->{'_subject_seq'} = $ss;
    $self->{'_homol_seq'} = $hs;

    $sframe = 0 unless $sframe;
    $qframe = 0 unless $qframe;
    $self->frame($qframe, $sframe);
    return $self;
}

=head1 Bio::Search::HSPI methods implemented

=head2 report_type

 Title   : report_type
 Usage   : my $r_type = $hsp->report_type
 Function: Obtain the report type for an HSP
 Returns : string
 Args    : none


=cut

sub report_type{
   my ($self) = @_;
   return $self->{'_reporttype'} ;
}

=head2 P

 Title   : P
 Usage   : my $pvalue = $hsp->P();
 Function: Returns the P-value for this HSP (same as evalue)
 Returns : float or exponential (2e-10)
 Args    : none

=head2 evalue

 Title   : evalue
 Usage   : my $pvalue = $hsp->evalue();
 Function: Returns the evalue value for this HSP
 Returns : float or exponential (2e-10)
 Args    : none

=cut

sub evalue {
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

=head2 percent_identity

 Title   : percent_identity (alias for percent_identical)

 Usage   : my $perc_id = $hsp->percent_identity();
 Function: Returns the Percent Identity (num positive / total HSP len)
           for this HSP
 Returns : Float in range 0 -> 100
 Args    : none

=cut

sub percent_identity { $_[0]->percent_identical }

=head2 percent_identical

 Title   : percent_identical
 Usage   : my $perc_id = $hsp->percent_identical();
 Function: Returns the Percent Identity (num positive / total HSP len)
           for this HSP
 Returns : Float in range 0 -> 100
 Args    : none

=cut

sub percent_identical {
   my ($self) = @_;
   return $self->{'_percent_id'};
}

=head2 positive

 Title    : positive
 Usage    : $hsp->positive();
 Function : returns the number of positive matches (symbols in the alignment
            with a positive score)
 Returns  : (int) number of positive matches in the alignment
 Args     : none

=cut

sub positive        {
    my ($self) = @_;
    return $self->{'_positive'};
}

=head2 gaps

 Title    : gaps
 Usage    : $hsp->gaps();
 Function : returns the number of gaps or 0 if none
 Returns  : (int) number of gaps or 0 if none
 Args     : none

=cut

sub gaps        {
    my ($self) = @_;
    return $self->{'_gaps'};
}

=head2 query_seq

 Title   : query_seq
 Usage   : my $qseq = $hsp->query_seq;
 Function: Retrieves the query sequence that is part of this HSP
 Returns : string
 Args    : none


=cut

sub query_seq{
   my ($self,@args) = @_;
   return $self->{'_query_seq'};
}

=head2 subject_seq

 Title   : subject_seq
 Usage   : my $sseq = $hsp->subject_seq;
 Function: Retrieves the subject sequence that is part of this HSP
 Returns : string
 Args    : none


=cut

sub subject_seq{
   my ($self,@args) = @_;
   return $self->{'_subject_seq'};
}

=head2 homology_seq


 Title   : homology_seq
 Usage   : my $homo_seq = $hsp->homology_seq;
 Function: Retrieves the homology sequence for this HSP
 Returns : string
 Args    : none

=cut

sub homology_seq{
   my ($self,@args) = @_;
   return $self->{'_homol_seq'};
}

=head2 hsp_length

 Title   : hsp_length
 Usage   : my $len = $hsp->hsp_length
 Function: Returns the aggregate length of the HSP
           (which may be greater than either subject or query )
 Returns : integer
 Args    : none

=cut

sub hsp_length{
   my ($self,@args) = @_;
   return $self->{'_hsp_length'};
}

=head2 frame

 Title   : frame
 Usage   : $hsp->frame($queryframe,$subjectframe)
 Function: Set the Frame for both query and subject and insure that
           they agree.
           This overrides the frame() method implementation in
           FeaturePair.
 Returns : array of query and subjects if return type wants an array
           or query frame if defined or subject frame
 Args    : none
 Note    : Frames are stored in the GFF way (0-2 +/-) not 1-3
           as they are in BLAST

=cut

sub frame {
    my ($self, $qframe, $sframe) = @_;
    if( defined $qframe ) {
	if( $qframe == 0 ) {
	    $qframe = 0;
	} elsif( $qframe !~ /^([+-])?([1-3])/ ) {
	    $self->warn("Specifying an invalid query frame ($qframe)");
	    $qframe = undef;
	} else {
	    my $dir = $1;
	    $dir = '+' unless defined $dir;
	    if( ($dir eq '-' && $self->query->strand >= 0) ||
		($dir eq '+' && $self->query->strand <= 0) ) {
		$self->warn("Query frame ($qframe) did not match strand of query (". $self->query->strand() . ")");
	    }
	    # Set frame to GFF [0-2] -
	    # what if someone tries to put in a GFF frame!
	    $qframe = $2 - 1;
	}
	$self->query->frame($qframe);
    }
    if( defined $sframe ) {
	  if( $sframe == 0 ) {
	    $sframe = 0;
	  } elsif( $sframe !~ /^([+-])?([1-3])/ ) {
	    $self->warn("Specifying an invalid subject frame ($sframe)");
	    $sframe = undef;
	  } else {
	      if( ($1 eq '-' && $self->subject->strand >= 0) ||
		  ($1 eq '+' && $self->subject->strand <= 0) )
	      {
		  $self->warn("Subject frame ($sframe) did not match strand of subject (". $self->subject->strand() . ")");
	      }

	      # Set frame to GFF [0-2]
	      $sframe = $2 - 1;
	  }
	  $self->subject->frame($sframe);
      }

    if (wantarray() &&
	$self->report_type eq 'TBLASTX')
    {
	return ($self->query->frame(), $self->subject->frame());
    } elsif (wantarray())  {
	($self->query->frame() &&
	 return ($self->query->frame(), undef)) ||
	     ($self->subject->frame() &&
	      return (undef, $self->subject->frame()));
    } else {
	($self->query->frame() &&
	 return $self->query->frame()) ||
	($self->subject->frame() &&
	 return $self->subject->frame());
    }
}

=head2 Bio::SeqFeature::SimilarityPair methods

=cut

=head2 significance

 Title   : significance
 Usage   : $evalue = $obj->significance();
           $obj->significance($evalue);
 Function:
 Returns :
 Args    :

=head2 bits

 Title   : bits
 Usage   : $bits = $obj->bits();
           $obj->bits($value);
 Function:
 Returns :
 Args    :

=head2 score

 Title   : score
 Usage   : $score = $obj->score();
           $obj->score($value);
 Function:
 Returns :

=head2 Bio::SeqFeature::FeaturePair methods

=cut

=head2 query

 Title   : query
 Usage   : my $query = $hsp->query;
 Function: Access to the SeqFeature::Similarity
           for the query sequence that makes up this HSP
 Returns : Bio::SeqFeature::Similarity
 Args    : none

=head2 subject

 Title   : subject
 Usage   : my $subject = $hsp->subject
 Function: Access to the SeqFeature::Similarity
           for the subject (hit) sequence that makes up this HSP
 Returns : Bio::SeqFeature::Similarity
 Args    : none

=head2 Bio::SeqFeature::FeaturePair methods

=cut

=head2 feature1

 Title   : feature1
 Usage   : $f = $featpair->feature1
           $featpair->feature1($feature)
 Function: Get/set for the query feature
 Returns : Bio::SeqFeatureI
 Args    : Bio::SeqFeatureI

=head2 feature2

 Title   : feature2
 Usage   : $f = $featpair->feature2
           $featpair->feature2($feature)
 Function: Get/set for the hit feature
 Returns : Bio::SeqFeatureI
 Args    : Bio::SeqFeatureI

=head2 hseqname

 Title   : hseqname
 Usage   : $featpair->hseqname($newval)
 Function: Get/set method for the name of
           feature2.
 Returns : value of $feature2->seqname
 Args    : newvalue (optional)

=head2 hstart

 Title   : hstart
 Usage   : $start = $featpair->hstart
           $featpair->hstart(20)
 Function: Get/set on the start coordinate of feature2
 Returns : integer
 Args    : none

=head2 hend

 Title   : hend
 Usage   : $end = $featpair->hend
           $featpair->hend($end)
 Function: get/set on the end coordinate of feature2
 Returns : integer
 Args    : none

=head2 hstrand

 Title   : hstrand
 Usage   : $strand = $feat->strand()
           $feat->strand($strand)
 Function: get/set on strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : none

=head2 hscore

 Title   : hscore
 Usage   : $score = $feat->score()
           $feat->score($score)
 Function: get/set on score information
 Returns : float
 Args    : none if get, the new value if set

=head2 hframe

 Title   : hframe
 Usage   : $frame = $feat->frame()
           $feat->frame($frame)
 Function: get/set on frame information
 Returns : 0,1,2
 Args    : none if get, the new value if set


=head2 hprimary_tag

 Title   : hprimary_tag
 Usage   : $ptag = $featpair->hprimary_tag
 Function: Get/set on the primary_tag of feature2
 Returns : 0,1,2
 Args    : none if get, the new value if set

=head2 hsource_tag

 Title   : hsource_tag
 Usage   : $tag = $feat->hsource_tag()
           $feat->source_tag('genscan');
 Function: Returns the source tag for a feature,
           eg, 'genscan'
 Returns : a string
 Args    : none

=head2 invert

 Title   : invert
 Usage   : $tag = $feat->invert
 Function: Swaps feature1 and feature2 around
 Returns : Nothing
 Args    : none

=head2 Bio::SeqFeatureI methods

=cut

=head2 start

 Title   : start
 Usage   : $start = $featpair->start
           $featpair->start(20)
 Function: Get/set on the start coordinate of feature1
 Returns : integer
 Args    : [optional] beginning of feature

=head2 end

 Title   : end
 Usage   : $end = $featpair->end
           $featpair->end($end)
 Function: get/set on the end coordinate of feature1
 Returns : integer
 Args    : [optional] ending point of feature

=head2 strand

 Title   : strand
 Usage   : $strand = $feat->strand()
           $feat->strand($strand)
 Function: get/set on strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : [optional] strand information to set

=head2 location

 Title   : location
 Usage   : $location = $featpair->location
           $featpair->location($location)
 Function: Get/set location object (using feature1)
 Returns : Bio::LocationI object
 Args    : [optional] LocationI to store

=head2 primary_tag

 Title   : primary_tag
 Usage   : $ptag = $featpair->primary_tag
 Function: get/set on the primary_tag of feature1
 Returns : 0,1,2
 Args    : none if get, the new value if set

=head2 source_tag

 Title   : source_tag
 Usage   : $tag = $feat->source_tag()
           $feat->source_tag('genscan');
 Function: Returns the source tag for a feature,
           eg, 'genscan'
 Returns : a string
 Args    : none

=head2 seqname

 Title   : seqname
 Usage   : $obj->seqname($newval)
 Function: There are many cases when you make a feature that you
           do know the sequence name, but do not know its actual
           sequence. This is an attribute such that you can store
           the seqname.

           This attribute should *not* be used in GFF dumping, as
           that should come from the collection in which the seq
           feature was found.
 Returns : value of seqname
 Args    : newvalue (optional)


1;
