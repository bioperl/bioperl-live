# $Id$
#
# BioPerl module for Bio::Search::HSP::GenericHSP
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::HSP::GenericHSP - A "Generic" implementation of a High Scoring Pair 

=head1 SYNOPSIS

    my $hsp = new Bio::Search::HSP::GenericHSP( -algorithm => 'blastp',
						-evalue    => '1e-30',
						);

    $r_type = $hsp->algorithm

    $pvalue = $hsp->p();

    $evalue = $hsp->evalue();

    $frac_id = $hsp->frac_identical( ['query'|'hit'|'total'] );

    $frac_cons = $hsp->frac_conserved( ['query'|'hit'|'total'] );

    $gaps = $hsp->gaps( ['query'|'hit'|'total'] );

    $qseq = $hsp->query_string;

    $hseq = $hsp->hit_string;

    $homo_string = $hsp->homology_string;

    $len = $hsp->length( ['query'|'hit'|'total'] );

    $len = $hsp->length( ['query'|'hit'|'total'] );
    
    $rank = $hsp->rank;


=head1 DESCRIPTION

This implementation is "Generic", meaning it is is suitable for
holding information about High Scoring pairs from most Search reports
such as BLAST and FastA.  Specialized objects can be derived from
this.

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
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich and Steve Chervitz

Email jason@bioperl.org
Email sac@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::HSP::GenericHSP;
use vars qw(@ISA $GAP_SYMBOL);
use strict;

use Bio::Root::Root;
use Bio::SeqFeature::Similarity;
use Bio::Search::HSP::HSPI;

@ISA = qw(Bio::Search::HSP::HSPI Bio::Root::Root );

BEGIN {
    $GAP_SYMBOL = '-';
}
=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Search::HSP::GenericHSP();
 Function: Builds a new Bio::Search::HSP::GenericHSP object 
 Returns : Bio::Search::HSP::GenericHSP
 Args    : -algorithm => algorithm used (BLASTP, TBLASTX, FASTX, etc)
           -evalue    => evalue
           -pvalue    => pvalue
           -bits      => bit value for HSP
           -score     => score value for HSP (typically z-score but depends on
					      analysis)
           -hsp_length=> Length of the HSP (including gaps)
           -identical => # of residues that that matched identically
           -conserved => # of residues that matched conservatively 
                           (only protein comparisions; 
			    conserved == identical in nucleotide comparisons)
           -hsp_gaps   => # of gaps in the HSP
           -query_gaps => # of gaps in the query in the alignment
           -hit_gaps   => # of gaps in the subject in the alignment    
           -query_name  => HSP Query sequence name (if available)
           -query_start => HSP Query start (in original query sequence coords)
           -query_end   => HSP Query end (in original query sequence coords)
           -hit_name    => HSP Hit sequence name (if available)
           -hit_start   => HSP Hit start (in original hit sequence coords)
           -hit_end     => HSP Hit end (in original hit sequence coords)
           -hit_length  => total length of the hit sequence
           -query_length=> total length of the query sequence
           -query_seq   => query sequence portion of the HSP
           -hit_seq     => hit sequence portion of the HSP
           -homology_seq=> homology sequence for the HSP
           -hit_frame   => hit frame (only if hit is translated protein)
           -query_frame => query frame (only if query is translated protein)
           -rank        => HSP rank

=cut

sub new {
    my($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
    my ($algo, $evalue, $pvalue, $identical, $conserved, 
	$gaps, $query_gaps, $hit_gaps,
	$hit_seq, $query_seq, $homology_seq,
	$hsp_len, $query_len,$hit_len,
	$hit_name,$query_name,$bits,$score,
	$hs,$he,$qs,$qe,
	$qframe,$hframe,
	$rank) = $self->_rearrange([qw(ALGORITHM
				       EVALUE
				       PVALUE
				       IDENTICAL
				       CONSERVED
				       HSP_GAPS
				       QUERY_GAPS
				       HIT_GAPS
				       HIT_SEQ
				       QUERY_SEQ
				       HOMOLOGY_SEQ
				       HSP_LENGTH
				       QUERY_LENGTH
				       HIT_LENGTH
				       HIT_NAME
				       QUERY_NAME
				       BITS
				       SCORE
				       HIT_START
				       HIT_END
				       QUERY_START
				       QUERY_END
				       QUERY_FRAME
				       HIT_FRAME
				       RANK )], @args);

    $algo = 'GENERIC' unless defined $algo;
    $self->algorithm($algo);

#    defined $evalue    && $self->evalue($evalue)
#    $hsp->significance is initialized by the 
#    the SimilarityPair object - let's only keep one
#    value, don't need 2 slots.

    defined $pvalue    && $self->pvalue($pvalue);
    defined $bits      && $self->bits($bits);
    defined $score     && $self->score($score);

    my ($queryfactor, $hitfactor) = (0,0);
    if( $algo eq 'TFASTN' || $algo eq 'TFASTY' || $algo eq 'TFASTXY' ||
	$algo eq 'TBLASTN' ) {
	$hitfactor = 1;	
    } elsif ($algo eq 'BLASTX' || 
	     $algo eq 'FASTX' || $algo eq 'FASTY' || $algo eq 'FASTXY'  ) {
	$queryfactor = 1;	
    } elsif ($algo eq 'TBLASTX' ||$algo eq 'TFASTX' ||
	     $algo eq 'TFASTXY' || $algo eq 'TFASTY' || 
	     $algo eq 'BLASTN' || 
	     $algo eq 'FASTN' || $algo eq 'WABA')  {
	$hitfactor = 1;
	$queryfactor = 1;
    } elsif( $algo eq 'RPSBLAST' ) {
	$queryfactor = $hitfactor = 0;
	$qframe = $hframe = 0;
    }
    # Store the aligned query as sequence feature
    my $strand;
    unless(  $qe && $qs ) { $self->throw("Did not specify a Query End or Query Begin @args"); }
    unless( $he && $hs ) { $self->throw("Did not specify a Hit End or Hit Begin"); }
    if ($qe > $qs) {  # normal query: start < end
	if ($queryfactor) { $strand = 1; } else { $strand = undef; }	
    } else { # reverse query (i dont know if this is possible, 
	     # but feel free to correct)
	if ($queryfactor) { $strand = -1; } else { $strand = undef; }
	($qs,$qe) = ($qe,$qs);
    
    }
    $self->query( new  Bio::SeqFeature::Similarity
		  ('-primary'  => $self->primary_tag,
		   '-start'    => $qs,
		   '-expect'   => $evalue,
		   '-bits'     => $bits,
		   '-end'      => $qe,
		   '-strand'   => $strand,
		   '-seq_id'   => $query_name,
		   '-seqlength'=> $query_len,
		   '-source'   => $algo,
		   ) );
    
    # to determine frame from something like FASTXY which doesn't
    # report the frame
    if( defined $strand && ! defined $qframe && $queryfactor ) {
	$qframe = ( $self->query->start % 3 ) * $strand;
    } elsif( ! defined $strand ) { 
	$qframe = 0;
    }
    # store the aligned subject as sequence feature
    if ($he > $hs) {		# normal subject
	if ($hitfactor) { $strand = 1; } else { $strand = undef; }
    } else {
	if ($hitfactor) { $strand = -1; } else { $strand = undef; }
	($hs,$he) = ( $he,$hs); # reverse subject: start bigger than end
    }

    $self->hit( Bio::SeqFeature::Similarity->new
		('-start'     => $hs,
		 '-end'       => $he,
		 '-strand'    => $strand,
		 '-expect'    => $evalue,
		 '-bits'      => $bits,
		 '-source'    => $algo,
		 '-seq_id'    => $hit_name,
		 '-seqlength' => $hit_len,
		 '-primary'   => $self->primary_tag ));
    
    if( defined $strand && ! defined $hframe && $hitfactor ) {
	$hframe = ( $hs % 3 ) * $strand;
    } elsif( ! defined $strand ) { 
	$hframe = 0;
    }

    $self->frame($qframe,$hframe);

    if( ! defined $query_len || ! defined $hit_len ) { 
	$self->throw("Must defined hit and query length");
    }

    if( ! defined $identical ) { 
	$self->warn("Did not defined the number of identical matches in the HSP assuming 0");
	$identical = 0;
    } 
    if( ! defined $conserved ) {
	$self->warn("Did not defined the number of conserved matches in the HSP assuming conserved == identical ($identical)") if( $algo !~ /(FAST|BLAST)N/i);
	$conserved = $identical;
    } 
    # protect for divide by zero if user does not specify 
    # hsp_len, query_len, or hit_len
    
    $self->num_identical($identical);
    $self->num_conserved($conserved);
    
    if( $hsp_len ) {
	$self->length('total', $hsp_len);
	$self->frac_identical( 'total', $identical / $self->length('total'));
	$self->frac_conserved( 'total', $conserved / $self->length('total'));
    }
    if( $hit_len ) {
#	$self->length('hit', $self->hit->length);
	$self->frac_identical( 'hit', $identical / $self->length('hit'));
	$self->frac_conserved( 'hit', $conserved / $self->length('hit'));
    }
    if( $query_len ) {
#	$self->length('query', $self->query->length);	
	$self->frac_identical( 'query', $identical / $self->length('query')) ;
	$self->frac_conserved( 'query', $conserved / $self->length('query'));
    }
    $self->query_string($query_seq);
    $self->hit_string($hit_seq);
    $self->homology_string($homology_seq);
    
    if( defined $query_gaps ) {
	$self->gaps('query', $query_gaps);
    } else {
	$self->gaps('query', scalar ( $query_seq =~ tr/\-//));
    } 
    if( defined $hit_gaps ) {
	$self->gaps('hit', $hit_gaps);
    } else {
	$self->gaps('hit', scalar ( $hit_seq =~ tr/\-//));
    }
    if(! defined $gaps ) {
	$gaps = $self->gaps("query") + $self->gaps("hit");
    } 
    $self->gaps('total', $gaps);

    $self->percent_identity($identical / $hsp_len ) if( $hsp_len > 0 );

    $rank && $self->rank($rank);
    return $self;
}



=head2 Bio::Search::HSP::HSPI methods

Implementation of Bio::Search::HSP::HSPI methods follow

=head2 algorithm

 Title   : algorithm
 Usage   : my $r_type = $hsp->algorithm
 Function: Obtain the name of the algorithm used to obtain the HSP
 Returns : string (e.g., BLASTP)
 Args    : [optional] scalar string to set value

=cut

sub algorithm{
    my ($self,$value) = @_;
    my $previous = $self->{'_algorithm'};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = '' unless defined $value;
	$self->{'_algorithm'} = $value;
    } 

    return $previous;   
}

=head2 pvalue

 Title   : pvalue
 Usage   : my $pvalue = $hsp->pvalue();
 Function: Returns the P-value for this HSP or undef 
 Returns : float or exponential (2e-10)
           P-value is not defined with NCBI Blast2 reports.
 Args    : [optional] numeric to set value

=cut

sub pvalue {
    my ($self,$value) = @_;
    my $previous = $self->{'_pvalue'};
    if( defined $value  ) { 	
	$self->{'_pvalue'} = $value;
    } 
    return $previous;   
}

=head2 evalue

 Title   : evalue
 Usage   : my $evalue = $hsp->evalue();
 Function: Returns the e-value for this HSP
 Returns : float or exponential (2e-10)
 Args    : [optional] numeric to set value

=cut

sub evalue { shift->significance(@_) }

=head2 frac_identical

 Title   : frac_identical
 Usage   : my $frac_id = $hsp->frac_identical( ['query'|'hit'|'total'] );
 Function: Returns the fraction of identitical positions for this HSP 
 Returns : Float in range 0.0 -> 1.0
 Args    : arg 1:  'query' = num identical / length of query seq (without gaps)
                   'hit'   = num identical / length of hit seq (without gaps)
                   'total' = num identical / length of alignment (with gaps)
                   default = 'total' 
           arg 2: [optional] frac identical value to set for the type requested

=cut

sub frac_identical {
   my ($self, $type,$value) = @_;

   $type = lc $type if defined $type;
   $type = 'total' if( ! defined $type ||
		       $type !~ /query|hit|total/); 
   my $previous = $self->{'_frac_identical'}->{$type};
   if( defined $value || ! defined $previous ) { 
       $value = $previous = '' unless defined $value;
       if( $type eq 'hit' || $type eq 'query' ) {
	   $self->$type()->frac_identical( $value);
       }
       $self->{'_frac_identical'}->{$type} = $value;
   } 
   return $previous;   

}

=head2 frac_conserved

 Title    : frac_conserved
 Usage    : my $frac_cons = $hsp->frac_conserved( ['query'|'hit'|'total'] );
 Function : Returns the fraction of conserved positions for this HSP.
            This is the fraction of symbols in the alignment with a 
            positive score.
 Returns : Float in range 0.0 -> 1.0
 Args    : arg 1: 'query' = num conserved / length of query seq (without gaps)
                  'hit'   = num conserved / length of hit seq (without gaps)
                  'total' = num conserved / length of alignment (with gaps)
                  default = 'total' 
           arg 2: [optional] frac conserved value to set for the type requested

=cut

sub frac_conserved {
    my ($self, $type,$value) = @_;
    $type = lc $type if defined $type;
    $type = 'total' if( ! defined $type ||
			$type !~ /query|hit|total/);
    my $previous = $self->{'_frac_conserved'}->{$type};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = '' unless defined $value;
	$self->{'_frac_conserved'}->{$type} = $value;
    } 
    return $previous;   
}

=head2 gaps

 Title    : gaps
 Usage    : my $gaps = $hsp->gaps( ['query'|'hit'|'total'] );
 Function : Get the number of gaps in the query, hit, or total alignment.
 Returns  : Integer, number of gaps or 0 if none
 Args     : arg 1: 'query' = num gaps in query seq
                   'hit'   = num gaps in hit seq
                   'total' = num gaps in whole alignment 
                   default = 'total' 
            arg 2: [optional] integer gap value to set for the type requested

=cut

sub gaps        {
    my ($self, $type,$value) = @_;
    $type = lc $type if defined $type;
    $type = 'total' if( ! defined $type ||
			$type !~ /query|hit|total/);
    my $previous = $self->{'_gaps'}->{$type};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = '' unless defined $value;
	$self->{'_gaps'}->{$type} = $value;
    }
    return $previous;
}

=head2 query_string

 Title   : query_string
 Usage   : my $qseq = $hsp->query_string;
 Function: Retrieves the query sequence of this HSP as a string
 Returns : string
 Args    : [optional] string to set for query sequence


=cut

sub query_string{
    my ($self,$value) = @_;
    my $previous = $self->{'_query_string'};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = '' unless defined $value;
	$self->{'_query_string'} = $value;
	# do some housekeeping so we know when to 
	# re-run _calculate_seq_positions
	$self->{'_sequenceschanged'} = 1;
    } 
    return $previous;   
}

=head2 hit_string

 Title   : hit_string
 Usage   : my $hseq = $hsp->hit_string;
 Function: Retrieves the hit sequence of this HSP as a string
 Returns : string
 Args    : [optional] string to set for hit sequence


=cut

sub hit_string{
    my ($self,$value) = @_;
    my $previous = $self->{'_hit_string'};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = '' unless defined $value;
	$self->{'_hit_string'} = $value;
	# do some housekeeping so we know when to 
	# re-run _calculate_seq_positions
	$self->{'_sequenceschanged'} = 1;
    } 
    return $previous;   
}

=head2 homology_string

 Title   : homology_string
 Usage   : my $homo_string = $hsp->homology_string;
 Function: Retrieves the homology sequence for this HSP as a string.
         : The homology sequence is the string of symbols in between the 
         : query and hit sequences in the alignment indicating the degree
         : of conservation (e.g., identical, similar, not similar).
 Returns : string
 Args    : [optional] string to set for homology sequence

=cut

sub homology_string{
    my ($self,$value) = @_;
    my $previous = $self->{'_homology_string'};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = '' unless defined $value;
	$self->{'_homology_string'} = $value;
	# do some housekeeping so we know when to 
	# re-run _calculate_seq_positions
	$self->{'_sequenceschanged'} = 1;
    } 
    return $previous;   
}

=head2 length

 Title    : length
 Usage    : my $len = $hsp->length( ['query'|'hit'|'total'] );
 Function : Returns the length of the query or hit in the alignment 
            (without gaps) 
            or the aggregate length of the HSP (including gaps;
            this may be greater than either hit or query )
 Returns  : integer
 Args     : arg 1: 'query' = length of query seq (without gaps)
                   'hit'   = length of hit seq (without gaps)
                   'total' = length of alignment (with gaps)
                   default = 'total' 
            arg 2: [optional] integer length value to set for specific type

=cut

sub length {

    my $self = shift;
    my $type = shift;

    $type = 'total' unless defined $type;
    $type = lc $type;

    if( $type =~ /^q/i ) {
	return $self->query()->length(shift);
    } elsif( $type =~ /^(hit|subject|sbjct)/ ) {
	return $self->hit()->length(shift);
    } else { 
	my $v = shift;
	if( defined $v ) { 
	    $self->{'_hsplength'} = $v;
	}
	return $self->{'_hsplength'};
   }
    return 0; # should never get here
}

=head2 hsp_length

 Title   : hsp_length
 Usage   : my $len = $hsp->hsp_length()
 Function: shortcut  length('hsp')
 Returns : floating point between 0 and 100 
 Args    : none

=cut

sub hsp_length { return shift->length('hsp', shift); }

=head2 percent_identity

 Title   : percent_identity
 Usage   : my $percentid = $hsp->percent_identity()
 Function: Returns the calculated percent identity for an HSP
 Returns : floating point between 0 and 100 
 Args    : none


=cut


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
 Note    : Frames are stored in the GFF way (0-2) not 1-3
           as they are in BLAST (negative frames are deduced by checking 
				 the strand of the query or hit)

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
	      my $dir = $1;
	      $dir = '+' unless defined $dir;
	      if( ($dir eq '-' && $self->hit->strand >= 0) ||
		  ($dir eq '+' && $self->hit->strand <= 0) )
	      {
		  $self->warn("Subject frame ($sframe) did not match strand of subject (". $self->hit->strand() . ")");
	      }

	      # Set frame to GFF [0-2]
	      $sframe = $2 - 1;
	  }
	  $self->hit->frame($sframe);
      }
    if (wantarray() &&
	$self->algorithm eq 'TBLASTX')
    {
	return ($self->query->frame(), $self->hit->frame());
    } elsif (wantarray())  {
	($self->query->frame() &&
	 return ($self->query->frame(), undef)) ||
	     ($self->hit->frame() &&
	      return (undef, $self->hit->frame()));
    } else {
	($self->query->frame() &&
	 return $self->query->frame()) ||
	($self->hit->frame() &&
	 return $self->hit->frame());
    }
}


=head2 get_aln

 Title   : get_aln
 Usage   : my $aln = $hsp->gel_aln
 Function: Returns a Bio::SimpleAlign representing the HSP alignment
 Returns : Bio::SimpleAlign
 Args    : none

=cut

sub get_aln {
    my ($self) = @_;
    require Bio::LocatableSeq;
    require Bio::SimpleAlign;
    my $aln = new Bio::SimpleAlign;
    my $hs = $self->hit_string();
    my $qs = $self->query_string();
    if( $self->algorithm  =~ /FAST/i ) {
	# fasta reports some extra 'regional' sequence information
	# we need to clear out first
	# this seemed a bit insane to me at first, but it appears to 
	# work --jason
	
	# we infer the end of the regional sequence where the first
	# non space is in the homology string
	# then we use the HSP->length to tell us how far to read
	# to cut off the end of the sequence

	# one possible problem is the sequence which 
	
	my ($start) = 0;
	if( $self->homology_string() =~ /^(\s+)/ ) {
	    $start = CORE::length($1);
	}
	$hs = substr($hs, $start,$self->length('total'));
	$qs = substr($qs, $start,$self->length('total'));
	foreach my $seq ( $qs,$hs)  {
	    foreach my $f ( '\\', '/', ' ') {
		my $index =  index($seq,$f);
		while( $index >=0 ) {
		    substr($hs,$index,1) = '';
		    substr($qs,$index,1) = '';
		    $index = index($seq,$f,$index+1);
		}
	    }
	}
    }

    my $seqonly = $qs;
    $seqonly =~ s/[\-\s]//g;
    my ($q_nm,$s_nm) = ($self->query->seq_id(),
			$self->hit->seq_id());
    unless( defined $q_nm && CORE::length ($q_nm) ) {
	$q_nm = 'query';
    }
    unless( defined $s_nm && CORE::length ($s_nm) ) {
	$s_nm = 'hit';
    }
    my $query = new Bio::LocatableSeq('-seq'   => $qs,
				      '-id'    => $q_nm,
				      '-start' => 1,
				      '-end' => CORE::length($seqonly),
				      );
    $seqonly = $hs;
    $seqonly =~ s/[\-\s]//g;
    my $hit =  new Bio::LocatableSeq('-seq'   => $hs,
				      '-id'    => $s_nm,
				      '-start' => 1,
				      '-end' => CORE::length($seqonly),
				      );
    $aln->add_seq($query);
    $aln->add_seq($hit);
    return $aln;
}

=head2 num_conserved

 Title   : num_conserved
 Usage   : $obj->num_conserved($newval)
 Function: returns the number of conserved residues in the alignment
 Returns : inetger
 Args    : integer (optional)


=cut

sub num_conserved{
   my ($self,$value) = @_;
   if( defined $value) {
       $self->{'num_conserved'} = $value;
   }
   return $self->{'num_conserved'};
}

=head2 num_identical

 Title   : num_identical
 Usage   : $obj->num_identical($newval)
 Function: returns the number of identical residues in the alignment
 Returns : integer
 Args    : integer (optional)


=cut

sub num_identical{
   my ($self,$value) = @_;
   if( defined $value) {
       $self->{'_num_identical'} = $value;
   }
   return $self->{'_num_identical'};
}

=head2 rank

 Usage     : $hsp->rank( [string] );
 Purpose   : Get the rank of the HSP within a given Blast hit.
 Example   : $rank = $hsp->rank;
 Returns   : Integer (1..n) corresponding to the order in which the HSP
             appears in the BLAST report.

=cut

sub rank { 
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'_rank'} = $value;
    }
    return $self->{'_rank'};
}


=head2 seq_inds

 Title   : seq_inds
 Purpose   : Get a list of residue positions (indices) for all identical 
           : or conserved residues in the query or sbjct sequence.
 Example   : @s_ind = $hsp->seq_inds('query', 'identical');
           : @h_ind = $hsp->seq_inds('hit', 'conserved');
           : @h_ind = $hsp->seq_inds('hit', 'conserved', 1);
 Returns   : List of integers 
           : May include ranges if collapse is true.
 Argument  : seq_type  = 'query' or 'hit' or 'sbjct'  (default = query)
           :  ('sbjct' is synonymous with 'hit') 
           : class     = 'identical' or 'conserved' or 'nomatch' or 'gap'
           :              (default = identical)
           :              (can be shortened to 'id' or 'cons')
           :              
           : collapse  = boolean, if true, consecutive positions are merged
           :             using a range notation, e.g., "1 2 3 4 5 7 9 10 11" 
           :             collapses to "1-5 7 9-11". This is useful for 
           :             consolidating long lists. Default = no collapse.
 Throws    : n/a.
 Comments  : 

See Also   : L<Bio::Search::BlastUtils::collapse_nums()|Bio::Search::BlastUtils>, L<Bio::Search::Hit::HitI::seq_inds()|Bio::Search::Hit::HitI>

=cut

sub seq_inds{
   my ($self, $seqType, $class, $collapse) = @_;

   # prepare the internal structures - this is cached so
   # if the strings have not changed we're okay
   $self->_calculate_seq_positions();

   $seqType  ||= 'query';
   $class ||= 'identical';
   $collapse ||= 0;
   $seqType = 'sbjct' if $seqType eq 'hit';
   my $t = lc(substr($seqType,0,1));
   if( $t eq 'q' ) {
       $seqType = 'query';
   } elsif ( $t eq 's' || $t eq 'h' ) {
       $seqType = 'sbjct';
   } else { 
       $self->warn("unknown seqtype $seqType using 'query'");
       $seqType = 'query';
   }

   $t = lc(substr($class,0,1));
   if( $t eq 'c' ) {
     $class = 'conserved';  
   } elsif( $t eq 'i' ) {
       $class = 'identical';
   } elsif( $t eq 'n' ) {
       $class = 'nomatch';
   } elsif( $t eq 'g' ) {
       $class = 'gap';
   } else { 
       $self->warn("unknown sequence class $class using 'identical'");
       $class = 'identical';
   }
   
   ## Sensitive to member name changes.
   $seqType  = "_\L$seqType\E";
   $class = "_\L$class\E";
   my @ary;
   if( $class eq '_gap' ) {
       # this means that we are remapping the gap length that is stored
       # in the hash (for example $self->{'_gapRes_query'} ) 
       # so we'll return an array which has the values of the position of the 
       # of the gap (the key in the hash) + the gap length (value in the
       # hash for this key - 1.

       @ary = map { $_ > 1 ?
			$_..($_ + $self->{"${class}Res$seqType"}->{$_} - 1) : 
			$_ }
              sort { $a <=> $b } keys %{ $self->{"${class}Res$seqType"}};
   } else {
       @ary = sort { $a <=> $b } keys %{ $self->{"${class}Res$seqType"}};
   }   
   require Bio::Search::BlastUtils if $collapse;
   
   return $collapse ? &Bio::Search::BlastUtils::collapse_nums(@ary) : @ary;
}


=head2 Inherited from Bio::SeqFeature::SimilarityPair

These methods come from Bio::SeqFeature::SimilarityPair

=head2 query

 Title   : query
 Usage   : my $query = $hsp->query
 Function: Returns a SeqFeature representing the query in the HSP
 Returns : Bio::SeqFeature::Similarity
 Args    : [optional] new value to set


=head2 hit

 Title   : hit
 Usage   : my $hit = $hsp->hit
 Function: Returns a SeqFeature representing the hit in the HSP
 Returns : Bio::SeqFeature::Similarity
 Args    : [optional] new value to set


=head2 significance

 Title   : significance
 Usage   : $evalue = $obj->significance();
           $obj->significance($evalue);
 Function: Get/Set the significance value
 Returns : numeric
 Args    : [optional] new value to set


=head2 score

 Title   : score
 Usage   : my $score = $hsp->score();
 Function: Returns the score for this HSP or undef 
 Returns : numeric           
 Args    : [optional] numeric to set value

=cut 

# overriding 

sub score {
    my ($self,$value) = @_;
    my $previous = $self->{'_score'};
    if( defined $value ) { 
        $self->{'_score'} = $value;
    } 
    return $previous;
}

=head2 bits

 Title   : bits
 Usage   : my $bits = $hsp->bits();
 Function: Returns the bit value for this HSP or undef 
 Returns : numeric
 Args    : none

=cut

# overriding 

sub bits {
    my ($self,$value) = @_;
    my $previous = $self->{'_bits'};
    if( defined $value ) { 
        $self->{'_bits'} = $value;
    } 
    return $previous;
}


=head2 strand

 Title   : strand
 Usage   : $hsp->strand('quer')
 Function: Retrieves the strand for the HSP component requested
 Returns : +1 or -1 (0 if unknown)
 Args    : 'hit' or 'subject' or 'sbjct' to retrieve the strand of the subject
           'query' to retrieve the query strand (default)

=cut

=head1 Private methods

=cut

=head2 _calculate_seq_positions

 Title   : _calculate_seq_positions
 Usage   : $self->_calculate_seq_positions
 Function:
 Returns : 
 Args    :


=cut

sub _calculate_seq_positions {
    my ($self,@args) = @_;
    return unless ( $self->{'_sequenceschanged'} );
    $self->{'_sequenceschanged'} = 0;
    my ($mchar, $schar, $qchar);
    my ($seqString, $qseq,$sseq) = ( $self->homology_string(),
				     $self->query_string(),
				     $self->hit_string() );

    # Using hashes to avoid saving duplicate residue numbers.
    my %identicalList_query = ();
    my %identicalList_sbjct = ();
    my %conservedList_query = ();
    my %conservedList_sbjct = ();
    
    my %gapList_query = ();
    my %gapList_sbjct = ();
    my %nomatchList_query = ();
    my %nomatchList_sbjct = ();

    my $qdir = $self->query->strand || 1;
    my $sdir = $self->hit->strand || 1;
    my $resCount_query = ($qdir >=0) ? $self->query->end : $self->query->start;
    my $resCount_sbjct = ($sdir >=0) ? $self->hit->end : $self->hit->start;
    
    my $prog = $self->algorithm;
    if( $prog  =~ /FAST/i ) {
	# fasta reports some extra 'regional' sequence information
	# we need to clear out first
	# this seemed a bit insane to me at first, but it appears to 
	# work --jason
	
	# we infer the end of the regional sequence where the first
	# non space is in the homology string
	# then we use the HSP->length to tell us how far to read
	# to cut off the end of the sequence

	# one possible problem is the sequence which 
	
	my ($start) = (0);
	if( $seqString =~ /^(\s+)/ ) {
	    $start = CORE::length($1);
	}

	$seqString = substr($seqString, $start,$self->length('total'));
	$qseq = substr($qseq, $start,$self->length('total'));
	$sseq = substr($sseq, $start,$self->length('total'));

	$qseq =~ s![\\\/]!!g;
	$sseq =~ s![\\\/]!!g;
    }
    if($prog eq 'TBLASTN' || $prog eq 'TFASTN' ) {
	$resCount_sbjct /= 3;
    } elsif($prog eq 'BLASTX' || $prog eq 'FASTX' || $prog eq 'FASTY' || 
	    $prog eq 'FASTXY' ) {
	$resCount_query /= 3;
    } elsif($prog eq 'TBLASTX' ||
	    $prog eq 'TFASTXY' || $prog eq 'TFASTY' || 
	    $prog eq 'TFASTX' ) {
	$resCount_query /= 3;
	$resCount_sbjct /= 3;
    }    
    while( $mchar = chop($seqString) ) {
	($qchar, $schar) = (chop($qseq), chop($sseq));
	if( $mchar eq '+' || $mchar eq '.' || $mchar eq ':' ) { 
	    $conservedList_query{ $resCount_query } = 1; 
	    $conservedList_sbjct{ $resCount_sbjct } = 1; 
	} elsif( $mchar ne ' ' ) { 
	    $identicalList_query{ $resCount_query } = 1; 
	    $identicalList_sbjct{ $resCount_sbjct } = 1;
	} elsif( $mchar eq ' ') { 
	    $nomatchList_query{ $resCount_query } = 1;
	    $nomatchList_sbjct{ $resCount_sbjct } = 1;
	}
	if( $qchar eq $GAP_SYMBOL ) {
	    $gapList_query{ $resCount_query } ++;
	} else { 	    
	    $resCount_query -= $qdir;
	}
	if( $schar eq $GAP_SYMBOL ) {
	    $gapList_sbjct{ $resCount_query } ++;
	} else {
	    $resCount_sbjct -=$sdir;
	}
    }
    $self->{'_identicalRes_query'} = \%identicalList_query;
    $self->{'_conservedRes_query'} = \%conservedList_query;
    $self->{'_nomatchRes_query'} = \%nomatchList_query;
    $self->{'_gapRes_query'} = \%gapList_query;

    $self->{'_identicalRes_sbjct'} = \%identicalList_sbjct;
    $self->{'_conservedRes_sbjct'} = \%conservedList_sbjct;
    $self->{'_nomatchRes_sbjct'} = \%nomatchList_sbjct;
    $self->{'_gapRes_sbjct'} = \%gapList_sbjct;
    return 1;
}

=head2 n

See documentation in L<Bio::Search::HSP::HSPI::n()|Bio::Search::HSP::HSPI>

=cut

#-----
sub n { 
    my $self = shift; 
    if(@_) { $self->{'_n'} = shift; }
    defined $self->{'_n'} ? $self->{'_n'} : '';
}

=head2 range

See documentation in L<Bio::Search::HSP::HSPI::range()|Bio::Search::HSP::HSPI>

=cut

#----------
sub range {
#----------
    my ($self, $seqType) = @_;

    $seqType ||= 'query';
    $seqType = 'sbjct' if $seqType eq 'hit';

    my ($start, $end);
    if( $seqType eq 'query' ) {
        $start = $self->query->start;
        $end = $self->query->end;
    }
    else {
        $start = $self->hit->start;
        $end = $self->hit->end;
    }
    return ($start, $end);
}


1;
