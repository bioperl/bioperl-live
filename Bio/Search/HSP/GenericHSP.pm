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

Bio::Search::HSP::GenericHSP - DESCRIPTION of Object

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


=head1 DESCRIPTION

Describe the object here

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
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::SeqFeature::Similarity;
use Bio::Search::HSP::HSPI;

@ISA = qw(Bio::Root::Root Bio::Search::HSP::HSPI);

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
	$qframe,$hframe) = $self->_rearrange([qw(ALGORITHM
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
						 HIT_FRAME)], @args);

    $algo = 'GENERIC' unless defined $algo;
    $self->algorithm($algo);

    defined $evalue    && $self->evalue($evalue);
    
    defined $pvalue    && $self->pvalue($pvalue);
    defined $bits      && $self->bits($bits);
    defined $score     && $self->score($score);

    my ($queryfactor, $hitfactor) = (0,0);
    if( $algo eq 'TFASTN' || $algo eq 'TFASTY' || $algo eq 'TFASTXY' ||
	$algo eq 'TBLASTN' ) {
	$hitfactor = 1;	
    } elsif ($algo eq 'FASTX' || $algo eq 'FASTY' || $algo eq 'FASTXY'  ) {
	$queryfactor = 1;	
    } elsif ($algo eq 'TBLASTX' ||$algo eq 'TFASTX' ||
	     $algo eq 'TFASTXY' || $algo eq 'TFASTY' || 
	     $algo eq 'BLASTN' || 
	     $algo eq 'FASTN' )  {
	$hitfactor = 1;
	$queryfactor = 1;
    } elsif( $algo eq 'RPSBLAST' ) {
	$queryfactor = $hitfactor = 0;
	$qframe = $hframe = 0;
    }

    # Store the aligned query as sequence feature
    my $strand;
    if( ! $qe || ! $qs ) { $self->throw("Did not specify a Query End or Query Begin"); }
    if ($qe > $qs) {  # normal query: start < end
	if ($queryfactor) { $strand = 1; } else { $strand = undef; }
	$self->query( Bio::SeqFeature::Similarity->new
		      ('-start' => $qs,
		       '-end'   => $qe,
		       '-strand'=> $strand,
		       '-source'=> $algo,
		      ) ) }
    else { # reverse query (i dont know if this is possible, 
	   # but feel free to correct)
	if ($queryfactor) { $strand = -1; } else { $strand = undef; }
	$self->query( Bio::SeqFeature::Similarity->new
		      ('-start' => $qe,
		       '-end'   => $qs,
		       '-strand'=> $strand,
		       '-source'=> $algo,
		      ) );
    }
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
	$self->hit( Bio::SeqFeature::Similarity->new
		    ('-start' => $hs,
		     '-end'   => $he,
		     '-strand'=> $strand,
		     '-source'=> $algo) ) }
    else {			# reverse subject: start bigger than end
	if ($hitfactor) { $strand = -1; } else { $strand = undef; }
	$self->hit( Bio::SeqFeature::Similarity->new
		    ('-start' => $he,
		     '-end'   => $hs,
		     '-strand'=> $strand,
		     '-source'=> $algo) );
    }

    if( defined $strand && ! defined $hframe && $hitfactor ) {
	$hframe = ( $self->hit->start % 3 ) * $strand;
    } elsif( ! defined $strand ) { 
	$hframe = 0;
    }

    $self->frame($qframe,$hframe);

    if( ! defined $query_len || ! defined $hit_len ) { 
	$self->throw("Must defined hit and query length");
    }

    $self->hit->seqname($hit_name);
    $self->query->seqname($query_name);

    $self->hit->seqlength($hit_len);
    $self->query->seqlength($query_len);

    if( ! defined $identical ) { 
	$self->warn("Did not defined the number of identical matches in the HSP assuming 0");
	$identical = 0;
    } 
    if( ! defined $conserved ) {
	$self->warn("Did not defined the number of conserved matches in the HSP assuming conserved == identical ($identical)") if( $algo !~ /T?(FAST|BLAST)N/i);
	$conserved = $identical;
    } 
    # protect for divide by zero if user does not specify 
    # hsp_len, query_len, or hit_len

    if( $hsp_len ) {
	$self->length('total', $hsp_len);
	$self->frac_identical( 'total', $identical / $self->length('total'));
	$self->frac_conserved( 'total', $conserved / $self->length('total'));
    }
    if( $hit_len ) {
	$self->length('hit', $hit_len);
	$self->frac_identical( 'hit', $identical / $self->length('hit'));
	$self->frac_conserved( 'hit', $conserved / $self->length('hit'));
    }
    if( $query_len ) {
	$self->length('query', $query_len);	
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
	$self->warn("Did not defined the number of gaps in the HSP calculating");
	$gaps = $self->gaps("query") + $self->gaps("hit");
    } 
    $self->gaps('total', $gaps);

    $self->percent_identity($identical / $hsp_len ) if( $hsp_len > 0 );


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

sub evalue {
    my ($self,$value) = @_;
    my $previous = $self->{'_evalue'};
    if( defined $value  ) { 	
	$self->{'_evalue'} = $value;
    } 
    return $previous;   
}

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

sub length{
    my ($self, $type,$value) = @_;
    $type = lc $type if defined $type;
    $type = 'total' if( ! defined $type ||
			$type !~ /query|hit|total/);
    my $previous = $self->{'_length'}->{$type};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = 0 unless defined $value;
	$self->{'_length'}->{$type} = $value;
    } 
    return $previous;
}

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
	      if( (defined $1 && $1 eq '-' && $self->hit->strand >= 0) ||
		  (defined $1 && $1 eq '+' && $self->hit->strand <= 0) )
	      {
		  $self->warn("Subject frame ($sframe) did not match strand of subject (". $self->hit->strand() . ")");
	      }

	      # Set frame to GFF [0-2]
	      $sframe = $2 - 1;
	  }
	  $self->hit->frame($sframe);
      }

    if (wantarray() &&
	$self->report_type eq 'TBLASTX')
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
    $hs =~ s/[\/\\]/\-/g;
    $qs =~ s/[\/\\]/\-/g;
    my $seqonly = $qs;
    $seqonly =~ s/\-//g;
    
    my $query = new Bio::LocatableSeq('-seq'   => $qs,
				      '-id'    => $self->query->seqname(),
				      '-start' => 1,
				      '-end' => CORE::length($seqonly),
				      );
    $seqonly = $hs;
    $seqonly =~ s/\-//g;
    
    my $hit =  new Bio::LocatableSeq('-seq'   => $hs,
				      '-id'    => $self->hit->seqname(),
				      '-start' => 1,
				      '-end' => CORE::length($seqonly),
				      );
    $aln->add_seq($query);
    $aln->add_seq($hit);
    return $aln;
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

1;
