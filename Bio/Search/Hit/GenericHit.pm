# $Id$
#
# BioPerl module for Bio::Search::Hit::GenericHit
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Hit::GenericHit - A generic implementation of the Bio::Search::Hit::HitI interface

=head1 SYNOPSIS

    use Bio::Search::Hit::GenericHit;
    my $hit = new Bio::Search::Hit::GenericHit(-algorithm => 'blastp');

    # more likely
    use Bio::SeachIO;
    my $parser = new Bio::SearchIO(-format => 'blast', -file => 'result.bls');

    my $result = $parser->next_result;
    my $hit    = $result->next_hit;


=head1 DESCRIPTION

This object handles the hit data from a Database Sequence Search such
as FASTA or BLAST.

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

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::Hit::GenericHit;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::Search::Hit::HitI;
require Bio::Search::BlastUtils;

@ISA = qw(Bio::Root::Root Bio::Search::Hit::HitI );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Search::Hit::GenericHit();
 Function: Builds a new Bio::Search::Hit::GenericHit object 
 Returns : Bio::Search::Hit::GenericHit
 Args    : -name         => Name of Hit (required)
           -description  => Description (optional)
           -accession    => Accession number (optional)
           -length       => Length of the Hit (optional)
           -score        => Raw Score for the Hit (optional)
           -significance => Significance value for the Hit (optional)
           -algorithm    => Algorithm used (BLASTP, FASTX, etc...)
           -hsps         => Array ref of HSPs for this Hit. 
           -iteration    => integer for the PSI-Blast iteration number

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($hsps, $name,$query_len,$desc, $acc, $length,
      $score,$algo,$signif,$bits,
      $iter,$rank) = $self->_rearrange([qw(HSPS 
				     NAME 
				     QUERY_LEN
				     DESCRIPTION
				     ACCESSION
				     LENGTH SCORE ALGORITHM 
				     SIGNIFICANCE BITS ITERATION
				     RANK )], @args);
  
  $self->{'_query_length'} = $query_len;

  if( ! defined $name ) { 
      $self->throw("Must have defined a valid name for Hit");
  } else { 
      $self->name($name);
  }  

  defined $acc    && $self->accession($acc);
  defined $desc   && $self->description($desc);
  defined $length && $self->length($length);
  defined $algo   && $self->algorithm($algo);
  defined $signif && $self->significance($signif);
  defined $score  && $self->raw_score($score);
  defined $bits   && $self->bits($bits);
  defined $iter   && $self->iteration($iter);
  defined $rank   && $self->rank($rank);

  $self->{'_iterator'} = 0;
  $self->{'_hsps'} = [];
  if( defined $hsps  ) {
      if( ref($hsps) !~ /array/i ) {
	  $self->warn("Did not specify a valid array ref for the param HSPS ($hsps)");
      } else {
	  while( @$hsps ) { 
	      $self->add_hsp(shift @$hsps );
	  }
      }
  } 
  return $self;
}

=head2 add_hsp

 Title   : add_hsp
 Usage   : $hit->add_hsp($hsp)
 Function: Add a HSP to the collection of HSPs for a Hit
 Returns : number of HSPs in the Hit
 Args    : Bio::Search::HSP::HSPI object


=cut

sub add_hsp {
   my ($self,$hsp) = @_;
   if( !defined $hsp || ! $hsp->isa('Bio::Search::HSP::HSPI') ) { 
       $self->warn("Must provide a valid Bio::Search::HSP::HSPI object to object: $self method: add_hsp value: $hsp");
       return undef;
   }
   push @{$self->{'_hsps'}}, $hsp;
   return scalar @{$self->{'_hsps'}};
}



=head2 Bio::Search::Hit::HitI methods

Implementation of Bio::Search::Hit::HitI methods

=head2 name

 Title   : name
 Usage   : $hit_name = $hit->name();
 Function: returns the name of the Hit sequence
 Returns : a scalar string
 Args    : [optional] scalar string to set the name

=cut

sub name {
    my ($self,$value) = @_;
    my $previous = $self->{'_name'};
    if( defined $value || ! defined $previous ) {
	$value = $previous = '' unless defined $value;
	$self->{'_name'} = $value;
    } 
    return $previous;
}

=head2 accession

 Title   : accession
 Usage   : $acc = $hit->accession();
 Function: Retrieve the accession (if available) for the hit
 Returns : a scalar string (empty string if not set)
 Args    : none

=cut

sub accession {
    my ($self,$value) = @_;
    my $previous = $self->{'_accession'};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = '' unless defined $value;
	$self->{'_accession'} = $value;
    } 
    return $previous;
}

=head2 description

 Title   : description
 Usage   : $desc = $hit->description();
 Function: Retrieve the description for the hit
 Returns : a scalar string
 Args    : [optional] scalar string to set the descrition

=cut

sub description {
    my ($self,$value) = @_;
    my $previous = $self->{'_description'};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = '' unless defined $value;
	$self->{'_description'} = $value;
    } 
    return $previous;
}

=head2 length

 Title   : length
 Usage   : my $len = $hit->length
 Function: Returns the length of the hit 
 Returns : integer
 Args    : [optional] integer to set the length

=cut

sub length {
    my ($self,$value) = @_;
    my $previous = $self->{'_length'};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = 0 unless defined $value;
	$self->{'_length'} = $value;
    } 
    return $previous;
}


=head2 algorithm

 Title   : algorithm
 Usage   : $alg = $hit->algorithm();
 Function: Gets the algorithm specification that was used to obtain the hit
           For BLAST, the algorithm denotes what type of sequence was aligned 
           against what (BLASTN: dna-dna, BLASTP prt-prt, BLASTX translated 
           dna-prt, TBLASTN prt-translated dna, TBLASTX translated 
           dna-translated dna).
 Returns : a scalar string 
 Args    : [optional] scalar string to set the algorithm

=cut

sub algorithm {
    my ($self,$value) = @_;
    my $previous = $self->{'_algorithm'};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = '' unless defined $value;
	$self->{'_algorithm'} = $value;
    } 
    return $previous;
}

=head2 raw_score

 Title   : raw_score
 Usage   : $score = $hit->raw_score();
 Function: Gets the "raw score" generated by the algorithm.  What
           this score is exactly will vary from algorithm to algorithm,
           returning undef if unavailable.
 Returns : a scalar value
 Args    : [optional] scalar value to set the raw score

=cut

sub raw_score {
    my ($self,$value) = @_;
    my $previous = $self->{'_score'};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = '' unless defined $value;
	$self->{'_score'} = $value;
    } 
    return $previous;
}

=head2 significance

 Title   : significance
 Usage   : $significance = $hit->significance();
 Function: Used to obtain the E or P value of a hit, i.e. the probability that
           this particular hit was obtained purely by random chance.  If
           information is not available (nor calculatable from other
           information sources), return undef.
 Returns : a scalar value or undef if unavailable
 Args    : [optional] scalar value to set the significance

=cut

sub significance {
    my ($self,$value) = @_;
    my $previous = $self->{'_significance'};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = '' unless defined $value;
	$self->{'_significance'} = $value;
    } 
    return $previous;
}

=head2 bits

 Usage     : $hit_object->bits();
 Purpose   : Gets the bit score of the best HSP for the current hit.
 Example   : $bits = $hit_object->bits();
 Returns   : Integer or undef if bit score is not set
 Argument  : n/a
 Comments  : For BLAST1, the non-bit score is listed in the summary line.

See Also   : L<score()|score>

=cut

#---------
sub bits { 
#---------
    my ($self) = @_; 
    
    my $bits = $self->{'_bits'};
    if( ! defined $bits ) {
	$bits = $self->{'_hsps'}->[0]->bits();
	$self->{'_bits'} = $bits;
    } 
    return $bits;
}

=head2 next_hsp

 Title    : next_hsp
 Usage    : while( $hsp = $obj->next_hsp()) { ... }
 Function : Returns the next available High Scoring Pair
 Example  : 
 Returns  : Bio::Search::HSP::HSPI object or null if finished
 Args     : none

=cut

sub next_hsp {
    my ($self) = @_;
    $self->{'_iterator'} = 0 unless defined $self->{'_iterator'};
    return undef if $self->{'_iterator'} > scalar @{$self->{'_hsps'}};
    return $self->{'_hsps'}->[$self->{'_iterator'}++];    
}


=head2 hsps

 Usage     : $hit_object->hsps();
 Purpose   : Get a list containing all HSP objects.
           : Get the numbers of HSPs for the current hit.
 Example   : @hsps = $hit_object->hsps();
           : $num  = $hit_object->hsps();  # alternatively, use num_hsps()
 Returns   : Array context : list of Bio::Search::HSP::BlastHSP.pm objects.
           : Scalar context: integer (number of HSPs).
           :                 (Equivalent to num_hsps()).
 Argument  : n/a. Relies on wantarray
 Throws    : Exception if the HSPs have not been collected.

See Also   : L<hsp()|hsp>, L<num_hsps()|num_hsps>

=cut

#---------
sub hsps {
#---------
   my $self = shift;
   
   if (not ref $self->{'_hsps'}) {
       $self->throw("Can't get HSPs: data not collected.");
   }
   
   return wantarray 
       #  returning list containing all HSPs.
       ? @{$self->{'_hsps'}}
   #  returning number of HSPs.
   : scalar(@{$self->{'_hsps'}});
}

=head2 num_hsps

 Usage     : $hit_object->num_hsps();
 Purpose   : Get the number of HSPs for the present Blast hit.
 Example   : $nhsps = $hit_object->num_hsps();
 Returns   : Integer
 Argument  : n/a
 Throws    : Exception if the HSPs have not been collected.

See Also   : L<hsps()|hsps>

=cut

#-------------
sub num_hsps {
    my $self = shift;
    
    if (not defined $self->{'_hsps'}) {
	$self->throw("Can't get HSPs: data not collected.");
    }

    return scalar(@{$self->{'_hsps'}});
}

=head2 rewind

 Title   : rewind
 Usage   : $hit->rewind;
 Function: Allow one to reset the HSP iteration to the beginning
           Since this is an in-memory implementation
 Returns : none
 Args    : none

=cut

sub rewind{
   my ($self) = @_;
   $self->{'_iterator'} = 0;
}

=head2 iteration

 Title   : iteration
 Usage   : $obj->iteration($newval)
 Function: PSI-BLAST iteration
 Returns : value of iteration
 Args    : newvalue (optional)


=cut

sub iteration{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_psiblast_iteration'} = $value;
    }
    return $self->{'_psiblast_iteration'};

}


=head2 ambiguous_aln

 Usage     : $ambig_code = $hit_object->ambiguous_aln();
 Purpose   : Sets/Gets ambiguity code data member.
 Example   : (see usage)
 Returns   : String = 'q', 's', 'qs', '-'
           :   'q'  = query sequence contains overlapping sub-sequences 
           :          while sbjct does not.
           :   's'  = sbjct sequence contains overlapping sub-sequences 
           :          while query does not.
           :   'qs' = query and sbjct sequence contains overlapping sub-sequences
           :          relative to each other.
           :   '-'  = query and sbjct sequence do not contains multiple domains 
           :          relative to each other OR both contain the same distribution
           :          of similar domains.
 Argument  : n/a
 Throws    : n/a
 Status    : Experimental

=cut

#--------------------
sub ambiguous_aln { 
#--------------------
    my $self = shift;
    if(@_) { $self->{'_ambiguous_aln'} = shift; }
    $self->{'_ambiguous_aln'} || '-';
}

=head2 overlap

 Usage     : $blast_object->overlap( [integer] );
 Purpose   : Gets/Sets the allowable amount overlap between different HSP sequences.
 Example   : $blast_object->overlap(5);
           : $overlap = $blast_object->overlap;
 Returns   : Integer.
 Argument  : integer.
 Throws    : n/a
 Status    : Experimental
 Comments  : Any two HSPs whose sequences overlap by less than or equal
           : to the overlap() number of resides will be considered separate HSPs
           : and will not get tiled by Bio::Search::BlastUtils::_adjust_contigs().

See Also   : L<Bio::Search::BlastUtils::_adjust_contigs()|Bio::Search::BlastUtils>, L<BUGS | BUGS>

=cut

#-------------
sub overlap { 
#-------------
    my $self = shift; 
    if(@_) { $self->{'_overlap'} = shift; }
    defined $self->{'_overlap'} ? $self->{'_overlap'} : 0;
}


=head2 n

 Usage     : $hit_object->n();
 Purpose   : Gets the N number for the current Blast hit.
           : This is the number of HSPs in the set which was ascribed
           : the lowest P-value (listed on the description line).
           : This number is not the same as the total number of HSPs.
           : To get the total number of HSPs, use num_hsps().
 Example   : $n = $hit_object->n();
 Returns   : Integer
 Argument  : n/a
 Throws    : Exception if HSPs have not been set (BLAST2 reports).
 Comments  : Note that the N parameter is not reported in gapped BLAST2.
           : Calling n() on such reports will result in a call to num_hsps().
           : The num_hsps() method will count the actual number of
           : HSPs in the alignment listing, which may exceed N in
           : some cases.

See Also   : L<num_hsps()|num_hsps>

=cut

#-----
sub n { 
#-----
    my $self = shift; 

    # The check for $self->{'_n'} is a remnant from the 'query' mode days
    # in which the sbjct object would collect data from the description 
    # line only.

    my ($n);
    if(not defined($self->{'_n'})) {
	$n = $self->hsp->n;
    } else {
	$n = $self->{'_n'}; 
    } 
    $n ||= $self->num_hsps;

    return $n;
}

=head2 p

 Usage     : $hit_object->p( [format] );
 Purpose   : Get the P-value for the best HSP of the given BLAST hit.
           : (Note that P-values are not provided with NCBI Blast2 reports).
 Example   : $p =  $sbjct->p;
           : $p =  $sbjct->p('exp');  # get exponent only.
           : ($num, $exp) =  $sbjct->p('parts');  # split sci notation into parts
 Returns   : Float or scientific notation number (the raw P-value, DEFAULT).
           : Integer if format == 'exp' (the magnitude of the base 10 exponent).
           : 2-element list (float, int) if format == 'parts' and P-value
           :                is in scientific notation (See Comments).
 Argument  : format: string of 'raw' | 'exp' | 'parts'
           :    'raw' returns value given in report. Default. (1.2e-34)
           :    'exp' returns exponent value only (34)
           :    'parts' returns the decimal and exponent as a 
           :            2-element list (1.2, -34) (See Comments).
 Throws    : Warns if no P-value is defined. Uses expect instead.
 Comments  : Using the 'parts' argument is not recommended since it will not
           : work as expected if the P-value is not in scientific notation.
           : That is, floats are not converted into sci notation before
           : splitting into parts.

See Also   : L<expect()|expect>, L<signif()|signif>, L<Bio::Search::BlastUtils::get_exponent()|Bio::Search::BlastUtils>

=cut

#--------
sub p { 
#--------
# Some duplication of logic for p(), expect() and signif() for the sake of performance.
    my ($self, $fmt) = @_;

    my $val = $self->{'_p'};

    # $val can be zero.
    if(not defined $val) {
        # P-value not defined, must be a NCBI Blast2 report.
        # Use expect instead.
        $self->warn( "P-value not defined. Using expect() instead.");
        $val = $self->{'_expect'};
    }

    return $val if not $fmt or $fmt =~ /^raw/i;
    ## Special formats: exponent-only or as list.
    return &Bio::Search::BlastUtils::get_exponent($val) if $fmt =~ /^exp/i;
    return (split (/eE/, $val)) if $fmt =~ /^parts/i;

    ## Default: return the raw P-value.
    return $val;
}

=head2 hsp

 Usage     : $hit_object->hsp( [string] );
 Purpose   : Get a single HSPI object for the present HitI object.
 Example   : $hspObj  = $hit_object->hsp;  # same as 'best'
           : $hspObj  = $hit_object->hsp('best');
           : $hspObj  = $hit_object->hsp('worst');
 Returns   : Object reference for a Bio::Search::HSP::BlastHSP.pm object.
 Argument  : String (or no argument).
           :   No argument (default) = highest scoring HSP (same as 'best').
           :   'best' or 'first' = highest scoring HSP.
           :   'worst' or 'last' = lowest scoring HSP.
 Throws    : Exception if the HSPs have not been collected.
           : Exception if an unrecognized argument is used.

See Also   : L<hsps()|hsps>, L<num_hsps>()

=cut

#----------
sub hsp {
#----------
    my( $self, $option ) = @_;
    $option ||= 'best';
    
    if (not ref $self->{'_hsps'}) {
	$self->throw("Can't get HSPs: data not collected.");
    }

    my @hsps = @{$self->{'_hsps'}};
    
    return $hsps[0]      if $option =~ /best|first|1/i;
    return $hsps[$#hsps] if $option =~ /worst|last/i;

    $self->throw("Can't get HSP for: $option\n" .
		 "Valid arguments: 'best', 'worst'");
}

=head2 logical_length

 Usage     : $hit_object->logical_length( [seq_type] );
           : (mostly intended for internal use).
 Purpose   : Get the logical length of the hit sequence.
           : If the Blast is a TBLASTN or TBLASTX, the returned length 
           : is the length of the would-be amino acid sequence (length/3).
           : For all other BLAST flavors, this function is the same as length().
 Example   : $len    = $hit_object->logical_length();
 Returns   : Integer 
 Argument  : seq_type = 'query' or 'hit' or 'sbjct' (default = 'query')
             ('sbjct' is synonymous with 'hit')
 Throws    : n/a
 Comments  : This is important for functions like frac_aligned_query()
           : which need to operate in amino acid coordinate space when dealing
           : with [T]BLAST[NX] type reports.

See Also   : L<length()|length>, L<frac_aligned_query()|frac_aligned_query>, L<frac_aligned_hit()|frac_aligned_hit>

=cut

#--------------------
sub logical_length {
#--------------------
    my $self = shift;
    my $seqType = shift || 'query';
    $seqType = 'sbjct' if $seqType eq 'hit';

    my $length;

    # For the sbjct, return logical sbjct length
    if( $seqType eq 'sbjct' ) {
	$length = $self->{'_logical_length'} || $self->{'_length'};
    }
    else {
        # Otherwise, return logical query length
        $length = $self->{'_query_length'};
	$self->throw("Must have defined query_len") unless ( $length );

        # Adjust length based on BLAST flavor.
        if($self->algorithm =~ /T?BLASTX/ ) {
            $length /= 3;
        }
    }
    return $length;
}

=head2 length_aln

 Usage     : $hit_object->length_aln( [seq_type] );
 Purpose   : Get the total length of the aligned region for query or sbjct seq.
           : This number will include all HSPs
 Example   : $len    = $hit_object->length_aln(); # default = query
           : $lenAln = $hit_object->length_aln('query');
 Returns   : Integer 
 Argument  : seq_Type = 'query' or 'hit' or 'sbjct' (Default = 'query')
             ('sbjct' is synonymous with 'hit')
 Throws    : Exception if the argument is not recognized.
 Comments  : This method will report the logical length of the alignment,
           : meaning that for TBLAST[NX] reports, the length is reported
           : using amino acid coordinate space (i.e., nucleotides / 3).
           : 
           : This method requires that all HSPs be tiled. If they have not
           : already been tiled, they will be tiled first automatically..
           : If you don't want the tiled data, iterate through each HSP
           : calling length() on each (use hsps() to get all HSPs).

See Also   : L<length()|length>, L<frac_aligned_query()|frac_aligned_query>, L<frac_aligned_hit()|frac_aligned_hit>, L<gaps()|gaps>, L<Bio::Search::BlastUtils::tile_hsps()|Bio::Search::BlastUtils>, L<Bio::Search::HSP::BlastHSP::length()|Bio::Search::HSP::BlastHSP>

=cut

#---------------'
sub length_aln {
#---------------
    my( $self, $seqType ) = @_;
    
    $seqType ||= 'query';
    $seqType = 'sbjct' if $seqType eq 'hit';

    Bio::Search::BlastUtils::tile_hsps($self) if not $self->{'_tile_hsps'};

    my $data = $self->{'_length_aln_'.$seqType};
    
    ## If we don't have data, figure out what went wrong.
    if(!$data) {
	$self->throw("Can't get length aln for sequence type \"$seqType\"" . 
		     "Valid types are 'query', 'hit', 'sbjct' ('sbjct' = 'hit')");
    }		
    $data;
}    

=head2 gaps

 Usage     : $hit_object->gaps( [seq_type] );
 Purpose   : Get the number of gaps in the aligned query, sbjct, or both sequences.
           : Data is summed across all HSPs.
 Example   : $qgaps = $hit_object->gaps('query');
           : $hgaps = $hit_object->gaps('hit');
           : $tgaps = $hit_object->gaps();    # default = total (query + hit)
 Returns   : scalar context: integer
           : array context without args: two-element list of integers  
           :    (queryGaps, sbjctGaps)
           : Array context can be forced by providing an argument of 'list' or 'array'.
           :
           : CAUTION: Calling this method within printf or sprintf is arrray context.
           : So this function may not give you what you expect. For example:
           :          printf "Total gaps: %d", $hit->gaps();
           : Actually returns a two-element array, so what gets printed 
           : is the number of gaps in the query, not the total
           :
 Argument  : seq_type: 'query' | 'hit' or 'sbjct' | 'total' | 'list'  (default = 'total')
             ('sbjct' is synonymous with 'hit')
 Throws    : n/a
 Comments  : If you need data for each HSP, use hsps() and then interate
           : through each HSP object.
           : This method requires that all HSPs be tiled. If they have not
           : already been tiled, they will be tiled first automatically..
           : Not relying on wantarray since that will fail in situations 
           : such as printf "%d", $hit->gaps() in which you might expect to 
           : be printing the total gaps, but evaluates to array context.

See Also   : L<length_aln()|length_aln>

=cut

#----------
sub gaps {
#----------
    my( $self, $seqType ) = @_;

    $seqType ||= (wantarray ? 'list' : 'total');
    $seqType = 'sbjct' if $seqType eq 'hit';

    Bio::Search::BlastUtils::tile_hsps($self) if not $self->{'_tile_hsps'};

    $seqType = lc($seqType);

    if($seqType =~ /list|array/i) {
	return ($self->{'_gaps_query'}, $self->{'_gaps_sbjct'});
    }

    if($seqType eq 'total') {
	return ($self->{'_gaps_query'} + $self->{'_gaps_sbjct'}) || 0;
    } else {
	return $self->{'_gaps_'.$seqType} || 0;
    }
}    



=head2 matches

 Usage     : $hit_object->matches( [class] );
 Purpose   : Get the total number of identical or conserved matches 
           : (or both) across all HSPs.
           : (Note: 'conservative' matches are indicated as 'positives' 
	   :         in the Blast report.)
 Example   : ($id,$cons) = $hit_object->matches(); # no argument
           : $id = $hit_object->matches('id');
           : $cons = $hit_object->matches('cons'); 
 Returns   : Integer or a 2-element array of integers 
 Argument  : class = 'id' | 'cons' OR none. 
           : If no argument is provided, both identical and conservative 
           : numbers are returned in a two element list.
           : (Other terms can be used to refer to the conservative
           :  matches, e.g., 'positive'. All that is checked is whether or
           :  not the supplied string starts with 'id'. If not, the 
           : conservative matches are returned.)
 Throws    : Exception if the requested data cannot be obtained.
 Comments  : If you need data for each HSP, use hsps() and then interate
           : through the HSP objects.
           : Does not rely on wantarray to return a list. Only checks for
           : the presence of an argument (no arg = return list).

See Also   : L<Bio::Search::HSP::BlastHSP::matches()|Bio::Search::HSP::BlastHSP>, L<hsps()|hsps>

=cut

#---------------
sub matches {
#---------------
    my( $self, $arg) = @_;
    my(@data,$data);

    if(!$arg) {
	@data = ($self->{'_totalIdentical'}, $self->{'_totalConserved'});

	return @data if @data;

    } else {

	if($arg =~ /^id/i) { 
	    $data = $self->{'_totalIdentical'};
	} else {
	    $data = $self->{'_totalConserved'};
	}
	return $data if $data;
    }
    
    ## Something went wrong if we make it to here.
    $self->throw("Can't get identical or conserved data: no data.");
}


=head2 start

 Usage     : $sbjct->start( [seq_type] );
 Purpose   : Gets the start coordinate for the query, sbjct, or both sequences
           : in the BlastHit object. If there is more than one HSP, the lowest start
           : value of all HSPs is returned.
 Example   : $qbeg = $sbjct->start('query');
           : $sbeg = $sbjct->start('hit');
           : ($qbeg, $sbeg) = $sbjct->start();
 Returns   : scalar context: integer 
           : array context without args: list of two integers (queryStart, sbjctStart)
           : Array context can be "induced" by providing an argument of 'list' or 'array'.
 Argument  : In scalar context: seq_type = 'query' or 'hit' or 'sbjct' (default = 'query')
             ('sbjct' is synonymous with 'hit')
 Throws    : n/a
 Comments  : This method requires that all HSPs be tiled. If there is more than one
           : HSP and they have not already been tiled, they will be tiled first automatically..
           : Remember that the start and end coordinates of all HSPs are 
           : normalized so that start < end. Strand information can be
           : obtained by calling $hit->strand().

See Also   : L<end()|end>, L<range()|range>, L<strand()|strand>, 
             L<Bio::Search::HSP::BlastHSP::start|Bio::Search::HSP::BlastHSP>

=cut

#----------
sub start {
#----------
    my ($self, $seqType) = @_;

    $seqType ||= (wantarray ? 'list' : 'query');
    $seqType = 'sbjct' if $seqType eq 'hit';

    # If there is only one HSP, defer this call to the solitary HSP.
    if($self->num_hsps == 1) {
	return $self->hsp->start($seqType);
    } else {
	Bio::Search::BlastUtils::tile_hsps($self) if not $self->{'_tile_hsps'};
	if($seqType =~ /list|array/i) {
	    return ($self->{'_queryStart'}, $self->{'_sbjctStart'});
	} else {
	    ## Sensitive to member name changes.
	    $seqType = "_\L$seqType\E";
	    return $self->{$seqType.'Start'};
	}
    }
}


=head2 end

 Usage     : $sbjct->end( [seq_type] );
 Purpose   : Gets the end coordinate for the query, sbjct, or both sequences
           : in the BlastHit object. If there is more than one HSP, the largest end
           : value of all HSPs is returned.
 Example   : $qend = $sbjct->end('query');
           : $send = $sbjct->end('hit');
           : ($qend, $send) = $sbjct->end();
 Returns   : scalar context: integer
           : array context without args: list of two integers (queryEnd, sbjctEnd)
           : Array context can be "induced" by providing an argument of 'list' or 'array'.
 Argument  : In scalar context: seq_type = 'query' or 'sbjct'
           :  (case insensitive). If not supplied, 'query' is used.
 Throws    : n/a
 Comments  : This method requires that all HSPs be tiled. If there is more than one
           : HSP and they have not already been tiled, they will be tiled first automatically..
           : Remember that the start and end coordinates of all HSPs are 
           : normalized so that start < end. Strand information can be
           : obtained by calling $hit->strand().

See Also   : L<start()|start>, L<range()|range>, L<strand()|strand>

=cut

#----------
sub end {
#----------
    my ($self, $seqType) = @_;

    $seqType ||= (wantarray ? 'list' : 'query');
    $seqType = 'sbjct' if $seqType eq 'hit';

    # If there is only one HSP, defer this call to the solitary HSP.
    if($self->num_hsps == 1) {
	return $self->hsp->end($seqType);
    } else {
	Bio::Search::BlastUtils::tile_hsps($self) if not $self->{'_tile_hsps'};
	if($seqType =~ /list|array/i) {
	    return ($self->{'_queryStop'}, $self->{'_sbjctStop'});
	} else {
	    ## Sensitive to member name changes.
	    $seqType = "_\L$seqType\E";
	    return $self->{$seqType.'Stop'};
	}
    }
}

=head2 range

 Usage     : $sbjct->range( [seq_type] );
 Purpose   : Gets the (start, end) coordinates for the query or sbjct sequence
           : in the HSP alignment.
 Example   : ($qbeg, $qend) = $sbjct->range('query');
           : ($sbeg, $send) = $sbjct->range('hit');
 Returns   : Two-element array of integers 
 Argument  : seq_type = string, 'query' or 'hit' or 'sbjct'  (default = 'query')
             ('sbjct' is synonymous with 'hit')
 Throws    : n/a

See Also   : L<start()|start>, L<end()|end>

=cut

#----------
sub range {
#----------
    my ($self, $seqType) = @_;
    $seqType ||= 'query';
    $seqType = 'sbjct' if $seqType eq 'hit';
    return ($self->start($seqType), $self->end($seqType));
}


=head2 frac_identical

 Usage     : $hit_object->frac_identical( [seq_type] );
 Purpose   : Get the overall fraction of identical positions across all HSPs.
           : The number refers to only the aligned regions and does not
           : account for unaligned regions in between the HSPs, if any.
 Example   : $frac_iden = $hit_object->frac_identical('query');
 Returns   : Float (2-decimal precision, e.g., 0.75).
 Argument  : seq_type: 'query' | 'hit' or 'sbjct' | 'total'
           : default = 'query' (but see comments below).
           : ('sbjct' is synonymous with 'hit')
 Throws    : n/a
 Comments  : Different versions of Blast report different values for the total
           : length of the alignment. This is the number reported in the
           : denominators in the stats section:
           : "Identical = 34/120 Positives = 67/120".
           : NCBI BLAST uses the total length of the alignment (with gaps)
           : WU-BLAST uses the length of the query sequence (without gaps).
           :
           : Therefore, when called with an argument of 'total',
           : this method will report different values depending on the
           : version of BLAST used. Total does NOT take into account HSP
           : tiling, so it should not be used.
           :
           : To get the fraction identical among only the aligned residues,
           : ignoring the gaps, call this method without an argument or 
           : with an argument of 'query' or 'hit'.
           :
           : If you need data for each HSP, use hsps() and then iterate
           : through the HSP objects.
           : This method requires that all HSPs be tiled. If they have not
           : already been tiled, they will be tiled first automatically.

See Also   : L<frac_conserved()|frac_conserved>, L<frac_aligned_query()|frac_aligned_query>, L<matches()|matches>, L<Bio::Search::BlastUtils::tile_hsps()|Bio::Search::BlastUtils>

=cut

#------------------
sub frac_identical {
#------------------
    my ($self, $seqType) = @_;
    $seqType ||= 'query';
    $seqType = 'sbjct' if $seqType eq 'hit';

    ## Sensitive to member name format.
    $seqType = lc($seqType);

    Bio::Search::BlastUtils::tile_hsps($self) if not $self->{'_tile_hsps'};

    my $ident = $self->{'_totalIdentical'};
    my $total = $self->{'_length_aln_'.$seqType};
    my $ratio = $ident / $total;
    my $ratio_rounded = sprintf( "%.3f", $ratio);

    # Round down iff normal rounding yields 1 (just like blast)
    $ratio_rounded = 0.999 if (($ratio_rounded == 1) && ($ratio < 1));
    return $ratio_rounded;
}



=head2 frac_conserved

 Usage     : $hit_object->frac_conserved( [seq_type] );
 Purpose   : Get the overall fraction of conserved positions across all HSPs.
           : The number refers to only the aligned regions and does not
           : account for unaligned regions in between the HSPs, if any.
 Example   : $frac_cons = $hit_object->frac_conserved('hit');
 Returns   : Float (2-decimal precision, e.g., 0.75).
 Argument  : seq_type: 'query' | 'hit' or 'sbjct' | 'total'
           : default = 'query' (but see comments below).
           : ('sbjct' is synonymous with 'hit')
 Throws    : n/a
 Comments  : Different versions of Blast report different values for the total
           : length of the alignment. This is the number reported in the
           : denominators in the stats section:
           : "Positives = 34/120 Positives = 67/120".
           : NCBI BLAST uses the total length of the alignment (with gaps)
           : WU-BLAST uses the length of the query sequence (without gaps).
           :
           : Therefore, when called with an argument of 'total',
           : this method will report different values depending on the
           : version of BLAST used. Total does NOT take into account HSP
           : tiling, so it should not be used.
           :
           : To get the fraction conserved among only the aligned residues,
           : ignoring the gaps, call this method without an argument or 
           : with an argument of 'query' or 'hit'.
           :
           : If you need data for each HSP, use hsps() and then interate
           : through the HSP objects.
           : This method requires that all HSPs be tiled. If they have not
           : already been tiled, they will be tiled first automatically.

See Also   : L<frac_identical()|frac_identical>, L<matches()|matches>, L<Bio::Search::BlastUtils::tile_hsps()|Bio::Search::BlastUtils>

=cut

#--------------------
sub frac_conserved {
#--------------------
    my ($self, $seqType) = @_;
    $seqType ||= 'query';
    $seqType = 'sbjct' if $seqType eq 'hit';

    ## Sensitive to member name format.
    $seqType = lc($seqType);

    Bio::Search::BlastUtils::tile_hsps($self) if not $self->{'_tile_hsps'};

    my $consv = $self->{'_totalConserved'};
    my $total = $self->{'_length_aln_'.$seqType};
    my $ratio = $consv / $total;
    my $ratio_rounded = sprintf( "%.3f", $ratio);

    # Round down iff normal rounding yields 1 (just like blast)
    $ratio_rounded = 0.999 if (($ratio_rounded == 1) && ($ratio < 1));
    return $ratio_rounded;
}




=head2 frac_aligned_query

 Usage     : $hit_object->frac_aligned_query();
 Purpose   : Get the fraction of the query sequence which has been aligned
           : across all HSPs (not including intervals between non-overlapping
           : HSPs).
 Example   : $frac_alnq = $hit_object->frac_aligned_query();
 Returns   : Float (2-decimal precision, e.g., 0.75).
 Argument  : n/a
 Throws    : n/a
 Comments  : If you need data for each HSP, use hsps() and then interate
           : through the HSP objects.
           : To compute the fraction aligned, the logical length of the query
           : sequence is used, meaning that for [T]BLASTX reports, the 
           : full length of the query sequence is converted into amino acids
           : by dividing by 3. This is necessary because of the way 
           : the lengths of aligned sequences are computed.
           : This method requires that all HSPs be tiled. If they have not
           : already been tiled, they will be tiled first automatically.

See Also   : L<frac_aligned_hit()|frac_aligned_hit>, L<logical_length()|logical_length>, L<length_aln()|length_aln>,  L<Bio::Search::BlastUtils::tile_hsps()|Bio::Search::BlastUtils>

=cut

#----------------------
sub frac_aligned_query {
#----------------------
    my $self = shift;

    Bio::Search::BlastUtils::tile_hsps($self) if not $self->{'_tile_hsps'};

    sprintf( "%.2f", $self->{'_length_aln_query'}/ 
	     $self->logical_length('query'));
}



=head2 frac_aligned_hit

 Usage     : $hit_object->frac_aligned_hit();
 Purpose   : Get the fraction of the hit (sbjct) sequence which has been aligned
           : across all HSPs (not including intervals between non-overlapping
           : HSPs).
 Example   : $frac_alnq = $hit_object->frac_aligned_hit();
 Returns   : Float (2-decimal precision, e.g., 0.75).
 Argument  : n/a
 Throws    : n/a
 Comments  : If you need data for each HSP, use hsps() and then interate
           : through the HSP objects.
           : To compute the fraction aligned, the logical length of the sbjct
           : sequence is used, meaning that for TBLAST[NX] reports, the 
           : full length of the sbjct sequence is converted into amino acids
           : by dividing by 3. This is necessary because of the way 
           : the lengths of aligned sequences are computed.
           : This method requires that all HSPs be tiled. If they have not
           : already been tiled, they will be tiled first automatically.

See Also   : L<frac_aligned_query()|frac_aligned_query>, L<matches()|matches>, , L<logical_length()|logical_length>, L<length_aln()|length_aln>,  L<Bio::Search::BlastUtils::tile_hsps()|Bio::Search::BlastUtils>

=cut

#--------------------
sub frac_aligned_hit {
#--------------------
    my $self = shift;

    Bio::Search::BlastUtils::tile_hsps($self) if not $self->{'_tile_hsps'};

    sprintf( "%.2f", $self->{'_length_aln_sbjct'}/$self->logical_length('sbjct'));
}


## These methods are being maintained for backward compatibility. 

=head2 frac_aligned_sbjct

Same as L<frac_aligned_hit()|frac_aligned_hit>

=cut

#----------------
sub frac_aligned_sbjct {  my $self=shift; $self->frac_aligned_hit(@_); }
#----------------

=head2 num_unaligned_sbjct

Same as L<num_unaligned_hit()|num_unaligned_hit>

=cut

#----------------
sub num_unaligned_sbjct {  my $self=shift; $self->num_unaligned_hit(@_); }
#----------------



=head2 num_unaligned_hit

 Usage     : $hit_object->num_unaligned_hit();
 Purpose   : Get the number of the unaligned residues in the hit sequence.
           : Sums across all all HSPs.
 Example   : $num_unaln = $hit_object->num_unaligned_hit();
 Returns   : Integer
 Argument  : n/a
 Throws    : n/a
 Comments  : See notes regarding logical lengths in the comments for frac_aligned_hit().
           : They apply here as well.
           : If you need data for each HSP, use hsps() and then interate
           : through the HSP objects.
           : This method requires that all HSPs be tiled. If they have not
           : already been tiled, they will be tiled first automatically..

See Also   : L<num_unaligned_query()|num_unaligned_query>,  L<Bio::Search::BlastUtils::tile_hsps()|Bio::Search::BlastUtils>, L<frac_aligned_hit()|frac_aligned_hit>

=cut

#---------------------
sub num_unaligned_hit {
#---------------------
    my $self = shift;

    Bio::Search::BlastUtils::tile_hsps($self) if not $self->{'_tile_hsps'};

    my $num = $self->logical_length('sbjct') - $self->{'_length_aln_sbjct'};
    ($num < 0 ? 0 : $num );
}


=head2 num_unaligned_query

 Usage     : $hit_object->num_unaligned_query();
 Purpose   : Get the number of the unaligned residues in the query sequence.
           : Sums across all all HSPs.
 Example   : $num_unaln = $hit_object->num_unaligned_query();
 Returns   : Integer
 Argument  : n/a
 Throws    : n/a
 Comments  : See notes regarding logical lengths in the comments for frac_aligned_query().
           : They apply here as well.
           : If you need data for each HSP, use hsps() and then interate
           : through the HSP objects.
           : This method requires that all HSPs be tiled. If they have not
           : already been tiled, they will be tiled first automatically..

See Also   : L<num_unaligned_hit()|num_unaligned_hit>, L<frac_aligned_query()|frac_aligned_query>,  L<Bio::Search::BlastUtils::tile_hsps()|Bio::Search::BlastUtils>

=cut

#-----------------------
sub num_unaligned_query {
#-----------------------
    my $self = shift;

    Bio::Search::BlastUtils::tile_hsps($self) if not $self->{'_tile_hsps'};

    my $num = $self->logical_length('query') - $self->{'_length_aln_query'};
    ($num < 0 ? 0 : $num );
}



=head2 seq_inds

 Usage     : $hit->seq_inds( seq_type, class, collapse );
 Purpose   : Get a list of residue positions (indices) across all HSPs
           : for identical or conserved residues in the query or sbjct sequence.
 Example   : @s_ind = $hit->seq_inds('query', 'identical');
           : @h_ind = $hit->seq_inds('hit', 'conserved');
           : @h_ind = $hit->seq_inds('hit', 'conserved', 1);
 Returns   : Array of integers 
           : May include ranges if collapse is non-zero.
 Argument  : [0] seq_type  = 'query' or 'hit' or 'sbjct'  (default = 'query')
           :                 ('sbjct' is synonymous with 'hit')
           : [1] class     = 'identical' or 'conserved' (default = 'identical')
           :              (can be shortened to 'id' or 'cons')
           :              (actually, anything not 'id' will evaluate to 'conserved').
           : [2] collapse  = boolean, if non-zero, consecutive positions are merged
           :             using a range notation, e.g., "1 2 3 4 5 7 9 10 11" 
           :             collapses to "1-5 7 9-11". This is useful for 
           :             consolidating long lists. Default = no collapse.
 Throws    : n/a.

See Also   : L<Bio::Search::HSP::BlastHSP::seq_inds()|Bio::Search::HSP::BlastHSP>

=cut

#-------------
sub seq_inds {
#-------------
    my ($self, $seqType, $class, $collapse) = @_;

    $seqType  ||= 'query';
    $class ||= 'identical';
    $collapse ||= 0;

    $seqType = 'sbjct' if $seqType eq 'hit';

    my (@inds, $hsp);
    foreach $hsp ($self->hsps) {
	# This will merge data for all HSPs together.
	push @inds, $hsp->seq_inds($seqType, $class);
    }
    
    # Need to remove duplicates and sort the merged positions.
    if(@inds) {
	my %tmp = map { $_, 1 } @inds;
	@inds = sort {$a <=> $b} keys %tmp;
    }

    $collapse ?  &Bio::Search::BlastUtils::collapse_nums(@inds) : @inds; 
}

=head2 strand

 Usage     : $sbjct->strand( [seq_type] );
 Purpose   : Gets the strand(s) for the query, sbjct, or both sequences
           : in the BlastHit object. 
 Example   : $qstrand = $sbjct->strand('query');
           : $sstrand = $sbjct->strand('hit');
           : ($qstrand, $sstrand) = $sbjct->strand();
 Returns   : scalar context: integer '1', '-1', or '-1/1'
           : if there are HSPs on both strands.
           : array context without args: list of two strings (queryStrand, sbjctStrand)
           : Array context can be "induced" by providing an argument of 'list' or 'array'.
 Argument  : In scalar context: seq_type = 'query' or 'hit' or 'sbjct' (default = 'query')
             ('sbjct' is synonymous with 'hit')
 Throws    : n/a

See Also   : B<Bio::Search::HSP::BlastHSP::strand>()

=cut

#----------
sub strand {
#----------
    my ($self, $seqType) = @_;

    $seqType ||= (wantarray ? 'list' : 'query');
    $seqType = 'sbjct' if $seqType eq 'hit';

    # If there is only one HSP, defer this call to the solitary HSP.
    if($self->num_hsps == 1) {
	return $self->hsp->strand($seqType);
    } else {
	# otherwise, iterate through all HSPs collecting strand info.
        my (%qstr, %hstr);
        foreach my $hsp( $self->hsps ) {
            my ( $q, $h ) = $hsp->strand();
            $qstr{ $q }++;
            $hstr{ $h }++;
        }
        my $qstr = join( '/', sort keys %qstr);
        my $hstr = join( '/', sort keys %hstr);
	if($seqType =~ /list|array/i) {
	    return ($qstr, $hstr);
	} elsif( $seqType eq 'query' ) {
	    return $qstr;
	} else {
	    return $hstr;
	}
    }
}

=head2 rank

 Title   : rank
 Usage   : $obj->rank($newval)
 Function: Get/Set the rank of this Hit in the Query search list
           i.e. this is the Nth hit for a specific query
 Returns : value of rank
 Args    : newvalue (optional)


=cut

sub rank{
    my $self = shift;
    return $self->{'_rank'} = shift if @_;
    return $self->{'_rank'} || 1;
}

1;
