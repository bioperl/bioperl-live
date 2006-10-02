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

    $r_type = $hsp->algorithm;

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

# TODO: Describe how to configure a SearchIO stream so that it generates
#       GenericHSP objects.


=head1 DESCRIPTION

This implementation is "Generic", meaning it is is suitable for
holding information about High Scoring pairs from most Search reports
such as BLAST and FastA.  Specialized objects can be derived from
this.

Unless you're writing a parser, you won't ever need to create a
GenericHSP or any other HSPI-implementing object. If you use
the SearchIO system, HSPI objects are created automatically from
a SearchIO stream which returns Bio::Search::Result::ResultI objects
and you get the HSPI objects via the ResultI API.

For documentation on what you can do with GenericHSP (and other HSPI
objects), please see the API documentation in
L<Bio::Search::HSP::HSPI|Bio::Search::HSP::HSPI>.

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

=head1 AUTHOR - Jason Stajich and Steve Chervitz

Email jason-at-bioperl.org
Email sac-at-bioperl.org

=head1 CONTRIBUTORS

Sendu Bala, bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Search::HSP::GenericHSP;
use vars qw($GAP_SYMBOL);
use strict;

use Bio::Root::Root;
use Bio::SeqFeature::Similarity;

use base qw(Bio::Search::HSP::HSPI);

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
           -percent_identity => (optional) percent identity
           -conserved => # of residues that matched conservatively
                           (only protein comparisions;
                            conserved == identical in nucleotide comparisons)
           -hsp_gaps   => # of gaps in the HSP
           -query_gaps => # of gaps in the query in the alignment
           -hit_gaps   => # of gaps in the subject in the alignment
           -query_name  => HSP Query sequence name (if available)
           -query_start => HSP Query start (in original query sequence coords)
           -query_end   => HSP Query end (in original query sequence coords)
           -query_length=> total length of the query sequence
           -query_seq   => query sequence portion of the HSP
           -query_desc  => textual description of the query
           -hit_name    => HSP Hit sequence name (if available)
           -hit_start   => HSP Hit start (in original hit sequence coords)
           -hit_end     => HSP Hit end (in original hit sequence coords)
           -hit_length  => total length of the hit sequence
           -hit_seq     => hit sequence portion of the HSP
           -hit_desc    => textual description of the hit
           -homology_seq=> homology sequence for the HSP
           -hit_frame   => hit frame (only if hit is translated protein)
           -query_frame => query frame (only if query is translated protein)
           -rank        => HSP rank
           -links       => HSP links information (WU-BLAST only)
           -hsp_group   => HSP Group informat (WU-BLAST only)
=cut

sub new {
    my($class,@args) = @_;

    # don't pass anything to SUPER; complex heirarchy results in lots of work
    # for nothing
    my $self = $class->SUPER::new();

    # for speed, don't use _rearrange and just store all input data directly
    # with no method calls and no work done. work can be carried
    # out just-in-time later if desired
    my %args = @args;
    while (my ($arg, $value) = each %args) {
        $arg =~ tr/a-z\055/A-Z/d;
        $self->{$arg} = $value;
    }
    my $bits = $self->{BITS};

    defined $self->{VERBOSE} && $self->verbose($self->{VERBOSE});

    $self->{ALGORITHM} ||= 'GENERIC';

    if (! defined $self->{QUERY_LENGTH} || ! defined $self->{HIT_LENGTH}) {
        $self->throw("Must define hit and query length");
    }

    $self->{'_sequenceschanged'} = 1;
    
    $self->{_finished_new} = 1;
    return $self;
}


sub _logical_length {
    my ($self, $type) = @_;
    my $algo = $self->algorithm;
    my $len = $self->length($type);
    my $logical = $len;
    if($algo =~ /^(PSI)?T(BLAST|FAST)[NY]/oi ) {
        $logical = $len/3 if $type =~ /sbjct|hit|tot/i;
    } elsif($algo =~ /^(BLAST|FAST)(X|Y|XY)/oi ) {
        $logical = $len/3 if $type =~ /query|tot/i;
    } elsif($algo =~ /^T(BLAST|FAST)(X|Y|XY)/oi ) {
        $logical = $len/3;
    }
    return $logical;
}

=head2 L<Bio::Search::HSP::HSPI> methods

Implementation of L<Bio::Search::HSP::HSPI> methods follow

=head2 algorithm

 Title   : algorithm
 Usage   : my $r_type = $hsp->algorithm
 Function: Obtain the name of the algorithm used to obtain the HSP
 Returns : string (e.g., BLASTP)
 Args    : [optional] scalar string to set value

=cut

sub algorithm{
    my ($self,$value) = @_;
    my $previous = $self->{'ALGORITHM'};
    if( defined $value || ! defined $previous ) {
        $value = $previous = '' unless defined $value;
        $self->{'ALGORITHM'} = $value;
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
    my $previous = $self->{'PVALUE'};
    if( defined $value  ) {
        $self->{'PVALUE'} = $value;
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

# Override significance to return the e-value or, if this is
# not defined (WU-BLAST), return the p-value.
sub significance {
    my $self = shift;
    my $signif = $self->query->significance(@_);
    return (defined $signif && $signif ne '') ? $signif : $self->pvalue(@_);
}


=head2 frac_identical

 Title   : frac_identical
 Usage   : my $frac_id = $hsp->frac_identical( ['query'|'hit'|'total'] );
 Function: Returns the fraction of identitical positions for this HSP
 Returns : Float in range 0.0 -> 1.0
 Args    : arg 1:  'query' = num identical / length of query seq (without gaps)
                   'hit'   = num identical / length of hit seq (without gaps)
                             synonyms: 'sbjct', 'subject'
                   'total' = num conserved / length of alignment (with gaps)
                             synonyms: 'hsp'
                   default = 'total'
           arg 2: [optional] frac identical value to set for the type requested

=cut

sub frac_identical {
   my ($self, $type,$value) = @_;

    unless ($self->{_did_prefrac}) {
        $self->_pre_frac;
    }

   $type = lc $type if defined $type;
   $type = 'hit' if( defined $type &&
		     $type =~ /subject|sbjct/);
   $type = 'total' if( ! defined $type || $type eq 'hsp' ||
                        $type !~ /query|hit|subject|sbjct|total/);
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
                             synonyms: 'sbjct', 'subject'
                  'total' = num conserved / length of alignment (with gaps)
                             synonyms: 'hsp'
                  default = 'total'
           arg 2: [optional] frac conserved value to set for the type requested

=cut

sub frac_conserved {
    my ($self, $type,$value) = @_;

    unless ($self->{_did_prefrac}) {
        $self->_pre_frac;
    }

    $type = lc $type if defined $type;
    $type = 'hit' if( defined $type && $type =~ /subject|sbjct/);
    $type = 'total' if( ! defined $type || $type eq 'hsp' ||
                        $type !~ /query|hit|subject|sbjct|total/);
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
                   'hit'   = num gaps in hit seq; synonyms: 'sbjct', 'subject'
                   'total' = num gaps in whole alignment;  synonyms: 'hsp'
                   default = 'total'
            arg 2: [optional] integer gap value to set for the type requested

=cut

sub gaps {
    my ($self, $type,$value) = @_;

    unless ($self->{_did_pregaps}) {
        $self->_pre_gaps;
    }

    $type = lc $type if defined $type;
    $type = 'total' if( ! defined $type || $type eq 'hsp' ||
                        $type !~ /query|hit|subject|sbjct|total/);
    $type = 'hit' if $type =~ /sbjct|subject/;
    my $previous = $self->{'_gaps'}->{$type};
    if( defined $value || ! defined $previous ) {
        $value = $previous = '' unless defined $value;
        $self->{'_gaps'}->{$type} = $value;
    }
    return $previous || 0;
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
    my $previous = $self->{QUERY_SEQ};
    if( defined $value || ! defined $previous ) {
        $value = $previous = '' unless defined $value;
        $self->{QUERY_SEQ} = $value;
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
    my $previous = $self->{HIT_SEQ};
    if( defined $value || ! defined $previous ) {
        $value = $previous = '' unless defined $value;
        $self->{HIT_SEQ} = $value;
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
    my $previous = $self->{HOMOLOGY_SEQ};
    if( defined $value || ! defined $previous ) {
        $value = $previous = '' unless defined $value;
        $self->{HOMOLOGY_SEQ} = $value;
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
                   'hit'   = length of hit seq (without gaps) (synonyms: sbjct, subject)
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
            $self->{HSP_LENGTH} = $v;
        }
        return $self->{HSP_LENGTH};
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

sub percent_identity {
    my $self = shift;

    unless ($self->{_did_prepi}) {
        $self->_pre_pi;
    }

    return $self->SUPER::percent_identity(@_);
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
 Note    : Frames are stored in the GFF way (0-2) not 1-3
           as they are in BLAST (negative frames are deduced by checking
                                 the strand of the query or hit)

=cut

sub frame {
    my $self = shift;

    unless (defined $self->{_did_preframe}) {
        $self->_pre_frame;
    }
    my $qframe = $self->{QUERY_FRAME};
    my $sframe = $self->{HIT_FRAME};

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
    if (wantarray() && $self->algorithm =~ /^T(BLAST|FAST)(X|Y|XY)/oi)
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
 Function: Returns a L<Bio::SimpleAlign> object representing the HSP alignment
 Returns : L<Bio::SimpleAlign>
 Args    : none

=cut

sub get_aln {
    my ($self) = @_;
    require Bio::LocatableSeq;
    require Bio::SimpleAlign;
    my $aln = new Bio::SimpleAlign;
    my $hs = $self->hit_string();
    my $qs = $self->query_string();
    # FASTA specific stuff moved to the FastaHSP object
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
                                      '-start' => $self->query->start,
                                      '-end'   => $self->query->end,
                                      );
    $seqonly = $hs;
    $seqonly =~ s/[\-\s]//g;
    my $hit =  new Bio::LocatableSeq('-seq'    => $hs,
                                      '-id'    => $s_nm,
                                      '-start' => $self->hit->start,
                                      '-end'   => $self->hit->end,
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

    unless ($self->{_did_presimilar}) {
        $self->_pre_similar_stats;
    }

    if (defined $value) {
        $self->{CONSERVED} = $value;
    }
    return $self->{CONSERVED};
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

   unless ($self->{_did_presimilar}) {
        $self->_pre_similar_stats;
    }

   if( defined $value) {
       $self->{IDENTICAL} = $value;
   }
   return $self->{IDENTICAL};
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
        $self->{RANK} = $value;
    }
    return $self->{RANK};
}


=head2 seq_inds

 Title   : seq_inds
 Purpose   : Get a list of residue positions (indices) for all identical
           : or conserved residues in the query or sbjct sequence.
 Example   : @s_ind = $hsp->seq_inds('query', 'identical');
           : @h_ind = $hsp->seq_inds('hit', 'conserved');
           : @h_ind = $hsp->seq_inds('hit', 'conserved-not-identical');
           : @h_ind = $hsp->seq_inds('hit', 'conserved', 1);
 Returns   : List of integers
           : May include ranges if collapse is true.
 Argument  : seq_type  = 'query' or 'hit' or 'sbjct'  (default = query)
           :  ('sbjct' is synonymous with 'hit')
           : class     = 'identical' or 'conserved' or 'nomatch' or 'gap'
           :              (default = identical)
           :              (can be shortened to 'id' or 'cons')
           :             or 'conserved-not-identical'
           : collapse  = boolean, if true, consecutive positions are merged
           :             using a range notation, e.g., "1 2 3 4 5 7 9 10 11"
           :             collapses to "1-5 7 9-11". This is useful for
           :             consolidating long lists. Default = no collapse.
 Throws    : n/a.
 Comments  :

See Also   : L<Bio::Search::SearchUtils::collapse_nums()|Bio::Search::SearchUtils>,
             L<Bio::Search::Hit::HitI::seq_inds()|Bio::Search::Hit::HitI>

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
       if( $class =~ /conserved\-not\-identical/ ) {
	   $class = 'conserved';
       } else {
	   $class = 'conservedall';
       }
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
   } elsif( $class eq '_conservedall' ) {
       @ary = sort { $a <=> $b }
       keys %{ $self->{"_conservedRes$seqType"}},
       keys %{ $self->{"_identicalRes$seqType"}},
   }  else {
       @ary = sort { $a <=> $b } keys %{ $self->{"${class}Res$seqType"}};
   }
   require Bio::Search::BlastUtils if $collapse;

   return $collapse ? &Bio::Search::SearchUtils::collapse_nums(@ary) : @ary;
}


=head2 Inherited from L<Bio::SeqFeature::SimilarityPair>

These methods come from L<Bio::SeqFeature::SimilarityPair>

=head2 query

 Title   : query
 Usage   : my $query = $hsp->query
 Function: Returns a SeqFeature representing the query in the HSP
 Returns : L<Bio::SeqFeature::Similarity>
 Args    : [optional] new value to set

=cut

sub query {
    my $self = shift;
    unless ($self->{_created_qff}) {
        $self->_query_seq_feature;
    }
    return $self->SUPER::query(@_);
}

sub feature1 {
    my $self = shift;
    if (! $self->{_finished_new} || $self->{_making_qff}) {
        return $self->{_sim1} if $self->{_sim1};
        $self->{_sim1} = Bio::SeqFeature::Similarity->new();
        return $self->{_sim1};
    }
    unless ($self->{_created_qff}) {
        $self->_query_seq_feature;
    }
    return $self->SUPER::feature1(@_);
}

=head2 hit

 Title   : hit
 Usage   : my $hit = $hsp->hit
 Function: Returns a SeqFeature representing the hit in the HSP
 Returns : L<Bio::SeqFeature::Similarity>
 Args    : [optional] new value to set

=cut

sub hit {
    my $self = shift;
    unless ($self->{_created_sff}) {
        $self->_subject_seq_feature;
    }
    return $self->SUPER::hit(@_);
}

sub feature2 {
    my $self = shift;
    if (! $self->{_finished_new} || $self->{_making_sff}) {
        return $self->{_sim2} if $self->{_sim2};
        $self->{_sim2} = Bio::SeqFeature::Similarity->new();
        return $self->{_sim2};
    }
    unless ($self->{_created_sff}) {
        $self->_subject_seq_feature;
    }
    return $self->SUPER::feature2(@_);
}

=head2 significance

 Title   : significance
 Usage   : $evalue = $obj->significance();
           $obj->significance($evalue);
 Function: Get/Set the significance value
 Returns : numeric
 Args    : [optional] new value to set


=head2 score

 Title   : score
 Usage   : $score = $obj->score();
           $obj->score($value);
 Function: Get/Set the score value
 Returns : numeric
 Args    : [optional] new value to set


=head2 bits

 Title   : bits
 Usage   : $bits = $obj->bits();
           $obj->bits($value);
 Function: Get/Set the bits value
 Returns : numeric
 Args    : [optional] new value to set


=head2 strand

 Title   : strand
 Usage   : $hsp->strand('query')
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
 Function: Internal function
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
    if( $prog  =~ /FAST|SSEARCH|SMITH-WATERMAN/i ) {
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

    if($prog =~ /^(PSI)?T(BLAST|FAST)N/oi ) {
	$resCount_sbjct = int($resCount_sbjct / 3);
    } elsif($prog =~ /^(BLAST|FAST)(X|Y|XY)/oi  ) {
	$resCount_query = int($resCount_query / 3);
    } elsif($prog =~ /^T(BLAST|FAST)(X|Y|XY)/oi ) {
	$resCount_query = int($resCount_query / 3);
	$resCount_sbjct = int($resCount_sbjct / 3);
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
    $self->{'_nomatchRes_query'}   = \%nomatchList_query;
    $self->{'_gapRes_query'}       = \%gapList_query;

    $self->{'_identicalRes_sbjct'} = \%identicalList_sbjct;
    $self->{'_conservedRes_sbjct'} = \%conservedList_sbjct;
    $self->{'_nomatchRes_sbjct'}   = \%nomatchList_sbjct;
    $self->{'_gapRes_sbjct'}       = \%gapList_sbjct;
    return 1;
}

=head2 n

See documentation in L<Bio::Search::HSP::HSPI::n()|Bio::Search::HSP::HSPI>

=cut

sub n {
    my $self = shift;
    if(@_) { $self->{'_n'} = shift; }
    defined $self->{'_n'} ? $self->{'_n'} : '';
}

=head2 range

See documentation in L<Bio::Search::HSP::HSPI::range()|Bio::Search::HSP::HSPI>

=cut

sub range {
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


=head2 links

 Title   : links
 Usage   : $obj->links($newval)
 Function: Get/Set the Links value (from WU-BLAST)
           Indicates the placement of the alignment in the group of HSPs
 Returns : Value of links
 Args    : On set, new value (a scalar or undef, optional)


=cut

sub links{
    my $self = shift;

    return $self->{LINKS} = shift if @_;
    return $self->{LINKS};
}

=head2 hsp_group

 Title   : hsp_group
 Usage   : $obj->hsp_group($newval)
 Function: Get/Set the Group value (from WU-BLAST)
           Indicates a grouping of HSPs
 Returns : Value of group
 Args    : On set, new value (a scalar or undef, optional)


=cut

sub hsp_group {
    my $self = shift;

    return $self->{HSP_GROUP} = shift if @_;
    return $self->{HSP_GROUP};
}

# The cigar string code is written by Juguang Xiao <juguang@fugu-sg.org>

=head1 Brief introduction on cigar string

NOTE: the concept is originally from EnsEMBL docs at
http://may2005.archive.ensembl.org/Docs/wiki/html/EnsemblDocs/CigarFormat.html
Please append to these docs if you have a better definition.

Sequence alignment hits can be stored in a database as ungapped alignments.
This imposes 2 major constraints on alignments:

a) alignments for a single hit record require multiple rows in the database,
and
b) it is not possible to accurately retrieve the exact original alignment.

Alternatively, sequence alignments can be stored as gapped alignments using
the CIGAR line format (where CIGAR stands for Concise Idiosyncratic Gapped
Alignment Report).

In the cigar line format alignments are stored as follows:

M: Match
D: Deletion
I: Insertion

An example of an alignment for a hypthetical protein match is shown below:


Query:   42 PGPAGLP----GSVGLQGPRGLRGPLP-GPLGPPL...

            PG    P    G     GP   R      PLGP

Sbjct: 1672 PGTP*TPLVPLGPWVPLGPSSPR--LPSGPLGPTD...


protein_align_feature table as the following cigar line:

7M4D12M2I2MD7M

=head2 cigar_string

  Name:     cigar_string
  Usage:    $cigar_string = $hsp->cigar_string
  Function: Generate and return cigar string for this HSP alignment
  Args:     No input needed
  Return:   a cigar string

=cut


sub cigar_string {
    my ($self, $arg) = @_;
    $self->warn("this is not a setter") if(defined $arg);

    unless(defined $self->{_cigar_string}){ # generate cigar string
        my $cigar_string = $self->generate_cigar_string($self->query_string, $self->hit_string);
        $self->{_cigar_string} = $cigar_string;
    } # end of unless

    return $self->{_cigar_string};
}

=head2 generate_cigar_string

  Name:     generate_cigar_string
  Usage:    my $cigar_string = Bio::Search::HSP::GenericHSP::generate_cigar_string ($qstr, $hstr);
  Function: generate cigar string from a simple sequence of alignment.
  Args:     the string of query and subject
  Return:   cigar string

=cut

sub generate_cigar_string {
    my ($self, $qstr, $hstr) = @_;
    my @qchars = split //, $qstr;
    my @hchars = split //, $hstr;

    unless(scalar(@qchars) == scalar(@hchars)){
        $self->throw("two sequences are not equal in lengths");
    }

    $self->{_count_for_cigar_string} = 0;
    $self->{_state_for_cigar_string} = 'M';

    my $cigar_string = '';
    for(my $i=0; $i <= $#qchars; $i++){
        my $qchar = $qchars[$i];
        my $hchar = $hchars[$i];
        if($qchar ne $GAP_SYMBOL && $hchar ne $GAP_SYMBOL){ # Match
            $cigar_string .= $self->_sub_cigar_string('M');
        }elsif($qchar eq $GAP_SYMBOL){ # Deletion
            $cigar_string .= $self->_sub_cigar_string('D');
        }elsif($hchar eq $GAP_SYMBOL){ # Insertion
            $cigar_string .= $self->_sub_cigar_string('I');
        }else{
            $self->throw("Impossible state that 2 gaps on each seq aligned");
        }
    }
    $cigar_string .= $self->_sub_cigar_string('X'); # not forget the tail.
    return $cigar_string;
}

# an internal method to help generate cigar string

sub _sub_cigar_string {
    my ($self, $new_state) = @_;

    my $sub_cigar_string = '';
    if($self->{_state_for_cigar_string} eq $new_state){
        $self->{_count_for_cigar_string} += 1; # Remain the state and increase the counter
    }else{
        $sub_cigar_string .= $self->{_count_for_cigar_string}
            unless $self->{_count_for_cigar_string} == 1;
        $sub_cigar_string .= $self->{_state_for_cigar_string};
        $self->{_count_for_cigar_string} = 1;
        $self->{_state_for_cigar_string} = $new_state;
    }
    return $sub_cigar_string;
}



# needed before seqfeatures can be made
sub _pre_seq_feature {
    my $self = shift;
    my $algo = $self->{ALGORITHM};

    my ($queryfactor, $hitfactor) = (0,0);
    if( $algo =~ /^(PSI)?T(BLAST|FAST|SW)[NY]/oi ) {
        $hitfactor = 1;
    }
    elsif ($algo =~ /^(FAST|BLAST)(X|Y|XY)/oi || $algo =~ /^P?GENEWISE/oi ) {
        $queryfactor = 1;
    }
    elsif ($algo =~ /^T(BLAST|FAST|SW)(X|Y|XY)/oi || $algo =~ /^(BLAST|FAST|SW)N/oi || $algo =~ /^WABA|AXT|BLAT|BLASTZ|PSL|MEGABLAST|EXONERATE|SW|SMITH\-WATERMAN|SIM4$/){
        $hitfactor = 1;
        $queryfactor = 1;
    }
    elsif ($algo =~ /^RPS-BLAST/) {
        $queryfactor = ($algo =~ /^RPS-BLAST\(BLASTX\)/) ? 1 : 0;
        $hitfactor = 0;
    }
    $self->{_query_factor} = $queryfactor;
    $self->{_hit_factor} = $hitfactor;
}

# make query seq feature
sub _query_seq_feature {
    my $self = shift;
    $self->{_making_qff} = 1;
    my $qs = $self->{QUERY_START};
    my $qe = $self->{QUERY_END};
    unless (defined $self->{_query_factor}) {
        $self->_pre_seq_feature;
    }
    my $queryfactor = $self->{_query_factor};

    unless( defined $qe && defined $qs ) { $self->throw("Did not specify a Query End or Query Begin"); }

    my $strand;
    if ($qe > $qs) {  # normal query: start < end
        if ($queryfactor) {
            $strand = 1;
        }
        else {
            $strand = undef;
        }
    }
    else {
        if ($queryfactor) {
            $strand = -1;
        }
        else {
            $strand = undef;
        }
        ($qs,$qe) = ($qe,$qs);
    }

    # Note: many of these data are not query- and hit-specific.
    # Only start, end, name, length are.
    # We could be more efficient by only storing this info once.
    # steve chervitz --- Sat Apr  5 00:55:07 2003

    my $sim1 = $self->{_sim1} || Bio::SeqFeature::Similarity->new();
    $sim1->start($qs);
    $sim1->end($qe);
    $sim1->significance($self->{EVALUE});
    $sim1->bits($self->{BITS});
    $sim1->score($self->{SCORE});
    $sim1->strand($strand);
    $sim1->seq_id($self->{QUERY_NAME});
    $sim1->seqlength($self->{QUERY_LENGTH});
    $sim1->source_tag($self->{ALGORITHM});
    $sim1->seqdesc($self->{QUERY_DESC});
    
    $self->SUPER::feature1($sim1);

    # to determine frame from something like FASTXY which doesn't
    # report the frame
    my $qframe = $self->{QUERY_FRAME};
    if (defined $strand && ! defined $qframe && $queryfactor) {
        $qframe = ( $self->query->start % 3 ) * $strand;
    }
    elsif (! defined $strand) {
        $qframe = 0;
    }
    $self->{QUERY_FRAME} = $qframe;

    $self->{_created_qff} = 1;
    $self->{_making_qff} = 0;
    $self->_pre_frame;
}

# make subject seq feature
sub _subject_seq_feature {
    my $self = shift;
    $self->{_making_sff} = 1;
    my $hs = $self->{HIT_START};
    my $he = $self->{HIT_END};
    unless (defined $self->{_hit_factor}) {
        $self->_pre_seq_feature;
    }
    my $hitfactor = $self->{_hit_factor};

    unless( defined $he && defined $hs ) { $self->throw("Did not specify a Hit End or Hit Begin"); }

    my $strand;
    if ($he > $hs) { # normal subject
        if ($hitfactor) {
            $strand = 1;
        }
        else {
            $strand = undef;
        }
    }
    else {
        if ($hitfactor) {
            $strand = -1;
        }
        else {
            $strand = undef;
        }
        ($hs,$he) = ( $he,$hs); # reverse subject: start bigger than end
    }

    my $sim2 = $self->{_sim2} || Bio::SeqFeature::Similarity->new();
    $sim2->start($hs);
    $sim2->end($he);
    $sim2->significance($self->{EVALUE});
    $sim2->bits($self->{BITS});
    $sim2->score($self->{SCORE});
    $sim2->strand($strand);
    $sim2->seq_id($self->{HIT_NAME});
    $sim2->seqlength($self->{HIT_LENGTH});
    $sim2->source_tag($self->{ALGORITHM});
    $sim2->seqdesc($self->{HIT_DESC});
    $self->SUPER::feature2($sim2);

    my $hframe = $self->{HIT_FRAME};
    if (defined $strand && ! defined $hframe && $hitfactor) {
        $hframe = ( $hs % 3 ) * $strand;
    }
    elsif (! defined $strand) {
        $hframe = 0;
    }
    $self->{HIT_FRAME} = $hframe;

    $self->{_created_sff} = 1;
    $self->{_making_sff} = 0;
    $self->_pre_frame;
}

# know the frame following seq feature creation
sub _pre_frame {
    my $self = shift;
    $self->{_created_qff} || $self->_query_seq_feature;
    $self->{_created_sff} || $self->_subject_seq_feature;
    $self->{_did_preframe} = 1;
    $self->frame;
}

# before calling the num_* methods
sub _pre_similar_stats {
    my $self = shift;
    my $identical = $self->{IDENTICAL};
    my $conserved = $self->{CONSERVED};
    my $percent_id = $self->{PERCENT_IDENTITY};

    if (! defined $identical) {
        if (! defined $percent_id) {
            $self->warn("Did not defined the number of identical matches or overall percent identity in the HSP assuming 0");
            $identical = 0;
        }
        else {
            $identical = int($percent_id * $self->{HSP_LENGTH});
        }
    }

    if (! defined $conserved) {
        $self->warn("Did not defined the number of conserved matches in the HSP assuming conserved == identical ($identical)")
            if( $self->{ALGORITHM} !~ /^((FAST|BLAST)N)|EXONERATE|SIM4|AXT|PSL|BLAT|BLASTZ|WABA/oi);
        $conserved = $identical;
    }
    $self->{IDENTICAL} = $identical;
    $self->{CONSERVED} = $conserved;
    $self->{_did_presimilar} = 1;
}

# before calling the frac_* methods
sub _pre_frac {
    my $self = shift;
    my $hsp_len = $self->{HSP_LENGTH};
    my $hit_len = $self->{HIT_LENGTH};
    my $query_len = $self->{QUERY_LENGTH};

    my $identical = $self->num_identical;
    my $conserved = $self->num_conserved;

    $self->{_did_prefrac} = 1;
    my $logical;
    if( $hsp_len ) {
        $self->length('total', $hsp_len);
        $logical = $self->_logical_length('total');
        $self->frac_identical( 'total', $identical / $hsp_len);
        $self->frac_conserved( 'total', $conserved / $hsp_len);
    }
    if( $hit_len ) {
        $logical = $self->_logical_length('hit');
        $self->frac_identical( 'hit', $identical / $logical);
        $self->frac_conserved( 'hit', $conserved / $logical);
    }
    if( $query_len ) {
        $logical = $self->_logical_length('query');
        $self->frac_identical( 'query', $identical / $logical) ;
        $self->frac_conserved( 'query', $conserved / $logical);
    }
}

# before calling gaps()
sub _pre_gaps {
    my $self = shift;
    my $query_gaps = $self->{QUERY_GAPS};
    my $query_seq = $self->{QUERY_SEQ};
    my $hit_gaps = $self->{HIT_GAPS};
    my $hit_seq = $self->{HIT_SEQ};
    my $gaps = $self->{HSP_GAPS};

    $self->{_did_pregaps} = 1; # well, we're in the process; avoid recursion
    if( defined $query_gaps ) {
        $self->gaps('query', $query_gaps);
    } elsif( defined $query_seq ) {
        $self->gaps('query', scalar ( $query_seq =~ tr/\-//));
    }
    if( defined $hit_gaps ) {
        $self->gaps('hit', $hit_gaps);
    } elsif( defined $hit_seq ) {
        $self->gaps('hit', scalar ( $hit_seq =~ tr/\-//));
    }
    if( ! defined $gaps ) {
        $gaps = $self->gaps("query") + $self->gaps("hit");
    }
    $self->gaps('total', $gaps);
}

# before percent_identity
sub _pre_pi {
    my $self = shift;
    $self->{_did_prepi} = 1;
    $self->percent_identity($self->{PERCENT_IDENTITY} || $self->frac_identical('total')*100) if( $self->{HSP_LENGTH} > 0 );
}

1;
