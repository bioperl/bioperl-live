#-----------------------------------------------------------------
#
# BioPerl module Bio::Search::Hit::HitI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Steve Chervitz <sac@bioperl.org>
#
# Originally created by Aaron Mackey <amackey@virginia.edu>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Hit::HitI - Interface for a hit in a similarity search result

=head1 SYNOPSIS

# Bio::Search::Hit::HitI objects should not be instantiated since this
# module defines a pure interface.

# Given an object that implements the Bio::Search::Hit::HitI  interface,
# you can do the following things with it:

    # Get a HitI object from a SearchIO stream:
    use Bio::SeachIO;
    my $searchio = Bio::SearchIO->new(-format => 'blast', -file => 'result.bls');
    my $result = $searchio->next_result;
    my $hit    = $result->next_hit;

    $hit_name = $hit->name();

    $desc = $hit->description();

    $len = $hit->length

    $alg = $hit->algorithm();

    $score = $hit->raw_score();

    $significance = $hit->significance();

    $rank = $hit->rank(); # the Nth hit for a specific query

    while( $hsp = $obj->next_hsp()) { ... } # process in iterator fashion

    for my $hsp ( $obj->hsps()()) { ... } # process in list fashion

=head1 DESCRIPTION

    Bio::Search::Hit::* objects are data structures that contain information
about specific hits obtained during a library search.  Some information will
be algorithm-specific, but others will be generally defined.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Aaron Mackey, Steve Chervitz

Email amackey@virginia.edu  (original author)
Email sac@bioperl.org

=head1 COPYRIGHT

Copyright (c) 1999-2001 Aaron Mackey, Steve Chervitz. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Search::Hit::HitI;


use strict;

use base qw(Bio::Root::RootI);


=head2 name

 Title   : name
 Usage   : $hit_name = $hit->name();
 Function: returns the name of the Hit sequence
 Returns : a scalar string
 Args    : none

The B<name> of a hit is unique within a Result or within an Iteration.

=cut

sub name {
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}

=head2 description

 Title   : description
 Usage   : $desc = $hit->description();
 Function: Retrieve the description for the hit
 Returns : a scalar string
 Args    : none

=cut

sub description {
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}


=head2 accession

 Title   : accession
 Usage   : $acc = $hit->accession();
 Function: Retrieve the accession (if available) for the hit
 Returns : a scalar string (empty string if not set)
 Args    : none

=cut

sub accession {
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}

=head2 locus

 Title   : locus
 Usage   : $acc = $hit->locus();
 Function: Retrieve the locus(if available) for the hit
 Returns : a scalar string (empty string if not set)
 Args    : none

=cut

sub locus {
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}

=head2 length

 Title   : length
 Usage   : my $len = $hit->length
 Function: Returns the length of the hit 
 Returns : integer
 Args    : none

=cut

sub length {
   my ($self,@args) = @_;
   $self->throw_not_implemented;
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
 Args    : none

=cut

sub algorithm {
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}

=head2 raw_score

 Title   : raw_score
 Usage   : $score = $hit->raw_score();
 Function: Gets the "raw score" generated by the algorithm.  What
           this score is exactly will vary from algorithm to algorithm,
           returning undef if unavailable.
 Returns : a scalar value
 Args    : none

=cut

sub raw_score {
    $_[0]->throw_not_implemented;
}

=head2 score

Equivalent to L<raw_score()|raw_score>

=cut

sub score { shift->raw_score(@_); }

=head2 significance

 Title   : significance
 Usage   : $significance = $hit->significance();
 Function: Used to obtain the E or P value of a hit, i.e. the probability that
           this particular hit was obtained purely by random chance.  If
           information is not available (nor calculatable from other
           information sources), return undef.
 Returns : a scalar value or undef if unavailable
 Args    : none

=cut

sub significance {
    $_[0]->throw_not_implemented;
}

=head2 bits

 Usage     : $hit_object->bits();
 Purpose   : Gets the bit score of the best HSP for the current hit.
 Example   : $bits = $hit_object->bits();
 Returns   : Integer or double for FASTA reports
 Argument  : n/a
 Comments  : For BLAST1, the non-bit score is listed in the summary line.

See Also   : L<score()|score>

=cut

#---------
sub bits { 
#---------
    $_[0]->throw_not_implemented();
}

=head2 next_hsp

 Title    : next_hsp
 Usage    : while( $hsp = $obj->next_hsp()) { ... }
 Function : Returns the next available High Scoring Pair
 Example  : 
 Returns  : L<Bio::Search::HSP::HSPI> object or null if finished
 Args     : none

=cut

sub next_hsp {
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}


=head2 hsps

 Usage     : $hit_object->hsps();
 Purpose   : Get a list containing all HSP objects.
           : Get the numbers of HSPs for the current hit.
 Example   : @hsps = $hit_object->hsps();
           : $num  = $hit_object->hsps();  # alternatively, use num_hsps()
 Returns   : Array context : list of L<Bio::Search::HSP::BlastHSP> objects.
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

    $self->throw_not_implemented();
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
#-------------
    shift->throw_not_implemented();
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

See Also   : L<Bio::Search::HSP::HSPI::seq_inds()|Bio::Search::HSP::HSPI>

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

=head2 rewind

 Title   : rewind
 Usage   : $hit->rewind;
 Function: Allow one to reset the HSP iterator to the beginning
           if possible
 Returns : none
 Args    : none

=cut

sub rewind{
   my ($self) = @_;
   $self->throw_not_implemented();
}


=head2 overlap

 Usage     : $hit_object->overlap( [integer] );
 Purpose   : Gets/Sets the allowable amount overlap between different HSP sequences.
 Example   : $hit_object->overlap(5);
           : $overlap = $hit_object->overlap;
 Returns   : Integer.
 Argument  : integer.
 Throws    : n/a
 Status    : Experimental
 Comments  : Any two HSPs whose sequences overlap by less than or equal
           : to the overlap() number of resides will be considered separate HSPs
           : and will not get tiled by L<Bio::Search::BlastUtils::_adjust_contigs()>.

See Also   : L<Bio::Search::BlastUtils::_adjust_contigs()|Bio::Search::BlastUtils>, L<BUGS | BUGS>

=cut

#-------------
sub overlap { shift->throw_not_implemented }


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
sub n { shift->throw_not_implemented }

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
sub p { shift->throw_not_implemented() }

=head2 hsp

 Usage     : $hit_object->hsp( [string] );
 Purpose   : Get a single HSPI object for the present HitI object.
 Example   : $hspObj  = $hit_object->hsp;  # same as 'best'
           : $hspObj  = $hit_object->hsp('best');
           : $hspObj  = $hit_object->hsp('worst');
 Returns   : Object reference for a L<Bio::Search::HSP::HSPI> object.
 Argument  : String (or no argument).
           :   No argument (default) = highest scoring HSP (same as 'best').
           :   'best' or 'first' = highest scoring HSP.
           :   'worst' or 'last' = lowest scoring HSP.
 Throws    : Exception if the HSPs have not been collected.
           : Exception if an unrecognized argument is used.

See Also   : L<hsps()|hsps>, L<num_hsps>()

=cut

#----------
sub hsp { shift->throw_not_implemented }

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
sub logical_length { shift->throw_not_implemented() }


=head2 rank

 Title   : rank
 Usage   : $obj->rank($newval)
 Function: Get/Set the rank of this Hit in the Query search list
           i.e. this is the Nth hit for a specific query
 Returns : value of rank
 Args    : newvalue (optional)


=cut

sub rank{
   my ($self,$value) = @_;
   $self->throw_not_implemented();
}

=head2 each_accession_number

 Title   : each_accession_number
 Usage   : $obj->each_accession_number
 Function: Get each accession number listed in the description of the hit.
           If there are no alternatives, then only the primary accession will 
           be given
 Returns : list of all accession numbers in the description
 Args    : none


=cut

sub each_accession_number{
   my ($self,$value) = @_;
   $self->throw_not_implemented();
}


=head2 tiled_hsps

 Usage     : $hit_object->tiled_hsps( [integer] );
 Purpose   : Gets/Sets an indicator for whether or not the HSPs in this Hit 
           : have been tiled.
           : Methods that rely on HSPs being tiled should check this
           : and then call SearchUtils::tile_hsps() if not.
 Example   : $hit_object->tiled_hsps(1);
           : if( $hit_object->tiled_hsps ) { # do something }
 Returns   : Boolean (1 or 0) 
 Argument  : integer (optional)
 Throws    : n/a

=cut

sub tiled_hsps { shift->throw_not_implemented }


=head2 strand

 Usage     : $sbjct->strand( [seq_type] );
 Purpose   : Gets the strand(s) for the query, sbjct, or both sequences
           : in the best HSP of the BlastHit object after HSP tiling.
           : Only valid for BLASTN, TBLASTX, BLASTX-query, TBLASTN-hit.
 Example   : $qstrand = $sbjct->strand('query');
           : $sstrand = $sbjct->strand('hit');
           : ($qstrand, $sstrand) = $sbjct->strand();
 Returns   : scalar context: integer '1', '-1', or '0'
           : array context without args: list of two strings (queryStrand, sbjctStrand)
           : Array context can be "induced" by providing an argument of 'list' or 'array'.
 Argument  : In scalar context: seq_type = 'query' or 'hit' or 'sbjct' (default = 'query')
             ('sbjct' is synonymous with 'hit')
 Throws    : n/a
 Comments  : This method requires that all HSPs be tiled. If they have not
           : already been tiled, they will be tiled first automatically..
           : If you don't want the tiled data, iterate through each HSP
           : calling strand() on each (use hsps() to get all HSPs).
           :
           : Formerly (prior to 10/21/02), this method would return the
           : string "-1/1" for hits with HSPs on both strands.
           : However, now that strand and frame is properly being accounted
           : for during HSP tiling, it makes more sense for strand()
           : to return the strand data for the best HSP after tiling.
           :
           : If you really want to know about hits on opposite strands,
           : you should be iterating through the HSPs using methods on the
           : HSP objects.
           :
           : A possible use case where knowing whether a hit has HSPs 
           : on both strands would be when filtering via SearchIO for hits with 
           : this property. However, in this case it would be better to have a
           : dedicated method such as $hit->hsps_on_both_strands(). Similarly
           : for frame. This could be provided if there is interest.

See Also   : L<Bio::Search::HSP::HSPI::strand>()

=cut

#---------'
sub strand { shift->throw_not_implemented }


=head2 frame

 Usage     : $hit_object->frame();
 Purpose   : Gets the reading frame for the best HSP after HSP tiling.
           : This is only valid for BLASTX and TBLASTN/X type reports.
 Example   : $frame = $hit_object->frame();
 Returns   : Integer (-2 .. +2)
 Argument  : n/a
 Throws    : Exception if HSPs have not been set.
 Comments  : This method requires that all HSPs be tiled. If they have not
           : already been tiled, they will be tiled first automatically..
           : If you don't want the tiled data, iterate through each HSP
           : calling frame() on each (use hsps() to get all HSPs).

See Also   : L<hsps()|hsps>

=cut

#---------'
sub frame { shift->throw_not_implemented }


=head2 matches

 Usage     : $hit_object->matches( [class] );
 Purpose   : Get the total number of identical or conserved matches 
           : (or both) across all HSPs.
           : (Note: 'conservative' matches are indicated as 'positives' 
           :         in BLAST reports.)
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
 Comments  : This method requires that all HSPs be tiled. If there is more than one
           : HSP and they have not already been tiled, they will be tiled first automatically..
           :
           : If you need data for each HSP, use hsps() and then interate
           : through the HSP objects.
           : Does not rely on wantarray to return a list. Only checks for
           : the presence of an argument (no arg = return list).

See Also   : L<Bio::Search::HSP::GenericHSP::matches()|Bio::Search::HSP::GenericHSP>, L<hsps()|hsps>

=cut

sub matches { shift->throw_not_implemented }


# aliasing for Steve's method names
sub hit_description { shift->description(@_) }
# aliasing for Steve's method names
sub hit_length { shift->length(@_) }


# sort method for HSPs

=head2 sort_hits

 Title		: sort_hsps
 Usage		: $result->sort_hsps(\&sort_function)
 Function	: Sorts the available HSP objects by a user-supplied function. Defaults to sort
                  by descending score.
 Returns	: n/a
 Args		: A coderef for the sort function.  See the documentation on the Perl sort()
                  function for guidelines on writing sort functions.  
 Note		: To access the special variables $a and $b used by the Perl sort() function 
                  the user function must access Bio::Search::Hit::HitI namespace. 
                  For example, use :
                  $hit->sort_hsps( sub{$Bio::Search::Result::HitI::a->length <=> 
					  $Bio::Search::Result::HitI::b->length});
                   NOT $hit->sort_hsps($a->length <=> $b->length);

=cut

sub sort_hsps {shift->throw_not_implemented }

=head2 _default sort_hsps

  Title	: _default_sort_hsps
  Usage	: Do not call directly.
  Function : Sort hsps in ascending order by evalue
  Args	: None
  Returns: 1 on success
  Note	: Used by $hit->sort_hsps()

=cut

sub _default_sort_hsps {
    $Bio::Search::Hit::HitI::a->evalue <=> 
	    $Bio::Search::Hit::HitI::a->evalue;
}

1;




