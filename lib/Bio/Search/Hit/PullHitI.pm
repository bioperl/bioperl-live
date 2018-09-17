#
# BioPerl module Bio::Search::Hit::PullHitI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Hit::PullHitI - Bio::Search::Hit::HitI interface for pull parsers.

=head1 SYNOPSIS

	# This is an interface and cannot be instantiated

    # typically one gets HitI objects from a SearchIO stream via a ResultI
    use Bio::SearchIO;
    my $parser = Bio::SearchIO->new(-format => 'hmmer_pull',
                                   -file => 't/data/hmmpfam.out');

    my $result = $parser->next_result;
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

This object handles the hit data from a database sequence search.

PullHitI is for fast implementations that only do parsing work on the hit
data when you actually request information by calling one of the HitI
methods.

Many methods of HitI are implemented in a way suitable for inheriting classes
that use Bio::PullParserI. It only really makes sense for PullHit modules to be
created by (and have as a -parent) PullResult modules.

In addition to the usual -chunk and -parent, -hit_data is all you should supply
when making a PullHit object. This will store that data and make it accessible
via _raw_hit_data, which you can access in your subclass. It would be best to
simply provide the data as the input -chunk instead, if the raw data is large
enough.

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

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 COPYRIGHT

Copyright (c) 2006 Sendu Bala. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Search::Hit::PullHitI;

use Bio::Search::SearchUtils;

use strict;

use base qw(Bio::PullParserI Bio::Search::Hit::HitI);

=head2 _setup

 Title   : _setup
 Usage   : $self->_setup(@args)
 Function: Implementers should call this to setup common fields and deal with
           common arguments to new().
 Returns : n/a
 Args    : @args received in new().

=cut

sub _setup {
    my ($self, @args) = @_;
    
    # fields most subclasses probably will want
    $self->_fields( { ( next_hsp => undef,
                        num_hsps => undef,
                        hsps => undef,
                        hit_start => undef,
                        query_start => undef,
                        hit_end => undef,
                        query_end => undef,
                        length => undef,
						name => undef ,
						accession => undef ) } );
    
    my ($parent, $chunk, $hit_data) = $self->_rearrange([qw(PARENT
                                                            CHUNK
                                                            HIT_DATA)], @args);
    $self->throw("Need -parent or -chunk to be defined") unless $parent || $chunk;
    
	$self->parent($parent) if $parent;
    
    if ($chunk) {
        my ($io, $start, $end) = (undef, 0, undef);
        if (ref($chunk) eq 'ARRAY') {
            ($io, $start, $end) = @{$chunk};
        }
        else {
            $io = $chunk;
        }
        $self->chunk($io, -start => $start, -end => $end);
    }
    
    $self->_raw_hit_data($hit_data) if $hit_data;
}

sub _raw_hit_data {
	my $self = shift;
	if (@_) {
		$self->{_raw_hit_data} = shift;
	}
	return $self->{_raw_hit_data};
}

#
# Some of these methods are written explitely to avoid HitI throwing not
# implemented; if it didn't do that then PullParserI AUTOLOAD would have
# cought them.
#

=head2 name

 Title   : name
 Usage   : $hit_name = $hit->name();
 Function: returns the name of the Hit sequence
 Returns : a scalar string
 Args    : none

The B<name> of a hit is unique within a Result or within an Iteration.

=cut

sub name {
    return shift->get_field('name');
}

=head2 description

 Title   : description
 Usage   : $desc = $hit->description();
 Function: Retrieve the description for the hit
 Returns : a scalar string
 Args    : none

=cut

sub description {
    return shift->get_field('description');
}

=head2 accession

 Title   : accession
 Usage   : $acc = $hit->accession();
 Function: Retrieve the accession (if available) for the hit
 Returns : a scalar string (empty string if not set)
 Args    : none

=cut

sub accession {
    return shift->get_field('accession');
}

=head2 locus

 Title   : locus
 Usage   : $acc = $hit->locus();
 Function: Retrieve the locus(if available) for the hit
 Returns : a scalar string (empty string if not set)
 Args    : none

=cut

sub locus {
    return shift->get_field('locus');
}

=head2 length

 Title   : length
 Usage   : my $len = $hit->length
 Function: Returns the length of the hit 
 Returns : integer
 Args    : none

=cut

sub length {
   return shift->get_field('length');
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
    return shift->get_field('algorithm');
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
    return shift->get_field('score');
}

=head2 score

Equivalent to L<raw_score()|raw_score>

=cut

sub score {
    return shift->get_field('score');
}

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
    return shift->get_field('significance');
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

sub bits {
    return shift->get_field('bits');
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
    return shift->get_field('next_hsp');
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

sub hsps {
    return shift->get_field('hsps');
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

sub num_hsps {
    return shift->get_field('num_hsps');
}

#
# HitI/ GenericHit methods that are unrelated to simply parsing information
# directly out of a file, but need more complex calculation; mostly not
# implemented here.
#

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
           : [1] class = 'identical' or 'conserved' or 'nomatch' or 'gap'
		   :         (default = 'identical')
           :         (can be shortened to 'id' or 'cons')
		   :         Note that 'conserved' includes identical unless you use
		   :         'conserved-not-identical'
           : [2] collapse = boolean, if non-zero, consecutive positions are
           :             merged using a range notation, e.g.,
           :             "1 2 3 4 5 7 9 10 11" collapses to "1-5 7 9-11". This
           :             is useful for  consolidating long lists. Default = no
           :             collapse.
 Throws    : n/a.

See Also   : L<Bio::Search::HSP::HSPI::seq_inds()|Bio::Search::HSP::HSPI>

=cut

sub seq_inds {
    my ($self, $seqType, $class, $collapse) = @_;
    
    $seqType  ||= 'query';
    $class ||= 'identical';
    $collapse ||= 0;
    
    $seqType = 'hit' if $seqType eq 'sbjct';
    
	my $storage_name = '_seq_inds_'.$seqType.'_'.$class;
	unless (defined $self->{$storage_name}) {
		my @inds;    
		foreach my $hsp ($self->hsps) {
			# This will merge data for all HSPs together.
			push @inds, $hsp->seq_inds($seqType, $class);
		}
		
		# Need to remove duplicates and sort the merged positions, unless gaps.
		if (@inds && $class ne 'gap') {
			my %tmp = map { $_, 1 } @inds;
			@inds = sort {$a <=> $b} keys %tmp;
		}
		
		$self->{$storage_name} = \@inds;
	}
	
	my @inds = @{$self->{$storage_name}};
    $collapse ? &Bio::Search::SearchUtils::collapse_nums(@inds) : @inds;
}

=head2 rewind

 Title   : rewind
 Usage   : $hit->rewind;
 Function: Allow one to reset the HSP iterator to the beginning if possible
 Returns : none
 Args    : none

=cut

sub rewind {
    shift->throw_not_implemented();
}

=head2 overlap

 Usage     : $hit_object->overlap( [integer] );
 Purpose   : Gets/Sets the allowable amount overlap between different HSP
             sequences.
 Example   : $hit_object->overlap(5);
           : $overlap = $hit_object->overlap;
 Returns   : Integer.
 Argument  : integer.
 Throws    : n/a
 Status    : Deprecated
 Comments  : This value isn't used for anything

=cut

sub overlap {
    my $self = shift;
    if (@_) { $self->{_overlap} = shift }
    return $self->{_overlap} || 0;
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

sub n {
    return shift->get_field('num_hsps');
}

=head2 p

 Usage     : $hit_object->p( [format] );
 Purpose   : Get the P-value for the best HSP of the given BLAST hit.
           : (Note that P-values are not provided with NCBI Blast2 reports).
 Example   : $p =  $sbjct->p;
           : $p =  $sbjct->p('exp');  # get exponent only.
           : ($num, $exp) =  $sbjct->p('parts'); # split sci notation into parts
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

See Also   : L<expect()|expect>, L<signif()|signif>,
             L<Bio::Search::BlastUtils::get_exponent()|Bio::Search::BlastUtils>

=cut

sub p {
    shift->throw_not_implemented;
}

=head2 hsp

 Usage     : $hit_object->hsp( [string] );
 Purpose   : Get a single HSPI object for the present HitI object.
 Example   : $hspObj  = $hit_object->hsp;  # same as 'best'
           : $hspObj  = $hit_object->hsp('best');
           : $hspObj  = $hit_object->hsp('worst');
 Returns   : Object reference for a L<Bio::Search::HSP::HSPI> object.
 Argument  : String (or no argument).
           :   No argument (default) = highest scoring HSP (same as 'best').
           :   'best'  = highest scoring HSP.
           :   'worst' = lowest scoring HSP.
 Throws    : Exception if an unrecognized argument is used.

See Also   : L<hsps()|hsps>, L<num_hsps>()

=cut

sub hsp {
    shift->throw_not_implemented;
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

See Also   : L<length()|length>, L<frac_aligned_query()|frac_aligned_query>,
             L<frac_aligned_hit()|frac_aligned_hit>

=cut

sub logical_length {
    my ($self, $type) = @_;
    $type ||= 'query';
    $type = lc($type);
	$type = 'hit' if $type eq 'sbjct';
    if ($type eq 'query') {
        return $self->get_field('query_length');
    }
    elsif ($type eq 'hit') {
        return $self->get_field('length');
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

sub rank {
    return shift->get_field('rank');
}

=head2 each_accession_number

 Title   : each_accession_number
 Usage   : $obj->each_accession_number
 Function: Get each accession number listed in the description of the hit.
           If there are no alternatives, then only the primary accession will 
           be given (if there is one).
 Returns : list of all accession numbers in the description
 Args    : none

=cut

sub each_accession_number {
    my $self = shift;
    my $accession = $self->get_field('accession') if $self->has_field('accession');
    my $desc = $self->get_field('description') if $self->has_field('description');
    return unless $accession || $desc;
    
    my @accnums;
    push (@accnums, $accession) if $accession;
    
    if (defined $desc) { 
        while ($desc =~ /(\b\S+\|\S*\|\S*\s?)/g) {
            my $id = $1;
            my $acc;
            if ($id =~ /(?:gb|emb|dbj|sp|pdb|bbs|ref|tp[gde])\|(.*)\|(?:.*)/) {
                ($acc) = split /\./, $1; 
            }
            elsif ($id =~ /(?:pir|prf|pat|gnl)\|(?:.*)\|(.*)/) {
                ($acc) = split /\./, $1;  
            }
            elsif ($id =~ /(?:gim|gi|bbm|bbs|lcl)\|(?:\d*)/) {
                $acc = $id;
            }
            elsif ($id =~ /(?:oth)\|(.*)\|(?:.*)\|(?:.*)/ ) { # discontinued...
                $acc = $1;
            }
            else {
                $acc = $id;
            }
            push(@accnums, $acc);
        }
    }
    return @accnums;
}

=head2 tiled_hsps

 Usage     : $hit_object->tiled_hsps( [integer] );
 Purpose   : Gets/Sets an indicator for whether or not the HSPs in this Hit 
           : have been tiled.
 Example   : $hit_object->tiled_hsps(1);
           : if( $hit_object->tiled_hsps ) { # do something }
 Returns   : Boolean (1 or 0) 
 Argument  : integer (optional)
 Throws    : n/a
 Status    : Deprecated
 Notes     : This value is not used for anything

=cut

sub tiled_hsps {
    my $self = shift;
    if (@_) { $self->{_hsps_are_tiled} = shift }
    return $self->{_hsps_are_tiled} || 0;
}

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

sub strand {
    shift->throw_not_implemented;
}

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

sub frame {
    shift->throw_not_implemented;
}

=head2 length_aln

 Usage     : $hit_object->length_aln( [seq_type] );
 Purpose   : Get the total length of the aligned region for query or sbjct seq.
           : This number will include all HSPs, and excludes gaps.
 Example   : $len    = $hit_object->length_aln(); # default = query
           : $lenAln = $hit_object->length_aln('query');
 Returns   : Integer 
 Argument  : seq_Type = 'query' or 'hit' or 'sbjct' (Default = 'query')
             ('sbjct' is synonymous with 'hit')
 Throws    : Exception if the argument is not recognized.
 Comments  : This method will report the logical length of the alignment,
           : meaning that for TBLAST[NX] reports, the length is reported
           : using amino acid coordinate space (i.e., nucleotides / 3).

See Also   : L<length()|length>, L<frac_aligned_query()|frac_aligned_query>,
             L<frac_aligned_hit()|frac_aligned_hit>, L<gaps()|gaps>,
             L<Bio::Search::SearchUtils::tile_hsps()|Bio::Search::SearchUtils>,
             L<Bio::Search::HSP::BlastHSP::length()|Bio::Search::HSP::BlastHSP>

=cut

sub length_aln {
    my ($self, $seqType) = @_;
    $seqType ||= 'query';
    $seqType = 'hit' if $seqType eq 'sbjct';
    
	my %non_gaps = map { $_, 1 } $self->seq_inds($seqType, 'conserved'),
							     $self->seq_inds($seqType, 'no_match');
	return scalar(keys %non_gaps);
}

=head2 gaps

 Usage     : $hit_object->gaps( [seq_type] );
 Purpose   : Get the number of gaps in the aligned query, hit, or both sequences.
           : Data is summed across all HSPs.
 Example   : $qgaps = $hit_object->gaps('query');
           : $hgaps = $hit_object->gaps('hit');
           : $tgaps = $hit_object->gaps();    # default = total (query + hit)
 Returns   : scalar context: integer
           : array context without args: two-element list of integers  
           :    (queryGaps, hitGaps)
           : Array context can be forced by providing an argument of 'list' or
		   : 'array'.
           :
           : CAUTION: Calling this method within printf or sprintf is arrray
		   : context.
           : So this function may not give you what you expect. For example:
           :          printf "Total gaps: %d", $hit->gaps();
           : Actually returns a two-element array, so what gets printed 
           : is the number of gaps in the query, not the total
           :
 Argument  : seq_type: 'query' | 'hit' or 'sbjct' | 'total' | 'list'
           : (default = 'total') ('sbjct' is synonymous with 'hit')
 Comments  : If you need data for each HSP, use hsps() and then interate
           : through each HSP object.

=cut

sub gaps {
    my ($self, $seqType) = @_;
    
    $seqType ||= (wantarray ? 'list' : 'total');
    $seqType = 'hit' if $seqType eq 'sbjct';
    
    if ($seqType =~ /list|array/i) {
        return (scalar($self->seq_inds('query', 'gap')), scalar($self->seq_inds('hit', 'gap')));
    }
    elsif ($seqType eq 'total') {
        return (scalar($self->seq_inds('query', 'gap')) + scalar($self->seq_inds('hit', 'gap'))) || 0;
    }
    else {
        return scalar($self->seq_inds($seqType, 'gap')) || 0;
    }
}

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
 Argument  : [0] class = 'id' | 'cons' OR none.
           : [1] seq_type  = 'query' or 'hit' or 'sbjct'  (default = 'query')
           :                 ('sbjct' is synonymous with 'hit')
           : If no argument is provided, both identical and conservative 
           : numbers are returned in a two element list.
           : (Other terms can be used to refer to the conservative
           :  matches, e.g., 'positive'. All that is checked is whether or
           :  not the supplied string starts with 'id'. If not, the 
           : conservative matches are returned.)

=cut

sub matches {
    my ($self, $class, $seqType) = @_;
    
	# no query/hit choice? The answer differs depending on sequence, since
	# hsps could overlap on one sequence and not the other. Added an option,
	# but otherwise will assume 'hit'
	$seqType ||= 'hit';
	$seqType = 'hit' if $seqType eq 'sbjct';
	
	unless (exists $self->{_id_matches}) {
		$self->{_id_matches}->{hit} = scalar($self->seq_inds('hit', 'identical'));
		$self->{_id_matches}->{query} = scalar($self->seq_inds('query', 'identical'));
	}
	unless (exists $self->{_con_matches}) {
		foreach my $type ('hit', 'query') {
			# 'conserved-not-identical' can give us 'identical' matches if hsps
			# overlapped so have to get the difference
			my %identicals = map { $_ => 1 } $self->seq_inds($type, 'identical');
			my @conserved = $self->seq_inds($type, 'conserved-not-identical');
			
			my $real_conserved;
			foreach (@conserved) {
				unless (exists $identicals{$_}) {
					$real_conserved++;
				}
			}
			$self->{_con_matches}->{$type} = $real_conserved;
		}
	}
	
	
    unless ($class) {
        return ($self->{_id_matches}->{$seqType}, $self->{_con_matches}->{$seqType});
    }
    else {
		if ($class =~ /^id/i) { 
            return $self->{_id_matches}->{$seqType};
        }
        else {
            return $self->{_con_matches}->{$seqType};
        }
    }
    return;
}

=head2 start

 Usage     : $sbjct->start( [seq_type] );
 Purpose   : Gets the start coordinate for the query, sbjct, or both sequences
           : in the object. If there is more than one HSP, the lowest start
           : value of all HSPs is returned.
 Example   : $qbeg = $sbjct->start('query');
           : $sbeg = $sbjct->start('hit');
           : ($qbeg, $sbeg) = $sbjct->start();
 Returns   : scalar context: integer 
           : array context without args: list of two integers (queryStart,
           : sbjctStart)
           : Array context can be "induced" by providing an argument of 'list'
           : or 'array'.
 Argument  : 'query' or 'hit' or 'sbjct' (default = 'query') ('sbjct' is
             synonymous with 'hit')

=cut

sub start {
    my ($self, $seqType) = @_;
    
    unless ($self->get_field('num_hsps')) {
        Bio::Search::SearchUtils::_warn_about_no_hsps($self);
        return;
    }
    
    $seqType ||= (wantarray ? 'list' : 'query');
    $seqType = 'hit' if $seqType eq 'sbjct';
    
    if ($seqType =~ /list|array/i) {
	    return ($self->get_field('query_start'), $self->get_field('hit_start'));
	}
    elsif ($seqType eq 'hit') {
        return $self->get_field('hit_start');
	}
    elsif ($seqType eq 'query') {
        return $self->get_field('query_start');
    }
    else {
        $self->throw("Unknown sequence type '$seqType'");
    }
}

=head2 end

 Usage     : $sbjct->end( [seq_type] );
 Purpose   : Gets the end coordinate for the query, sbjct, or both sequences
           : in the object. If there is more than one HSP, the largest end
           : value of all HSPs is returned.
 Example   : $qend = $sbjct->end('query');
           : $send = $sbjct->end('hit');
           : ($qend, $send) = $sbjct->end();
 Returns   : scalar context: integer
           : array context without args: list of two integers 
           : (queryEnd, sbjctEnd)
           : Array context can be "induced" by providing an argument 
           : of 'list' or 'array'.
 Argument  : 'query' or 'hit' or 'sbjct' (default = 'query') ('sbjct' is
             synonymous with 'hit')

=cut

sub end {
    my ($self, $seqType) = @_;
    
    unless ($self->get_field('num_hsps')) {
        Bio::Search::SearchUtils::_warn_about_no_hsps($self);
        return;
    }
    
    $seqType ||= (wantarray ? 'list' : 'query');
    $seqType = 'hit' if $seqType eq 'sbjct';
    
    if ($seqType =~ /list|array/i) {
	    return ($self->get_field('query_end'), $self->get_field('hit_end'));
	}
    elsif ($seqType eq 'hit') {
        return $self->get_field('hit_end');
	}
    elsif ($seqType eq 'query') {
        return $self->get_field('query_end');
    }
    else {
        $self->throw("Unknown sequence type '$seqType'");
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

sub range {
    my ($self, $seqType) = @_;
    $seqType ||= 'query';
    $seqType = 'hit' if $seqType eq 'sbjct';
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

=cut

sub frac_identical {
    my ($self, $seqType) = @_;
    $seqType ||= 'query';
    $seqType = lc($seqType);
    $seqType = 'hit' if $seqType eq 'sbjct';
    
    my $ident = $self->matches('id', $seqType);
    my $total = $self->length_aln($seqType);
    my $ratio = $ident / $total;
    my $ratio_rounded = sprintf( "%.3f", $ratio);
    
    # Round down if normal rounding yields 1 (just like blast)
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

=cut

sub frac_conserved {
    my ($self, $seqType) = @_;
    $seqType ||= 'query';
    $seqType = lc($seqType);
    $seqType = 'hit' if $seqType eq 'sbjct';
    
    my $consv = $self->matches('cons');
    my $total = $self->length_aln($seqType);
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
 Argument  : none
 Comments  : If you need data for each HSP, use hsps() and then interate
           : through the HSP objects.

=cut

sub frac_aligned_query {
    my $self = shift;
    return sprintf("%.2f", $self->length_aln('query') / $self->logical_length('query'));
}

=head2 frac_aligned_hit

 Usage     : $hit_object->frac_aligned_hit();
 Purpose   : Get the fraction of the hit (sbjct) sequence which has been aligned
           : across all HSPs (not including intervals between non-overlapping
           : HSPs).
 Example   : $frac_alnq = $hit_object->frac_aligned_hit();
 Returns   : Float (2-decimal precision, e.g., 0.75).
 Argument  : none
 Comments  : If you need data for each HSP, use hsps() and then interate
           : through the HSP objects.

=cut

sub frac_aligned_hit {
    my $self = shift;
    return sprintf( "%.2f", $self->length_aln('sbjct') / $self->logical_length('sbjct'));
}

=head2 num_unaligned_hit

 Usage     : $hit_object->num_unaligned_hit();
 Purpose   : Get the number of the unaligned residues in the hit sequence.
           : Sums across all all HSPs.
 Example   : $num_unaln = $hit_object->num_unaligned_hit();
 Returns   : Integer
 Argument  : none
 Comments  : If you need data for each HSP, use hsps() and then interate
           : through the HSP objects.

=cut

sub num_unaligned_hit {
    my $self = shift;
    # why does this method even exist?!
    return $self->gaps('hit');
}

=head2 num_unaligned_query

 Usage     : $hit_object->num_unaligned_query();
 Purpose   : Get the number of the unaligned residues in the query sequence.
           : Sums across all all HSPs.
 Example   : $num_unaln = $hit_object->num_unaligned_query();
 Returns   : Integer
 Argument  : none
 Comments  : If you need data for each HSP, use hsps() and then interate
           : through the HSP objects.

=cut

sub num_unaligned_query {
    my $self = shift;
	# why does this method even exist?!
    return $self->gaps('query');
}

# aliasing for Steve's method names
*hit_description = \&description;
*hit_length = \&length;

1;
