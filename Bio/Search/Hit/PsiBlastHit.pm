#-----------------------------------------------------------------
#
# BioPerl module Bio::Search::Hit::PsiBlastHit
#
# (This module was originally called Bio::Tools::Blast::Sbjct)
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Steve Chervitz <sac@bioperl.org>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

## POD Documentation:

=head1 NAME

Bio::Search::Hit::PsiBlastHit - Bioperl BLAST Hit object

=head1 SYNOPSIS

See L<Bio::Search::Result::BlastResult>.

=head1 DESCRIPTION

The Bio::Search::Hit::PsiBlastHit.pm module encapsulates data and
methods for manipulating "hits" from a BLAST report. A BLAST hit is a
collection of HSPs along with other metadata such as sequence name and
score information. Hit objects are accessed via
L<Bio::Search::Result::BlastResult> objects after parsing a BLAST
report using the L<Bio::SearchIO> system.

In Blast lingo, the "sbjct" sequences are all the sequences in a
target database which were compared against a "query" sequence.  The
terms "sbjct" and "hit" will be used interchangeably in this module.
All methods that take 'sbjct' as an argument also support 'hit' as a
synonym.

This module supports BLAST versions 1.x and 2.x, gapped and ungapped,
and PSI-BLAST.

The construction of PsiBlastHit objects is performed by
Bio::SearchIO::blast::PsiBlastHitFactory in a process that is
orchestrated by the Blast parser (L<Bio::SearchIO::blast>).
The resulting PsiBlastHits are then accessed via
L<Bio::Search::Result::BlastResult>). Therefore, you do not need to
use L<Bio::Search::Hit::PsiBlastHit>) directly. If you need to
construct PsiBlastHits directly, see the C<new()> function for details.

For L<Bio::SearchIO> BLAST parsing usage examples, see the
C<examples/search-blast> directory of the Bioperl distribution.


=head2 HSP Tiling and Ambiguous Alignments

If a Blast hit has more than one HSP, the Bio::Search::Hit::PsiBlastHit.pm
object has the ability to merge overlapping HSPs into contiguous
blocks. This permits the PsiBlastHit object to sum data across all HSPs
without counting data in the overlapping regions multiple times, which
would happen if data from each overlapping HSP are simply summed.  HSP
tiling is performed automatically when methods of the PsiBlastHit object
that rely on tiled data are invoked. These include
L<frac_identical()|frac_identical>, L<frac_conserved()|frac_conserved>, L<gaps()|gaps>,
L<frac_aligned_query()|frac_aligned_query>, L<frac_aligned_hit()|frac_aligned_hit>,
L<num_unaligned_query()|num_unaligned_query>, L<num_unaligned_hit()|num_unaligned_hit>.

It also permits the assessment of an "ambiguous alignment" if the
query (or sbjct) sequences from different HSPs overlap
(see L<ambiguous_aln()|ambiguous_aln>). The existence
of an overlap could indicate a biologically interesting region in the
sequence, such as a repeated domain.  The PsiBlastHit object uses the
C<-OVERLAP> parameter to determine when two sequences overlap; if this is
set to 2 -- the default -- then any two sbjct or query HSP sequences
must overlap by more than two residues to get merged into the same
contig and counted as an overlap. See the L<BUGS | BUGS> section below for
"issues" with HSP tiling.


The results of the HSP tiling is reported with the following ambiguity codes:

   'q' = Query sequence contains multiple sub-sequences matching
         a single region in the sbjct sequence.

   's' = Subject (PsiBlastHit) sequence contains multiple sub-sequences matching
         a single region in the query sequence.

   'qs' = Both query and sbjct sequences contain more than one
          sub-sequence with similarity to the other sequence.


For addition information about ambiguous BLAST alignments, see
L<Bio::Search::BlastUtils::tile_hsps()|Bio::Search::BlastUtils>

=head1 DEPENDENCIES

Bio::Search::Hit::PsiBlastHit.pm is a concrete class that inherits from
L<Bio::Root::Root> and L<Bio::Search::Hit::HitI>.  and relies on
L<Bio::Search::HSP::BlastHSP>.


=head1 BUGS

One consequence of the HSP tiling is that methods that rely on HSP
tiling such as L<frac_identical()|frac_identical>, L<frac_conserved()|frac_conserved>, L<gaps()|gaps>
etc. may report misleading numbers when C<-OVERLAP> is set to a large
number.  For example, say we have two HSPs and the query sequence tile
as follows:

            1      8             22      30        40             60
 Full seq:  ------------------------------------------------------------
                    *  ** *   **
 HSP1:             ---------------                    (6 identical matches)
                              **   **  **
 HSP2:                        -------------           (6 identical matches)


If C<-OVERLAP> is set to some number over 4, HSP1 and HSP2 will not be
tiled into a single contig and their numbers of identical matches will
be added, giving a total of 12, not 10 if they had be combined into
one contig. This can lead to number greater than 1.0 for methods
L<frac_identical()|frac_identical> and L<frac_conserved()|frac_conserved>. This is less of an issue
with gapped Blast since it tends to combine HSPs that would be listed
separately without gapping.  (Fractions E<gt>1.0 can be viewed as a
signal for an interesting alignment that warrants further inspection,
thus turning this bug into a feature :-).

Using large values for C<-OVERLAP> can lead to incorrect numbers
reported by methods that rely on HSP tiling but can be useful if you
care more about detecting ambiguous alignments.  Setting C<-OVERLAP>
to zero will lead to the most accurate numbers for the
tiling-dependent methods but will be useless for detecting overlapping
HSPs since all HSPs will appear to overlap.


=head1 SEE ALSO

 Bio::Search::HSP::BlastHSP.pm         - Blast HSP object.
 Bio::Search::Result::BlastResult.pm   - Blast Result object.
 Bio::Search::Hit::HitI.pm             - Interface implemented by PsiBlastHit.pm
 Bio::Root::Root.pm                    - Base class for PsiBlastHit.pm

Links:

 http://bio.perl.org/                       - Bioperl Project Homepage


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules.  Send your comments and suggestions preferably to one
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
the bugs and their resolution. Bug reports can be submitted via the
web:

    https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Steve Chervitz E<lt>sac@bioperl.orgE<gt>

See L<the FEEDBACK section | FEEDBACK> for where to send bug reports and comments.

=head1 ACKNOWLEDGEMENTS

This software was originally developed in the Department of Genetics
at Stanford University. I would also like to acknowledge my
colleagues at Affymetrix for useful feedback.

=head1 COPYRIGHT

Copyright (c) 1996-2001 Steve Chervitz. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Search::Hit::PsiBlastHit;

use strict;
use Bio::Search::BlastUtils;
use vars qw(%SUMMARY_OFFSET);

use overload
    '""' => \&to_string;

use base qw(Bio::Root::Root Bio::Search::Hit::HitI);


=head2 new

 Usage     : $hit = Bio::Search::Hit::PsiBlastHit->new( %named_params );
           : Bio::Search::Hit::PsiBlastHit.pm objects are constructed
           : automatically by Bio::SearchIO::PsiBlastHitFactory.pm,
           : so there is no need for direct instantiation.
 Purpose   : Constructs a new PsiBlastHit object and Initializes key variables
           : for the hit.
 Returns   : A Bio::Search::Hit::PsiBlastHit object
 Argument  : Named Parameters:
           : Parameter keys are case-insensitive.
           :     -RAW_DATA   => array reference holding raw BLAST report data
           :                    for a single hit. This includes all lines
           :                    within the HSP alignment listing section of a
           :                    traditional BLAST or PSI-BLAST (non-XML) report,
           :                    starting at (or just after) the leading '>'.
           :         -HOLD_RAW_DATA => boolean, should -RAW_DATA be saved within the object.
           :         -QUERY_LEN  => Length of the query sequence
           :         -ITERATION  => integer (PSI-BLAST iteration number in which hit was found)
           :         -OVERLAP    => integer (maximum overlap between adjacent
           :                    HSPs when tiling)
           :         -PROGRAM    => string (type of Blast: BLASTP, BLASTN, etc)
           :         -SIGNIF     => significance
           :         -IS_PVAL    => boolean, true if -SIGNIF contains a P-value
           :         -SCORE      => raw BLAST score
           :         -FOUND_AGAIN   => boolean, true if this was a hit from the
           :                       section of a PSI-BLAST with iteration > 1
           :                       containing sequences that were also found
           :                       in iteration 1.
 Comments  : This object accepts raw Blast report data not because it
           : is required for parsing, but in order to retrieve it
           : (only available if -HOLD_RAW_DATA is set to true).

See Also   : L<Bio::Search::BlastUtils::tile_hsps()|Bio::Search::BlastUtils>, L<Bio::Root::Root::new()|Bio::Root::Root>

=cut

#-------------------
sub new {
#-------------------
    my ($class, @args ) = @_;
    my $self = $class->SUPER::new( @args );

    my ($raw_data, $signif, $is_pval, $hold_raw);

    ($self->{'_blast_program'}, $self->{'_query_length'}, $raw_data, $hold_raw,
     $self->{'_overlap'}, $self->{'_iteration'}, $signif, $is_pval,
     $self->{'_score'}, $self->{'_found_again'} ) =
       $self->_rearrange( [qw(PROGRAM
                              QUERY_LEN
                              RAW_DATA
                              HOLD_RAW_DATA
                              OVERLAP
                              ITERATION
                              SIGNIF
                              IS_PVAL
                              SCORE
                              FOUND_AGAIN )], @args );

    # TODO: Handle this in parser. Just pass in name parameter.
    $self->_set_id( $raw_data->[0] );

    if($is_pval) {
        $self->{'_p'} = $signif;
    } else {
        $self->{'_expect'} = $signif;
    }

    if( $hold_raw ) {
        $self->{'_hit_data'} = $raw_data;
    }

    return $self;
}

sub DESTROY {
    my $self=shift;
    #print STDERR "-->DESTROYING $self\n";
}


#=================================================
# Begin Bio::Search::Hit::HitI implementation
#=================================================

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

#----------------
sub algorithm {
#----------------
    my ($self,@args) = @_;
    return $self->{'_blast_program'};
}

=head2 name

 Usage     : $hit->name([string]);
 Purpose   : Set/Get a string to identify the hit.
 Example   : $name = $hit->name;
           : $hit->name('M81707');
 Returns   : String consisting of the hit's name or undef if not set.
 Comments  : The name is parsed out of the "Query=" line as the first chunk of
             non-whitespace text. If you want the rest of the line, use
             $hit->description().

See Also: L<accession()|accession>

=cut

#'

#----------------
sub name {
#----------------
    my $self = shift;
    if (@_) {
        my $name = shift;
        $name =~ s/^\s+|(\s+|,)$//g;
        $self->{'_name'} = $name;
    }
    return $self->{'_name'};
}

=head2 description

 Usage     : $hit_object->description( [integer] );
 Purpose   : Set/Get a description string for the hit.
             This is parsed out of the "Query=" line as everything after
             the first chunk of non-whitespace text. Use $hit->name()
             to get the first chunk (the ID of the sequence).
 Example   : $description = $hit->description;
           : $desc_60char = $hit->description(60);
 Argument  : Integer (optional) indicating the desired length of the
           : description string to be returned.
 Returns   : String consisting of the hit's description or undef if not set.

=cut

#'

#----------------
sub description {
#----------------
    my( $self, $len ) = @_;
    $len = (defined $len) ? $len : (CORE::length $self->{'_description'});
    return substr( $self->{'_description'}, 0 ,$len );
}

=head2 accession

 Title   : accession
 Usage   : $acc = $hit->accession();
 Function: Retrieve the accession (if available) for the hit
 Returns : a scalar string (empty string if not set)
 Args    : none
 Comments: Accession numbers are extracted based on the assumption that they
           are delimited by | characters (NCBI-style). If this is not the case,
           use the name() method and parse it as necessary.

See Also: L<name()|name>

=cut

#--------------------
sub accession {
#--------------------
    my $self = shift;
    if(@_) { $self->{'_accession'} = shift; }
    $self->{'_accession'} || '';
}

=head2 raw_score

 Usage     : $hit_object->raw_score();
 Purpose   : Gets the BLAST score of the best HSP for the current Blast hit.
 Example   : $score = $hit_object->raw_score();
 Returns   : Integer
 Argument  : n/a
 Throws    : n/a

See Also   : L<bits()|bits>

=cut

#----------
sub raw_score {
#----------
    my $self = shift;

    # The check for $self->{'_score'} is a remnant from the 'query' mode days
    # in which the sbjct object would collect data from the description line only.

    my ($score);
    if(not defined($self->{'_score'})) {
        $score = $self->hsp->score;
    } else {
        $score = $self->{'_score'};
    }
    return $score;
}


=head2 length

 Usage     : $hit_object->length();
 Purpose   : Get the total length of the hit sequence.
 Example   : $len = $hit_object->length();
 Returns   : Integer
 Argument  : n/a
 Throws    : n/a
 Comments  : Developer note: when using the built-in length function within
           : this module, call it as CORE::length().

See Also   : L<logical_length()|logical_length>,  L<length_aln()|length_aln>

=cut

#-----------
sub length {
#-----------
    my $self = shift;
    return $self->{'_length'};
}

=head2 significance

Equivalent to L<signif()|signif>

=cut

#----------------
sub significance { shift->signif( @_ ); }
#----------------


=head2 next_hsp

 Title    : next_hsp
 Usage    : $hsp = $obj->next_hsp();
 Function : returns the next available High Scoring Pair object
 Example  :
 Returns  : Bio::Search::HSP::BlastHSP or undef if finished
 Args     : none

=cut

#----------------
sub next_hsp {
#----------------
    my $self = shift;

    unless($self->{'_hsp_queue_started'}) {
        $self->{'_hsp_queue'} = [$self->hsps()];
        $self->{'_hsp_queue_started'} = 1;
    }
    pop @{$self->{'_hsp_queue'}};
}

#=================================================
# End Bio::Search::Hit::HitI implementation
#=================================================


# Providing a more explicit method for getting name of hit
# (corresponds with column name in HitTableWriter)
#----------------
sub hit_name {
#----------------
    my $self = shift;
    $self->name( @_ );
}

# Older method Delegates to description()
#----------------
sub desc {
#----------------
    my $self = shift;
    return $self->description( @_ );
}

# Providing a more explicit method for getting description of hit
# (corresponds with column name in HitTableWriter)
#----------------
sub hit_description {
#----------------
    my $self = shift;
    return $self->description( @_ );
}

=head2 score

Equivalent to L<raw_score()|raw_score>

=cut

#----------------
sub score { shift->raw_score( @_ ); }
#----------------


=head2 hit_length

Equivalent to L<length()|length>

=cut

# Providing a more explicit method for getting length of hit
#----------------
sub hit_length { shift->length( @_ ); }
#----------------


=head2 signif

 Usage     : $hit_object->signif( [format] );
 Purpose   : Get the P or Expect value for the best HSP of the given BLAST hit.
           : The value returned is the one which is reported in the description
           : section of the Blast report. For Blast1 and WU-Blast2, this
           : is a P-value, for Blast2, it is an Expect value.
 Example   : $obj->signif()        # returns 1.3e-34
           : $obj->signif('exp')   # returns -34
           : $obj->signif('parts') # returns (1.3, -34)
 Returns   : Float or scientific notation number (the raw P/Expect value, DEFAULT).
           : Integer if format == 'exp' (the magnitude of the base 10 exponent).
           : 2-element list (float, int) if format == 'parts' and P/Expect value
           :                is in scientific notation (see Comments).
 Argument  : format: string of 'raw' | 'exp' | 'parts'
           :    'raw' returns value given in report. Default. (1.2e-34)
           :    'exp' returns exponent value only (34)
           :    'parts' returns the decimal and exponent as a
           :            2-element list (1.2, -34)  (see Comments).
 Throws    : n/a
 Comments  : The signif() method provides a way to deal with the fact that
           : Blast1 and Blast2 formats (and WU- vs. NCBI-BLAST) differ in
           : what is reported in the description lines of each hit in the
           : Blast report. The signif() method frees any client code from
           : having to know if this is a P-value or an Expect value,
           : making it easier to write code that can process both
           : Blast1 and Blast2 reports. This is not necessarily a good thing,
           : since one should always know when one is working with P-values or
           : Expect values (hence the deprecated status).
           : Use of expect() is recommended since all hits will have an Expect value.
           :
           : Using the 'parts' argument is not recommended since it will not
           : work as expected if the expect value is not in scientific notation.
           : That is, floats are not converted into sci notation before
           : splitting into parts.

See Also   : L<p()|p>, L<expect()|expect>, L<Bio::Search::BlastUtils::get_exponent()|Bio::Search::BlastUtils>

=cut

#-------------
sub signif {
#-------------
# Some duplication of logic for p(), expect() and signif() for the sake of performance.
    my ($self, $fmt) = @_;

    my $val = defined($self->{'_p'}) ? $self->{'_p'} : $self->{'_expect'};

    # $val can be zero.
    defined($val) or $self->throw("Can't get P- or Expect value: HSPs may not have been set.");

    return $val if not $fmt or $fmt =~ /^raw/i;
    ## Special formats: exponent-only or as list.
    return &Bio::Search::BlastUtils::get_exponent($val) if $fmt =~ /^exp/i;
    return (split (/eE/, $val)) if $fmt =~ /^parts/i;

    ## Default: return the raw P/Expect-value.
    return $val;
}

#----------------
sub raw_hit_data {
#----------------
    my $self = shift;
    my $data = '>';
    # Need to add blank lines where we've removed them.
    foreach( @{$self->{'_hit_data'}} ) {
        if( $_ eq 'end') {
            $data .= "\n";
        }
        else {
            $data .= /^\s*(Score|Query)/ ? "\n$_" : $_;
        }
    }
    return $data;
}


#=head2 _set_length
#
# Usage     : $hit_object->_set_length( "233" );
# Purpose   : Set the total length of the hit sequence.
# Example   : $hit_object->_set_length( $len );
# Returns   : n/a
# Argument  : Integer (only when setting). Any commas will be stripped out.
# Throws    : n/a
#
#=cut

#-----------
sub _set_length {
#-----------
    my ($self, $len) = @_;
    $len =~ s/,//g; # get rid of commas
    $self->{'_length'} = $len;
}

#=head2 _set_description
#
# Usage     : Private method; called automatically during construction
# Purpose   : Sets the description of the hit sequence.
#            : For sequence without descriptions, does not set any description.
# Argument  : Array containing description (multiple lines).
# Comments  : Processes the supplied description:
#                1. Join all lines into one string.
#                2. Remove sequence id at the beginning of description.
#                3. Removes junk charactes at begin and end of description.
#
#=cut

#--------------
sub _set_description {
#--------------
    my( $self, @desc ) = @_;
    my( $desc);

#    print STDERR "PsiBlastHit: RAW DESC:\n@desc\n";

    $desc = join(" ", @desc);

    my $name = $self->name;

    if($desc) {
        $desc =~ s/^\s*\S+\s+//; # remove the sequence ID(s)
                                 # This won't work if there's no description.
        $desc =~ s/^\s*$name//;  # ...but this should.
        $desc =~ s/^[\s!]+//;
        $desc =~ s/ \d+$//;
        $desc =~ s/\.+$//;
        $self->{'_description'} = $desc;
    }

#    print STDERR "PsiBlastHit: _set_description =  $desc\n";
}

=head2 to_string

 Title   : to_string
 Usage   : print $hit->to_string;
 Function: Returns a string representation for the Blast Hit.
           Primarily intended for debugging purposes.
 Example : see usage
 Returns : A string of the form:
           [PsiBlastHit] <name> <description>
           e.g.:
           [PsiBlastHit] emb|Z46660|SC9725 S.cerevisiae chromosome XIII cosmid
 Args    : None

=cut

#----------------
sub to_string {
#----------------
    my $self = shift;
    return "[PsiBlastHit] " . $self->name . " " . $self->description;
}


#=head2 _set_id
#
# Usage     : Private method; automatically called by new()
# Purpose   : Sets the name of the PsiBlastHit sequence from the BLAST summary line.
#           : The identifier is assumed to be the first
#           : chunk of non-whitespace characters in the description line
#           : Does not assume any semantics in the structure of the identifier
#           : (Formerly, this method attempted to extract database name from
#           : the seq identifiers, but this was prone to break).
# Returns   : n/a
# Argument  : String containing description line of the hit from Blast report
#           : or first line of an alignment section (with or without the leading '>').
# Throws    : Warning if cannot locate sequence ID.
#
#See Also   : L<new()|new>, L<accession()|accession>
#
#=cut

#---------------
sub _set_id {
#---------------
    my( $self, $desc ) = @_;

    # New strategy: Assume only that the ID is the first white space
    # delimited chunk. Not attempting to extract accession & database name.
    # Clients will have to interpret it as necessary.
    if($desc =~ /^>?(\S+)\s*(.*)/) {
        my ($name, $desc) = ($1, $2);
        $self->name($name);
        $self->{'_description'} = $desc;
        # Note that this description comes from the summary section of the
        # BLAST report and so may be truncated. The full description will be
        # set from the alignment section. We're setting description here in case
        # the alignment section isn't being parsed.

        # Assuming accession is delimited with | symbols (NCBI-style)
        my @pieces = split(/\|/,$name);
        my $acc = pop @pieces;
        $self->accession( $acc );
    }
    else {
        $self->warn("Can't locate sequence identifier in summary line.", "Line = $desc");
        $desc = 'Unknown sequence ID' if not $desc;
        $self->name($desc);
    }
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

See Also   : L<Bio::Search::BlastUtils::tile_hsps>, L<HSP Tiling and Ambiguous Alignments>

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






=head2 bits

 Usage     : $hit_object->bits();
 Purpose   : Gets the BLAST bit score of the best HSP for the current Blast hit.
 Example   : $bits = $hit_object->bits();
 Returns   : Integer
 Argument  : n/a
 Throws    : Exception if bit score is not set.
 Comments  : For BLAST1, the non-bit score is listed in the summary line.

See Also   : L<score()|score>

=cut

#---------
sub bits {
#---------
    my $self = shift;

    # The check for $self->{'_bits'} is a remnant from the 'query' mode days
    # in which the sbjct object would collect data from the description line only.

    my ($bits);
    if(not defined($self->{'_bits'})) {
        $bits = $self->hsp->bits;
    } else {
        $bits = $self->{'_bits'};
    }
    return $bits;
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
    # in which the sbjct object would collect data from the description line only.

    my ($n);
    if(not defined($self->{'_n'})) {
        $n = $self->hsp->n;
    } else {
        $n = $self->{'_n'};
    }
    $n ||= $self->num_hsps;

    return $n;
}



=head2 frame

 Usage     : $hit_object->frame();
 Purpose   : Gets the reading frame for the best HSP after HSP tiling.
           : This is only valid for BLASTX and TBLASTN/X reports.
 Example   : $frame = $hit_object->frame();
 Returns   : Integer (-2 .. +2)
 Argument  : n/a
 Throws    : Exception if HSPs have not been set (BLAST2 reports).
 Comments  : This method requires that all HSPs be tiled. If they have not
           : already been tiled, they will be tiled first automatically..
           : If you don't want the tiled data, iterate through each HSP
           : calling frame() on each (use hsps() to get all HSPs).

See Also   : L<hsps()|hsps>

=cut

#----------'
sub frame {
#----------
    my $self = shift;

    Bio::Search::BlastUtils::tile_hsps($self) if not $self->{'_tile_hsps'};

    # The check for $self->{'_frame'} is a remnant from the 'query' mode days
    # in which the sbjct object would collect data from the description line only.

    my ($frame);
    if(not defined($self->{'_frame'})) {
        $frame = $self->hsp->frame('hit');
    } else {
        $frame = $self->{'_frame'};
    }
    return $frame;
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



=head2 expect

 Usage     : $hit_object->expect( [format] );
 Purpose   : Get the Expect value for the best HSP of the given BLAST hit.
 Example   : $e =  $sbjct->expect;
           : $e =  $sbjct->expect('exp');  # get exponent only.
           : ($num, $exp) = $sbjct->expect('parts');  # split sci notation into parts
 Returns   : Float or scientific notation number (the raw expect value, DEFAULT).
           : Integer if format == 'exp' (the magnitude of the base 10 exponent).
           : 2-element list (float, int) if format == 'parts' and Expect
           :                is in scientific notation (see Comments).
 Argument  : format: string of 'raw' | 'exp' | 'parts'
           :    'raw' returns value given in report. Default. (1.2e-34)
           :    'exp' returns exponent value only (34)
           :    'parts' returns the decimal and exponent as a
           :            2-element list (1.2, -34)  (see Comments).
 Throws    : Exception if the Expect value is not defined.
 Comments  : Using the 'parts' argument is not recommended since it will not
           : work as expected if the expect value is not in scientific notation.
           : That is, floats are not converted into sci notation before
           : splitting into parts.

See Also   : L<p()|p>, L<signif()|signif>, L<Bio::Search::BlastUtils::get_exponent()|Bio::Search::BlastUtils>

=cut

#-----------
sub expect {
#-----------
# Some duplication of logic for p(), expect() and signif() for the sake of performance.
    my ($self, $fmt) = @_;

    my $val;

    # For Blast reports that list the P value on the description line,
    # getting the expect value requires fully parsing the HSP data.
    # For NCBI blast, there's no problem.
    if(not defined($self->{'_expect'})) {
        if( defined $self->{'_hsps'}) {
            $self->{'_expect'} = $val = $self->hsp->expect;
        } else {
            # If _expect is not set and _hsps are not set,
            # then this must be a P-value-based report that was
            # run without setting the HSPs (shallow parsing).
            $self->throw("Can't get expect value. HSPs have not been set.");
        }
    } else {
        $val = $self->{'_expect'};
    }

    # $val can be zero.
    defined($val) or $self->throw("Can't get Expect value.");

    return $val if not $fmt or $fmt =~ /^raw/i;
    ## Special formats: exponent-only or as list.
    return &Bio::Search::BlastUtils::get_exponent($val) if $fmt =~ /^exp/i;
    return (split (/eE/, $val)) if $fmt =~ /^parts/i;

    ## Default: return the raw Expect-value.
    return $val;
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



=head2 hsp

 Usage     : $hit_object->hsp( [string] );
 Purpose   : Get a single BlastHSP.pm object for the present PsiBlastHit.pm object.
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
    my $self = shift;

    if (not defined $self->{'_hsps'}) {
        $self->throw("Can't get HSPs: data not collected.");
    }

    return scalar(@{$self->{'_hsps'}});
}



=head2 logical_length

 Usage     : $hit_object->logical_length( [seq_type] );
           : (mostly intended for internal use).
 Purpose   : Get the logical length of the hit sequence.
           : For query sequence of BLASTX and TBLASTX reports and the hit
           : sequence of TBLASTN and TBLASTX reports, the returned length
           : is the length of the would-be amino acid sequence (length/3).
           : For all other BLAST flavors, this function is the same as length().
 Example   : $len = $hit_object->logical_length();
 Returns   : Integer
 Argument  : seq_type = 'query' or 'hit' or 'sbjct' (default = 'query')
             ('sbjct' is synonymous with 'hit')
 Throws    : n/a
 Comments  : This is important for functions like frac_aligned_query()
           : which need to operate in amino acid coordinate space when dealing
           : with T?BLASTX type reports.

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

        # Adjust length based on BLAST flavor.
        if($self->{'_blast_program'} =~ /T?BLASTX/ ) {
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
           : in the PsiBlastHit object. If there is more than one HSP, the lowest start
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

See Also   : L<end()|end>, L<range()|range>, L<strand()|strand>, L<HSP Tiling and Ambiguous Alignments>, L<Bio::Search::HSP::BlastHSP::start|Bio::Search::HSP::BlastHSP>

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
           : in the PsiBlastHit object. If there is more than one HSP, the largest end
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

See Also   : L<start()|start>, L<range()|range>, L<strand()|strand>, L<HSP Tiling and Ambiguous Alignments>, L<Bio::Search::HSP::BlastHSP::end|Bio::Search::HSP::BlastHSP>

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

    sprintf( "%.2f", $self->{'_totalIdentical'}/$self->{'_length_aln_'.$seqType});
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

    sprintf( "%.2f", $self->{'_totalConserved'}/$self->{'_length_aln_'.$seqType});
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

    sprintf( "%.2f", $self->{'_length_aln_query'}/$self->logical_length('query'));
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
 Comments  : Note that HSPs are not tiled for this. This could be a problem
           : for hits containing mutually exclusive HSPs.
           : TODO: Consider tiling and then reporting seq_inds for the
           : best HSP contig.

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


=head2 iteration

 Usage     : $sbjct->iteration( );
 Purpose   : Gets the iteration number in which the Hit was found.
 Example   : $iteration_num = $sbjct->iteration();
 Returns   : Integer greater than or equal to 1
             Non-PSI-BLAST reports will report iteration as 1, but this number
             is only meaningful for PSI-BLAST reports.
 Argument  : none
 Throws    : none

See Also   : L<found_again()|found_again>

=cut

#----------------
sub iteration { shift->{'_iteration'} }
#----------------


=head2 found_again

 Usage     : $sbjct->found_again;
 Purpose   : Gets a boolean indicator whether or not the hit has
             been found in a previous iteration.
             This is only applicable to PSI-BLAST reports.

              This method indicates if the hit was reported in the
              "Sequences used in model and found again" section of the
              PSI-BLAST report or if it was reported in the
              "Sequences not found previously or not previously below threshold"
              section of the PSI-BLAST report. Only for hits in iteration > 1.

 Example   : if( $sbjct->found_again()) { ... };
 Returns   : Boolean (1 or 0) for PSI-BLAST report iterations greater than 1.
             Returns undef for PSI-BLAST report iteration 1 and non PSI_BLAST
             reports.
 Argument  : none
 Throws    : none

See Also   : L<found_again()|found_again>

=cut

#----------------
sub found_again { shift->{'_found_again'} }
#----------------


=head2 strand

 Usage     : $sbjct->strand( [seq_type] );
 Purpose   : Gets the strand(s) for the query, sbjct, or both sequences
           : in the best HSP of the PsiBlastHit object after HSP tiling.
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

See Also   : L<Bio::Search::HSP::BlastHSP::strand>()

=cut

#----------'
sub strand {
#----------
    my ($self, $seqType) = @_;

    Bio::Search::BlastUtils::tile_hsps($self) if not $self->{'_tile_hsps'};

    $seqType ||= (wantarray ? 'list' : 'query');
    $seqType = 'sbjct' if $seqType eq 'hit';

    my ($qstr, $hstr);
    # If there is only one HSP, defer this call to the solitary HSP.
    if($self->num_hsps == 1) {
        return $self->hsp->strand($seqType);
    }
    elsif( defined $self->{'_qstrand'}) {
        # Get the data computed during hsp tiling.
        $qstr = $self->{'_qstrand'};
        $hstr = $self->{'_sstrand'};
    }
    else {
        # otherwise, iterate through all HSPs collecting strand info.
        # This will return the string "-1/1" if there are HSPs on different strands.
        # NOTE: This was the pre-10/21/02 procedure which will no longer be used,
        # (unless the above elsif{} is commented out).
        my (%qstr, %hstr);
        foreach my $hsp( $self->hsps ) {
            my ( $q, $h ) = $hsp->strand();
            $qstr{ $q }++;
            $hstr{ $h }++;
        }
        $qstr = join( '/', sort keys %qstr);
        $hstr = join( '/', sort keys %hstr);
    }

    if($seqType =~ /list|array/i) {
        return ($qstr, $hstr);
    } elsif( $seqType eq 'query' ) {
        return $qstr;
    } else {
        return $hstr;
    }
}


1;
__END__

#####################################################################################
#                                END OF CLASS                                       #
#####################################################################################


=head1 FOR DEVELOPERS ONLY

=head2 Data Members

Information about the various data members of this module is provided for those
wishing to modify or understand the code. Two things to bear in mind:

=over 4

=item 1 Do NOT rely on these in any code outside of this module.

All data members are prefixed with an underscore to signify that they are private.
Always use accessor methods. If the accessor doesn't exist or is inadequate,
create or modify an accessor (and let me know, too!). (An exception to this might
be for BlastHSP.pm which is more tightly coupled to PsiBlastHit.pm and
may access PsiBlastHit data members directly for efficiency purposes, but probably
should not).

=item 2 This documentation may be incomplete and out of date.

It is easy for these data member descriptions to become obsolete as
this module is still evolving. Always double check this info and search
for members not described here.

=back

An instance of Bio::Search::Hit::PsiBlastHit.pm is a blessed reference to a hash containing
all or some of the following fields:

 FIELD           VALUE
 --------------------------------------------------------------
 _hsps          : Array ref for a list of Bio::Search::HSP::BlastHSP.pm objects.
                :
 _db            : Database identifier from the summary line.
                :
 _desc          : Description data for the hit from the summary line.
                :
 _length        : Total length of the hit sequence.
                :
 _score         : BLAST score.
                :
 _bits          : BLAST score (in bits). Matrix-independent.
                :
 _p             : BLAST P value. Obtained from summary section. (Blast1/WU-Blast only)
                :
 _expect        : BLAST Expect value. Obtained from summary section.
                :
 _n             : BLAST N value (number of HSPs) (Blast1/WU-Blast2 only)
                :
 _frame         : Reading frame for TBLASTN and TBLASTX analyses.
                :
 _totalIdentical: Total number of identical aligned monomers.
                :
 _totalConserved: Total number of conserved aligned monomers (a.k.a. "positives").
                :
 _overlap       : Maximum number of overlapping residues between adjacent HSPs
                : before considering the alignment to be ambiguous.
                :
 _ambiguous_aln : Boolean. True if the alignment of all HSPs is ambiguous.
                :
 _length_aln_query : Length of the aligned region of the query sequence.
                   :
 _length_aln_sbjct : Length of the aligned region of the sbjct sequence.


=cut

1;
