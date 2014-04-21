#
# BioPerl module for Bio::LocatableSeq
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::LocatableSeq - A Bio::PrimarySeq object with start/end points on it
that can be projected into a MSA or have coordinates relative to
another seq.

=head1 SYNOPSIS

    use Bio::LocatableSeq;
    my $seq = Bio::LocatableSeq->new(-seq => "CAGT-GGT",
                    -id  => "seq1",
                    -start => 1,
                    -end   => 7);

    # a normal sequence object
    $locseq->seq();
    $locseq->id();

    # has start,end points
    $locseq->start();
    $locseq->end();

    # inherits off RangeI, so range operations possible

=head1 DESCRIPTION

The LocatableSeq sequence object was developed mainly because the SimpleAlign
object requires this functionality, and in the rewrite of the Sequence object we
had to decide what to do with this.

It is, to be honest, not well integrated with the rest of bioperl. For example,
the trunc() function does not return a LocatableSeq object, as some might have
thought. Also, the sequence is not a Bio::SeqI, so the location is simply
inherited from Bio::RangeI and is not stored in a Bio::Location. 

There are all sorts of nasty gotcha's about interactions between coordinate
systems when these sort of objects are used. Some mapping now occurs to deal
with HSP data, however it can probably be integrated in better and most methods
do not implement it correctly yet. Also, several PrimarySeqI methods (subseq(),
trunc(), etc.) do not behave as expected and must be used with care. Due to this,
LocatableSeq functionality is to be refactored in a future BioPerl release.
However, for alignment functionality it works adequately for the time being.

If you do not need alignment functionality, L<Bio::SeqfeatureI>-implementing
modules may be a suitable alternative to L<Bio::LocatableSeq>. For example,
L<Bio::SeqFeature::Generic> and L<Bio::SeqFeature::Lite> provide methods to
attach a sequence to a specific region of a parent sequence and to set other
useful attributes.

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

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut



package Bio::LocatableSeq;
use strict;

use Bio::Location::Simple;
use Bio::Location::Fuzzy;
use vars qw($GAP_SYMBOLS $OTHER_SYMBOLS $FRAMESHIFT_SYMBOLS $RESIDUE_SYMBOLS $MATCHPATTERN);

# The following global variables contain symbols used to represent gaps,
# frameshifts, residues, and other valid symbols. These are set at compile-time;
# expect scoping errors when using 'local' and resetting $MATCHPATTERN (see
# LocatableSeq.t)

$GAP_SYMBOLS = '\-\.=~';
$FRAMESHIFT_SYMBOLS = '\\\/';
$OTHER_SYMBOLS = '\?';
$RESIDUE_SYMBOLS = '0-9A-Za-z\*';
$MATCHPATTERN = $RESIDUE_SYMBOLS.$GAP_SYMBOLS.$FRAMESHIFT_SYMBOLS.$OTHER_SYMBOLS;

use base qw(Bio::PrimarySeq Bio::RangeI);


sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($start,$end,$strand, $mapping, $fs, $nse) =
    $self->_rearrange( [qw(START
                        END
                        STRAND
                        MAPPING
                        FRAMESHIFTS
                        FORCE_NSE
                        )],
               @args);
    
    $mapping ||= [1,1];
    $self->mapping($mapping);
    $nse || 0;
    $self->force_nse($nse);
    defined $fs    && $self->frameshifts($fs);
    defined $start && $self->start($start);
    defined $end   && $self->end($end);
    defined $strand && $self->strand($strand);

    return $self; # success - we hope!
}


=head2 start

 Title   : start
 Usage   : $obj->start($newval)
 Function: Get/set the 1-based start position of this sequence in the original
           sequence. '0' means before the original sequence starts.
 Returns : value of start
 Args    : newvalue (optional)

=cut

sub start {
    my $self = shift;
    if( @_ ) {
        my $value = shift;
        $self->{'start'} = $value;
    }
    return $self->{'start'} if defined $self->{'start'};
    return 1                if $self->seq;
    return;
}


=head2 end

 Title   : end
 Usage   : $obj->end($newval)
 Function: Get/set the 1-based end position of this sequence in the original
           sequence. '0' means before the original sequence starts.
 Returns : value of end
 Args    : newvalue (optional)
 Note    : although this is a get/set, it checks passed values against the
           calculated end point ( derived from the sequence and based on
           $GAP_SYMBOLS and possible frameshifts() ).  If there is no match,
           it will warn and set the proper value.  Probably best used for
           debugging proper sequence calculations.

=cut

sub end {
    my $self = shift;
    if( @_ ) {
        my $value = shift;
        my $st = $self->start;
        # start of 0 usually means the sequence is all gaps but maps to
        # other sequences in an alignment
        if ($self->seq && $st != 0 ) {
            my $len = $self->_ungapped_len;
            my $calend = $st + $len - 1;
            my $id = $self->id || 'unknown';
            if ($calend != $value) {
                $self->warn("In sequence $id residue count gives end value ".
                "$calend.  \nOverriding value [$value] with value $calend for ".
                "Bio::LocatableSeq::end().\n".$self->seq);
                $value = $calend;
            }
        }
        $self->{'end'} = $value;
    }

    if (defined $self->{'end'}) {
        return $self->{'end'}
    } elsif ( my $len = $self->_ungapped_len) {
        return $len + $self->start - 1;
    } else {
        return;
    }
}


# changed 08.10.26 to return ungapped length, not the calculated end
# of the sequence
sub _ungapped_len {
    my $self = shift;
    return unless my $string = $self->seq;
    my ($map_res, $map_coord) = $self->mapping;
    my $offset = 0;
    if (my %data = $self->frameshifts) {
        map {$offset += $_} values %data;
    }
    $string =~ s{[$GAP_SYMBOLS$FRAMESHIFT_SYMBOLS]+}{}g;
    return CORE::length($string)/($map_res/$map_coord) + $offset/($map_coord/$map_res);
}

#sub length {
#    my $self = shift;
#    return unless my $string = $self->seq;
#    $string =~ s{[$GAP_SYMBOLS$FRAMESHIFT_SYMBOLS]+}{}g;
#    return CORE::length($string);
#}


=head2 strand

 Title   : strand
 Usage   : $obj->strand($newval)
 Function: return or set the strandedness
 Returns : the value of the strandedness (-1, 0 or 1)
 Args    : the value of the strandedness (-1, 0 or 1)

=cut

sub strand {
   my $self = shift;
   if( @_ ) {
        my $value = shift;
        $self->{'strand'} = $value;
    }
    return $self->{'strand'};
}


=head2 mapping

 Title   : mapping
 Usage   : $obj->mapping($newval)
 Function: return or set the mapping indices (indicates # symbols/positions in
           the source string mapping to # of coordinate positions)
 Returns : two-element array (# symbols => # coordinate pos)
 Args    : two elements (# symbols => # coordinate pos); this can also be
           passed in as an array reference of the two elements (as might be
           passed upon Bio::LocatableSeq instantiation, for instance).

=cut

sub mapping {
    my $self = shift;
    if( @_ ) {
        my @mapping = (ref($_[0]) eq 'ARRAY') ? @{$_[0]} : @_;
        $self->throw("Must pass two values (# residues mapped to # positions)")
            if @mapping != 2;
        if ((grep {$_ != 1 && $_ != 3} @mapping) || ($mapping[0] == 3 && $mapping[1] == 3)) {
            $self->throw("Mapping values other than 1 or 3 are not currently supported")
        }
        $self->{'_mapping'} = \@mapping;
    }
    $self->throw('Mapping for LocatableSeq not set') if !exists $self->{'_mapping'};
    return @{ $self->{'_mapping'} };
}


=head2 frameshifts

 Title   : frameshifts
 Usage   : $obj->frameshifts($newval)
 Function: get/set the frameshift hash, which contains sequence positions as
           keys and the shift (-2, -1, 1, 2) as the value
 Returns : hash
 Args    : hash or hash reference

=cut

sub frameshifts {
    my $self = shift;
    if( @_ ) {
        if (ref $_[0] eq 'HASH') {
            $self->{_frameshifts} = $_[0];
        } else {
            # assume this is a full list to be converted to a hash
            $self->{_frameshifts} = \%{@_} # coerce into hash ref
        }
    }
    (defined $self->{_frameshifts} && ref $self->{_frameshifts} eq 'HASH') ?
        return %{$self->{_frameshifts}} : return ();
}


=head2 get_nse

 Title   : get_nse
 Usage   :
 Function: read-only name of form id/start-end
 Example :
 Returns :
 Args    :

=cut

sub get_nse {
   my ($self,$char1,$char2) = @_;

   $char1 ||= "/";
   $char2 ||= "-";
   
   my ($id, $st, $end, $strand)  = ($self->id(), $self->start(),
                                    $self->end(), $self->strand || 0);
   
   if ($self->force_nse) {
        $id  ||= '';
        $st  ||= 0;
        $end ||= 0;
   }
   
   $self->throw("Attribute id not set") unless defined($id);
   $self->throw("Attribute start not set") unless defined($st);
   $self->throw("Attribute end not set") unless defined($end);
   
   if ($strand && $strand == -1) {
      ($st, $end) = ($end, $st);
   }
   
   #Stockholm Rfam includes version if present so it is optional
   my $v = $self->version ? '.'.$self->version : '';
   return join('',$id, $v, $char1, $st, $char2, $end);
}


=head2 force_nse

 Title   : force_nse
 Usage   : $ls->force_nse()
 Function: Boolean which forces get_nse() to build an NSE, regardless
           of whether id(), start(), or end() is set
 Returns : Boolean value
 Args    : (optional) Boolean (1 or 0)
 Note    : This will convert any passed value evaluating as TRUE/FALSE to 1/0
           respectively

=cut

sub force_nse {
    my ($self, $flag) = @_;
    if (defined $flag) {
        $flag ? (return $self->{'_force_nse'} = 1) : (return $self->{'_force_nse'} = 0);
    }
    return $self->{'_force_nse'};
}


=head2 num_gaps

 Title   : num_gaps
 Usage   :$self->num_gaps('.')
 Function:Gets number of gaps in the sequence. The count excludes
           leading or trailing gap characters.

           Valid bioperl sequence characters are [A-Za-z\-\.\*]. Of
           these, '.' and '-' are counted as gap characters unless an
           optional argument specifies one of them.

 Returns : number of internal gaps in the sequence.
 Args    : a gap character (optional)
 Status  : Stable
 Note    : replaces no_gaps

=cut

sub num_gaps {
    my ($self,$char) = @_;
    my ($seq, $count) = (undef, 0);

    # default gap characters
    $char ||= $GAP_SYMBOLS;

    $self->warn("I hope you know what you are doing setting gap to [$char]")
        unless $char =~ /[$GAP_SYMBOLS]/;

    $seq = $self->seq;
    return 0 unless $seq; # empty sequence does not have gaps

    $seq =~ s/^([$char]+)//;
    $seq =~ s/([$char]+)$//;
    while ( $seq =~ /[$char]+/g ) {
        $count++;
    }

    return $count;
}


=head2 column_from_residue_number

 Title   : column_from_residue_number
 Usage   : $col = $seq->column_from_residue_number($resnumber)
 Function:

           This function gives the position in the alignment
           (i.e. column number) of the given residue number in the
           sequence. For example, for the sequence

         Seq1/91-97 AC..DEF.GH

           column_from_residue_number(94) returns 6.

           An exception is thrown if the residue number would lie
           outside the length of the aligment
           (e.g. column_from_residue_number( "Seq2", 22 )

 Returns : A column number for the position of the
           given residue in the given sequence (1 = first column)
 Args    : A residue number in the whole sequence (not just that
           segment of it in the alignment)

=cut

sub column_from_residue_number {
    my ($self, $resnumber) = @_;

    $self->throw("Residue number has to be a positive integer, not [$resnumber]")
    unless $resnumber =~ /^\d+$/ and $resnumber > 0;

    if ($resnumber >= $self->start() and $resnumber <= $self->end()) {
        my @chunks;
        my $column_incr;
        my $current_column;
        my $current_residue = $self->start - 1;
        my $seq = $self->seq;
        my $strand = $self->strand || 0;

        if ($strand == -1) {
           #@chunks = reverse $seq =~ m/[^\.\-]+|[\.\-]+/go;
            @chunks = reverse $seq =~ m/[$RESIDUE_SYMBOLS]+|[$GAP_SYMBOLS]+/go;
            $column_incr = -1;
            $current_column = (CORE::length $seq) + 1;
        }
        else {
            #@chunks = $seq =~ m/[^\.\-]+|[\.\-]+/go;
            @chunks = $seq =~ m/[$RESIDUE_SYMBOLS]+|[$GAP_SYMBOLS]+/go;
            $column_incr = 1;
            $current_column = 0;
        }

        while (my $chunk = shift @chunks) {
            #if ($chunk =~ m|^[\.\-]|o) {
            if ($chunk =~ m|^[$GAP_SYMBOLS]|o) {
                $current_column += $column_incr * CORE::length($chunk);
            }
            else {
                if ($current_residue + CORE::length($chunk) < $resnumber) {
                    $current_column += $column_incr * CORE::length($chunk);
                    $current_residue += CORE::length($chunk);
                }
                else {
                    if ($strand == -1) {
                        $current_column -= $resnumber - $current_residue;
                    }
                    else {
                        $current_column += $resnumber - $current_residue;
                    }
                    return $current_column;
                }
            }
        }
    }

    $self->throw("Could not find residue number $resnumber");

}


=head2 location_from_column

 Title   : location_from_column
 Usage   : $loc = $ali->location_from_column($column_number)
 Function:

           This function gives the residue number for a given position
           in the alignment (i.e. column number) of the given. Gaps
           complicate this process and force the output to be a
           L<Bio::Location::Simple> where values can be undefined. 
           For example, for the sequence:

         Seq/91-96 .AC..DEF.G.

           location_from_column( 3 ) position 92
           location_from_column( 4 ) position 92^93
           location_from_column( 9 ) position 95^96
           location_from_column( 1 ) position undef

           An exact position returns a Bio::Location::Simple object
           where where location_type() returns 'EXACT', if a position
           is between bases location_type() returns 'IN-BETWEEN'.
           Column before the first residue returns undef. Note that if
           the position is after the last residue in the alignment,
           that there is no guarantee that the original sequence has
           residues after that position.

           An exception is thrown if the column number is not within
           the sequence.

 Returns : Bio::Location::Simple or undef
 Args    : A column number
 Throws  : If column is not within the sequence

See L<Bio::Location::Simple> for more.

=cut

sub location_from_column {
    my ($self, $column) = @_;

    $self->throw("Column number has to be a positive integer, not [$column]")
    unless $column =~ /^\d+$/ and $column > 0;
    $self->throw("Column number [$column] is larger than".
         " sequence length [". $self->length. "]")
    unless $column <= $self->length;

    my ($loc);
    my $s = $self->subseq(1,$column);
    $s =~ s/[^a-zA-Z\*]//g;

    my $pos = CORE::length $s;

    my $start = $self->start || 0 ;
    my $strand = $self->strand() || 1;
    my $relative_pos = ($strand == -1)
        ? ($self->end - $pos + 1)
    : ($pos + $start - 1);
    if ($self->subseq($column, $column) =~ /[a-zA-Z\*]/ ) {
    $loc = Bio::Location::Simple->new
        (-start  => $relative_pos,
         -end    => $relative_pos,
         -strand => 1,
         );
    } elsif ($pos == 0 and $self->start == 1) {
    } else {
      my ($start,$end) = ($relative_pos, $relative_pos + $strand);
      if ($strand == -1) {
    ($start,$end) = ($end,$start);
      }
    $loc = Bio::Location::Simple->new
        (-start         => $start,
         -end           => $end,
         -strand        => 1,
         -location_type => 'IN-BETWEEN'
         );
    }
    return $loc;
}


=head2 revcom

 Title   : revcom
 Usage   : $rev = $seq->revcom()
 Function: Produces a new Bio::LocatableSeq object which
           has the reversed complement of the sequence. For protein
           sequences this throws an exception of "Sequence is a
           protein. Cannot revcom"

 Returns : A new Bio::LocatableSeq object
 Args    : none

=cut

sub revcom {
    my ($self) = @_;
    # since we don't know whether sequences without 1 => 1 correlation can be
    # revcom'd, kick back
    if (grep {$_ != 1} $self->mapping) {
        $self->warn('revcom() not supported for sequences with mapped values of > 1');
        return;
    }
    my $new = $self->SUPER::revcom;
    $new->strand($self->strand * -1) if $self->strand;
    $new->start($self->start) if $self->start;
    $new->end($self->end) if $self->end;
    return $new;
}


=head2 trunc

 Title   : trunc
 Usage   : $subseq = $myseq->trunc(10,100);
 Function: Provides a truncation of a sequence,
 Returns : a fresh Bio::PrimarySeqI implementing object
 Args    : Two integers denoting first and last columns of the
           sequence to be included into sub-sequence.

=cut

sub trunc {
    my ($self, $start, $end) = @_;
    my $new = $self->SUPER::trunc($start, $end);
    $new->strand($self->strand);

    # end will be automatically calculated
    $start = $end if $self->strand && $self->strand == -1;

    $start = $self->location_from_column($start);
    $start ? ($start = $start->end) : ($start = 1);
    $new->start($start) if $start;

    return $new;
}


=head2 validate_seq

 Title   : validate_seq
 Usage   : if(! $seqobj->validate_seq($seq_str) ) {
                print "sequence $seq_str is not valid for an object of
                alphabet ",$seqobj->alphabet, "\n";
           }
 Function: Test that the given sequence is valid, i.e. contains only valid
           characters. The allowed characters are all letters (A-Z) and '-','.',
           '*','?','=' and '~'. Spaces are not valid. Note that this
           implementation does not take alphabet() into account.
 Returns : 1 if the supplied sequence string is valid, 0 otherwise.
 Args    : - Sequence string to be validated
           - Boolean to throw an error if the sequence is invalid

=cut

sub validate_seq {
    my ($self, $seqstr, $throw) = @_;
    $seqstr = '' if not defined $seqstr;
    $throw  = 0  if not defined $throw ; # 0 for backward compatiblity
    if ( (CORE::length $seqstr > 0         ) &&
         ($seqstr !~ /^([$MATCHPATTERN]+)$/) ) {
        if ($throw) {
            $self->throw("Failed validation of sequence '".(defined($self->id) ||
            '[unidentified sequence]')."'. Invalid characters were: " .
            join('',($seqstr =~ /([^$MATCHPATTERN]+)/g)));
        }
        return 0;
    }
    return 1;
}


################## DEPRECATED METHODS ##################


=head2 no_gap

 Title     : no_gaps
 Usage     : $self->no_gaps('.')
 Function  : Gets number of gaps in the sequence. The count excludes
             leading or trailing gap characters.

             Valid bioperl sequence characters are [A-Za-z\-\.\*]. Of
             these, '.' and '-' are counted as gap characters unless an
             optional argument specifies one of them.

 Returns   : number of internal gaps in the sequence.
 Args      : a gap character (optional)
 Status    : Deprecated (in favor of num_gaps()) 

=cut

sub no_gaps {
    my $self = shift;
    $self->deprecated( -warn_version  => 1.0069,
                       -throw_version => 1.0075,
                       -message => 'Use of method no_gaps() is deprecated, use num_gaps() instead' );
    return $self->num_gaps(@_);
}


=head2 no_sequences

 Title     : no_sequences
 Usage     : $gaps = $seq->no_sequences
 Function  : number of sequence in the sequence alignment
 Returns   : integer
 Argument  :
 Status    : Deprecated (in favor of num_sequences())

=cut

sub no_sequences {
    my $self = shift;
    $self->deprecated( -warn_version  => 1.0069,
                       -throw_version => 1.0075,
                       -message => 'Use of method no_sequences() is deprecated, use num_sequences() instead' );
    return $self->num_sequences(@_);
}

1;
