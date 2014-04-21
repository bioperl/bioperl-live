#
# BioPerl module for Bio::Seq::EncodedSeq
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Aaron Mackey
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::EncodedSeq - subtype of L<Bio::LocatableSeq|Bio::LocatableSeq> to store DNA that encodes a protein

=head1 SYNOPSIS

  $obj = Bio::Seq::EncodedSeq->new( -seq      => $dna,
                                    -encoding => "CCCCCCCIIIIICCCCC",
                                    -start    => 1,
                                    -strand   => 1,
                                    -length   => 17 );

  # splice out (and possibly revcomp) the coding sequence
  $cds = obj->cds;

  # obtain the protein translation of the sequence
  $prot = $obj->translate;

  # other access/inspection routines as with Bio::LocatableSeq and
  # Bio::SeqI; note that coordinates are relative only to the DNA
  # sequence, not it's implicit encoded protein sequence.

=head1 DESCRIPTION

Bio::Seq::EncodedSeq is a L<Bio::LocatableSeq|Bio::LocatableSeq>
object that holds a DNA sequence as well as information about the
coding potential of that DNA sequence.  It is meant to be useful in an
alignment context, where the DNA may contain frameshifts, gaps and/or
introns, or in describing the transcript of a gene.  An EncodedSeq
provides the ability to access the "spliced" coding sequence, meaning
that all introns and gaps are removed, and any frameshifts are
adjusted to provide a "clean" CDS.

In order to make simultaneous use of either the DNA or the implicit
encoded protein sequence coordinates, please see
L<Bio::Coordinate::EncodingPair>.

=head1 ENCODING

We use the term "encoding" here to refer to the series of symbols that
we use to identify which residues of a DNA sequence are protein-coding
(i.e. part of a codon), intronic, part of a 5' or 3', frameshift
"mutations", etc.  From this information, a Bio::Seq::EncodedSeq is
able to "figure out" its translational CDS.  There are two sets of
coding characters, one termed "implicit" and one termed "explicit".

The "implicit" encoding is a bit simpler than the "explicit" encoding:
'C' is used for any nucleotide that's part of a codon, 'U' for any
UTR, etc.  The full list is shown below:

 Code  Meaning
 ----  -------
  C    coding
  I    intronic
  U    untranslated
  G    gapped (for use in alignments)
  F    a "forward", +1 frameshift
  B    a "backward", -1 frameshift

The "explicit" encoding is just an expansion of the "implicit"
encoding, to denote phase:

 Code  Meaning
 ----  -------
  C    coding, 1st codon position
  D    coding, 2nd codon position
  E    coding, 3rd codon position

  I    intronic, phase 0 (relative to intron begin)
  J    intronic, phase 1
  K    intronic, phase 2

  U    untranslated 3'UTR
  V    untranslated 5'UTR

  G    gapped (for use in alignments)
  F    a "forward", +1 frameshift
  B    a "backward", -1 frameshift

Note that the explicit coding is meant to provide easy access to
position/phase specific nucleotides:

  $obj = Bio::Seq::EncodedSeq->new(-seq => "ACAATCAGACTACG...",
                                   -encoding => "CCCCCCIII..."
                                  );

  # fetch arrays of nucleotides at each codon position:
  my @pos1 = $obj->dnaseq(encoding => 'C', explicit => 1);
  my @pos2 = $obj->dnaseq(encoding => 'D');
  my @pos3 = $obj->dnaseq(encoding => 'E');

  # fetch arrays of "3-1" codon dinucleotides, useful for genomic
  # signature analyses without compounding influences of codon bias:
  my @pairs = $obj->dnaseq(encoding => 'EC');

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

=head1 AUTHOR - Aaron Mackey

Email amackey@virginia.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut



package Bio::Seq::EncodedSeq;

use strict;

use base qw(Bio::LocatableSeq);


=head2 new

 Title   : new
 Usage   : $obj = Bio::Seq::EncodedSeq->new(-seq      => "AGTACGTGTCATG",
                                            -encoding => "CCCCCCFCCCCCC",
                                            -id       => "myseq",
                                            -start    => 1,
                                            -end      => 13,
                                            -strand   => 1
                                      );
 Function: creates a new Bio::Seq::EncodedSeq object from a supplied DNA
           sequence
 Returns : a new Bio::Seq::EncodedSeq object

 Args    : seq      - primary nucleotide sequence used to encode the
                      protein; note that any positions involved in a
                      gap ('G') or backward frameshift ('B') should
                      have one or more gap characters; if the encoding
                      specifies G or B, but no (or not enough) gap
                      characters exist, *they'll be added*; similary,
                      if there are gap characters without a
                      corresponding G or B encoding, G's will be
                      inserted into the encoding.  This allows some
                      flexibility in specifying your sequence and
                      coding without having to calculate a lot of the
                      encoding for yourself.

           encoding - a string of characters (see Encoding Table)
                      describing backwards frameshifts implied by the
                      encoding but not present in the sequence will be
                      added (as '-'s) to the sequence.  If not
                      supplied, it will be assumed that all positions
                      are coding (C).  Encoding may include either
                      implicit phase encoding characters (i.e. "CCC")
                      and/or explicit encoding characters (i.e. "CDE").
                      Additionally, prefixed numbers may be used to
                      denote repetition (i.e. "27C3I28C").

                      Alternatively, encoding may be a hashref
                      datastructure, with encoding characters as keys
                      and Bio::LocationI objects (or arrayrefs of
                      Bio::LocationI objects) as values, e.g.:

                      { C => [ Bio::Location::Simple->new(1,9),
                               Bio::Location::Simple->new(11,13) ],
                        F => Bio::Location::Simple->new(10,10),
                      } # same as "CCCCCCCCCFCCC"

                      Note that if the location ranges overlap, the
                      behavior of the encoding will be undefined
                      (well, it will be defined, but only according to
                      the order in which the hash keys are read, which
                      is basically undefined ... just don't do that).

           id, start, end, strand - as with Bio::LocatableSeq; note
                      that the coordinates are relative to the
                      encoding DNA sequence, not the implicit protein
                      sequence.  If strand is reversed, then the
                      encoding is assumed to be relative to the
                      reverse strand as well.

=cut

sub new {
    my ($self, @args) = @_;
    $self = $self->SUPER::new(@args, -alphabet => 'dna');
    my ($enc) = $self->_rearrange([qw(ENCODING)], @args);
    # set the real encoding:
    if ($enc) {
        $self->encoding($enc);
    } else {
        $self->_recheck_encoding;
    }
    return $self;
}


=head2 encoding

 Title   : encoding
 Usage   : $obj->encoding("CCCCCC");
           $obj->encoding( -encoding => { I => $location } );
           $enc = $obj->encoding(-explicit => 1);
           $enc = $obj->encoding("CCCCCC", -explicit => 1);
           $enc = $obj->encoding(-location => $location,
                                 -explicit => 1,
                                 -absolute => 1 );
 Function: get/set the objects encoding, either globally or by location(s).
 Returns : the (possibly new) encoding string.
 Args    : encoding - see the encoding argument to the new() function.

           explicit - whether or not to return explicit phase
                      information in the coding (i.e. "CCC" becomes
                      "CDE", "III" becomes "IJK", etc); defaults to 0.

           location - optional; location to get/set the encoding.
                      Defaults to the entire sequence.

           absolute - whether or not the locational elements (either
                      in the encoding hashref or the location
                      argument) are relative to the absolute start/end
                      of the Bio::LocatableSeq, or to the internal,
                      relative coordinate system (beginning at 1);
                      defaults to 0 (i.e. relative)

=cut

sub encoding {
    my ($self, @args) = @_;
    my ($enc, $loc, $exp, $abs) = $self->_rearrange([qw(ENCODING LOCATION EXPLICIT ABSOLUTE)], @args);

    if (!$enc) {
        # do nothing; _recheck_encoding will fix for us, if necessary
    } elsif (ref $enc eq 'HASH') {
        $self->throw( -class => 'Bio::Root::NotImplemented',
                      -text  => "Hashref functionality not yet implemented;\nplease email me if you really need this.");
        #TODO: finish all this
        while (my ($char, $locs) = each %$enc) {
            if (ref $locs eq 'ARRAY') {
            } elsif (UNIVERSAL::isa($locs, "Bio::LocationI")) {
            } else {
                $self->throw("Only a scalar or a ref to a hash; not a ref to a @{{lc ref $enc}}");
            }
        }
    } elsif (! ref $enc) {
        $enc = uc $enc;
        $exp = 1 if (!defined $exp && $enc =~ m/[DEJKV]/o);

        if ($enc =~ m/\d/o) { # numerically "enhanced" encoding
            my $numenc = $enc;
            $enc = '';
            while ($numenc =~ m/\G(\d*)([CDEIJKUVGFB])/g) {
                my ($num, $char) = ($1, $2);
                $num = 1 unless $num;
                $enc .= $char x $num;
            }
        }

        if (defined $exp && $exp == 0 && $enc =~ m/([^CIUGFB])/) {
            $self->throw("Unrecognized character '$1' in implicit encoding");
        } elsif ($enc =~ m/[^CDEIJKUVGFB]/) {
            $self->throw("Unrecognized character '$1' in explicit encoding");
        }

        if ($loc) { # a global location, over which to apply the specified encoding.

            # balk if too many non-gap characters present in encoding:
            my ($ct) = $enc =~ tr/GB/GB/;
            $ct = length($enc) - $ct;
            $self->throw("Location length doesn't match number of coding chars in encoding!")
                if ($loc->location_type eq 'EXACT' &&  $loc->length != $ct);

            my $start = $loc->start;
            my $end = $loc->end;

            # strip any encoding that hangs off the ends of the sequence:
            if ($start < $self->start) {
                my $diff = $self->start - $start;
                $start = $self->start;
                $enc = substr($enc, $diff);
            }
            if ($end > $self->end) {
                my $diff = $end - $self->end;
                $end = $self->end;
                $enc = substr($enc, -$diff);
            }

            my $currenc = $self->{_encoding};
            my $currseq = $self->seq;

            my ($spanstart, $spanend) = ($self->column_from_residue_number($start),
                                         $self->column_from_residue_number($end) );

            if ($currseq) {
                # strip any gaps in sequence spanned by this location:
                ($spanstart, $spanend) = ($spanend, $spanstart) if $self->strand < 0;
                my ($before, $in, $after) = $currseq =~ m/(.{@{[ $spanstart - ($loc->location_type eq 'IN-BETWEEN' ? 0 : 1) ]}})
                                                          (.{@{[ $spanend - $spanstart + ($loc->location_type eq 'IN-BETWEEN' ? -1 : 1) ]}})
                                                          (.*)
                                                         /x;
                $in ||= '';
                $in =~ s/[\.\-]+//g;
                $currseq = ($before||'') . $in. ($after||'');
                # change seq without changing the alphabet
                $self->seq($currseq,$self->alphabet());
            }

            $currenc = reverse $currenc if $self->strand < 0;
            substr($currenc, $spanstart, $spanend - $spanstart + ($loc->location_type eq 'IN-BETWEEN' ? -1 : 1),
                   $self->strand >= 0 ? $enc : reverse $enc);
            $currenc = reverse $currenc if $self->strand < 0;

            $self->{_encoding} = $currenc;
            $self->_recheck_encoding;

            $currenc = $self->{_encoding};
            $currenc = reverse $currenc if $self->strand < 0;
            $enc = substr($currenc, $spanstart, length $enc);
            $enc = reverse $enc if $self->strand < 0;

            return $exp ? $enc: $self->_convert2implicit($enc);

        } else {
            # presume a global redefinition; strip any current gap
            # characters in the sequence so they don't corrupt the
            # encoding
            my $dna = $self->seq;
            $dna =~ s/[\.\-]//g;
            $self->seq($dna, 'dna');
            $self->{_encoding} = $enc;
        }
    } else {
        $self->throw("Only a scalar or a ref to a hash; not a ref to a @{{lc ref $enc}}");
    }

    $self->_recheck_encoding();

    return $exp ? $self->{_encoding} : $self->_convert2implicit($self->{_encoding});
}


sub _convert2implicit {
    my ($self, $enc) = @_;

    $enc =~ s/[DE]/C/g;
    $enc =~ s/[JK]/I/g;
    $enc =~ s/V/U/g;

    return $enc;
}


sub _recheck_encoding {

    my $self = shift;

    my @enc = split //, ($self->{_encoding} || '');

    my @nt = split(//, $self->SUPER::seq);
    @nt = reverse @nt if $self->strand && $self->strand < 0;

    # make sure an encoding exists!
    @enc = ('C') x scalar grep { !/[\.\-]/o } @nt
        unless @enc;

    # check for gaps to be truly present in the sequence
    # and vice versa
    my $i;
    for ($i = 0 ; $i < @nt && $i < @enc ; $i++) {
        if ($nt[$i] =~ /[\.\-]/o && $enc[$i] !~ m/[GB]/o) {
            splice(@enc, $i, 0, 'G');
        } elsif ($nt[$i] !~ /[\.\-]/o && $enc[$i] =~ m/[GB]/o) {
            splice(@nt, $i, 0, '-');
        }
    }
    if ($i < @enc) {
        # extra encoding; presumably all gaps?
        for (  ; $i < @enc ; $i++) {
            if ($enc[$i] =~ m/[GB]/o) {
                push @nt, '-';
            } else {
                $self->throw("Extraneous encoding info: " . join('', @enc[$i..$#enc]));
            }
        }
    } elsif ($i < @nt) {
        for (  ; $i < @nt ; $i++) {
            if ($nt[$i] =~ m/[\.\-]/o) {
                push @enc, 'G';
            } else {
                push @enc, 'C';
            }
        }
    }

    my @cde_array = qw(C D E);
    my @ijk_array = qw(I J K);
    # convert any leftover implicit coding into explicit coding
    my ($Cct, $Ict, $Uct, $Vct, $Vwarned) = (0, 0, 0, 0);
    for ($i = 0 ; $i < @enc ; $i++) {
        if ($enc[$i] =~ m/[CDE]/o) {
            my  $temp_index = $Cct %3;
            $enc[$i] = $cde_array[$temp_index];
            $Cct++; $Ict = 0; $Uct = 1;
            $self->warn("3' untranslated encoding (V) seen prior to other coding symbols")
                if ($Vct && !$Vwarned++);
        } elsif ($enc[$i] =~ m/[IJK]/o) {
            $enc[$i] = $ijk_array[$Ict % 3];
            $Ict++; $Uct = 1;
            $self->warn("3' untranslated encoding (V) seen before other coding symbols")
                if ($Vct && !$Vwarned++);
        } elsif ($enc[$i] =~ m/[UV]/o) {
            if ($Uct == 1) {
                $enc[$i] = 'V';
                $Vct = 1;
            }
        } elsif ($enc[$i] eq 'B') {
            $Cct++; $Ict++
        } elsif ($enc[$i] eq 'G') {
            # gap; leave alone
        }
    }

    @nt = reverse @nt if $self->strand && $self->strand < 0;
    $self->seq(join('', @nt), 'dna');

    $self->{_encoding} = join '', @enc;
}


=head2 cds

 Title   : cds
 Usage   : $cds = $obj->cds(-nogaps => 1);
 Function: obtain the "spliced" DNA sequence, by removing any
           nucleotides that participate in an UTR, forward frameshift
           or intron, and replacing any unknown nucleotide implied by
           a backward frameshift or gap with N's.
 Returns : a Bio::Seq::EncodedSeq object, with an encoding consisting only
           of "CCCC..".
 Args    : nogaps - strip any gap characters (resulting from 'G' or 'B'
           encodings), rather than replacing them with N's.

=cut

sub cds {
    my ($self, @args) = @_;

    my ($nogaps, $loc) = $self->_rearrange([qw(NOGAPS LOCATION)], @args);
    $nogaps = 0 unless defined $nogaps;

    my @nt  = split //, $self->strand < 0 ? $self->revcom->seq : $self->seq;
    my @enc = split //, $self->_convert2implicit($self->{_encoding});

    my ($start, $end) = (0, scalar @nt);

    if ($loc) {
        $start = $loc->start;
        $start++ if $loc->location_type eq 'IN-BETWEEN';
        $start = $self->column_from_residue_number($start);
        $start--;

        $end = $loc->end;
        $end = $self->column_from_residue_number($end);

        ($start, $end) = ($end, $start) if $self->strand < 0;
        $start--;
    }

    for (my $i = $start ; $i < $end ; $i++) {
        if ($enc[$i] eq 'I' || $enc[$i] eq 'U' || $enc[$i] eq 'F') {
            # remove introns, untranslated and forward frameshift nucleotides
            $nt[$i] = undef;
        } elsif ($enc[$i] eq 'G' || $enc[$i] eq 'B') {
            # replace gaps and backward frameshifts with N's, unless asked not to.
            $nt[$i] = $nogaps ? undef : 'N';
        }
    }

    return ($self->can_call_new ? ref($self) : __PACKAGE__)->new(
        -seq      => join('', grep { defined } @nt[$start..--$end]),
        -start    => $self->start,
        -end      => $self->end,
        -strand   => 1,
        -alphabet => 'dna' );
}


=head2 translate

 Title   : translate
 Usage   : $prot = $obj->translate(@args);
 Function: obtain the protein sequence encoded by the underlying DNA
           sequence; same as $obj->cds()->translate(@args).
 Returns : a Bio::PrimarySeq object.
 Args    : same as the translate() function of Bio::PrimarySeqI

=cut

sub translate { shift->cds(-nogaps => 1, @_)->SUPER::translate(@_) };


=head2 protseq

 Title   : seq
 Usage   : $protseq = $obj->protseq();
 Function: obtain the raw protein sequence encoded by the underlying
           DNA sequence; This is the same as calling
           $obj->translate()->seq();
 Returns : a string of single-letter amino acid codes
 Args :    same as the seq() function of Bio::PrimarySeq; note that this
           function may not be used to set the protein sequence; see
           the dnaseq() function for that.

=cut

sub protseq { shift->cds(-nogaps => 1, @_)->SUPER::translate(@_)->seq };


=head2 dnaseq

 Title   : dnaseq
 Usage   : $dnaseq = $obj->dnaseq();
           $obj->dnaseq("ACGTGTCGT", "CCCCCCCCC");
           $obj->dnaseq(-seq      => "ATG",
                        -encoding => "CCC",
                        -location => $loc );
           @introns = $obj->$dnaseq(-encoding => 'I')
 Function: get/set the underlying DNA sequence; will overwrite any
           current DNA and/or encoding information present.
 Returns : a string of single-letter nucleotide codes, including any
           gaps implied by the encoding.
 Args    : seq      - the DNA sequence to be used as a replacement
           encoding - the encoding of the DNA sequence (see the new()
                      constructor); defaults to all 'C' if setting a
                      new DNA sequence.  If no new DNA sequence is
                      being provided, then the encoding is used as a
                      "filter" for which to return fragments of
                      non-overlapping DNA that match the encoding.
           location - optional, the location of the DNA sequence to
                      get/set; defaults to the entire sequence.

=cut

sub dnaseq {
    my ($self, @args) = @_;
    my ($seq, $enc, $loc) = $self->_rearrange([qw(DNASEQ ENCODING LOCATION)], @args);
    return $self;
}


# need to overload this so that we truncate both the seq and the encoding!
sub trunc {
    my ($self, $start, $end) = @_;
    my $new = $self->SUPER::trunc($start, $end);
    $start--;
    my $enc = $self->{_encoding};
    $enc = reverse $enc if $self->strand < 0;
    $enc = substr($enc, $start, $end - $start);
    $enc = reverse $enc if $self->strand < 0;
    $new->encoding($enc);
    return $new;
}

1;
