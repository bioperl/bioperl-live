#
# BioPerl module for Bio::PrimarySeqI
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

Bio::PrimarySeqI - Interface definition for a Bio::PrimarySeq

=head1 SYNOPSIS

    # Bio::PrimarySeqI is the interface class for sequences.
    # If you are a newcomer to bioperl, you might want to start with
    # Bio::Seq documentation.

    # Test if this is a seq object
    $obj->isa("Bio::PrimarySeqI") ||
      $obj->throw("$obj does not implement the Bio::PrimarySeqI interface");

    # Accessors
    $string    = $obj->seq();
    $substring = $obj->subseq(12,50);
    $display   = $obj->display_id();       # for human display
    $id        = $obj->primary_id();       # unique id for this object,
                                           # implementation defined
    $unique_key= $obj->accession_number(); # unique biological id


    # Object manipulation
    eval {
       $rev = $obj->revcom();
    };
    if( $@ ) {
       $obj->throw( "Could not reverse complement. ".
                    "Probably not DNA. Actual exception\n$@\n" );
    }

    $trunc = $obj->trunc(12,50);
    # $rev and $trunc are Bio::PrimarySeqI compliant objects


=head1 DESCRIPTION

This object defines an abstract interface to basic sequence
information - for most users of the package the documentation (and
methods) in this class are not useful - this is a developers-only
class which defines what methods have to be implmented by other Perl
objects to comply to the Bio::PrimarySeqI interface. Go "perldoc
Bio::Seq" or "man Bio::Seq" for more information on the main class for
sequences.

PrimarySeq is an object just for the sequence and its name(s), nothing
more. Seq is the larger object complete with features. There is a pure
perl implementation of this in L<Bio::PrimarySeq>. If you just want to
use L<Bio::PrimarySeq> objects, then please read that module first. This
module defines the interface, and is of more interest to people who
want to wrap their own Perl Objects/RDBs/FileSystems etc in way that
they "are" bioperl sequence objects, even though it is not using Perl
to store the sequence etc.

This interface defines what bioperl considers necessary to "be" a
sequence, without providing an implementation of this, an
implementation is provided in L<Bio::PrimarySeq>. If you want to provide
a Bio::PrimarySeq-compliant object which in fact wraps another
object/database/out-of-perl experience, then this is the correct thing
to wrap, generally by providing a wrapper class which would inherit
from your object and this Bio::PrimarySeqI interface. The wrapper class
then would have methods lists in the "Implementation Specific
Functions" which would provide these methods for your object.

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

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


package Bio::PrimarySeqI;
use strict;
use Bio::Tools::CodonTable;

use base qw(Bio::Root::RootI);


=head1 Implementation-specific Functions

These functions are the ones that a specific implementation must
define.

=head2 seq

 Title   : seq
 Usage   : $string = $obj->seq()
 Function: Returns the sequence as a string of letters. The
           case of the letters is left up to the implementer.
           Suggested cases are upper case for proteins and lower case for
           DNA sequence (IUPAC standard), but implementations are suggested to
           keep an open mind about case (some users... want mixed case!)
 Returns : A scalar
 Status  : Virtual

=cut

sub seq {
    my ($self) = @_;
    $self->throw_not_implemented();
}


=head2 subseq

 Title   : subseq
 Usage   : $substring = $obj->subseq(10,40);
 Function: Returns the subseq from start to end, where the first base
           is 1 and the number is inclusive, i.e. 1-2 are the first two
           bases of the sequence.

           Start cannot be larger than end but can be equal.

 Returns : A string
 Args    :
 Status  : Virtual

=cut

sub subseq {
    my ($self) = @_;
    $self->throw_not_implemented();
}


=head2 display_id

 Title   : display_id
 Usage   : $id_string = $obj->display_id();
 Function: Returns the display id, also known as the common name of the Sequence
           object.

           The semantics of this is that it is the most likely string
           to be used as an identifier of the sequence, and likely to
           have "human" readability.  The id is equivalent to the ID
           field of the GenBank/EMBL databanks and the id field of the
           Swissprot/sptrembl database. In fasta format, the >(\S+) is
           presumed to be the id, though some people overload the id
           to embed other information. Bioperl does not use any
           embedded information in the ID field, and people are
           encouraged to use other mechanisms (accession field for
           example, or extending the sequence object) to solve this.

           Notice that $seq->id() maps to this function, mainly for
           legacy/convenience reasons.
 Returns : A string
 Args    : None
 Status  : Virtual

=cut

sub display_id {
    my ($self) = @_;
    $self->throw_not_implemented();
}


=head2 accession_number

 Title   : accession_number
 Usage   : $unique_biological_key = $obj->accession_number;
 Function: Returns the unique biological id for a sequence, commonly
           called the accession_number. For sequences from established
           databases, the implementors should try to use the correct
           accession number. Notice that primary_id() provides the
           unique id for the implemetation, allowing multiple objects
           to have the same accession number in a particular implementation.

           For sequences with no accession number, this method should return
           "unknown".
 Returns : A string
 Args    : None
 Status  : Virtual

=cut

sub accession_number {
    my ($self,@args) = @_;
    $self->throw_not_implemented();
}


=head2 primary_id

 Title   : primary_id
 Usage   : $unique_implementation_key = $obj->primary_id;
 Function: Returns the unique id for this object in this
           implementation. This allows implementations to manage their
           own object ids in a way the implementaiton can control
           clients can expect one id to map to one object.

           For sequences with no accession number, this method should
           return a stringified memory location.

 Returns : A string
 Args    : None
 Status  : Virtual

=cut

sub primary_id {
    my ($self,@args) = @_;
    $self->throw_not_implemented();
}


=head2 can_call_new

 Title   : can_call_new
 Usage   : if( $obj->can_call_new ) {
             $newobj = $obj->new( %param );
         }
 Function: Can_call_new returns 1 or 0 depending
           on whether an implementation allows new
           constructor to be called. If a new constructor
           is allowed, then it should take the followed hashed
           constructor list.

           $myobject->new( -seq => $sequence_as_string,
                           -display_id  => $id
                           -accession_number => $accession
                           -alphabet => 'dna',
                           );
 Returns : 1 or 0
 Args    :


=cut

sub can_call_new {
    my ($self,@args) = @_;
    # we default to 0 here
    return 0;
}


=head2 alphabet

 Title   : alphabet
 Usage   : if( $obj->alphabet eq 'dna' ) { /Do Something/ }
 Function: Returns the type of sequence being one of
           'dna', 'rna' or 'protein'. This is case sensitive.

           This is not called "type" because this would cause
           upgrade problems from the 0.5 and earlier Seq objects.

 Returns : A string either 'dna','rna','protein'. NB - the object must
           make a call of the alphabet, if there is no alphabet specified it
           has to guess.
 Args    : None
 Status  : Virtual

=cut

sub alphabet {
    my ( $self ) = @_;
    $self->throw_not_implemented();
}


=head2 moltype

 Title   : moltype
 Usage   : Deprecated. Use alphabet() instead.

=cut

sub moltype {
    my ($self,@args) = @_;
    $self->warn("moltype: pre v1.0 method. Calling alphabet() instead...");
    return $self->alphabet(@args);
}


=head1 Implementation-optional Functions

The following functions rely on the above functions. An
implementing class does not need to provide these functions, as they
will be provided by this class, but is free to override these
functions.

The revcom(), trunc(), and translate() methods create new sequence
objects. They will call new() on the class of the sequence object
instance passed as argument, unless can_call_new() returns FALSE. In
the latter case a Bio::PrimarySeq object will be created. Implementors
which really want to control how objects are created (eg, for object
persistence over a database, or objects in a CORBA framework), they
are encouraged to override these methods

=head2 revcom

 Title   : revcom
 Usage   : $rev = $seq->revcom()
 Function: Produces a new Bio::PrimarySeqI implementing object which
           is the reversed complement of the sequence. For protein
           sequences this throws an exception of "Sequence is a
           protein. Cannot revcom".

           The id is the same id as the original sequence, and the
           accession number is also indentical. If someone wants to
           track that this sequence has be reversed, it needs to
           define its own extensions.

           To do an inplace edit of an object you can go:

           $seq = $seq->revcom();

           This of course, causes Perl to handle the garbage
           collection of the old object, but it is roughly speaking as
           efficient as an inplace edit.

 Returns : A new (fresh) Bio::PrimarySeqI object
 Args    : None


=cut

sub revcom {
    my ($self) = @_;

    # Create a new fresh object if $self is 'Bio::Seq::LargePrimarySeq'
    # or 'Bio::Seq::LargeSeq', if not take advantage of
    # Bio::Root::clone to get an object copy
    my $out;
    if (   $self->isa('Bio::Seq::LargePrimarySeq')
        or $self->isa('Bio::Seq::LargeSeq')
        ) {
        my ($seqclass, $opts) = $self->_setup_class;
        $out = $seqclass->new(
            -seq              => $self->_revcom_from_string($self->seq, $self->alphabet),
            -is_circular      => $self->is_circular,
            -display_id       => $self->display_id,
            -accession_number => $self->accession_number,
            -alphabet         => $self->alphabet,
            -desc             => $self->desc,
            -verbose          => $self->verbose,
            %$opts,
        );
    } else {
        $out = $self->clone;
        $out->seq( $out->_revcom_from_string($out->seq, $out->alphabet) );
    }
    return $out;
}


sub _revcom_from_string {
    my ($self, $string, $alphabet) = @_;

    # Check that reverse-complementing makes sense
    if( $alphabet eq 'protein' ) {
        $self->throw("Sequence is a protein. Cannot revcom.");
    }
    if( $alphabet ne 'dna' && $alphabet ne 'rna' ) {
        my $msg = "Sequence is not dna or rna, but [$alphabet]. Attempting to revcom, ".
                  "but unsure if this is right.";
        if( $self->can('warn') ) {
            $self->warn($msg);
        } else {
            warn("[$self] $msg");
        }
    }

    # If sequence is RNA, map to DNA (then map back later)
    if( $alphabet eq 'rna' ) {
        $string =~ tr/uU/tT/;
    }

    # Reverse-complement now
    $string =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    $string = CORE::reverse $string;

    # Map back RNA to DNA
    if( $alphabet eq 'rna' ) {
        $string =~ tr/tT/uU/;
    }

    return $string;
}


=head2 trunc

 Title   : trunc
 Usage   : $subseq = $myseq->trunc(10,100);
 Function: Provides a truncation of a sequence.
 Returns : A fresh Bio::PrimarySeqI implementing object.
 Args    : Two integers denoting first and last base of the sub-sequence.


=cut

sub trunc {
    my ($self,$start,$end) = @_;

    my $str;
    if( defined $start && ref($start) &&
        $start->isa('Bio::LocationI') ) {
        $str = $self->subseq($start); # start is a location actually
    } elsif( !$end ) {
        $self->throw("trunc start,end -- there was no end for $start");
    } elsif( $end < $start ) {
        my $msg = "start [$start] is greater than end [$end]. \n".
                  "If you want to truncated and reverse complement, \n".
                  "you must call trunc followed by revcom. Sorry.";
        $self->throw($msg);
    } else {
        $str = $self->subseq($start,$end);
    }

    # Create a new fresh object if $self is 'Bio::Seq::LargePrimarySeq'
    # or 'Bio::Seq::LargeSeq', if not take advantage of
    # Bio::Root::clone to get an object copy
    my $out;
    if (   $self->isa('Bio::Seq::LargePrimarySeq')
        or $self->isa('Bio::Seq::LargeSeq')
        ) {
        my ($seqclass, $opts) = $self->_setup_class;
        $out = $seqclass->new(
            -seq              => $str,
            -is_circular      => $self->is_circular,
            -display_id       => $self->display_id,
            -accession_number => $self->accession_number,
            -alphabet         => $self->alphabet,
            -desc             => $self->desc,
            -verbose          => $self->verbose,
            %$opts,
        );
    } else {
        $out = $self->clone;
        $out->seq($str);
    }
    return $out;
}


=head2 translate

 Title   : translate
 Usage   : $protein_seq_obj = $dna_seq_obj->translate

           Or if you expect a complete coding sequence (CDS) translation,
           with initiator at the beginning and terminator at the end:

           $protein_seq_obj = $cds_seq_obj->translate(-complete => 1);

           Or if you want translate() to find the first initiation
           codon and return the corresponding protein:

           $protein_seq_obj = $cds_seq_obj->translate(-orf => 1);

 Function: Provides the translation of the DNA sequence using full
           IUPAC ambiguities in DNA/RNA and amino acid codes.

           The complete CDS translation is identical to EMBL/TREMBL
           database translation. Note that the trailing terminator
           character is removed before returning the translated protein
           object.

           Note: if you set $dna_seq_obj->verbose(1) you will get a
           warning if the first codon is not a valid initiator.

 Returns : A Bio::PrimarySeqI implementing object
 Args    : -terminator
               character for terminator, default '*'
           -unknown
               character for unknown, default 'X'
           -frame
               positive integer frame shift (in bases), default 0
           -codontable_id
               integer codon table id, default 1
           -complete
               boolean, if true, complete CDS is expected. default false
           -complete_codons
               boolean, if true, codons which are incomplete are translated if a
               suitable amino acid is found. For instance, if the incomplete
               codon is 'GG', the completed codon is 'GGN', which is glycine
               (G). Defaults to 'false'; setting '-complete' also makes this
               true.
           -throw
               boolean, throw exception if ORF not complete, default false
           -orf
               if 'longest', find longest ORF.  other true value, find
               first ORF.  default 0
           -codontable
               optional L<Bio::Tools::CodonTable> object to use for
               translation
           -start
               optional three-character string to force as initiation
               codon (e.g. 'atg'). If unset, start codons are
               determined by the CodonTable.  Case insensitive.
           -offset
               optional positive integer offset for fuzzy locations.
               if set, must be either 1, 2, or 3

=head3 Notes

The -start argument only applies when -orf is set to 1. By default all
initiation codons found in the given codon table are used but when
"start" is set to some codon this codon will be used exclusively as
the initiation codon. Note that the default codon table (NCBI
"Standard") has 3 initiation codons!

By default translate() translates termination codons to the some
character (default is *), both internal and trailing codons. Setting
"-complete" to 1 tells translate() to remove the trailing character.

-offset is used for seqfeatures which contain the the \codon_start tag
and can be set to 1, 2, or 3.  This is the offset by which the
sequence translation starts relative to the first base of the feature

For details on codon tables used by translate() see L<Bio::Tools::CodonTable>.

Deprecated argument set (v. 1.5.1 and prior versions) where each argument is an
element in an array:

  1: character for terminator (optional), defaults to '*'.
  2: character for unknown amino acid (optional), defaults to 'X'.
  3: frame (optional), valid values are 0, 1, 2, defaults to 0.
  4: codon table id (optional), defaults to 1.
  5: complete coding sequence expected, defaults to 0 (false).
  6: boolean, throw exception if not complete coding sequence
     (true), defaults to warning (false)
  7: codontable, a custom Bio::Tools::CodonTable object (optional).

=cut

sub translate {
    my ($self,@args) = @_;
    my ($terminator, $unknown, $frame, $codonTableId, $complete,
        $complete_codons, $throw, $codonTable, $orf, $start_codon, $offset);

    ## new API with named parameters, post 1.5.1
    if ($args[0] && $args[0] =~ /^-[A-Z]+/i) {
        ($terminator, $unknown, $frame, $codonTableId, $complete,
        $complete_codons, $throw,$codonTable, $orf, $start_codon, $offset) =
            $self->_rearrange([qw(TERMINATOR
                                               UNKNOWN
                                               FRAME
                                               CODONTABLE_ID
                                               COMPLETE
                                               COMPLETE_CODONS
                                               THROW
                                               CODONTABLE
                                               ORF
                                               START
                                               OFFSET)], @args);
    ## old API, 1.5.1 and preceding versions
    } else {
        ($terminator, $unknown, $frame, $codonTableId,
         $complete, $throw, $codonTable, $offset) = @args;
    }
    
    ## Initialize termination codon, unknown codon, codon table id, frame
    $terminator = '*'    unless (defined($terminator) and $terminator ne '');
    $unknown = "X"       unless (defined($unknown) and $unknown ne '');
    $frame = 0           unless (defined($frame) and $frame ne '');
    $codonTableId = 1    unless (defined($codonTableId) and $codonTableId ne '');
    $complete_codons ||= $complete || 0;
    
    ## Get a CodonTable, error if custom CodonTable is invalid
    if ($codonTable) {
        $self->throw("Need a Bio::Tools::CodonTable object, not ". $codonTable)
            unless $codonTable->isa('Bio::Tools::CodonTable');
    } else {
        
        # shouldn't this be cached?  Seems wasteful to have a new instance
        # every time...
        $codonTable = Bio::Tools::CodonTable->new( -id => $codonTableId);
    }

    ## Error if alphabet is "protein"
    $self->throw("Can't translate an amino acid sequence.") if
        ($self->alphabet =~ /protein/i);

    ## Error if -start parameter isn't a valid codon
    if ($start_codon) {
        $self->throw("Invalid start codon: $start_codon.") if
            ( $start_codon !~ /^[A-Z]{3}$/i );
    }

    my $seq;
    if ($offset) {
        $self->throw("Offset must be 1, 2, or 3.") if
            ( $offset !~ /^[123]$/ );
        my ($start, $end) = ($offset, $self->length);
        ($seq) = $self->subseq($start, $end);
    } else {
        ($seq) = $self->seq();
    }

         ## ignore frame if an ORF is supposed to be found
    if ( $orf ) {
        my ($orf_region) = $self->_find_orfs_nucleotide( $seq, $codonTable, $start_codon, $orf eq 'longest' ? 0 : 'first_only' );
        $seq = $self->_orf_sequence( $seq, $orf_region );
    } else {
        ## use frame, error if frame is not 0, 1 or 2
        $self->throw("Valid values for frame are 0, 1, or 2, not $frame.")
            unless ($frame == 0 or $frame == 1 or $frame == 2);
        $seq = substr($seq,$frame);
    }

    ## Translate it
    my $output = $codonTable->translate($seq, $complete_codons);
    # Use user-input terminator/unknown
    $output =~ s/\*/$terminator/g;
    $output =~ s/X/$unknown/g;

    ## Only if we are expecting to translate a complete coding region
    if ($complete) {
        my $id = $self->display_id;
        # remove the terminator character
        if( substr($output,-1,1) eq $terminator ) {
            chop $output;
        } else {
            $throw && $self->throw("Seq [$id]: Not using a valid terminator codon!");
            $self->warn("Seq [$id]: Not using a valid terminator codon!");
        }
        # test if there are terminator characters inside the protein sequence!
        if ($output =~ /\Q$terminator\E/) {
            $id ||= '';
            $throw && $self->throw("Seq [$id]: Terminator codon inside CDS!");
            $self->warn("Seq [$id]: Terminator codon inside CDS!");
        }
        # if the initiator codon is not ATG, the amino acid needs to be changed to M
        if ( substr($output,0,1) ne 'M' ) {
            if ($codonTable->is_start_codon(substr($seq, 0, 3)) ) {
                $output = 'M'. substr($output,1);
            } elsif ($throw) {
                $self->throw("Seq [$id]: Not using a valid initiator codon!");
            } else {
                $self->warn("Seq [$id]: Not using a valid initiator codon!");
            }
        }
    }

    # Create a new fresh object if $self is 'Bio::Seq::LargePrimarySeq'
    # or 'Bio::Seq::LargeSeq', if not take advantage of
    # Bio::Root::clone to get an object copy
    my $out;
    if (   $self->isa('Bio::Seq::LargePrimarySeq')
        or $self->isa('Bio::Seq::LargeSeq')
        ) {
        my ($seqclass, $opts) = $self->_setup_class;
        $out = $seqclass->new(
            -seq              => $output,
            -is_circular      => $self->is_circular,
            -display_id       => $self->display_id,
            -accession_number => $self->accession_number,
            -alphabet         => 'protein',
            -desc             => $self->desc,
            -verbose          => $self->verbose,
            %$opts,
        );
    } else {
        $out = $self->clone;
        $out->seq($output);
        $out->alphabet('protein');
    }
    return $out;
}


=head2 transcribe()

 Title   : transcribe
 Usage   : $xseq = $seq->transcribe
 Function: Convert base T to base U
 Returns : PrimarySeqI object of alphabet 'rna' or
           undef if $seq->alphabet ne 'dna'
 Args    :

=cut

sub transcribe {
    my $self = shift;
    return unless $self->alphabet eq 'dna';
    my $s = $self->seq;
    $s =~ tr/tT/uU/;
    my $desc = $self->desc || '';

    # Create a new fresh object if $self is 'Bio::Seq::LargePrimarySeq'
    # or 'Bio::Seq::LargeSeq', if not take advantage of
    # Bio::Root::clone to get an object copy
    my $out;
    if (   $self->isa('Bio::Seq::LargePrimarySeq')
        or $self->isa('Bio::Seq::LargeSeq')
        ) {
        my ($seqclass, $opts) = $self->_setup_class;
        $out = $seqclass->new(
            -seq              => $s,
            -is_circular      => $self->is_circular,
            -display_id       => $self->display_id,
            -accession_number => $self->accession_number,
            -alphabet         => 'rna',
            -desc             => "${desc}[TRANSCRIBED]",
            -verbose          => $self->verbose,
            %$opts,
        );
    } else {
        $out = $self->clone;
        $out->seq($s);
        $out->alphabet('rna');
        $out->desc($desc . "[TRANSCRIBED]");
    }
    return $out;
}


=head2 rev_transcribe()

 Title   : rev_transcribe
 Usage   : $rtseq = $seq->rev_transcribe
 Function: Convert base U to base T
 Returns : PrimarySeqI object of alphabet 'dna' or
           undef if $seq->alphabet ne 'rna'
 Args    :

=cut

sub rev_transcribe {
    my $self = shift;
    return unless $self->alphabet eq 'rna';
    my $s = $self->seq;
    $s =~ tr/uU/tT/;
    my $desc = $self->desc || '';

    # Create a new fresh object if $self is 'Bio::Seq::LargePrimarySeq'
    # or 'Bio::Seq::LargeSeq', if not take advantage of
    # Bio::Root::clone to get an object copy
    my $out;
    if (   $self->isa('Bio::Seq::LargePrimarySeq')
        or $self->isa('Bio::Seq::LargeSeq')
        ) {
        my ($seqclass, $opts) = $self->_setup_class;
        $out = $seqclass->new(
            -seq              => $s,
            -is_circular      => $self->is_circular,
            -display_id       => $self->display_id,
            -accession_number => $self->accession_number,
            -alphabet         => 'dna',
            -desc             => $self->desc . "[REVERSE TRANSCRIBED]",
            -verbose          => $self->verbose,
            %$opts,
        );
    } else {
        $out = $self->clone;
        $out->seq($s);
        $out->alphabet('dna');
        $out->desc($desc . "[REVERSE TRANSCRIBED]");
    }
    return $out;
}


=head2 id

 Title   : id
 Usage   : $id = $seq->id()
 Function: ID of the sequence. This should normally be (and actually is in
           the implementation provided here) just a synonym for display_id().
 Returns : A string.
 Args    :

=cut

sub id {
   my ($self)= @_;
   return $self->display_id();
}


=head2 length

 Title   : length
 Usage   : $len = $seq->length()
 Function:
 Returns : Integer representing the length of the sequence.
 Args    :

=cut

sub length {
    my ($self)= @_;
    $self->throw_not_implemented();
}


=head2 desc

 Title   : desc
 Usage   : $seq->desc($newval);
           $description = $seq->desc();
 Function: Get/set description text for a seq object
 Returns : Value of desc
 Args    : newvalue (optional)

=cut

sub desc {
   shift->throw_not_implemented();
}


=head2 is_circular

 Title   : is_circular
 Usage   : if( $obj->is_circular) { # Do something }
 Function: Returns true if the molecule is circular
 Returns : Boolean value
 Args    : none

=cut

sub is_circular {
    shift->throw_not_implemented;
}


=head1 Private functions

These are some private functions for the PrimarySeqI interface. You do not
need to implement these functions

=head2 _find_orfs_nucleotide

 Title   : _find_orfs_nucleotide
 Usage   :
 Function: Finds ORF starting at 1st initiation codon in nucleotide sequence.
           The ORF is not required to have a termination codon.
 Example :
 Returns : a list of string coordinates of ORF locations (0-based half-open),
           sorted descending by length (so that the longest is first)
           as: [ start, end, frame, length ], [ start, end, frame, length ], ...
 Args    : Nucleotide sequence,
           CodonTable object,
           (optional) alternative initiation codon (e.g. 'ATA'),
           (optional) boolean that, if true, stops after finding the
                      first available ORF

=cut

sub _find_orfs_nucleotide {
    my ( $self, $sequence, $codon_table, $start_codon, $first_only ) = @_;
    $sequence    = uc $sequence;
    $start_codon = uc $start_codon if $start_codon;

    my $is_start = $start_codon
        ? sub { shift eq $start_codon }
        : sub { $codon_table->is_start_codon( shift ) };

    # stores the begin index of the currently-running ORF in each
    # reading frame
    my @current_orf_start = (-1,-1,-1);

    #< stores coordinates of longest observed orf (so far) in each
    #  reading frame
    my @orfs;

    # go through each base of the sequence, and each reading frame for each base
    my $seqlen = CORE::length $sequence;
    for( my $j = 0; $j <= $seqlen-3; $j++ ) {
        my $frame = $j % 3;

        my $this_codon = substr( $sequence, $j, 3 );

        # if in an orf and this is either a stop codon or the last in-frame codon in the string
        if ( $current_orf_start[$frame] >= 0 ) {
            if ( $codon_table->is_ter_codon( $this_codon ) ||( my $is_last_codon_in_frame = ($j >= $seqlen-5)) ) {
                # record ORF start, end (half-open), length, and frame
                my @this_orf = ( $current_orf_start[$frame], $j+3, undef, $frame );
                my $this_orf_length = $this_orf[2] = ( $this_orf[1] - $this_orf[0] );

                $self->warn( "Translating partial ORF "
                                 .$self->_truncate_seq( $self->_orf_sequence( $sequence, \@this_orf ))
                                 .' from end of nucleotide sequence'
                            )
                    if $first_only && $is_last_codon_in_frame;

                return \@this_orf if $first_only;
                push @orfs, \@this_orf;
                $current_orf_start[$frame] = -1;
            }
        }
        # if this is a start codon
        elsif ( $is_start->($this_codon) ) {
            $current_orf_start[$frame] = $j;
        }
    }

    return sort { $b->[2] <=> $a->[2] } @orfs;
}


sub _truncate_seq {
    my ($self, $seq) = @_;
    return CORE::length($seq) > 200 ? substr($seq,0,50).'...'.substr($seq,-50) : $seq;
}


sub _orf_sequence {
    my ($self, $seq, $orf ) = @_;
    return '' unless $orf;
    return substr( $seq, $orf->[0], $orf->[2] )
}


=head2 _attempt_to_load_Seq

 Title   : _attempt_to_load_Seq
 Usage   :
 Function:
 Example :
 Returns :
 Args    :

=cut

sub _attempt_to_load_Seq {
    my ($self) = @_;

    if( $main::{'Bio::PrimarySeq'} ) {
        return 1;
    } else {
        eval {
            require Bio::PrimarySeq;
        };
        if( $@ ) {
            my $text = "Bio::PrimarySeq could not be loaded for [$self]\n".
                "This indicates that you are using Bio::PrimarySeqI ".
                "without Bio::PrimarySeq loaded or without providing a ".
                "complete implementation.\nThe most likely problem is that there ".
                "has been a misconfiguration of the bioperl environment\n".
                "Actual exception:\n\n";
            $self->throw("$text$@\n");
            return 0;
        }
        return 1;
    }
}


sub _setup_class {
    # Return name of class and setup some default parameters
    my ($self) = @_;
    my $seqclass;
    if ($self->can_call_new()) {
        $seqclass = ref($self);
    } else {
        $seqclass = 'Bio::PrimarySeq';
        $self->_attempt_to_load_Seq();
    }
    my %opts;
    if ($seqclass eq 'Bio::PrimarySeq') {
        # Since sequence is in a Seq object, it has already been validated.
        # We do not need to validate its trunc(), revcom(), etc
        $opts{ -direct } = 1;
    }
    return $seqclass, \%opts;
}


1;
