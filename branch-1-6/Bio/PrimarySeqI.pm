# $Id$
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
	   $obj->throw("Could not reverse complement. ".
		    "Probably not DNA. Actual exception\n$@\n");
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

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::PrimarySeqI;
use strict;
use Bio::Tools::CodonTable;

use base qw(Bio::Root::RootI);

=head1 Implementation Specific Functions

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

sub subseq{
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

sub can_call_new{
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

sub alphabet{
    my ( $self ) = @_;
    $self->throw_not_implemented();
}

=head2 moltype

 Title   : moltype
 Usage   : Deprecated. Use alphabet() instead.

=cut

sub moltype{
   my ($self,@args) = @_;

   $self->warn("moltype: pre v1.0 method. Calling alphabet() instead...");
   $self->alphabet(@args);
}


=head1 Optional Implementation Functions

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
           define its own extensionsj.

           To do an inplace edit of an object you can go:

           $seq = $seq->revcom();

           This of course, causes Perl to handle the garbage
           collection of the old object, but it is roughly speaking as
           efficient as an inplace edit.

 Returns : A new (fresh) Bio::PrimarySeqI object
 Args    : None


=cut

sub revcom{
   my ($self) = @_;


   # check the type is good first.
   my $t = $self->alphabet;

   if( $t eq 'protein' ) {
       $self->throw("Sequence is a protein. Cannot revcom");
   }

   if( $t ne 'dna' && $t ne 'rna' ) {
       if( $self->can('warn') ) {
	   $self->warn("Sequence is not dna or rna, but [$t]. ".
		       "Attempting to revcom, but unsure if this is right");
       } else {
	   warn("[$self] Sequence is not dna or rna, but [$t]. ".
		"Attempting to revcom, but unsure if this is right");
       }
   }

   # yank out the sequence string

   my $str = $self->seq();

   # if is RNA - map to DNA then map back

   if( $t eq 'rna' ) {
       $str =~ tr/uU/tT/;
   }

   # revcom etc...

   $str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
   my $revseq = CORE::reverse $str;

   if( $t eq 'rna' ) {
       $revseq =~ tr/tT/uU/;
   }

   my $seqclass;
   if($self->can_call_new()) {
       $seqclass = ref($self);
   } else {
       $seqclass = 'Bio::PrimarySeq';
       $self->_attempt_to_load_Seq();
   }
   my $out = $seqclass->new( '-seq' => $revseq,
			     '-is_circular'  => $self->is_circular,
			     '-display_id'  => $self->display_id,
			     '-accession_number' => $self->accession_number,
			     '-alphabet' => $self->alphabet,
			     '-desc' => $self->desc(),
                             '-verbose' => $self->verbose
			     );
   return $out;

}

=head2 trunc

 Title   : trunc
 Usage   : $subseq = $myseq->trunc(10,100);
 Function: Provides a truncation of a sequence.
 Returns : A fresh Bio::PrimarySeqI implementing object.
 Args    : Two integers denoting first and last base of the sub-sequence.


=cut

sub trunc{
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

   my $seqclass;
   if($self->can_call_new()) {
       $seqclass = ref($self);
   } else {
       $seqclass = 'Bio::PrimarySeq';
       $self->_attempt_to_load_Seq();
   }

   my $out = $seqclass->new( '-seq' => $str,
			     '-display_id'  => $self->display_id,
			     '-accession_number' => $self->accession_number,
			     '-alphabet' => $self->alphabet,
			     '-desc' => $self->desc(),
                             '-verbose' => $self->verbose
			     );
   return $out;
}


=head2 translate

 Title   : translate
 Usage   : $protein_seq_obj = $dna_seq_obj->translate

           Or if you expect a complete coding sequence (CDS) translation,
           with inititator at the beginning and terminator at the end:

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
 Args    : -terminator    - character for terminator        default is *
           -unknown       - character for unknown           default is X
           -frame         - frame                           default is 0
           -codontable_id - codon table id                  default is 1
           -complete      - complete CDS expected           default is 0
           -throw         - throw exception if not complete default is 0
           -orf           - find 1st ORF                    default is 0
           -start         - alternative initiation codon
           -codontable    - Bio::Tools::CodonTable object
		   -offset        - offset for fuzzy locations      default is 0

 Notes   : The -start argument only applies when -orf is set to 1. By default
           all initiation codons found in the given codon table are used
           but when "start" is set to some codon this codon will be used
           exclusively as the initiation codon. Note that the default codon
           table (NCBI "Standard") has 3 initiation codons!

           By default translate() translates termination codons to
           the some character (default is *), both internal and trailing
           codons. Setting "-complete" to 1 tells translate() to remove
           the trailing character.

		   -offset is used for seqfeatures which contain the the \codon_start
		   tag and can be set to 1, 2, or 3.  This is the offset by which the
		   sequence translation starts relative to the first base of the
		   feature

For details on codon tables used by translate() see L<Bio::Tools::CodonTable>.

           Deprecated argument set (v. 1.5.1 and prior versions)
           where each argument is an element in an array:

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
	 my ($terminator, $unknown, $frame, $codonTableId, $complete, $throw,
		  $codonTable, $orf, $start_codon, $offset);

	 ## new API with named parameters, post 1.5.1
	 if ($args[0] && $args[0] =~ /^-[A-Z]+/i) {
		 ($terminator, $unknown, $frame, $codonTableId, $complete, $throw,
		  $codonTable, $orf, $start_codon, $offset) =
			 $self->_rearrange([qw(TERMINATOR
								UNKNOWN
								FRAME
								CODONTABLE_ID
								COMPLETE
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

    ## Get a CodonTable, error if custom CodonTable is invalid
    if ($codonTable) {
		 $self->throw("Need a Bio::Tools::CodonTable object, not ". $codonTable)
			unless $codonTable->isa('Bio::Tools::CodonTable');
    } else {
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
	 if ($orf) {
		 $seq = $self->_find_orf($seq,$codonTable,$start_codon);
	 } else {
	 ## use frame, error if frame is not 0, 1 or 2
		 $self->throw("Valid values for frame are 0, 1, or 2, not $frame.")
			unless ($frame == 0 or $frame == 1 or $frame == 2);
		 $seq = substr($seq,$frame);
    }

    ## Translate it
    my $output = $codonTable->translate($seq);
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
		 if ($output =~ /\*/) {
			 $throw && $self->throw("Seq [$id]: Terminator codon inside CDS!");
			 $self->warn("Seq [$id]: Terminator codon inside CDS!");
		 }
		 # if the initiator codon is not ATG, the amino acid needs to be changed to M
		 if ( substr($output,0,1) ne 'M' ) {
			 if ($codonTable->is_start_codon(substr($seq, 0, 3)) ) {
				 $output = 'M'. substr($output,1);
			 }	elsif ($throw) {
				 $self->throw("Seq [$id]: Not using a valid initiator codon!");
			 } else {
				 $self->warn("Seq [$id]: Not using a valid initiator codon!");
			 }
		 }
    }

    my $seqclass;
    if ($self->can_call_new()) {
		 $seqclass = ref($self);
    } else {
		 $seqclass = 'Bio::PrimarySeq';
		 $self->_attempt_to_load_Seq();
    }
    my $out = $seqclass->new( '-seq' => $output,
										'-display_id'  => $self->display_id,
										'-accession_number' => $self->accession_number,
										# is there anything wrong with retaining the
										# description?
										'-desc' => $self->desc(),
										'-alphabet' => 'protein',
                              '-verbose' => $self->verbose
			      );
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

sub  id {
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

sub  length {
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
 Usage   : if( $obj->is_circular) { /Do Something/ }
 Function: Returns true if the molecule is circular
 Returns : Boolean value
 Args    : none

=cut

sub is_circular{
    shift->throw_not_implemented;
}

=head1 Private functions

These are some private functions for the PrimarySeqI interface. You do not
need to implement these functions

=head2 _find_orf

 Title   : _find_orf
 Usage   :
 Function: Finds ORF starting at 1st initiation codon in nucleotide sequence.
           The ORF is not required to have a termination codon.
 Example :
 Returns : A nucleotide sequence or nothing, if no initiation codon is found.
 Args    : Nucleotide sequence, CodonTable object, alternative initiation
           codon (optional).

=cut

sub _find_orf {
	my ($self,$sequence,$codonTable,$start_codon) = @_;

	# find initiation codon and remove leading sequence
	while ($sequence) {
		my $codon = substr($sequence,0,3);
		if ($start_codon) {
			last if ( $codon =~ /$start_codon/i );
		} else {
			last if ($codonTable->is_start_codon($codon));
		}
		$sequence = substr($sequence,1);
	}
	return unless $sequence;

	# find termination codon and remove trailing sequence
	my $len = CORE::length($sequence);
	my $offset = 3;
	while ($offset < $len) {
		my $codon = substr($sequence,$offset,3);
		if ( $codonTable->is_ter_codon($codon) ){
			$sequence = substr($sequence, 0, $offset + 3);
			return $sequence;
		}
		$offset += 3;
	}
	$self->warn("No termination codon found, will translate - sequence:\n$sequence");
	$sequence;
}

=head2 _attempt_to_load_Seq

 Title   : _attempt_to_load_Seq
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


=cut

sub _attempt_to_load_Seq{
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

1;
