# $Id$
#
# BioPerl module for Bio::PrimarySeqI
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::PrimarySeqI [Developers] - Interface definition for a Bio::PrimarySeq

=head1 SYNOPSIS


    # Bio::PrimarySeqI is the interface class for sequences.

    # If you are a newcomer to bioperl, you should
    # start with Bio::Seq documentation. This
    # documentation is mainly for developers using
    # Bioperl.

    # to test this is a seq object

    $obj->isa("Bio::PrimarySeqI") ||
      $obj->throw("$obj does not implement the Bio::PrimarySeqI interface");

    # accessors

    $string    = $obj->seq();
    $substring = $obj->subseq(12,50);
    $display   = $obj->display_id(); # for human display
    $id        = $obj->primary_id(); # unique id for this object,
                                     # implementation defined
    $unique_key= $obj->accession_number();
                       # unique biological id

    # object manipulation

    eval {
	$rev    = $obj->revcom();
    };
    if( $@ ) {
	$obj->throw(-class => 'Bio::Root::Exception',
                    -text => "Could not reverse complement. ".
		             "Probably not DNA. Actual exception\n$@\n",
                    -value => $@);
    }

    $trunc = $obj->trunc(12,50);

    # $rev and $trunc are Bio::PrimarySeqI compliant objects


=head1 DESCRIPTION

This object defines an abstract interface to basic sequence
information - for most users of the package the documentation (and
methods) in this class are not useful - this is a developers only
class which defines what methods have to be implmented by other Perl
objects to comply to the Bio::PrimarySeqI interface. Go "perldoc
Bio::Seq" or "man Bio::Seq" for more information on the main class for
sequences.


PrimarySeq is an object just for the sequence and its name(s), nothing
more. Seq is the larger object complete with features. There is a pure
perl implementation of this in Bio::PrimarySeq. If you just want to
use Bio::PrimarySeq objects, then please read that module first. This
module defines the interface, and is of more interest to people who
want to wrap their own Perl Objects/RDBs/FileSystems etc in way that
they "are" bioperl sequence objects, even though it is not using Perl
to store the sequence etc.


This interface defines what bioperl consideres necessary to "be" a
sequence, without providing an implementation of this. (An
implementation is provided in Bio::PrimarySeq). If you want to provide
a Bio::PrimarySeq 'compliant' object which in fact wraps another
object/database/out-of-perl experience, then this is the correct thing
to wrap, generally by providing a wrapper class which would inheriet
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

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Ewan Birney

Email birney@sanger.ac.uk

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
 Usage   : $string    = $obj->seq()
 Function: Returns the sequence as a string of letters. The
           case of the letters is left up to the implementer.
           Suggested cases are upper case for proteins and lower case for
           DNA sequence (IUPAC standard),
           but implementations are suggested to keep an open mind about
           case (some users... want mixed case!)
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
 Function: returns the subseq from start to end, where the first base
           is 1 and the number is inclusive, ie 1-2 are the first two
           bases of the sequence

           Start cannot be larger than end but can be equal

 Returns : a string
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
 Function: returns the display id, aka the common name of the Sequence object.

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
           legacy/convience issues
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

           [Note this method name is likely to change in 1.3]

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
 Function: can_call_new returns 1 or 0 depending
           on whether an implementation allows new
           constructor to be called. If a new constructor
           is allowed, then it should take the followed hashed
           constructor list.

           $myobject->new( -seq => $sequence_as_string,
			   -display_id  => $id
			   -accession_number => $accession
			   -alphabet => 'dna',
			   );
 Example :
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

           This is not called <type> because this would cause
           upgrade problems from the 0.5 and earlier Seq objects.

 Returns : a string either 'dna','rna','protein'. NB - the object must
           make a call of the type - if there is no type specified it
           has to guess.
 Args    : none
 Status  : Virtual


=cut

sub alphabet{
    my ( $self ) = @_;
    $self->throw_not_implemented();
}

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

All of revcom(), trunc(), and translate() create new sequence
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
           protein. Cannot revcom"

           The id is the same id as the original sequence, and the
           accession number is also indentical. If someone wants to
           track that this sequence has be reversed, it needs to
           define its own extensions

           To do an inplace edit of an object you can go:

           $seq = $seq->revcom();

           This of course, causes Perl to handle the garbage
           collection of the old object, but it is roughly speaking as
           efficient as an inplace edit.

 Returns : A new (fresh) Bio::PrimarySeqI object
 Args    : none


=cut

sub revcom{
   my ($self) = @_;


   # check the type is good first.
   my $t = $self->alphabet;

   if( $t eq 'protein' ) {
       $self->throw(-class => 'Bio::Root::Exception',
		    -text  => "Sequence is a protein. Cannot revcom");
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
			     '-display_id'  => $self->display_id,
			     '-accession_number' => $self->accession_number,
			     '-alphabet' => $self->alphabet,
			     '-desc' => $self->desc()
			     );
   return $out;

}

=head2 trunc

 Title   : trunc
 Usage   : $subseq = $myseq->trunc(10,100);
 Function: Provides a truncation of a sequence,

 Example :
 Returns : a fresh Bio::PrimarySeqI implementing object
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
			     '-desc' => $self->desc()
			     );
   return $out;
}


=head2 translate

 Title   : translate
 Usage   : $protein_seq_obj = $dna_seq_obj->translate

 Function:

           Provides the translation of the DNA sequence using full
           IUPAC ambiguities in DNA/RNA and amino acid codes.

           The full CDS translation is identical to EMBL/TREMBL
           database translation. Note that the trailing terminator
           character is removed before returning the translation
           object.

           Note: if you set $dna_seq_obj->verbose(1) you will get a
           warning if the first codon is not a valid initiator.


 Returns : A Bio::PrimarySeqI implementing object
 Args    : character for terminator (optional) defaults to '*'
           character for unknown amino acid (optional) defaults to 'X'
           frame (optional) valid values 0, 1, 2, defaults to 0
           codon table id (optional) defaults to 1
           complete coding sequence expected, defaults to 0 (false)
           boolean, throw exception if not complete CDS (true) or defaults to warning (false)

=cut

sub translate {
    my($self) = shift;
    my($stop, $unknown, $frame, $tableid, $fullCDS, $throw) = @_;
    my($i, $len, $output) = (0,0,'');
    my($codon)   = "";
    my $aa;

    ## User can pass in symbol for stop and unknown codons
    unless(defined($stop) and $stop ne '')    { $stop = "*"; }
    unless(defined($unknown) and $unknown ne '') { $unknown = "X"; }
    unless(defined($frame) and $frame ne '') { $frame = 0; }

    ## the codon table ID
    unless(defined($tableid) and $tableid ne '')    { $tableid = 1; }

    ##Error if monomer is "Amino"
    if ($self->alphabet eq 'protein') {
	$self->throw(-class => 'Bio::Root::Exception',
		     -text => "Can't translate an amino acid sequence.")
    }

    ##Error if frame is not 0, 1 or 2
    unless ($frame == 0 or $frame == 1 or $frame == 2) {
	$self->throw(-class => 'Bio::Root::BadParameter',
		     -text   => "Valid values for frame are 0, 1, 2, not [$frame].",
		     -value => $frame);
    }

    #warns if ID is invalid
    my $codonTable = Bio::Tools::CodonTable->new( -id => $tableid);

    my ($seq) = $self->seq();

    # deal with frame offset.
    if( $frame ) {
	$seq = substr ($seq,$frame);
    }

    # Translate it
    $output = $codonTable->translate($seq);
    # Use user-input stop/unknown
    $output =~ s/\*/$stop/g;
    $output =~ s/X/$unknown/g;
	
    # only if we are expecting to translate a complete coding region
    if ($fullCDS) {
	my $id = $self->display_id;
	#remove the stop character
	if( substr($output,-1,1) eq $stop ) {
	    chop $output;
	} else {
	    $throw && $self->throw(-class => 'Bio::Root::Exception',
				   -text  => "Seq [$id]: Not using a valid terminator codon!: ". substr($output,-1,1),
				   -value => substr($output,-1,1));
	    $self->warn("Seq [$id]: Not using a valid terminator codon!: ". substr($output,-1,1));
	}
	# test if there are terminator characters inside the protein sequence!
	if ($output =~ /\*/) {
	    $throw && $self->throw(-class => 'Bio::Root::Exception',
				   -text  => "Seq [$id]: Terminator codon inside CDS!");
	    $self->warn("Seq [$id]: Terminator codon inside CDS!");
	}
	# if the initiator codon is not ATG, the amino acid needs to changed into M
	if ( substr($output,0,1) ne 'M' ) {
	    if ($codonTable->is_start_codon(substr($seq, 0, 3)) ) {
		$output = 'M'. substr($output,1);
	    }
	    elsif ($throw) {
		$self->throw(-class => 'Bio::Root::Exception',
			     -text => "Seq [$id]: Not using a valid initiator codon!: ". substr($seq, 0, 3),
			     -value => substr($seq, 0, 3));
	    } else {
		$self->warn("Seq [$id]: Not using a valid initiator codon!");
	    }
	}
    }

    my $seqclass;
    if($self->can_call_new()) {
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
			      '-alphabet' => 'protein'
			    );
    return $out;

}

=head2 id

 Title   : id
 Usage   : $id = $seq->id()
 Function: ID of the sequence. This should normally be (and actually is in
           the implementation provided here) just a synonym for display_id().
 Example :
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
 Example :
 Returns : integer representing the length of the sequence.
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
 Example :
 Returns : value of desc
 Args    : newvalue (optional)


=cut

sub desc {
   my ($self,$value) = @_;
   $self->throw_not_implemented();
}


=head2 is_circular

 Title   : is_circular
 Usage   : if( $obj->is_circular) { /Do Something/ }
 Function: Returns true if the molecule is circular
 Returns : Boolean value
 Args    : none

=cut

sub is_circular{
    my ($self,$value) = @_;
    if (defined $value) {
	$self->{'_is_circular'} = $value;
    }
    return $self->{'_is_circular'};
}

=head1 Private functions

These are some private functions for the PrimarySeqI interface. You do not
need to implement these functions

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
	   $self->throw(-class => 'Bio::Root::Exception',
			-text => "$text$@\n",
			-value => $@);
	   return 0;
       }
       return 1;
   }

}

1;
