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

Bio::PrimarySeqI - Interface definition for a Bio::PrimarySeq

=head1 SYNOPSIS

    # get a Bio::PrimarySeqI compliant object somehow

    # to test this is a seq object

    $obj->isa("Bio::PrimarySeqI") || $obj->throw("$obj does not implement the Bio::PrimarySeqI interface");

    # accessors

    $string    = $obj->seq();
    $substring = $obj->subseq(12,50);
    $display   = $obj->display_id(); # for human display
    $id        = $obj->primary_id(); # unique id for this object, implementation defined
    $unique_key= $obj->accession_number();
                       # unique biological id

    # object manipulation

    eval {
	$rev    = $obj->revcom();
    };
    if( $@ ) {
	$obj->throw("Could not reverse complement. Probably not DNA. Actual exception\n$@\n");
    }

    $trunc = $obj->trunc(12,50);

    # $rev and $trunc are Bio::PrimarySeqI compliant objects



=head1 DESCRIPTION

This object defines an abstract interface to basic sequence
information. PrimarySeq is an object just for the sequence and its
name(s), nothing more. Seq is the larger object complete with
features. There is a pure perl implementation of this in
Bio::PrimarySeq. If you just want to use Bio::PrimarySeq objects, then
please read that module first. This module defines the interface, and
is of more interest to people who want to wrap their own Perl
Objects/RDBs/FileSystems etc in way that they "are" bioperl sequence
objects, even though it is not using Perl to store the sequence etc.


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

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@sanger.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::PrimarySeqI;
use vars qw(@ISA);
use strict;
use Carp;
use Bio::Tools::CodonTable;

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

   if( $self->can('throw') ) {
       $self->throw("Bio::PrimarySeqI definition of seq - implementing class did not provide this method");
   } else {
       confess("Bio::PrimarySeqI definition of seq - implementing class did not provide this method");
   }
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

   if( $self->can('throw') ) {
       $self->throw("Bio::PrimarySeqI definition of subseq - implementing class did not provide this method");
   } else {
       confess("Bio::PrimarySeqI definition of subseq - implementing class did not provide this method");
   }

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

   if( $self->can('throw') ) {
       $self->throw("Bio::PrimarySeqI definition of id - implementing class did not provide this method");
   } else {
       confess("Bio::PrimarySeqI definition of id - implementing class did not provide this method");
   }

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

   if( $self->can('throw') ) {
       $self->throw("Bio::PrimarySeqI definition of seq - implementing class did not provide this method");
   } else {
       confess("Bio::PrimarySeqI definition of seq - implementing class did not provide this method");
   }

}



=head2 primary_id

 Title   : primary_id
 Usage   : $unique_implementation_key = $obj->primary_id;
 Function: Returns the unique id for this object in this
           implementation. This allows implementations to manage
           their own object ids in a way the implementaiton can control
           clients can expect one id to map to one object.

           For sequences with no accession number, this method should return
           a stringified memory location.
 Returns : A string
 Args    : None
 Status  : Virtual


=cut

sub primary_id {
   my ($self,@args) = @_;

   if( $self->can('throw') ) {
       $self->throw("Bio::PrimarySeqI definition of seq - implementing class did not provide this method");
   } else {
       confess("Bio::PrimarySeqI definition of seq - implementing class did not provide this method");
   }

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
			   -moltype => 'dna',
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

=head2 moltype

 Title   : moltype
 Usage   : if( $obj->moltype eq 'dna' ) { /Do Something/ }
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

sub moltype{
   my ($self,@args) = @_;

   if( $self->can('throw') ) {
       $self->throw("Bio::PrimarySeqI definition of seq - implementing class did not provide this method");
   } else {
       confess("Bio::PrimarySeqI definition of seq - implementing class did not provide this method");
   }


}

=head1 Optional Implementation Functions

The following functions rely on the above functions. A implementing
class does not need to provide these functions, as they will be
provided by this class, but is free to override these functions.

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
           sequences this throws an exception of "Sequence is a protein. Cannot revcom"

           The id is the same id as the orginal sequence, and the accession number
           is also indentical. If someone wants to track that this sequence has be
           reversed, it needs to define its own extensions

           To do an inplace edit of an object you can go:

           $seq = $seq->revcom();

           This of course, causes Perl to handle the garbage collection of the old
           object, but it is roughly speaking as efficient as an inplace edit.

 Returns : A new (fresh) Bio::PrimarySeqI object
 Args    : none


=cut

sub revcom{
   my ($self) = @_;


   # check the type is good first.
   my $t = $self->moltype;

   if( $t eq 'protein' ) {
       if( $self->can('throw') ) {
	   $self->throw("Sequence is a protein. Cannot revcom");
       } else {
	   confess("[$self] Sequence is a protein. Cannot revcom");
       }
   }

   if( $t ne 'dna' && $t ne 'rna' ) {
       if( $self->can('warn') ) {
	   $self->warn("Sequence is not dna or rna, but [$t]. Attempting to revcom, but unsure if this is right");
       } else {
	   warn("[$self] Sequence is not dna or rna, but [$t]. Attempting to revcom, but unsure if this is right");
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
			     '-moltype' => $self->moltype,
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

   if( !$end ) {
       if( $self->can('throw')  ) {
	   $self->throw("trunc start,end");
       } else {
	   confess("[$self] trunc start,end");
       }
   }

   if( $end < $start ) {
       if( $self->can('throw')  ) {
	   $self->throw("$end is smaller than $start. if you want to truncated and reverse complement, you must call trunc followed by revcom. Sorry.");
       } else {
	   confess("[$self] $end is smaller than $start. If you want to truncated and reverse complement, you must call trunc followed by revcom. Sorry.");
       }
   }

   my $str = $self->subseq($start,$end);

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
			     '-moltype' => $self->moltype,
			     '-desc' => $self->desc()
			     );
   return $out;
}


=head2 translate

 Title   : translate
 Usage   : $protein_seq_obj = $dna_seq_obj->translate
           #if full CDS expected:
           $protein_seq_obj = $cds_seq_obj->translate(undef,undef,undef,undef,1);
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
    $self->throw("Can't translate an amino acid sequence.") if
	($self->moltype eq 'protein');
    
    ##Error if frame is not 0, 1 or 2
    $self->throw("Valid values for frame are 0, 1, 2, not [$frame].") unless
	($frame == 0 or $frame == 1 or $frame == 2);
    
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
	    $throw && $self->throw("Seq [$id]: Not using a valid terminator codon!");
	    $self->warn("Seq [$id]: Not using a valid terminator codon!");
	}
	# test if there are terminator characters inside the protein sequence!
	if ($output =~ /\*/) {
	    $throw && $self->throw("Seq [$id]: Terminator codon inside CDS!");
	    $self->warn("Seq [$id]: Terminator codon inside CDS!");
	}
	# if the initiator codon is not ATG, the amino acid needs to changed into M
	if ( substr($output,0,1) ne 'M' ) {
	    if ($codonTable->is_start_codon(substr($seq, 0, 3)) ) {
		$output = 'M'. substr($output,1);
	    }
	    elsif ($throw) {
		$self->warn("Seq [$id]: Not using a valid initiator codon!");
	    } else {
		$self->throw("Seq [$id]: Not using a valid initiator codon!");
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
			      '-moltype' => 'protein'
			      );
    return $out;
    
}

sub translate_old {
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
  $self->throw("Can't translate an amino acid sequence.") if
      ($self->moltype eq 'protein');

  ##Error if frame is not 0, 1 or 2
  $self->throw("Valid values for frame are 0, 1, 2, not [$frame].") unless
      ($frame == 0 or $frame == 1 or $frame == 2);

  #thows a warning if ID is invalid
  my $codonTable = Bio::Tools::CodonTable->new( -id => $tableid);

  my ($seq) = $self->seq();

  # deal with frame offset.
  if( $frame ) {
      $seq = substr ($seq,$frame);
  }
  $seq .= 'n';
  my $length = (length $seq) - 2;
  for ($i = 0 ; $i < $length ; $i += 3)  {
      my $codon = substr($seq, $i, 3);
      my $aa = $codonTable->translate_old($codon);
      if ($aa eq '*') {
   	   $output .= $stop;
      }
      elsif ($aa eq 'X') {
   	   $output .= $unknown;
      }
      else {
	  $output .= $aa ;
      }
  }
  # only if we are expecting to translate a complete coding region
  if ($fullCDS) {
      my $id = $self->display_id;
      #remove the stop character
      if( substr($output,-1,1) eq $stop ) {
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
      # if the initiator codon is not ATG, the amino acid needs to changed into M
      if ( substr($output,0,1) ne 'M' ) {
	  if ($codonTable->is_start_codon(substr($seq, 0, 3)) ) {
	      $output = 'M'. substr($output,1);
	  }
	  elsif ($throw) {
	      $self->warn("Seq [$id]: Not using a valid initiator codon!");
	  } else {
	      $self->throw("Seq [$id]: Not using a valid initiator codon!");
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
			    '-moltype' => 'protein'
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

   if( $self->can('throw') ) {
       $self->throw("Bio::PrimarySeqI definition of length - implementing class did not provide this method");
   } else {
       confess("Bio::PrimarySeqI definition of length - implementing class did not provide this method");
   }
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
   if( $self->can('warn') ) {
       $self->warn("Bio::PrimarySeqI definition of desc - implementing class did not provide this method");
   } else {
       warn("Bio::PrimarySeqI definition of desc - implementing class did not provide this method");
   }
   return '';
}

#  These methods are here for backward compatibility with the old, 0.5
#  Seq objects. They all throw warnings that someone is using a
#  deprecated method, and may eventually be removed completely from
#  this object. However, they are important to ease the transition from
#  the old system.

sub str{
   my ($self,$start,$end) = @_;

   # we assumme anyone using this is using vanilla bioperl object
   my ($p,$f,$l) = caller;
   $self->warn("$f:$l Seq::str - deprecated method. You should use \$obj->seq in preference");

   if( defined $end ) {
       return $self->subseq($start,$end);
   } else {
       return $self->seq();
   }
}

sub ary{
   my ($self,$start,$end) = @_;

   # we assumme anyone using this is using vanilla bioperl object
   my ($p,$f,$l) = caller;
   $self->warn("$f:$l Seq::ary - deprecated method. You should use \$obj->seq in preference, followed by your split to an array");

   my $str;
   if( defined $end ) {
       $str = $self->subseq($start,$end);
   } else {
      $str = $self->seq();
   }

   return split(//,$str);
}

sub getseq{
   my ($self,@args) = @_;

   if( wantarray ) {
       return $self->ary(@args);
   } else {
       return $self->str(@args);
   }
}

sub setseq {
   my ($self,$seq) = @_;

   # we assumme anyone using this is using vanilla bioperl object
   my ($p,$f,$l) = caller;
   $self->warn("$f:$l Seq::setseq - deprecated method. You should use \$obj->seq in preference, followed by your split to an array");

   return $self->seq($seq);
}

sub type{
   my ($self) = @_;

   # we assumme anyone using this is using vanilla bioperl object
   my ($p,$f,$l) = caller;
   $self->warn("$f:$l Seq::type - deprecated method. You should use \$obj->moltype in preference (notice that moltype returns lowercase strings)");

   my $t = $self->moltype;
   $t eq "dna" && return "DNA";
   $t eq "rna" && return "RNA";
   $t eq "protein" && return "PROTEIN";
   return "UNKNOWN";
}

sub seq_len {
    my $self = shift;
   # we assumme anyone using this is using vanilla bioperl object
    my ($p,$f,$l) = caller;
    $self->warn("$f:$l Seq::seq_len - deprecated method. You should use \$obj->length in preference");
    return $self->length();
}

sub out_fasta{
   my ($self,@args) = @_;

   my ($p,$f,$l) = caller;
   $self->warn("$f:$l Seq::out_fasta - deprecated method. You should use the SeqIO package in preference");

   my $str = $self->seq;
   $str =~ tr/a-z/A-Z/;
   $str=~ s/(.{1,60})/$1\n/g;
#  return ">". $self->id(). " ".$self->desc()."\n";
   return ">". $self->id(). "\n" . $str . "\n";  #ps 3/25/00
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
	   if( $self->can('throw') ) {
	       $self->throw("Bio::PrimarySeq could not be loaded for $self\nThis indicates that you are usnig Bio::PrimarySeqI without Bio::PrimarySeq loaded and without providing a complete solution\nThe most likely problem is that there has been a misconfiguration of the bioperl environment\nActual exception\n\n$@\n");
	   } else {
	       confess("Bio::PrimarySeq could not be loaded for $self\nThis indicates that you are usnig Bio::PrimarySeqI without Bio::PrimarySeq loaded and without providing a complete solution\nThe most likely problem is that there has been a misconfiguration of the bioperl environment\nActual exception\n\n$@\n");
	   }
	   return 0;
       }
       return 1;
   }

}




1;


