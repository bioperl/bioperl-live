
#
# bioperl module for Bio::NewSeq
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq - Bioperl lightweight Sequence Object

=head1 SYNOPSIS

  # The Bio::SeqIO for file reading, Bio::DB::GenBank for
  # database reading
  use Bio::Seq;
  use Bio::SeqIO;
  use Bio::DB::GenBank;

  #make from memory
  $seqobj = Bio::Seq->new ( -seq => 'ATGGGGTGGGCGGTGGGTGGTTTG',
			 -id  => 'GeneFragment-12',
			 -accession => 'X78121',
			 -moltype => 'dna'
			 );
  
  # read from file
  $inputstream = Bio::SeqIO->new(-file => "myseq.fa",-format => 'Fasta');
  $seqobj = $inputstream->next_seq();

  # get from database
  $db = Bio::DB::GenBank->new();
  $seqobj = $db->get_Seq_by_acc('X78121');

  # to get out parts of the sequence.

  print "Sequence ", $seqobj->id(), " with accession ", $seqobj->accession, " and desc ", $seqobj->desc, "\n";

  $string  = $seqobj->seq();
  $string2 = $seqobj->subseq(1,40);

  

=head1 DESCRIPTION

Seq is a lightweight Sequence object, storing little more than the
sequence, its name, a computer useful unique name. It does not contain
sequence features or other information.  To have a sequence with
sequence features you should use the AnnSeq object (AnnotatedSequence)
which uses this object.

Sequence objects are defined by the Bio::SeqI interface, and this
object is a pure Perl implementation of the interface (if that's
gibberish to you, don't worry. The take home message is that this
object is the bioperl default sequence object, but other people can
use their own objects as sequences if they so wish). If you are
interested in wrapping your own objects as compliant Bioperl sequence
objects, then you should read the Bio::SeqI documentation

The documenation of this object is a merge of the Bio::Seq and
Bio::SeqI documentation.  This allows all the methods which you can
call on sequence objects here.

=head1 Reimplementation

The Sequence object was completely rewritten for the 0.6 series. This
was because the old Sequence object was becoming heavily bloated and
difficult to maintain. There are some key changes from the old object
to the new object, but basically, everything should work with the new
object with a minimal number of changes.

The key change is that the format IO has been removed from this object
and moved to the Bio::SeqIO system, which provides a much better way
to encapsulate the sequence format reading. Please read the SeqIO
documentation, but the take home message is that lines like

    # old style reading from files
    $seq = Bio::Seq->new( -file => "myfile");

Becomes

    # new style reading from files.
    $inputstream = Bio::SeqIO->new( -file => "myfile", -format => 'Fasta');
    $seqobj = $inputstream->next_seq();

For writing files, a similar system is used

     # old style writing to files
     print OUTPUT $seq->layout_fasta;

     # new style writing to files
     $outputstream = Bio::SeqIO->new( -fh => \*OUTPUT, -format => 'Fasta');
     $outputstream->write_seq($seqobj);

=head2 Deprecated methods

A number of methods which were present in the old 0.04/0.05 series
have been deprecated.  Most of these methods work as before, but
provide a warning that someone has called a deprecated method.

=over 

=item getseq - use seq/subseq instead

=item str - use seq/subseq instead

=item ary - use seq/subseq with your own split afterwards

=item type - use moltype, but notice that moltype returns different values (lowercase)

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.  Bug reports can be submitted via
 email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@sanger.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::NewSeq;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;
use Bio::SeqI;


@ISA = qw(Bio::Root::Object Bio::SeqI);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my($seq,$id,$acc,$desc,$moltype) =
      $self->_rearrange([qw(SEQ
			    ID
			    ACCESSION
			    DESC
			    MOLTYPE
			    )],
			@args);

  my $make = $self->SUPER::_initialize(@args);
  $self->seq($seq);
  $self->id($id);
  $self->accession($acc);
  $self->desc($desc);
  $self->moltype($moltype);

# set stuff in self from @args
  return $make; # success - we hope!
}

=head2 seq

 Title   : seq
 Usage   : $string    = $obj->seq()
 Function: Returns the sequence as a string of letters. The
           case of the letters is left up to the implementer.
           Suggested cases are upper case for proteins and lower case for
           DNA sequence (IUPAC standard), but you should not rely on this
 Returns : A scalar

=cut

sub seq{
   my ($obj,$value) = @_;
   if( defined $value) {
       if( $value !~ /^[A-Za-z]+$/ ) {
	   $obj->throw("Attempting to set the sequence to [$value] which does not look healthy");
       }
       $obj->{'seq'} = $value;
    }
    return $obj->{'seq'};

}

=head2 subseq

 Title   : subseq
 Usage   : $substring = $obj->subseq(10,40);
 Function: returns the subseq from start to end, where the first base
           is 1 and the number is inclusive, ie 1-2 are the first two
           bases of the sequence
 Returns : a string
 Args    :


=cut

sub subseq{
   my ($self,$start,$end) = @_;

   if( $start > $end ){
       $self->throw("in subseq, start [$start] has to be greater than end [$end]");
   }

   # remove one from start, and then length is end-start

   $start--;

   return substr $self->seq, $start, ($end-$start);

}

=head2 id

 Title   : id
 Usage   : $id_string = $obj->id();
 Function: returns the id, aka the common name of the Sequence object.
           
         The semantics of this is that it is the most likely string to be
         used as an identifier of the sequence, and likely to have "human" readability.
         The id is equivalent to the ID field of the GenBank/EMBL databanks and
         the id field of the Swissprot/sptrembl database. In fasta format, the >(\S+)
         is presumed to be the id, though some people overload the id to embed other 
         information. Bioperl does not use any embedded information in the ID field,
         and people are encouraged to use other mechanisms (accession field for example,
         or extending the sequence object) to solve this. 

 Returns : A string
 Args    : None


=cut

sub id{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'id'} = $value;
    }
    return $obj->{'id'};

}


=head2 accession

 Title   : accession
 Usage   : $unique_key = $obj->accession;
 Function: Returns a computer defined unique key for this sequence. The
           accession field name was chosen as in bioinformatics accession
           numbers are a known concept and well understood.

           The aim of this field is that it provides a unique
           placeholder of where this sequence came from. This allows
           the rest of the bioperl system to retrieve additional
           information on the sequence that it might want to
           store. Effectively it is the computer's unique id.

           The semantics for this should follow the URL type system
           (or interoperable name service from the OMG standard) being

           /xxx/yyy/zzz/unique_key.<version_number>

           any of these portions can be omitted except for unique key.

           The unique key indicates the actual unique key for this
           sequence object.  For sequences from the public databases
           these should be the accession numbers from
           GenBank/EMBL/DDBJ and accession numbers from
           swissprot/sptrembl.  For internal databases, the natural
           thing is to use your own database's primary key for the
           sequence.

           The version number is optional, and indicates a version
           which changes on the semantics for the underlying
           database. The only semantics which is enforced about
           version numbers is that higher numbers are more up to date
           than smaller numbers.

           The information before the unique_key is also optional but
           indicates the database (also called 'context') of the
           key. For example, 

           /mycompany/dna/DID138338

           would have a unique id of DID138338 and a context of /mycompany/dna

           For public databases, we suggest the following contexts are used:

           /bio/dna/X12671 - accession number X12671 from GenBank/EMBL/DDBJ
           /bio/pid/g496898 - protein pid number
           /bio/protein/P09651 - accession number from swissprot/sptrembl
           /bio/pdb/1HA1   - PDB identifier for protein structures

           For sequences with no context, ie accession numbers by
           themselves the correct interpretation is that for dna
           sequences these are /bio/dna/accession_number and for
           protein sequences these are /bio/protein/accession_number.

          Sequence From Files:

          A number of times in bioinformatics, one doesn't have a
          database but rather a sequence from a file. Here there is
          some ambiguity of what happens to this field. If the
          sequence file contains an accession number field then one
          should use that as providing the accession number
          information, probably interpreting it as one of the
          "standard" contexts above.  However, a different view is to
          claim that the accession number should indicate the file
          this was made from. For file formats that have no accession
          number field (eg, plain FASTA format, with no overloading of
          the ID line, or raw format, or PIR format), this provides a
          mechanism for identifying the sequence. The proposal is to
          extend the context now into a full URL, including the
          filename, with the "unique_id" now becoming the byte offset
          into the file for this sequence. To make this concept
          useful, the format of the file also needs to be encoded, so
          that this context can be used.  The proposal is that a
          ::<format-string> is placed after the machine specification
          of the URL. For example:

          file://localhost::EMBL/nfs/data/roa1.dat/556760

          would indicate a EMBL formatted file found on the nfs system 
          with byte offset 556760.


          Sequence From Raw Memory:

          Sequences created inside bioperl with no given accession number,
          should have the following accession field:

          /bioperl/<scalar_memory_location>
  
 Returns : A string
 Args    : None

=cut

sub accession{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'accession'} = $value;
    }
    return $obj->{'accession'};

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


=cut

sub moltype{
   my ($obj,$value) = @_;
   if( defined $value) {
       if( $value ne 'dna' && $value ne 'rna' && $value ne 'protein' ) {
	   $obj->throw("Molecular type $value is not a valid type (dna,rna or protein, lowercase)");
       }
       $obj->{'moltype'} = $value;
   }
   return $obj->{'moltype'};

}

=head2 desc

 Title   : desc
 Usage   : $obj->desc($newval)
 Function: 
 Example : 
 Returns : value of desc
 Args    : newvalue (optional)


=cut

sub desc{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'desc'} = $value;
    }
    return $obj->{'desc'};

}

=head2 can_call_new

 Title   : can_call_new
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub can_call_new{
   my ($self) = @_;

   return 1;

}

=head1 Methods Inherieted from Bio::SeqI

These methods are available on Bio::Seq, although they are actually
implemented on Bio::SeqI

=head2 revcom

 Title   : revcom
 Usage   : $rev = $seq->revcom()
 Function: Produces a new Bio::SeqI implementing object which
           is the reversed complement of the sequence. For protein
           sequences this throws an exception of 
           "Sequence is a protein. Cannot revcom"

           The id is the same id as the orginal sequence, and the
           accession number is also indentical. If someone wants to
           track that this sequence has be reversed, it needs to
           define its own extensions

           To do an inplace edit of an object you can go:
   
           $seqobj = $seqobj->revcom();

           This of course, causes Perl to handle the garbage
           collection of the old object, but it is roughly speaking as
           efficient as an inplace edit.

 Returns : A new (fresh) Bio::SeqI object
 Args    : none

=cut

=head2 trunc

 Title   : trunc
 Usage   : $subseq = $myseq->trunc(10,100);
 Function: Provides a truncation of a sequence,
           
 Example :
 Returns : a fresh Bio::SeqI implementing object
 Args    :


=cut


1;






