
#
# BioPerl module for Bio::SeqI
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqI - Interface definition for a Bio::SeqI

=head1 SYNOPSIS

    # get a Bio::SeqI compliant object somehow

    # to test this is a seq object

    $obj->isa("Bio::SeqI") || $obj->throw("$obj does not implement the Bio::SeqI interface");

    # accessors

    $string = $obj->seq();
    $substring = $obj->seq(12,50);
    $id        = $obj->id();
    $unique_key= $obj->accession(); # unique idenitifer assign to the sequence by the system
    
    # object manipulation
   
    eval {
	$rev    = $obj->revcom();
    };
    if( $@ ) {
	$obj->throw("Could not reverse complement. Probably not DNA. Actual exception\n$@\n");
    }

    $trunc = $obj->trunc(12,50);
    
    # $rev and $trunc are Bio::SeqI compliant objects



=head1 DESCRIPTION

 This object defines an abstract interface to sequences. There is a
pure perl implementation of this in Bio::Seq. If you just want to use
Bio::Seq objects, then please read that module first. This module
defines the interface, and is of more interest to people who want to 
wrap their own Perl Objects/RDBs/FileSystems etc in way that they "are"
bioperl sequence objects, even though it is not using Perl to store the
sequence etc.

This interface defines what bioperl consideres necessary to "be" a sequence,
without providing an implementation of this. (An implementation is provided in
Bio::Seq). If you want to provide a Bio::Seq 'compliant' object which in fact
wraps another object/database/out-of-perl experience, then this is the correct
thing to wrap.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@sanger.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqI;
use vars qw($AUTOLOAD @ISA);
use strict;
use Carp;

=head1 Implementation Specific Functions

These functions are the ones that a specific implementation must
define.

=head2 seq

 Title   : seq
 Usage   : $string    = $obj->seq()
           $substring = $obj->seq(10,40)
 Function: Returns the sequence as a string of letters. The
           case of the letters is left up to the implementer.
           Suggested cases are upper case at all times (IUPAC standard),
           but implementations are suggested to keep an open mind about
           case (some users... want lower case!)
 Returns : A scalar
 Args    : Optional start,end paramters. This start,end is in absolute
           coordinates.


=cut

sub seq{
   my ($self) = @_;

   if( $self->can('throw') ) {
       $self->throw("Bio::SeqI definition of seq - implementing class did not provide this method");
   } else {
       confess("Bio::SeqI definition of seq - implementing class did not provide this method");
   }
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
         and people are encouraged to use other mechanisms (accession field for example)
         to solve this. 

 Returns : A string
 Args    : None


=cut

sub id{
   my ($self) = @_;

   if( $self->can('throw') ) {
       $self->throw("Bio::SeqI definition of id - implementing class did not provide this method");
   } else {
       confess("Bio::SeqI definition of id - implementing class did not provide this method");
   }

}


=head2 accession

 Title   : accession
 Usage   : $unique_key = $obj->accession;
 Function: Returns a computer defined unique key for this sequence. The
           accession field name was chosen as in bioinformatics accession
           numbers are a known concept and well understood.

           The aim of this field is that it provides a unique placeholder
           of where this sequence came from. This allows the rest of the
           bioperl system to retrieve additional information on the sequence
           that it might want to store. Effectively it is the computer's unique id.

           The semantics for this should follow the URL type system (or interoperable
           name service from the OMG standard) being
   
           /xxx/yyy/zzz/unique_key.<version_number>

           any of these portions can be omitted except for unique key.

           The unique key indicates the actual unique key for this sequence object.
           For sequences from the public databases these should be the accession
           numbers from GenBank/EMBL/DDBJ and accession numbers from swissprot/sptrembl.
           For internal databases, the natural thing is to use your own database's
           primary key for the sequence.

           The version number is optional, and indicates a version which changes on
           the semantics for the underlying database. The only semantics which is
           enforced about version numbers is that higher numbers are more up to date
           than smaller numbers. 

           The information before the unique_key is also optional but indicates
           the database (also called 'context') of the key. For example,
 
           /mycompany/dna/DID138338 

           would have a unique id of DID138338 and a context of /mycompany/dna

           For public databases, we suggest the following contexts are used:

           /bio/dna/X12671 - accession number X12671 from GenBank/EMBL/DDBJ
           /bio/pid/g496898 - protein pid number
           /bio/protein/P09651 - accession number from swissprot/sptrembl
           /bio/pdb/1HA1   - PDB identifier for protein structures

           For sequences with no context, ie accession numbers by themselves
           the correct interpretation is that for dna sequences these are
           /bio/dna/accession_number and for protein sequences these are
           /bio/protein/accession_number.

          Sequence From Files:

          A number of times in bioinformatics, one doesn't have a database but rather
          a sequence from a file. Here there is some ambiguity of what happens to this
          field. If the sequence file contains an accession number field then
          one should use that as providing the accession number information,
          probably interpreting it as one of the "standard" contexts above.
          However, a different view is to claim that the accession number
          should indicate the file this was made from. For file formats that
          have no accession number field (eg, plain FASTA format, with no
          overloading of the ID line, or raw format, or PIR format), this
          provides a mechanism for identifying the sequence. The proposal is
          to extend the context now into a full URL, including the filename, 
          with the "unique_id" now becoming the byte offset into the file
          for this sequence. To make this concept useful, the format of the
          file also needs to be encoded, so that this context can be used. 
          The proposal is that a ::<format-string> is placed after the
          machine specification of the URL. For example:

          file://localhost::EMBL/nfs/data/roa1.dat/556760

          would indicate a EMBL formatted file found on the nfs system 
          with byte offset 556760.


          Sequence From Raw Memory:

          Sequences created inside bioperl should have the following accession field:

          /bioperl/<scalar_memory_location>
  
 Returns : A string
 Args    : None


=cut

sub accession{
   my ($self,@args) = @_;

   if( $self->can('throw') ) {
       $self->throw("Bio::SeqI definition of seq - implementing class did not provide this method");
   } else {
       confess("Bio::SeqI definition of seq - implementing class did not provide this method");
   }

}


=head2 can_call_new

 Title   : can_call_new
 Usage   : if( $obj->can_call_new ) {
             $newobj = $obj->new( %param );
	 }
 Function: can_call_new returns 1 or 0 depending
           on whether an implementation allows new
           constructor 
 Example :
 Returns : 
 Args    :


=cut

sub can_call_new{
   my ($self,@args) = @_;


}


1;


