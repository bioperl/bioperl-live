# $Id$
#
# BioPerl module for Bio::Seq
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq - Sequence object, with features

=head1 SYNOPSIS

    $seqio  = Bio::SeqIO->new( '-format' => 'embl' , -file => 'myfile.dat');
    $seqobj = $seqio->next_seq();

    # features must implement Bio::SeqFeatureI

    @features = $seqobj->top_SeqFeatures(); # just top level
    @features = $seqobj->all_SeqFeatures(); # descend into sub features

    $seq      = $seqobj->seq(); # actual sequence as a string
    $seqstr   = $seqobj->subseq(10,50);
    $ann      = $seqobj->annotation(); # annotation object

=head1 DESCRIPTION

A Seq object is a sequence with sequence features placed on it. The
Seq object contains a PrimarySeq object for the actual sequence and
also implements its interface.

In Bioperl we have 3 main players that people are going to use frequently

  Bio::PrimarySeq  - just the sequence and its names, nothing else.
  Bio::SeqFeatureI - a location on a sequence, potentially with a sequence
                     and annotation.
  Bio::Seq         - A sequence and a collection of sequence features
                     (an aggregate) with its own annotation.

Although Bioperl is not tied heavily to file formats these distinctions do
map to file formats sensibly and for some bioinformaticians this might help

  Bio::PrimarySeq  - Fasta file of a sequence
  Bio::SeqFeatureI - A single entry in an EMBL/GenBank/DDBJ feature table
  Bio::Seq         - A single EMBL/GenBank/DDBJ entry

By having this split we avoid a lot of nasty circular references
(sequence features can hold a reference to a sequence without the sequence
holding a reference to the sequence feature).

Ian Korf really helped in the design of the Seq and SeqFeature system.

=head1 EXAMPLES

A simple and fundamental block of code

  use Bio::SeqIO;

  my $seqIOobj = Bio::SeqIO->new(-file=>"1.fa"); # create a SeqIO object
  my $seqobj = $seqIOobj->next_seq;              # get a Seq object

With the Seq object in hand one has access to a powerful set of Bioperl
methods and Bioperl objects. This next script will take a file of sequences
in EMBL format and create a file of the reverse-complemented sequences
in Fasta format using Seq objects. It also prints out details about the
exons it finds as sequence features in Genbank Flat File format.

  use Bio::Seq;
  use Bio::SeqIO;

  $seqin = Bio::SeqIO->new( -format => 'EMBL' , -file => 'myfile.dat');
  $seqout= Bio::SeqIO->new( -format => 'Fasta', -file => '>output.fa');

  while((my $seqobj = $seqin->next_seq())) {
      print "Seen sequence ",$seqobj->display_id,", start of seq ",
            substr($seqobj->seq,1,10),"\n";
      if( $seqobj->moltype eq 'dna') {
	    $rev = $seqobj->revcom;
	    $id  = $seqobj->display_id();
            $id  = "$id.rev";
            $rev->display_id($id);
            $seqout->write_seq($rev);
      }

      foreach $feat ( $seqobj->top_SeqFeatures() ) {
           if( $feat->primary_tag eq 'exon' ) {
              print STDOUT "Location ",$feat->start,":",
                    $feat->end," GFF[",$feat->gff_string,"]\n";
	   }
      }
  }

Let's examine the script. The lines below import the Bioperl modules.
Seq is the main Bioperl sequence object and SeqIO is the Bioperl support
for reading sequences from files and to files

  use Bio::Seq;
  use Bio::SeqIO;

These two lines create two SeqIO streams: one for reading in sequences
and one for outputting sequences:

  $seqin = Bio::SeqIO->new( -format => 'EMBL' , -file => 'myfile.dat');
  $seqout= Bio::SeqIO->new( -format => 'Fasta', -file => '>output.fa');

Notice that in the "$seqout" case there is a greater-than sign,
indicating the file is being opened for writing.

Using the

  '-argument' => value

syntax is common in Bioperl. The file argument is like an argument
to open() . You can also pass in filehandles or FileHandle objects by
using the -fh argument (see L<Bio::SeqIO> documentation for details).
Many formats in Bioperl are handled, including Fasta, EMBL, GenBank,
Swissprot (swiss), PIR, and GCG.

  $seqin = Bio::SeqIO->new( -format => 'EMBL' , -file => 'myfile.dat');
  $seqout= Bio::SeqIO->new( -format => 'Fasta', -file => '>output.fa');

This is the main loop which will loop progressively through sequences
in a file, and each call to $seqio-E<gt>next_seq() provides a new Seq
object from the file:

  while((my $seqobj = $seqio->next_seq())) {

This print line below accesses fields in the Seq object directly. The
$seqobj-E<gt>display_id is the way to access the display_id attribute
of the Seq object. The $seqobj-E<gt>seq method gets the actual
sequence out as string. Then you can do manipulation of this if
you want to (there are however easy ways of doing truncation,
reverse-complement and translation).

  print "Seen sequence ",$seqobj->display_id,", start of seq ",
               substr($seqobj->seq,1,10),"\n";

Bioperl has to guess the type of the sequence, being either 'dna',
'rna', or 'protein'. The moltype attribute is one of these three
possibilities.

  if( $seqobj->moltype eq 'dna') {

The $seqobj-E<gt>revcom method provides the reverse complement of the Seq
object as another Seq object. Thus, the $rev variable is a reference to
another Seq object. For example, one could repeat the above print line
for this Seq object (putting $rev in place of $seqobj). In this
case we are going to output the object into the file stream we built
earlier on.

  $rev = $seqobj->revcom;

When we output it, we want the id of the outputted object
to be changed to "$id.rev", ie, with .rev on the end of the name. The
following lines retrieve the id of the sequence object, add .rev
to this and then set the display_id of the rev sequence object to
this. Notice that to set the display_id attribute you just need
call the same method, display_id(), with the new value as an argument.
Getting and setting values with the same method is common in Bioperl.

  $id  = $seqobj->display_id();
  $id  = "$id.rev";
  $rev->display_id($id);

The write_seq method on the SeqIO output object, $seqout, writes the
$rev object to the filestream we built at the top of the script.
The filestream knows that it is outputting in fasta format, and
so it provides fasta output.

  $seqout->write_seq($rev);

This block of code loops over sequence features in the sequence
object, trying to find ones who have been tagged as 'exon'.
Features have start and end attributes and can be outputted
in Genbank Flat File format, GFF, a standarized format for sequence
features.

  foreach $feat ( $seqobj->top_SeqFeatures() ) {
      if( $feat->primary_tag eq 'exon' ) {
          print STDOUT "Location ",$feat->start,":",
             $feat->end," GFF[",$feat->gff_string,"]\n";
      }
  }

The code above shows how a few Bio::Seq methods suffice to read, parse,
reformat and analyze sequences from a file. A full list of methods
available to Bio::Seq objects is shown below. Bear in mind that some of
these methods come from PrimarySeq objects, which are simpler
than Seq objects, stripped of features (see L<Bio::PrimarySeq> for
more information).

  # these methods return strings, and accept strings in some cases:

  $seqobj->seq();              # string of sequence
  $seqobj->subseq(5,10);       # part of the sequence as a string
  $seqobj->accession_number(); # when there, the accession number
  $seqobj->moltype();          # one of 'dna','rna',or 'protein'
  $seqobj->seq_version()       # when there, the version
  $seqobj->keywords();         # when there, the Keywords line
  $seqobj->length()            # length
  $seqobj->desc();             # description
  $seqobj->primary_id();       # a unique id for this sequence regardless
                               # of its display_id or accession number
  $seqobj->display_id();       # the human readable id of the sequence

Some of these values map to fields in common formats. For example, The
display_id() method returns the LOCUS name of a Genbank entry,
the (\S+) following the E<gt> character in a Fasta file, the ID from
a SwissProt file, and so on. The desc() method will return the DEFINITION
line of a Genbank file, the description following the display_id in a
Fasta file, and the DE field in a SwissProt file.

  # the following methods return new Seq objects, but
  # do not transfer features across to the new object:

  $seqobj->trunc(5,10)  # truncation from 5 to 10 as new object
  $seqobj->revcom       # reverse complements sequence
  $seqobj->translate    # translation of the sequence

  # if new() can be called this method returns 1, else 0

  $seqobj->can_call_new

  # the following method determines if the given string will be accepted
  # by the seq() method - if the string is acceptable then validate()
  # returns 1, or 0 if not

  $seqobj->validate_seq($string)

  # the following method returns or accepts a Species object:

  $seqobj->species();

Please see L<Bio::Species> for more information on this object.

  # the following method returns or accepts an Annotation object
  # which in turn allows access to Annotation::Reference
  # and Annotation::Comment objects:

  $seqobj->annotation();

These annotations typically refer to entire sequences, unlike features.
See L<Bio::Annotation>, L<Bio::Annotation::Reference>, and
L<Bio::Annotation::Comment> for details.

It is also important to be able to describe defined portions of a
sequence. The combination of some description and the corresponding
sub-sequence is called a feature - an exon and its coordinates within
a gene is an example of a feature, or a domain within a protein.

  # the following methods return an array of SeqFeatureI objects:

  $seqobj->top_SeqFeatures # The 'top level' sequence features
  $seqobj->all_SeqFeatures # All sequence features, including sub-seq
                           # features, such as features in an exon

  # to find out the number of features use:

  $seqobj->feature_count

Here are just some of the methods available to SeqFeatureI objects:

  # these methods return numbers:

  $feat->start          # start position (1 is the first base)
  $feat->end            # end position (2 is the second base)
  $feat->strand         # 1 means forward, -1 reverse, 0 not relevant

  # these methods return or accept strings:

  $feat->primary_tag    # the name of the sequence feature, eg
                        # 'exon', 'glycoslyation site', 'TM domain'
  $feat->source_tag     # where the feature comes from, eg, 'EMBL_GenBank',
                        # or 'BLAST'

  # this method returns the more austere PrimarySeq object, not a
  # Seq object - the main difference is that PrimarySeq objects do not
  # themselves contain sequence features

  $feat->seq            # the sequence between start,end on the
                        # correct strand of the sequence

See L<Bio::PrimarySeq> for more details on PrimarySeq objects.

  # useful methods for feature comparisons, for start/end points

  $feat->overlaps($other)  # do $feat and $other overlap?
  $feat->contains($other)  # is $other completely within $feat?
  $feat->equals($other)    # do $feat and $other completely agree?

  # one can also add features

  $seqobj->add_SeqFeature($feat)     # returns 1 if successful
  $seqobj->add_SeqFeature(@features) # returns 1 if successful

  # sub features. For complex join() statements, the feature
  # is one sequence feature with many sub SeqFeatures

  $feat->sub_SeqFeature  # returns array of sub seq features

Please see L<Bio::SeqFeatureI> and L<Bio::SeqFeature::Generic>,
for more information on sequence features.


=head1 FEEDBACK


=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bio.perl.org/MailList.html  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bioperl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney, inspired by Ian Korf objects


Email birney@sanger.ac.uk

Describe contact details here

=head1 APPENDIX


The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a "_".

=cut

#'
# Let the code begin...


package Bio::Seq;
use vars qw(@ISA $VERSION);
use strict;


# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Root;
use Bio::SeqI;
use Bio::RangeI;
use Bio::Annotation::Collection;
use Bio::PrimarySeq;

$VERSION = '1.0';
@ISA = qw(Bio::Root::Root Bio::SeqI Bio::RangeI);

=head2 new

 Title   : new
 Usage   : $seq = Bio::Seq->new( -seq => 'ATGGGGGTGGTGGTACCCT',
                                 -id  => 'human_id',
				 -accession_number => 'AL000012',
			       );

 Function: Returns a new Seq object from
           basic constructors, being a string for the sequence
           and strings for id and accession_number
 Returns : a new Bio::Seq object

=cut

sub new {
    # standard new call..
    my($caller,@args) = @_;
    my $self = $caller->SUPER::new(@args);
    # this is way too sneaky probably. We delegate the construction of
    # the Seq object onto PrimarySeq and then pop primary_seq into
    # our primary_seq slot

    my $pseq = Bio::PrimarySeq->new(@args);
    $self->{'_as_feat'} = [];

    my $ann = new Bio::Annotation::Collection;
    $self->annotation($ann);
    $self->primary_seq($pseq);

    return $self;
}

=head1 PrimarySeq interface


The PrimarySeq interface provides the basic sequence getting
and setting methods for on all sequences.

These methods implement the Bio::PrimarySeq interface by delegating
to the primary_seq inside the object. This means that you
can use a Seq object wherever there is a PrimarySeq, and
of course, you are free to use these functions anyway.

=cut

=head2 seq

 Title   : seq
 Usage   : $string = $obj->seq()
 Function: Returns the sequence as a string of letters. The
           case of the letters is left up to the implementer.
           Suggested cases are upper case for proteins and lower case for
           DNA sequence (IUPAC standard),
           but implementations are suggested to keep an open mind about
           case (some users... want mixed case!)
 Returns : A scalar
 Args    : None

=cut

sub seq {
   my ($self,@args) = @_;

   return $self->primary_seq()->seq(@args);
}

=head2 validate_seq

 Title   : validate_seq
 Usage   : if(! $seq->validate_seq($seq_str) ) {
                print "sequence $seq_str is not valid for an object of type ",
		      ref($seq), "\n";
	   }
 Function: Validates a given sequence string. A validating sequence string
           must be accepted by seq(). A string that does not validate will
           lead to an exception if passed to seq().

           The implementation provided here does not take alphabet() into
           account. Allowed are all letters (A-Z) and '-','.', and '*'.

 Example :
 Returns : 1 if the supplied sequence string is valid for the object, and
           0 otherwise.
 Args    : The sequence string to be validated.


=cut

sub validate_seq {
   my ($self,@args) = @_;

   return $self->primary_seq()->validate_seq(@args);
}

=head2 length

 Title   : length
 Usage   : $len = $seq->length()
 Function:
 Example :
 Returns : Integer representing the length of the sequence.
 Args    : None

=cut

sub length {
   my ($self) = @_;
   return $self->primary_seq()->length();
}

=head2 start

 Title   : start
 Usage   : $start = $seq->start()
 Function:
 Example :
 Returns : Integer representing the start of the sequence.
 Args    : None

=cut

sub start {
  my $self = shift;
  return 1;
}

=head2 end

 Title   : start
 Usage   : $start = $seq->end()
 Function:
 Example :
 Returns : Integer representing the end of the sequence.
 Args    : None

=cut

sub end {
  my $self = shift;
  return $self->length;
}

=head2 strand

 Title   : strand
 Usage   : $stand = $seq->strand()
 Function:
 Example :
 Returns : Integer representing the strand of the sequence [0,1,-1]
 Args    : None

=cut

sub strand {
  my $self = shift;
  return 0;
}

=head2 subseq

 Title   : subseq
 Usage   : $substring = $obj->subseq(10,40);
 Function: Returns the subseq from start to end, where the first base
           is 1 and the number is inclusive, ie 1-2 are the first two
           bases of the sequence

           Start cannot be larger than end but can be equal

 Returns : A string
 Args    : 2 integers


=cut

sub subseq {
   my ($self,$s,$e) = @_;
   return $self->primary_seq()->subseq($s,$e);
}

=head2 display_id

 Title   : display_id
 Usage   : $id = $obj->display_id or $obj->display_id($newid);
 Function: Gets or sets the display id, also known as the common name of
           the Seq object.

           The semantics of this is that it is the most likely string
           to be used as an identifier of the sequence, and likely to
           have "human" readability.  The id is equivalent to the LOCUS
           field of the GenBank/EMBL databanks and the ID field of the
           Swissprot/sptrembl database. In fasta format, the >(\S+) is
           presumed to be the id, though some people overload the id
           to embed other information. Bioperl does not use any
           embedded information in the ID field, and people are
           encouraged to use other mechanisms (accession field for
           example, or extending the sequence object) to solve this.

           Notice that $seq->id() maps to this function, mainly for
           legacy/convenience issues.
 Returns : A string
 Args    : None or a new id


=cut

sub display_id {
   my ($self,$value) = @_;
   if( defined $value ) {
       return $self->primary_seq->display_id($value);
   }
   return $self->primary_seq->display_id();
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

           Can also be used to set the accession number.
 Example : $key = $seq->accession_number or $seq->accession_number($key)
 Returns : A string
 Args    : None or an accession number


=cut

sub accession_number {
   my ($self,$value) = @_;
   if( defined $value ) {
       return $self->primary_seq->accession_number($value);
   }
   return $self->primary_seq->accession_number();
}

sub accession {
    my ($self,$value) = @_;

    return $self->accession_number($value);
}

=head2 desc

 Title   : desc
 Usage   : $seqobj->desc($string) or $seqobj->desc()
 Function: Sets or gets the description of the sequence
 Example :
 Returns : The description
 Args    : The description or none


=cut

sub desc {
   my ($self,$value) = @_;

   if( defined $value ) {
       return $self->primary_seq->desc($value);
   }
   return $self->primary_seq->desc();
}




=head2 primary_id

 Title   : primary_id
 Usage   : $unique_implementation_key = $obj->primary_id;
 Function: Returns the unique id for this object in this
           implementation. This allows implementations to manage
           their own object ids in a way the implementation can control
           clients can expect one id to map to one object. 

           For sequences with no natural id, this method should return
           a stringified memory location.

           Can also be used to set the primary_id.

           Also notice that this method is not delegated to the
           internal Bio::PrimarySeq object
 Example : $id = $seq->primary_id or $seq->primary_id($id)
 Returns : A string
 Args    : None or an id


=cut

sub primary_id {
   my ($obj,$value) = @_;

   if( defined $value) {
      $obj->{'primary_id'} = $value;
    }
   if( ! exists $obj->{'primary_id'} ) {
       return "$obj";
   }
   return $obj->{'primary_id'};
}

=head2 can_call_new

 Title   : can_call_new
 Usage   : if ( $obj->can_call_new ) {
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
 Args    : None


=cut

sub can_call_new {
   my ($self) = @_;

   return 1;
}

=head2 alphabet

 Title   : alphabet
 Usage   : if ( $obj->alphabet eq 'dna' ) { /Do Something/ }
 Function: Returns the type of sequence being one of
           'dna', 'rna' or 'protein'. This is case sensitive.

           This is not called <type> because this would cause
           upgrade problems from the 0.5 and earlier Seq objects.

 Returns : A string either 'dna','rna','protein'. NB - the object must
           make a call of the type - if there is no type specified it
           has to guess.
 Args    : None


=cut

sub alphabet {
   my ($self,$value) = @_;
   if( defined $value ) {
       return $self->primary_seq->alphabet($value);
   }
   return $self->primary_seq->alphabet();
}

=head1 Methods provided in the Bio::PrimarySeqI interface


These methods are inherited from the PrimarySeq interface
and work as one expects, building new Bio::Seq objects
or other information as expected. See L<Bio::PrimarySeq>
for more information.

Sequence Features are B<not> transfered to the new objects.
This is possibly a mistake. Anyone who feels the urge in
dealing with this is welcome to give it a go.

=head2 revcom

 Title   : revcom
 Usage   : $rev = $seq->revcom()
 Function: Produces a new Bio::Seq object which
           is the reversed complement of the sequence. For protein
           sequences this throws an exception of "Sequence is a protein.
           Cannot revcom"

           The id is the same id as the original sequence, and the
           accession number is also identical. If someone wants to track
           that this sequence has be reversed, it needs to define its own
           extensions

           To do an in-place edit of an object you can go:

           $seq = $seq->revcom();

           This of course, causes Perl to handle the garbage collection of
           the old object, but it is roughly speaking as efficient as an
           in-place edit.

 Returns : A new (fresh) Bio::Seq object
 Args    : None


=cut

=head2 trunc

 Title   : trunc
 Usage   : $subseq = $myseq->trunc(10,100);
 Function: Provides a truncation of a sequence

 Example :
 Returns : A fresh Seq object
 Args    : A Seq object


=cut

=head2 id

 Title   : id
 Usage   : $id = $seq->id()
 Function: This is mapped on display_id
 Example :
 Returns :
 Args    :


=cut

sub  id {
   my ($self,$value) = @_;

   if( defined $value ) {
	return $self->display_id($value);
   }
   return $self->display_id();
}


=head1 Seq only methods


These methods are specific to the Bio::Seq object, and not
found on the Bio::PrimarySeq object

=head2 primary_seq

 Title   : primary_seq
 Usage   : $seq->primary_seq or $seq->primary_seq($newval)
 Function: Get or set a PrimarySeq object
 Example :
 Returns : PrimarySeq object
 Args    : None or PrimarySeq object


=cut

sub primary_seq {
   my ($obj,$value) = @_;

   if( defined $value) {
       if( ! ref $value || ! $value->isa('Bio::PrimarySeqI') ) {
	   $obj->throw("$value is not a Bio::PrimarySeq compliant object");
       }

       $obj->{'primary_seq'} = $value;
       # descend down over all seqfeature objects, seeing whether they
       # want an attached seq.

       foreach my $sf ( $obj->top_SeqFeatures() ) {
	   if( $sf->can("attach_seq") ) {
	       $sf->attach_seq($value);
	   } else {
	       $obj->warn("In Seq primary_seq, a sequence feature cannot attach seq. Bugger");
	   }
       }

   }
   return $obj->{'primary_seq'};

}


=head2 add_SeqFeature

 Title   : add_SeqFeature
 Usage   : $seq->add_SeqFeature($feat);
           $seq->add_SeqFeature(@feat);
 Function: Adds the given feature object (or each of an array of feature
           objects to the feature array of this
           sequence. The object passed is required to implement the
           Bio::SeqFeatureI interface.
 Example :
 Returns : 1 on success
 Args    : A Bio::SeqFeatureI implementing object, or an array of such objects.


=cut

sub add_SeqFeature {
   my ($self,@feat) = @_;
   my ($fseq,$aseq);


   foreach my $feat ( @feat ) {
       if( !$feat->isa("Bio::SeqFeatureI") ) {
	   $self->throw("$feat is not a SeqFeatureI and that's what we expect...");
       }

       if( $feat->can("entire_seq") ) {
	   $fseq = $feat->entire_seq;
	   $aseq = $self->primary_seq;

	   if( defined $aseq ) {
	       if( defined $fseq ) {
		   unless ($aseq == $fseq) {
		       $self->warn("$feat has an attached sequence which is not in this annseq. I worry about this");
		   }
	       } else {
		   if( $feat->can("attach_seq") ) {
		       # attach it
		       $feat->attach_seq($aseq);
		   }
	       }
	   } # end of if aseq
       } # end of if the feat can entire_seq

       push(@{$self->{'_as_feat'}},$feat);
   }
   return 1;
}

=head2 flush_SeqFeatures

 Title   : flush_SeqFeatures
 Usage   : $seq->flush_SeqFeatures();
 Function: Flushes all attached SeqFeatureI objects. To remove individual
           feature objects, first obtain all using all_SeqFeatures(), then
           flush and re-add those you want to keep.
 Example :
 Returns : 1 on success
 Args    : None


=cut

sub flush_SeqFeatures {
   my ($self) = @_;

   $self->{'_as_feat'} = [];
   return 1;
}

# in all other modules we use the object in the singular --
# lack of consistency sucks
sub flush_SeqFeature {
    return shift()->flush_SeqFeatures();
}

=head2 top_SeqFeatures

 Title   : top_SeqFeatures
 Usage   : @feat_ary = $seq->top_SeqFeatures();
 Function: Returns the array of top-level features for this sequence object.
           Features which are not top-level are subfeatures of one or more
           of the returned feature objects, which means that you must
           traverse the subfeature arrays of each top-level feature object
           in order to traverse all features associated with this sequence.

           Use all_SeqFeatures() if you want the feature tree flattened into
           one single array.
 Example :
 Returns : An array of Bio::SeqFeatureI implementing objects.
 Args    : None


=cut

sub top_SeqFeatures {
   my ($self) = @_;

   return @{$self->{'_as_feat'}};
}

=head2 all_SeqFeatures

 Title   : all_SeqFeatures
 Usage   : @feat_ary = $seq->all_SeqFeatures();
 Function: Returns the tree of feature objects attached to this sequence
           object flattened into one single array. Top-level features will
           still contain their subfeature-arrays, which means that you
           will encounter subfeatures twice if you traverse the subfeature
           tree of the returned objects.

           Use top_SeqFeatures() if you want the array to contain only the
           top-level features.
 Example :
 Returns : An array of Bio::SeqFeatureI implementing objects.
 Args    : None


=cut

sub all_SeqFeatures {
   my ($self) = @_;
   my (@array);
   foreach my $feat ( $self->top_SeqFeatures() ){
       push(@array,$feat);
       &_retrieve_subSeqFeature(\@array,$feat);
   }

   return @array;
}

=head2 feature_count

 Title   : feature_count
 Usage   : $seq->feature_count()
 Function: Return the number of SeqFeatures attached to a sequence
 Example :
 Returns : Number of SeqFeatures
 Args    : None


=cut

sub feature_count {
    my ($self) = @_;

    if (defined($self->{'_as_feat'})) {
	return ($#{$self->{'_as_feat'}} + 1);
    } else {
	return 0;
    }
}

sub _retrieve_subSeqFeature {
    my ($arrayref,$feat) = @_;

    foreach my $sub ( $feat->sub_SeqFeature() ) {
	push(@$arrayref,$sub);
	&_retrieve_subSeqFeature($arrayref,$sub);
    }

}


=head2 species

 Title   : species
 Usage   : $species = $seq->species() or $seq->species($species)
 Function: Gets or sets the species
 Example :
 Returns : Bio::Species object
 Args    : None or Bio::Species object

See L<Bio::Species> for more information

=cut

sub species {
    my ($self, $species) = @_;
    if ($species) {
        $self->{'species'} = $species;
    } else {
        return $self->{'species'};
    }
}


=head2 annotation

 Title   : annotation
 Usage   : $ann = $seq->annotation or $seq->annotation($annotation)
 Function: Gets or sets the annotation
 Example :
 Returns : Bio::Annotation object
 Args    : None or Bio::Annotation object

See L<Bio::Annotation> for more information

=cut

sub annotation {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'annotation'} = $value;
    }
    return $obj->{'annotation'};

}


# keep AUTOLOAD happy
sub DESTROY {
    my ($self) = @_;
}

1;
