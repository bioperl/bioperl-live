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

    $seqio  = Bio::SeqIO->new ( '-format' => 'embl' , -file => 'myfile.dat');
    $seqobj = $seqio->next_seq();

    # features must implement Bio::SeqFeatureI

    @features = $seqobj->top_SeqFeatures(); # just top level
    @features = $seqobj->all_SeqFeatures(); # descend into sub features

    $seq      = $seqobj->seq(); # actual sequence as a string
    $seqstr   = $seqobj->subseq(10,50);
    $ann      = $seqobj->annotation(); # annotation object

=head1 DESCRIPTION

A Seq object is a sequence with sequence features placed on them. The
Seq object contains a PrimarySeq object for the actual sequence and
also implements its interface.

In bioperl we have 3 main players that people are going to use

  Bio::PrimarySeq - just the sequence and its names, nothing else.
  Bio::SeqFeatureI - a location on a sequence, potentially with a sequence.
                    and annotation
  Bio::Seq        - A sequence and a collection of seqfeatures (an aggregate) with
                    its own annotation.

Although bioperl is not tied to file formats heavily, these distinctions do
map to file formats sensibly and for some bioinformaticians this might help
you:

  Bio::PrimarySeq - Fasta file of a sequence
  Bio::SeqFeatureI - A single entry in an EMBL/GenBank/DDBJ feature table
  Bio::Seq        - A single EMBL/GenBank/DDBJ entry

By having this split we avoid alot of nasty ciricular references
(seqfeatures can hold a reference to a sequence without the sequence
holding a reference to the seqfeature).

Ian Korf really helped in the design of the Seq and SeqFeature system.

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

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq;
use vars qw(@ISA);
use strict;
use Bio::SeqI;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::RootI;
use Bio::Annotation;
use Bio::PrimarySeq;

@ISA = qw(Bio::Root::RootI Bio::SeqI);


=head2 new

 Title   : new
 Usage   : $seq    = Bio::Seq->new( -seq => 'ATGGGGGTGGTGGTACCCT',
                                    -id  => 'human_id',
				    -accession_number => 'AL000012',
				   );

 Function: Returns a new seq object from
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

    my $ann = new Bio::Annotation;
    $self->annotation($ann);
    $self->primary_seq($pseq);

    return $self;
}

=head1 PrimarySeq interface

The primaryseq interface is the basic sequence getting
and setting methods found on all sequences.

These methods implement the PrimarySeq interface by delegating
to the primary_seq inside the object. This means that you
can use a Seq object wherever there is a PrimarySeq, and
of course, you are free to use these functions anyway.

=cut

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

           The implementation provided here does not take moltype() into
           account. Allowed are all letters (A-Z) and '-','.', and '*'.

 Example :
 Returns : TRUE if the supplied sequence string is valid for the object, and
           FALSE otherwise.
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
 Returns : integer representing the length of the sequence.
 Args    :

=cut

sub length {
   my ($self) = @_;
   return $self->primary_seq()->length();
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


=cut

sub subseq {
   my ($self,$s,$e) = @_;
   return $self->primary_seq()->subseq($s,$e);
}

=head2 display_id

 Title   : display_id
 Usage   : $id_string = $obj->display_id($newid);
 Function: returns or sets the display id, aka the common name of the
           Sequence object.

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
 Args    : newid (optional)


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
 Returns : A string
 Args    : None


=cut

sub accession_number {
   my ($self,$value) = @_;
   if( defined $value ) {
       return $self->primary_seq->accession_number($value);
   }
   return $self->primary_seq->accession_number();
}


=head2 desc

 Title   : desc
 Usage   : $seqobj->desc()
 Function: Sets/Gets the description of the sequnce
 Example :
 Returns :
 Args    :


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
           their own object ids in a way the implementaiton can control
           clients can expect one id to map to one object.

           For sequences with no natural id, this method should return
           a stringified memory location.

           Also notice that this method is B<not> delegated to the
           internal PrimarySeq object
 Returns : A string
 Args    : None


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

sub can_call_new {
   my ($self) = @_;

   return 1;
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

sub moltype {
   my ($self,$value) = @_;
   if( defined $value ) {
       return $self->primary_seq->moltype($value);
   }
   return $self->primary_seq->moltype();
}

=head1 Methods provided in the Bio::PrimarySeqI interface

These methods are inherited from the PrimarySeq interface
and work as one expects, building new Bio::Seq objects
or other information as expected.

Sequence Features are B<not> transfered to the new objects.
This is possibly a mistake. Anyone who feels the urge in
dealing with this is welcome to give it a go.

=head2 revcom

 Title   : revcom
 Usage   : $rev = $seq->revcom()
 Function: Produces a new Bio::Seq object which
           is the reversed complement of the sequence. For protein
           sequences this throws an exception of "Sequence is a protein. Cannot revcom"

           The id is the same id as the orginal sequence, and the accession number
           is also indentical. If someone wants to track that this sequence has be
           reversed, it needs to define its own extensions

           To do an inplace edit of an object you can go:

           $seq = $seq->revcom();

           This of course, causes Perl to handle the garbage collection of the old
           object, but it is roughly speaking as efficient as an inplace edit.

 Returns : A new (fresh) Bio::Seq object
 Args    : none


=cut

=head2 trunc

 Title   : trunc
 Usage   : $subseq = $myseq->trunc(10,100);
 Function: Provides a truncation of a sequence,

 Example :
 Returns : a fresh Bio::Seq object
 Args    :


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
 Usage   : $obj->primary_seq($newval)
 Function:
 Example :
 Returns : value of primary_seq
 Args    : newvalue (optional)


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

=head2 annotation

 Title   : annotation
 Usage   : $obj->annotation($seq_obj)
 Function:
 Example :
 Returns : value of annotation
 Args    : newvalue (optional)


=cut

sub annotation {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'annotation'} = $value;
    }
    return $obj->{'annotation'};

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
 Returns : TRUE on success
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
 Returns : TRUE on success
 Args    : 


=cut

sub flush_SeqFeatures {
   my ($self) = @_;

   $self->{'_as_feat'} = [];
   return 1;
}

# in all other modules we use the object in the singular -- lack of consistency
# sucks
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
 Args    :


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
 Args    :


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
 Returns : number of SeqFeatures
 Args    : none


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
 Usage   :
 Function: Gets or sets the species
 Example : $species = $self->species();
 Returns : Bio::Species object
 Args    : Bio::Species object or none;


=cut

sub species {
    my ($self, $species) = @_;
    if ($species) {
        $self->{'species'} = $species;
    } else {
        return $self->{'species'};
    }
}

# keep AUTOLOAD happy
sub DESTROY {
    my ($self) = @_;
}

1;
