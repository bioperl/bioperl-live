
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

    $seqio  = Bio::SeqIO->new ( '-format' => 'Fasta' , -file => 'myfile.fasta');
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

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  bioperl-guts-l@bioperl.org            - Technically-oriented discussion
  http://bio.perl.org/MailList.html  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

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

use Bio::Root::Object;
use Bio::Annotation;
use Bio::PrimarySeq;

@ISA = qw(Bio::Root::Object Bio::SeqI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  my($ann);
  my $make = $self->SUPER::_initialize(@args);


  # this is way too sneaky probably. We delegate the construction of
  # the Seq object onto PrimarySeq and then pop primary_seq into
  # our primary_seq slot

  my $pseq = Bio::PrimarySeq->new(@args);
  $self->{'_as_feat'} = [];
  $self->{'date'} = [];
  $self->{'secondary_accession'} = [];
  
  $ann = new Bio::Annotation;
  $self->annotation($ann);
  $self->primary_seq($pseq);

# set stuff in self from @args
 return $make; # success - we hope!
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
   my ($self,$value) = @_;
   if( defined $value ) {
       return $self->primary_seq()->seq($value);
   }
   return $self->primary_seq()->seq();
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

These methods are inherieted from the PrimarySeq interface
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


=head1 Methods for Backward Compatibility

These methods are here for backward compatibility with the old, 0.5
Seq objects. They all throw warnings that someone is using a 
deprecated method, and may eventually be removed completely from
this object. However, they are important to ease the transition from
the old system.

=head2 str

 Title   : str
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

=head2 ary

 Title   : ary
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

=head2 getseq

 Title   : getseq
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

=head2 type

 Title   : type
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


=cut

=head1 Seq only methods

These methods are specific to the Bio::Seq object, and not
found on the Bio::PrimarySeq object

=head2 primary_seq

 Title   : seq
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
 Usage   : $annseq->add_SeqFeature($feat);
 Function: Adds t
 Example :
 Returns : TRUE on success
 Args    :


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

=head2 top_SeqFeatures

 Title   : top_SeqFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub top_SeqFeatures {
   my ($self) = @_;
   my ($p,$f,$l) = caller;

   return @{$self->{'_as_feat'}};
}

=head2 all_SeqFeatures

 Title   : all_SeqFeatures
 Usage   :
 Function:
 Example :
 Returns : 
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

=head2 fetch_SeqFeatures

 Title   : fetch_SeqFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_SeqFeatures {
   my ($self,@args) = @_;

   $self->throw("Not implemented yet");
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
        return $self->{'species'}
    }
}

=head1 EMBL/GenBank/DDBJ methods

These methods are here to support the EMBL/GenBank/DDBJ format.
The problem is that these formats require a certain amount
of additional information (eg, what division they are from), but
to make bioperl slavishly involved with this is just a bad idea.

If you want to use these methods, B<please> preface them with
a $as->can('method-name'). If this fails, then do something
sensible. This means that we do not have to think about
being in lock-step with EMBL/GenBank/DDBJ but can still support
all the information that is required.

=head2 division

 Title   : division
 Usage   : $obj->division($newval)
 Function: 
 Returns : value of division
 Args    : newvalue (optional)


=cut

sub division {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'division'} = $value;
    }
    return $obj->{'division'};

}

=head2 molecule

 Title   : molecule
 Usage   : $obj->molecule($newval)
 Function: 
 Returns : type of molecule (DNA, mRNA)
 Args    : newvalue (optional)


=cut

sub molecule {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'molecule'} = $value;
    }
    return $obj->{'molecule'};

}

=head2 add_date

 Title   : add_date
 Usage   : $self->add_date($ref)
 Function: adds a date
 Example :
 Returns : 
 Args    :


=cut

sub add_date {
   my ($self) = shift;
   foreach my $dt ( @_ ) {
       push(@{$self->{'date'}},$dt);
   }
}

=head2 each_date

 Title   : each_date
 Usage   : foreach $dt ( $self->each_date() )
 Function: gets an array of dates
 Example :
 Returns : 
 Args    :


=cut

sub each_date {
   my ($self) = @_;
   return @{$self->{'date'}}; 
}

=head2 accession

 Title   : accession
 Usage   : $obj->accession($newval)
 Function: Whilst the underlying sequence object does not 
           have an accession, so we need one here. Won't stay
           when we do the reimplementation.
 Example : 
 Returns : value of accession
 Args    : newvalue (optional)


=cut
#'

sub accession {
   my ($obj,$value) = @_;
   return $obj->accession_number($value);
}

=head2 add_secondary_accession

 Title   : add_secondary_accession
 Usage   : $self->add_domment($ref)
 Function: adds a secondary_accession
 Example :
 Returns : 
 Args    :


=cut

sub add_secondary_accession {
   my ($self) = shift;
   foreach my $dt ( @_ ) {
       push(@{$self->{'secondary_accession'}},$dt);
   }
}

=head2 each_secondary_accession

 Title   : each_secondary_accession
 Usage   : foreach $dt ( $self->each_secondary_accession() )
 Function: gets an array of secondary_accessions
 Example :
 Returns : 
 Args    :


=cut

sub each_secondary_accession {
   my ($self) = @_;
   return @{$self->{'secondary_accession'}}; 
}

=head2 sv

 Title   : sv
 Usage   : $obj->sv($newval)
 Function: 
 Returns : value of sv
 Args    : newvalue (optional)


=cut

sub sv {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'sv'} = $value;
    }
    return $obj->{'sv'};

}

=head2 keywords

 Title   : keywords
 Usage   : $obj->keywords($newval)
 Function: 
 Returns : value of keywords
 Args    : newvalue (optional)


=cut

sub keywords {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'keywords'} = $value;
    }
    return $obj->{'keywords'};

}

1;
