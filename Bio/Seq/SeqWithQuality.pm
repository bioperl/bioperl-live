#
# BioPerl module for Bio::Seq::QualI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chad Matsalla <bioinformatics@dieselwurks.com
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::SeqWithQuality - Bioperl object packaging a sequence with its quality.
Deprecated class, use Bio::Seq::Quality instead!

=head1 SYNOPSIS

    use Bio::PrimarySeq;
    use Bio::Seq::PrimaryQual;
    use Bio::Seq::SeqWithQuality;

    # make from memory
    my $qual = Bio::Seq::SeqWithQuality->new
    ( -qual => '10 20 30 40 50 50 20 10',
      -seq => 'ATCGATCG',
      -id  => 'human_id',
      -accession_number => 'AL000012',
    );

    # make from objects
    # first, make a PrimarySeq object
    my $seqobj = Bio::PrimarySeq->new
    ( -seq => 'atcgatcg',
      -id  => 'GeneFragment-12',
      -accession_number => 'X78121',
      -alphabet => 'dna'
    );

    # now make a PrimaryQual object
    my $qualobj = Bio::Seq::PrimaryQual->new
    ( -qual => '10 20 30 40 50 50 20 10',
      -id  => 'GeneFragment-12',
      -accession_number => 'X78121',
      -alphabet => 'dna'
    );

    # now make the SeqWithQuality object
    my $swqobj = Bio::Seq::SeqWithQuality->new
    ( -seq  => $seqobj,
      -qual => $qualobj
    );
    # done!

    $swqobj->id(); # the id of the SeqWithQuality object
    # may not match the the id of the sequence or
    # of the quality (check the pod, luke)
    $swqobj->seq(); # the sequence of the SeqWithQuality object
    $swqobj->qual(); # the quality of the SeqWithQuality object

    # to get out parts of the sequence.

    print "Sequence ", $seqobj->id(), " with accession ",
    $seqobj->accession, " and desc ", $seqobj->desc, "\n";

    $string2 = $seqobj->subseq(1,40);

=head1 DESCRIPTION

This object stores base quality values together with the sequence string.

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

=head1 AUTHOR - Chad Matsalla

Email bioinformatics@dieselwurks.com

=head1 CONTRIBUTORS 

Jason Stajich, jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


package Bio::Seq::SeqWithQuality;


use strict;
use Bio::PrimarySeq;
use Bio::Seq::PrimaryQual;

use base qw(Bio::Root::Root Bio::PrimarySeqI Bio::Seq::QualI);

=head2 new()

 Title   : new()
 Usage   : $qual = Bio::Seq::SeqWithQuality ->new
             ( -qual => '10 20 30 40 50 50 20 10',
               -seq => 'ATCGATCG',
               -id  => 'human_id',
               -accession_number => 'AL000012',
               -trace_indices    => '0 5 10 15 20 25 30 35'
             );
 Function: Returns a new Bio::Seq::SeqWithQual object from basic
           constructors.
 Returns : a new Bio::Seq::PrimaryQual object
 Args    : -qual can be a quality string (see Bio::Seq::PrimaryQual for more
            information on this) or a reference to a Bio::Seq::PrimaryQual
            object.
           -seq can be a sequence string (see Bio::PrimarySeq for more
            information on this) or a reference to a Bio::PrimaryQual object.
           -seq, -id, -accession_number, -primary_id, -desc, -id behave like
            this:
             1. if they are provided on construction of the
                 Bio::Seq::SeqWithQuality they will be set as the descriptors for
                 the object unless changed by one of the following mechanisms:
              a) $obj->set_common_descriptors() is used and both the -seq and
                 the -qual object have the same descriptors. These common
                 descriptors will then become the descriptors for the
                 Bio::Seq::SeqWithQual object.
              b) the descriptors are manually set using the seq(), id(),
                 desc(), or accession_number(), primary_id(),
             2. if no descriptors are provided, the new() constructor will see
                if the descriptor used in the PrimarySeq and in the
                PrimaryQual objects match. If they do, they will become
                the descriptors for the SeqWithQuality object.
             To eliminate ambiguity, I strongly suggest you set the
             descriptors manually on construction of the object. Really.
           -trace_indices : a space_delimited list of trace indices
            (where would the peaks be drawn if this list of qualities
            was to be plotted?)

=cut

sub new {
   my ($class, @args) = @_;
   my $self = $class->SUPER::new(@args);
   # default: turn OFF the warnings
   $self->{supress_warnings} = 1;
   my($qual,$seq,$id,$acc,$pid,$desc,$given_id,$alphabet,$trace_indices) =
   $self->_rearrange([qw( QUAL SEQ DISPLAY_ID ACCESSION_NUMBER PRIMARY_ID DESC
     ID ALPHABET TRACE_INDICES )], @args);
   # Deal with the ID
   if ( defined $id && defined $given_id ) {
      if( $id ne $given_id ) {
         $self->throw("Provided both id and display_id constructor functions. [$id] [$given_id]");
      }
   }
   if( defined $given_id ) {
      $self->display_id($given_id);
      $id = $given_id;
   }
   # Import sequence first
   if (!$seq) {
      my $id;
      unless ($self->{supress_warnings} == 1) {
         $self->warn("You did not provide sequence information during the ".
           "construction of a Bio::Seq::SeqWithQuality object. Sequence ".
           "components for this object will be empty.");
      }
      if (!$alphabet) {
         $self->throw("If you want me to create a PrimarySeq object for your ".
           "empty sequence <boggle> you must specify a -alphabet to satisfy ".
           "the constructor requirements for a Bio::PrimarySeq object with no ".
           "sequence. Read the POD for it, luke.");
      }
      $self->{seq_ref} = Bio::PrimarySeq->new( -seq              =>  "",
                                               -accession_number =>  $acc,
                                               -primary_id       =>  $pid,
                                               -desc             =>  $desc,
                                               -display_id       =>  $id,
                                               -alphabet         =>  $alphabet );
   } elsif ($seq->isa('Bio::PrimarySeqI') || $seq->isa('Bio::SeqI')) {
      $self->{seq_ref} = $seq;
   } elsif (ref($seq)) {
      $self->throw("You passed a seq argument into a SeqWithQUality object and".
        " it was a reference ($seq) which did not inherit from Bio::SeqI or ".
        "Bio::PrimarySeqI. I don't know what to do with this!");
   } else {
      my $seqobj = Bio::PrimarySeq->new( -seq              => $seq,
                                         -accession_number => $acc,
                                         -primary_id       => $pid,
                                         -desc             => $desc,
                                         -display_id       => $id   );
      $self->{seq_ref} = $seqobj;
   }
   # Then import the quality scores
   if (!defined($qual)) {
      $self->{qual_ref} = Bio::Seq::PrimaryQual->new( -qual             => "",
                                                      -accession_number => $acc,
                                                      -primary_id       => $pid,
                                                      -desc             => $desc,
                                                      -display_id       => $id, );
   } elsif (ref($qual) eq "Bio::Seq::PrimaryQual") {
      $self->{qual_ref} = $qual;
   } else {
      my $qualobj = Bio::Seq::PrimaryQual->new( -qual             => $qual,
                                                -accession_number => $acc,
                                                -primary_id       => $pid,
                                                -desc             => $desc,
                                                -display_id       => $id,
                                                -trace_indices    => $trace_indices );
      $self->{qual_ref} = $qualobj;
   }
   # Now try to set the descriptors for this object
   $self->_set_descriptors($qual,$seq,$id,$acc,$pid,$desc,$given_id,$alphabet);
   $self->length();
   $self->deprecated("deprecated class - use Bio::Seq::Quality instead");

   return $self;
}

=head2 _common_id()

 Title   : _common_id()
 Usage   : $common_id = $self->_common_id();
 Function: Compare the display_id of {qual_ref} and {seq_ref}.
 Returns : Nothing if they don't match. If they do return {seq_ref}->display_id()
 Args    : None.

=cut

#'
sub _common_id {
   my $self = shift;
   return if (!$self->{seq_ref} || !$self->{qual_ref});
   my $sid = $self->{seq_ref}->display_id();
   return if (!$sid);
   return if (!$self->{qual_ref}->display_id());
   return $sid if ($sid eq $self->{qual_ref}->display_id());
   # should this become a warning?
   # print("ids $sid and $self->{qual_ref}->display_id() do not match. Bummer.\n");
}

=head2 _common_display_id()

 Title   : _common_id()
 Usage   : $common_id = $self->_common_display_id();
 Function: Compare the display_id of {qual_ref} and {seq_ref}.
 Returns : Nothing if they don't match. If they do return {seq_ref}->display_id()
 Args    : None.

=cut

#'
sub _common_display_id {
   my $self = shift;
   $self->common_id();
}

=head2 _common_accession_number()

 Title   : _common_accession_number()
 Usage   : $common_id = $self->_common_accession_number();
 Function: Compare the accession_number() of {qual_ref} and {seq_ref}.
 Returns : Nothing if they don't match. If they do return {seq_ref}->accession_number()
 Args    : None.

=cut

#'
sub _common_accession_number {
   my $self = shift;
   return if ($self->{seq_ref} || $self->{qual_ref});
   my $acc = $self->{seq_ref}->accession_number();
   # if (!$acc) { print("the seqref has no acc.\n"); }
   return if (!$acc);
   # if ($acc eq $self->{qual_ref}->accession_number()) { print("$acc matches ".$self->{qual_ref}->accession_number()."\n"); }
   return $acc if ($acc eq $self->{qual_ref}->accession_number());
   # should this become a warning?
   # print("accession numbers $acc and $self->{qual_ref}->accession_number() do not match. Bummer.\n");
}

=head2 _common_primary_id()

 Title   : _common_primary_id()
 Usage   : $common_primard_id = $self->_common_primary_id();
 Function: Compare the primary_id of {qual_ref} and {seq_ref}.
 Returns : Nothing if they don't match. If they do return {seq_ref}->primary_id()
 Args    : None.

=cut

#'
sub _common_primary_id {
   my $self = shift;
   return if ($self->{seq_ref} || $self->{qual_ref});
   my $pid = $self->{seq_ref}->primary_id();
   return if (!$pid);
   return $pid if ($pid eq $self->{qual_ref}->primary_id());
   # should this become a warning?
   # print("primary_ids $pid and $self->{qual_ref}->primary_id() do not match. Bummer.\n");
}

=head2 _common_desc()

 Title   : _common_desc()
 Usage   : $common_desc = $self->_common_desc();
 Function: Compare the desc of {qual_ref} and {seq_ref}.
 Returns : Nothing if they don't match. If they do return {seq_ref}->desc()
 Args    : None.

=cut

#'
sub _common_desc {
   my $self = shift;
   return if ($self->{seq_ref} || $self->{qual_ref});
   my $des = $self->{seq_ref}->desc();
   return if (!$des);
   return $des if ($des eq $self->{qual_ref}->desc());
   # should this become a warning?
   # print("descriptions $des and $self->{qual_ref}->desc() do not match. Bummer.\n");
}

=head2 set_common_descriptors()

 Title   : set_common_descriptors()
 Usage   : $self->set_common_descriptors();
 Function: Compare the descriptors (id,accession_number,display_id,
           primary_id, desc) for the PrimarySeq and PrimaryQual objects
           within the SeqWithQuality object. If they match, make that
           descriptor the descriptor for the SeqWithQuality object.
 Returns : Nothing.
 Args    : None.

=cut

sub set_common_descriptors {
   my $self = shift;
   return if ($self->{seq_ref} || $self->{qual_ref});
   &_common_id();
   &_common_display_id();
   &_common_accession_number();
   &_common_primary_id();
   &_common_desc();
}

=head2 alphabet()

 Title   : alphabet();
 Usage   : $molecule_type = $obj->alphabet();
 Function: Get the molecule type from the PrimarySeq object.
 Returns : What what PrimarySeq says the type of the sequence is.
 Args    : None.

=cut

sub alphabet {
   my $self = shift;
   return $self->{seq_ref}->alphabet();
}

=head2 display_id()

 Title   : display_id()
 Usage   : $id_string = $obj->display_id();
 Function: Returns the display id, aka the common name of the Quality object.
           The semantics of this is that it is the most likely string to be
           used as an identifier of the quality sequence, and likely to have
           "human" readability.  The id is equivalent to the ID field of the
           GenBank/EMBL databanks and the id field of the Swissprot/sptrembl
           database. In fasta format, the >(\S+) is presumed to be the id,
           though some people overload the id to embed other information.
           Bioperl does not use any embedded information in the ID field,
           and people are encouraged to use other mechanisms (accession
           field for example, or extending the sequence object) to solve
           this. Notice that $seq->id() maps to this function, mainly for
           legacy/convience issues.
           This method sets the display_id for the SeqWithQuality object.
 Returns : A string
 Args    : If a scalar is provided, it is set as the new display_id for
           the SeqWithQuality object.
 Status  : Virtual

=cut

sub display_id {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'display_id'} = $value;
    }
    return $obj->{'display_id'};

}

=head2 accession_number()

 Title   : accession_number()
 Usage   : $unique_biological_key = $obj->accession_number();
 Function: Returns the unique biological id for a sequence, commonly
           called the accession_number. For sequences from established
           databases, the implementors should try to use the correct
           accession number. Notice that primary_id() provides the unique id
           for the implemetation, allowing multiple objects to have the same
           accession number in a particular implementation. For sequences
           with no accession number, this method should return "unknown".
           This method sets the accession_number for the SeqWithQuality
           object. 
 Returns : A string (the value of accession_number)
 Args    : If a scalar is provided, it is set as the new accession_number
           for the SeqWithQuality object.
 Status  : Virtual

=cut

sub accession_number {
    my( $obj, $acc ) = @_;

    if (defined $acc) {
        $obj->{'accession_number'} = $acc;
    } else {
        $acc = $obj->{'accession_number'};
        $acc = 'unknown' unless defined $acc;
    }
    return $acc;
}

=head2 primary_id()

 Title   : primary_id()
 Usage   : $unique_implementation_key = $obj->primary_id();
 Function: Returns the unique id for this object in this implementation.
           This allows implementations to manage their own object ids in a
           way the implementaiton can control clients can expect one id to
           map to one object. For sequences with no accession number, this
           method should return a stringified memory location.
           This method sets the primary_id for the SeqWithQuality object.
 Returns : A string. (the value of primary_id)
 Args    : If a scalar is provided, it is set as the new primary_id for
           the SeqWithQuality object.

=cut

sub primary_id {
   my ($obj,$value) = @_;
   if ($value) {
      $obj->{'primary_id'} = $value;
    }
   return $obj->{'primary_id'};

}

=head2 desc()

 Title   : desc()
 Usage   : $qual->desc($newval); _or_ 
           $description = $qual->desc();
 Function: Get/set description text for this SeqWithQuality object.
 Returns : A string. (the value of desc)
 Args    : If a scalar is provided, it is set as the new desc for the
           SeqWithQuality object.

=cut

sub desc {
   # a mechanism to set the disc for the SeqWithQuality object.
   # probably will be used most often by set_common_features()
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'desc'} = $value;
   }
   return $obj->{'desc'};
}

=head2 id()

 Title   : id()
 Usage   : $id = $qual->id();
 Function: Return the ID of the quality. This should normally be (and
           actually is in the implementation provided here) just a synonym
           for display_id().
 Returns : A string. (the value of id)
 Args    : If a scalar is provided, it is set as the new id for the
           SeqWithQuality object.

=cut

sub id {
   my ($self,$value) = @_;
   if (!$self) { $self->throw("no value for self in $value"); }
   if( defined $value ) {
      return $self->display_id($value);
   }
   return $self->display_id();
}

=head2 seq

 Title   : seq()
 Usage   : $string    = $obj->seq(); _or_ $obj->seq("atctatcatca");
 Function: Returns the sequence that is contained in the imbedded in the
           PrimarySeq object within the SeqWithQuality object
 Returns : A scalar (the seq() value for the imbedded PrimarySeq object.)
 Args    : If a scalar is provided, the SeqWithQuality object will
           attempt to set that as the sequence for the imbedded PrimarySeq
           object. Otherwise, the value of seq() for the PrimarySeq object
           is returned.
 Notes   : This is probably not a good idea because you then should call
           length() to make sure that the sequence and quality are of the
           same length. Even then, how can you make sure that this sequence
           belongs with that quality? I provided this to give you rope to
           hang yourself with. Tie it to a strong device and use a good
           knot.

=cut

sub seq {
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{seq_ref}->seq($value);
      $self->length();
   }
   return $self->{seq_ref}->seq();
}

=head2 qual()

 Title   : qual()
 Usage   : @quality_values  = @{$obj->qual()}; _or_
           $obj->qual("10 10 20 40 50");
 Function: Returns the quality as imbedded in the PrimaryQual object
           within the SeqWithQuality object.
 Returns : A reference to an array containing the quality values in the 
           PrimaryQual object.
 Args    : If a scalar is provided, the SeqWithQuality object will
           attempt to set that as the quality for the imbedded PrimaryQual
           object. Otherwise, the value of qual() for the PrimaryQual
           object is returned.
 Notes   : This is probably not a good idea because you then should call
           length() to make sure that the sequence and quality are of the
           same length. Even then, how can you make sure that this sequence
           belongs with that quality? I provided this to give you a strong
           board with which to flagellate yourself. 

=cut

sub qual {
   my ($self,$value) = @_;

   if( defined $value) {
      $self->{qual_ref}->qual($value);
      # update the lengths
      $self->length();
   }
   return $self->{qual_ref}->qual();
}



=head2 trace_indices()

 Title   : trace_indices()
 Usage   : @trace_indice_values  = @{$obj->trace_indices()}; _or_
           $obj->trace_indices("10 10 20 40 50");
 Function: Returns the trace_indices as imbedded in the Primaryqual object
           within the SeqWithQualiity object.
 Returns : A reference to an array containing the trace_indice values in the 
           PrimaryQual object.
 Args    : If a scalar is provided, the SeqWithuQuality object will
           attempt to set that as the trace_indices for the imbedded PrimaryQual
           object. Otherwise, the value of trace_indices() for the PrimaryQual
           object is returned.
 Notes   : This is probably not a good idea because you then should call
           length() to make sure that the sequence and trace_indices are of the
           same length. Even then, how can you make sure that this sequence
           belongs with that trace_indicex? I provided this to give you a strong
           board with which to flagellate yourself. 

=cut

sub trace_indices {
   my ($self,$value) = @_;

   if( defined $value) {
      $self->{qual_ref}->trace_indices($value);
      # update the lengths
      $self->length();
   }
   return $self->{qual_ref}->trace_indices();
}




=head2 length()

 Title   : length()
 Usage   : $length = $seqWqual->length();
 Function: Get the length of the SeqWithQuality sequence/quality.
 Returns : Returns the length of the sequence and quality if they are
           both the same. Returns "DIFFERENT" if they differ.
 Args    : None.

=cut

sub length {
   my $self = shift;
   if (!$self->{seq_ref}) {
      unless ($self->{supress_warnings} == 1) {
         $self->warn("Can't find {seq_ref} here in length().");
      }
      return;
   }
   if (!$self->{qual_ref}) {
      unless ($self->{supress_warnings} == 1) {
         $self->warn("Can't find {qual_ref} here in length().");
      }
      return;
   }
   my $seql = $self->{seq_ref}->length();

   if ($seql != $self->{qual_ref}->length()) {
      unless ($self->{supress_warnings} == 1) {
         $self->warn("Sequence length (".$seql.") is different from quality ".
           "length (".$self->{qual_ref}->length().") in the SeqWithQuality ".
           "object. This can only lead to problems later.");
      }
      $self->{'length'} = "DIFFERENT";
   } else {
      $self->{'length'} = $seql;
   }
   return $self->{'length'};
}


=head2 qual_obj

 Title   : qual_obj($different_obj)
 Usage   : $qualobj = $seqWqual->qual_obj(); _or_
           $qualobj = $seqWqual->qual_obj($ref_to_primaryqual_obj);
 Function: Get the PrimaryQual object that is imbedded in the
           SeqWithQuality object or if a reference to a PrimaryQual object
           is provided, set this as the PrimaryQual object imbedded in the
           SeqWithQuality object.
 Returns : A reference to a Bio::Seq::SeqWithQuality object.

=cut

sub qual_obj {
   my ($self,$value) = @_;
   if (defined($value)) {
      if (ref($value) eq "Bio::Seq::PrimaryQual") {
         $self->{qual_ref} = $value;
         $self->debug("You successfully changed the PrimaryQual object within ".
           "a SeqWithQuality object. ID's for the SeqWithQuality object may ".
           "now not be what you expect. Use something like ".
           "set_common_descriptors() to fix them if you care,");
      } else {
         $self->debug("You tried to change the PrimaryQual object within a ".
           "SeqWithQuality object but you passed a reference to an object that".
           " was not a Bio::Seq::PrimaryQual object. Thus your change failed. ".
           "Sorry.\n");	    
      }
   }
   return $self->{qual_ref};
}


=head2 seq_obj

 Title   : seq_obj()
 Usage   : $seqobj = $seqWqual->qual_obj(); _or_
           $seqobj = $seqWqual->seq_obj($ref_to_primary_seq_obj);
 Function: Get the PrimarySeq object that is imbedded in the
           SeqWithQuality object or if a reference to a PrimarySeq object is
           provided, set this as the PrimarySeq object imbedded in the
           SeqWithQuality object.
 Returns : A reference to a Bio::PrimarySeq object.

=cut

sub seq_obj {
   my ($self,$value) = @_;
   if( defined $value) {
      if (ref($value) eq "Bio::PrimarySeq") {
         $self->debug("You successfully changed the PrimarySeq object within a".
           " SeqWithQuality object. ID's for the SeqWithQuality object may now".
           " not be what you expect. Use something like ".
           "set_common_descriptors() to fix them if you care,");
      } else {
         $self->debug("You tried to change the PrimarySeq object within a ".
           "SeqWithQuality object but you passed a reference to an object that".
           " was not a Bio::PrimarySeq object. Thus your change failed. Sorry.\n");
      }
   }
   return $self->{seq_ref};
}

=head2 _set_descriptors

 Title   : _set_descriptors()
 Usage   : $seqWqual->_qual_obj($qual,$seq,$id,$acc,$pid,$desc,$given_id,
           $alphabet);
 Function: Set the descriptors for the SeqWithQuality object. Try to
           match the descriptors in the PrimarySeq object and in the
           PrimaryQual object if descriptors were not provided with
           construction.
 Returns : Nothing.
 Args    : $qual,$seq,$id,$acc,$pid,$desc,$given_id,$alphabet as found
           in the new() method.
 Notes   : Really only intended to be called by the new() method. If
           you want to invoke a similar function try set_common_descriptors().

=cut


sub _set_descriptors {
   my ($self,$qual,$seq,$id,$acc,$pid,$desc,$given_id,$alphabet) = @_;
   my ($c_id,$c_acc,$c_pid,$c_desc);
   if (!$self->display_id()) {
      if ($c_id = $self->_common_id() ) { $self->display_id($c_id); }
      else {
         if ($self->{seq_ref}) {
            # print("Using seq_ref to set id to ".$self->{seq_ref}->display_id()."\n");
            # ::dumpValue($self->{seq_ref});
            $self->display_id($self->{seq_ref}->id());
         } elsif ($self->{qual_ref}) {
            $self->display_id($self->{qual_ref}->id());
         }
      }
   }
   if ($acc) { $self->accession_number($acc); }
   elsif ($c_acc = $self->_common_accession_number() ) { $self->accession_number($c_acc); }
   if ($pid) { $self->primary_id($pid); }
   elsif ($c_pid = $self->_common_primary_id() ) { $self->primary_id($c_pid); }
   if ($desc) { $self->desc($desc); }
   elsif ($c_desc = $self->_common_desc() ) { $self->desc($c_desc); }
}

=head2 subseq($start,$end)

 Title   : subseq($start,$end)
 Usage   : $subsequence = $obj->subseq($start,$end);
 Function: Returns the subseq from start to end, where the first base
           is 1 and the number is inclusive, ie 1-2 are the first two
           bases of the sequence.
 Returns : A string.
 Args    : Two positions.

=cut

sub subseq {
    my ($self,@args) = @_;
    # does a single value work?
    return $self->{seq_ref}->subseq(@args);
}

=head2 baseat($position)

 Title   : baseat($position)
 Usage   : $base_at_position_6 = $obj->baseat("6");
 Function: Returns a single base at the given position, where the first
           base is 1 and the number is inclusive, ie 1-2 are the first two
           bases of the sequence.
 Returns : A scalar.
 Args    : A position.

=cut

sub baseat {
    my ($self,$val) = @_;
    return $self->{seq_ref}->subseq($val,$val);
}

=head2 subqual($start,$end)

 Title   : subqual($start,$end)
 Usage   : @qualities = @{$obj->subqual(10,20);
 Function: returns the quality values from $start to $end, where the
           first value is 1 and the number is inclusive, ie 1-2 are the
           first two bases of the sequence. Start cannot be larger than
           end but can be equal.
 Returns : A reference to an array.
 Args    : a start position and an end position

=cut

sub subqual {
    my ($self,@args) = @_;
    return $self->{qual_ref}->subqual(@args);
}

=head2 qualat($position)

 Title   : qualat($position)
 Usage   : $quality = $obj->qualat(10);
 Function: Return the quality value at the given location, where the
           first value is 1 and the number is inclusive, ie 1-2 are the
           first two bases of the sequence. Start cannot be larger than
           end but can be equal.
 Returns : A scalar.
 Args    : A position.

=cut

sub qualat {
    my ($self,$val) = @_;
    return $self->{qual_ref}->qualat($val);
}

=head2 sub_trace_index($start,$end)

 Title   : sub_trace_index($start,$end)
 Usage   : @trace_indices = @{$obj->sub_trace_index(10,20);
 Function: returns the trace index values from $start to $end, where the
           first value is 1 and the number is inclusive, ie 1-2 are the
           first two bases of the sequence. Start cannot be larger than
           end but can be e_trace_index.
 Returns : A reference to an array.
 Args    : a start position and an end position

=cut

sub sub_trace_index {
    my ($self,@args) = @_;
    return $self->{qual_ref}->sub_trace_index(@args);
}

=head2 trace_index_at($position)

 Title   : trace_index_at($position)
 Usage   : $trace_index = $obj->trace_index_at(10);
 Function: Return the trace_index value at the given location, where the
           first value is 1 and the number is inclusive, ie 1-2 are the
           first two bases of the sequence. Start cannot be larger than
           end but can be etrace_index_.
 Returns : A scalar.
 Args    : A position.

=cut

sub trace_index_at {
    my ($self,$val) = @_;
    return $self->{qual_ref}->trace_index_at($val);
}

=head2 to_string()

 Title   : to_string()
 Usage   : $quality = $obj->to_string();
 Function: Return a textual representation of what the object contains.
           For this module, this function will return:
                qual
                seq
                display_id
                accession_number
                primary_id
                desc
                id
                length_sequence
                length_quality
 Returns : A scalar.
 Args    : None.

=cut

sub to_string {
        my ($self,$out,$result) = shift;
        $out = "qual: ".join(',',@{$self->qual()})."\n";
        foreach (qw(seq display_id accession_number primary_id desc id)) {
                $result = $self->$_();
                if (!$result) { $result = "<unset>"; }
                $out .= "$_: $result\n";
        }
        return $out;
}
1;
