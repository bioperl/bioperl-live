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

Bio::Seq::QualI - Interface definition for a Bio::Seq::Qual

=head1 SYNOPSIS

    # get a Bio::Seq::Qual compliant object somehow

    # to test this is a seq object

    $obj->isa("Bio::Seq::QualI") 
       || $obj->throw("$obj does not implement the Bio::Seq::QualI interface");

    # accessors

    $string    = $obj->qual();
    $substring = $obj->subqual(12,50);
    $display   = $obj->display_id(); # for human display
    $id        = $obj->primary_id(); # unique id for this object,
                                     # implementation defined
    $unique_key= $obj->accession_number();
                       # unique biological id



=head1 DESCRIPTION

This object defines an abstract interface to basic quality
information. PrimaryQual is an object just for the quality and its
name(s), nothing more. There is a pure perl implementation of this in
Bio::Seq::PrimaryQual. If you just want to use Bio::Seq::PrimaryQual
objects, then please read that module first. This module defines the
interface, and is of more interest to people who want to wrap their
own Perl Objects/RDBs/FileSystems etc in way that they "are" bioperl
quality objects, even though it is not using Perl to store the
sequence etc.

This interface defines what bioperl consideres necessary to "be" a
sequence of qualities, without providing an implementation of
this. (An implementation is provided in Bio::Seq::PrimaryQual). If you
want to provide a Bio::Seq::PrimaryQual 'compliant' object which in
fact wraps another object/database/out-of-perl experience, then this
is the correct thing to wrap, generally by providing a wrapper class
which would inherit from your object and this Bio::Seq::QualI
interface. The wrapper class then would have methods lists in the
"Implementation Specific Functions" which would provide these methods
for your object.


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

This module is heavily based on Bio::Seq::PrimarySeq and is modeled after
or outright copies sections of it. Thanks Ewan!

Email bioinformatics@dieselwurks.com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::QualI;
use strict;
use Carp;

use base qw(Bio::Root::RootI);

=head1 Implementation Specific Functions

These functions are the ones that a specific implementation must
define.

=head2 qual()

 Title   : qual()
 Usage   : @quality_values  = @{$obj->qual()};
 Function: Returns the quality as a reference to an array containing the
	quality values. The individual elements of the quality array are
	not validated and can be any numeric value.
 Returns : A reference to an array.
 Status  : 

=cut

sub qual {
   my ($self) = @_;
   if( $self->can('throw') ) {
       $self->throw("Bio::Seq::QualI definition of qual - implementing class did not provide this method");
   } else {
       confess("Bio::Seq::QualI definition of qual - implementing class did not provide this method");
   }
}

=head2 subqual($start,$end)

 Title   : subqual($start,$end)
 Usage   : @subset_of_quality_values = @{$obj->subseq(10,40)};
 Function: returns the quality values from $start to $end, where the
	first value is 1 and the number is inclusive, ie 1-2 are the first
	two bases of the sequence. Start cannot be larger than end but can
	be equal.
 Returns : A reference to an array.
 Args    : a start position and an end position


=cut

sub subqual {
   my ($self) = @_;

   if( $self->can('throw') ) {
       $self->throw("Bio::Seq::QualI definition of subqual - implementing class did not provide this method");
   } else {
       confess("Bio::Seq::QualI definition of subqual - implementing class did not provide this method");
   }

}

=head2 display_id()

 Title   : display_id()
 Usage   : $id_string = $obj->display_id() _or_
	   $id_string = $obj->display_id($new_display_id);
 Function: Returns the display id, aka the common name of the Quality 
	object.
	The semantics of this is that it is the most likely string to be
	used as an identifier of the quality sequence, and likely to have
	"human" readability.  The id is equivalent to the ID field of the
	GenBank/EMBL databanks and the id field of the Swissprot/sptrembl
	database. In fasta format, the >(\S+) is presumed to be the id,
	though some people overload the id to embed other information.
	Bioperl does not use any embedded information in the ID field,
	and people are encouraged to use other mechanisms (accession field
	for example, or extending the sequence object) to solve this.
	Notice that $seq->id() maps to this function, mainly for 
	legacy/convience issues
 Returns : A string
 Args    : If an arg is provided, it will replace the existing display_id
	in the object.


=cut

sub display_id {
   my ($self) = @_;

   if( $self->can('throw') ) {
       $self->throw("Bio::Seq::QualI definition of id - implementing class did not provide this method");
   } else {
       confess("Bio::Seq::QualI definition of id - implementing class did not provide this method");
   }

}


=head2 accession_number()

 Title   : accession_number()
 Usage   : $unique_biological_key = $obj->accession_number(); _or_
	   $unique_biological_key = $obj->accession_number($new_acc_num);
 Function: Returns the unique biological id for a sequence, commonly 
	called the accession_number. For sequences from established 
	databases, the implementors should try to use the correct 
	accession number. Notice that primary_id() provides the unique id 
	for the implemetation, allowing multiple objects to have the same
	accession number in a particular implementation. For sequences
	with no accession number, this method should return "unknown".
 Returns : A string.
 Args    : If an arg is provided, it will replace the existing
	accession_number in the object.

=cut

sub accession_number {
   my ($self,@args) = @_;

   if( $self->can('throw') ) {
       $self->throw("Bio::Seq::QualI definition of seq - implementing class did not provide this method");
   } else {
       confess("Bio::Seq::QualI definition of seq - implementing class did not provide this method");
   }

}



=head2 primary_id()

 Title   : primary_id()
 Usage   : $unique_implementation_key = $obj->primary_id(); _or_
	   $unique_implementation_key = $obj->primary_id($new_prim_id);
 Function: Returns the unique id for this object in this implementation.
	This allows implementations to manage their own object ids in a
	way the implementaiton can control clients can expect one id to
	map to one object. For sequences with no accession number, this
	method should return a stringified memory location.
 Returns : A string
 Args    : If an arg is provided, it will replace the existing
	primary_id in the object.

=cut

sub primary_id {
   my ($self,@args) = @_;

   if( $self->can('throw') ) {
       $self->throw("Bio::Seq::QualI definition of qual - implementing class did not provide this method");
   } else {
       confess("Bio::Seq::QualI definition of qual - implementing class did not provide this method");
   }

}


=head2 can_call_new()

 Title   : can_call_new()
 Usage   : if( $obj->can_call_new ) {
             $newobj = $obj->new( %param );
	 }
 Function: can_call_new returns 1 or 0 depending on whether an
	implementation allows new constructor to be called. If a new
	constructor is allowed, then it should take the followed hashed
	constructor list.
           $myobject->new( -qual => $quality_as_string,
			   -display_id  => $id,
			   -accession_number => $accession,
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

=head2 qualat($position)

 Title   : qualat($position)
 Usage   : $quality = $obj->qualat(10);
 Function: Return the quality value at the given location, where the
        first value is 1 and the number is inclusive, ie 1-2 are the first
        two bases of the sequence. Start cannot be larger than end but can
        be equal.
 Returns : A scalar.
 Args    : A position.

=cut

sub qualat {
   my ($self,$value) = @_;
   if( $self->can('warn') ) {
       $self->warn("Bio::Seq::QualI definition of qualat - implementing class did not provide this method");
   } else {
       warn("Bio::Seq::QualI definition of qualat - implementing class did not provide this method");
   }
   return '';
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
 Usage   : @rev = @{$qual->revcom()};
 Function: Produces a new Bio::Seq::QualI implementing object which
	is reversed from the original quality array.
	The id is the same id as the orginal sequence, and the accession number
	is also indentical. If someone wants to track that this sequence has
	been reversed, it needs to define its own extensions

	To do an inplace edit of an object you can go:

	$qual = $qual->revcom();

	This of course, causes Perl to handle the garbage collection of the old
	object, but it is roughly speaking as efficient as an inplace edit.
 Returns : A new (fresh) Bio::Seq::PrimaryQualI object
 Args    : none

=cut

sub revcom{
   my ($self) = @_;
		# this is the cleanest way
	my @qualities = @{$self->qual()};
	my @reversed_qualities = reverse(@qualities);
   my $seqclass;
   if($self->can_call_new()) {
       $seqclass = ref($self);
   } else {
       $seqclass = 'Bio::Seq::PrimaryQual';
		# Wassat?
		# $self->_attempt_to_load_Seq();
   }
	# the \@reverse_qualities thing works simply because I will it to work.
   my $out = $seqclass->new( '-qual' => \@reversed_qualities,
			     '-display_id'  => $self->display_id,
			     '-accession_number' => $self->accession_number,
			     '-desc' => $self->desc()
			     );
   return $out;
}

=head2 trunc()

 Title   : trunc
 Usage   : $subseq = $myseq->trunc(10,100);
 Function: Provides a truncation of a sequence,
 Returns : a fresh Bio::Seq::QualI implementing object
 Args    : Two integers denoting first and last base of the sub-sequence.


=cut

sub trunc {
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

   my $r_qual = $self->subqual($start,$end);

   my $seqclass;
   if($self->can_call_new()) {
       $seqclass = ref($self);
   } else {
       $seqclass = 'Bio::Seq::PrimaryQual';
		# wassat?
		# $self->_attempt_to_load_Seq();
   }
   my $out = $seqclass->new( '-qual' => $r_qual,
			     '-display_id'  => $self->display_id,
			     '-accession_number' => $self->accession_number,
			     '-desc' => $self->desc()
			     );
   return $out;
}


=head2 translate()

 Title   : translate()
 Usage   : $protein_seq_obj = $dna_seq_obj->translate
           #if full CDS expected:
           $protein_seq_obj = $cds_seq_obj->translate(undef,undef,undef,undef,1);
 Function: Completely useless in this interface.
 Returns : Nothing.
 Args    : Nothing.

=cut


sub translate {
	return 0;
}


=head2 id()

 Title   : id()
 Usage   : $id = $qual->id()
 Function: ID of the quality. This should normally be (and actually is in
           the implementation provided here) just a synonym for display_id().
 Example :
 Returns : A string.
 Args    :


=cut

sub  id {
   my ($self)= @_;
   return $self->display_id();
}

=head2 length()

 Title   : length()
 Usage   : $length = $qual->length();
 Function: Return the length of the array holding the quality values.
        Under most circumstances, this should match the number of quality
        values but no validation is done when the PrimaryQual object is
        constructed and non-digits could be put into this array. Is this a
        bug? Just enough rope...
 Returns : A scalar (the number of elements in the quality array).
 Args    : None.

=cut

sub length {
   my ($self)= @_;
   if( $self->can('throw') ) {
       $self->throw("Bio::Seq::QualI definition of length - implementing class did not provide this method");
   } else {
       confess("Bio::Seq::QualI definition of length - implementing class did not provide this method");
   }
}


=head2 desc()

 Title   : desc()
 Usage   : $qual->desc($newval);
           $description = $seq->desc();
 Function: Get/set description text for a qual object
 Example :
 Returns : value of desc
 Args    : newvalue (optional)

=cut

sub desc {
   my ($self,$value) = @_;
   if( $self->can('warn') ) {
       $self->warn("Bio::Seq::QualI definition of desc - implementing class did not provide this method");
   } else {
       warn("Bio::Seq::QualI definition of desc - implementing class did not provide this method");
   }
   return '';
}

#  These methods are here for backward compatibility with the old, 0.5
#  Seq objects. They all throw warnings that someone is using a
#  deprecated method, and may eventually be removed completely from
#  this object. However, they are important to ease the transition from
#  the old system.

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

   if( $main::{'Bio::Seq::PrimaryQual'} ) {
       return 1;
   } else {
       eval {
	   require Bio::Seq::PrimaryQual;
       };
       if( $@ ) {
	   if( $self->can('throw') ) {
	       $self->throw("Bio::Seq::PrimaryQual could not be loaded for $self\nThis indicates that you are using Bio::Seq::PrimaryQualI without Bio::Seq::PrimaryQual loaded and without providing a complete solution\nThe most likely problem is that there has been a misconfiguration of the bioperl environment\nActual exception\n\n$@\n");
	   } else {
	       confess("Bio::Seq::PrimarySeq could not be loaded for $self\nThis indicates that you are usnig Bio::Seq::PrimaryQualI without Bio::Seq::PrimaryQual loaded and without providing a complete solution\nThe most likely problem is that there has been a misconfiguration of the bioperl environment\nActual exception\n\n$@\n");
	   }
	   return 0;
       }
       return 1;
   }

}


=head2 qualtype()

 Title   : qualtype()
 Usage   : if( $obj->qualtype eq 'phd' ) { /Do Something/ }
 Function: At this time, this function is not used for 
	Bio::Seq::PrimaryQual objects. In fact, now it is a month later and
	I just completed the Bio::Seq::SeqWithQuality object and this is
	definitely deprecated.
 Returns : Nothing. (not implemented)
 Args    : none
 Status  : Virtual


=cut

sub qualtype {
   my ($self,@args) = @_;
   if( $self->can('throw') ) {
	# $self->throw("Bio::Seq::QualI definition of qual - implementing class did not provide this method");
       $self->throw("qualtypetype is not used with quality objects.");
   } else {
	# confess("Bio::Seq::QualI definition of qual - implementing class did not provide this method");
	confess("qualtype is not used with quality objects.");
   }


}




1;
