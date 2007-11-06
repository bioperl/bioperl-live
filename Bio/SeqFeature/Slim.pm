# $Id$
#
# BioPerl module for Bio::SeqFeature::Slim
#
# Cared for by Jason Stajich <jason_AT_bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Slim - A very lightweight Bio::SeqFeatureI implementation

=head1 SYNOPSIS

use Bio::SeqFeature::Slim;

=head1 DESCRIPTION

Describe the object here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich

Email jason_AT_bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Slim;
use strict;
use base 'Bio::SeqFeatureI';

use constant {
    SEQ_ID      => 0,
    SOURCE      => 1,
    PRIMARY     => 2,
    START       => 3,
    STOP        => 4,
    SCORE       => 5,
    STRAND      => 6,
    FRAME       => 7,
    TAGS        => 8,
    NAME        => 9,
    PARENT      => 10,
    GFF_TYPE    => 11,
    SUBFEATURES => 12,
    SEQ_OBJ     => 13,
};


=head2 new

 Title   : new
 Usage   : my $obj = new Bio::SeqFeature::Slim();
 Function: Builds a new Bio::SeqFeature::Slim object 
 Returns : an instance of Bio::SeqFeature::Slim
 Args    :


=cut

sub new {
  my($class) = shift;
  my ($start, $end, $strand, $primary_tag, $source_tag, $primary, 
      $source, $frame, $score, $tag, $gff_string, $gff1_string,
      $seqname, $seqid, $annot, $location,$display_name) =
	  Bio::Root::RootI->_rearrange([qw
					(START
					 END
					 STRAND
					 PRIMARY_TAG
					 SOURCE_TAG
					 PRIMARY
					 SOURCE
					 FRAME
					 SCORE
					 TAG
					 GFF_STRING
					 GFF1_STRING
					 SEQNAME
					 SEQ_ID
					 ANNOTATION
					 LOCATION
					 DISPLAY_NAME
					)], @_);
  if( defined $primary_tag && defined $primary ) {
      Bio::Root::RootI->warn("Both primary and primary_tag are defined, only use one");
  } 
  if( defined $source_tag && defined $source ) {
      Bio::Root::RootI->warn("Both source and source_tag are defined, only use one");
  } 

  $primary_tag = $primary if defined $primary && ! defined $primary_tag;
  $source_tag  = $source  if defined $source && ! defined $source_tag;
  my $self = bless [$seqid,        #0
		    $source_tag,   #1
		    $primary_tag,  #2
		    $start,        #3
		    $end,          #4
		    $score,        #5
		    $strand,       #6
		    $frame,        #7
		    {},            #8 tags
		    $display_name, #9 display name
		    undef,         #10 parent
		    undef,         #11 gff_type
		    [],            #12 seqfeatures
		    undef,         #13 seqobj
      ], $class;

  $tag            && do {
      foreach my $t ( keys %$tag ) {
	  $self->add_tag_value($t, UNIVERSAL::isa($tag->{$t}, "ARRAY") ? 
			       @{$tag->{$t}} : $tag->{$t});
      }
  };
  return $self;
}

=head1 Bio::SeqFeatureI specific methods

New method interfaces.

=cut

=head2 get_SeqFeatures

 Title   : get_SeqFeatures
 Usage   : @feats = $feat->get_SeqFeatures();
 Function: Returns an array of sub Sequence Features
 Returns : An array
 Args    : none

=cut

sub get_SeqFeatures{
   return @{shift->[SUBFEATURES] || []};
}

=head2 display_name

 Title   : display_name
 Usage   : $name = $feat->display_name()
 Function: Returns the human-readable name of the feature for displays.
 Returns : a string
 Args    : none

=cut

sub display_name {
    my ($self) = shift;
    if( @_) {
	($self->[NAME]) = shift @_;
    }
    return $self->[NAME];
}

=head2 primary_tag

 Title   : primary_tag
 Usage   : $tag = $feat->primary_tag()
 Function: Returns the primary tag for a feature,
           eg 'exon'
 Returns : a string
 Args    : none


=cut

sub primary_tag{
   my ($self) = shift;
    if( @_) {
	($self->[PRIMARY]) = shift @_;
    }
    return $self->[PRIMARY];
}

=head2 source_tag

 Title   : source_tag
 Usage   : $tag = $feat->source_tag()
 Function: Returns the source tag for a feature,
           eg, 'genscan'
 Returns : a string
 Args    : none


=cut

sub source_tag{
   my ($self) = shift;
    if( @_) {
	($self->[SOURCE]) = shift @_;
    }
    return $self->[SOURCE];
}

=head2 frame

 Title   : frame
 Usage   : $frame = $feat->frame()
 Function: Returns the frame for a feature,
           eg, '1'
 Returns : '.', 0,1,2
 Args    : none


=cut

sub frame{
   my ($self) = shift;
    if( @_) {
	($self->[FRAME]) = shift @_;
    }
    return $self->[FRAME];
}

=head2 has_tag

 Title   : has_tag
 Usage   : $tag_exists = $self->has_tag('some_tag')
 Function: 
 Returns : TRUE if the specified tag exists, and FALSE otherwise
 Args    :


=cut

sub has_tag{
   my ($self,$tag) = @_;
   return unless defined $tag;   
   return exists($self->[TAGS]->{$tag});
}

=head2 get_tag_values

 Title   : get_tag_values
 Usage   : @values = $self->get_tag_values('some_tag')
 Function: 
 Returns : An array comprising the values of the specified tag.
 Args    : a string

throws an exception if there is no such tag

=cut

sub get_tag_values {
   my ($self,$tag) = @_;
   return unless defined $tag; 
   return $self->[TAGS]->{$tag};
}

=head2 get_tagset_values

 Title   : get_tagset_values
 Usage   : @values = $self->get_tagset_values(qw(label transcript_id product))
 Function: 
 Returns : An array comprising the values of the specified tags, in order of tags
 Args    : An array of strings

does NOT throw an exception if none of the tags are not present

this method is useful for getting a human-readable label for a
SeqFeatureI; not all tags can be assumed to be present, so a list of
possible tags in preferential order is provided

=cut

# interface + abstract method
sub get_tagset_values {
    my ($self, @args) = @_;
    my @vals = ();
    foreach my $arg (@args) {
        if ($self->has_tag($arg)) {
            push(@vals, $self->get_tag_values($arg));
        }
    }
    return @vals;
}

=head2 get_all_tags

 Title   : get_all_tags
 Usage   : @tags = $feat->get_all_tags()
 Function: gives all tags for this feature
 Returns : an array of strings
 Args    : none


=cut

sub get_all_tags{
    my ($self) = shift;
   return keys %{$self->[TAGS] || {}};
}

=head2 attach_seq

 Title   : attach_seq
 Usage   : $sf->attach_seq($seq)
 Function: Attaches a Bio::Seq object to this feature. This
           Bio::Seq object is for the *entire* sequence: ie
           from 1 to 10000

           Note that it is not guaranteed that if you obtain a feature from
           an object in bioperl, it will have a sequence attached. Also,
           implementors of this interface can choose to provide an empty
           implementation of this method. I.e., there is also no guarantee
           that if you do attach a sequence, seq() or entire_seq() will not
           return undef.

           The reason that this method is here on the interface is to enable
           you to call it on every SeqFeatureI compliant object, and
           that it will be implemented in a useful way and set to a useful
           value for the great majority of use cases. Implementors who choose
           to ignore the call are encouraged to specifically state this in
           their documentation.

 Example :
 Returns : TRUE on success
 Args    : a Bio::PrimarySeqI compliant object


=cut

sub attach_seq {
    my ($self) = shift;
    if(@_) {
	$self->[SEQ_OBJ] = shift @_;
	return 1 if defined $self->[SEQ_OBJ];
    }
    return 0;
}

=head2 seq

 Title   : seq
 Usage   : $tseq = $sf->seq()
 Function: returns the truncated sequence (if there is a sequence attached)
           for this feature
 Example :
 Returns : sub seq (a Bio::PrimarySeqI compliant object) on attached sequence
           bounded by start & end, or undef if there is no sequence attached
 Args    : none


=cut

sub seq {
    my ($self) = shift;
    if(defined $self->[SEQ_OBJ] ) {
	if( ! ref($self->[SEQ_OBJ]) ||
	    ! $self->[SEQ_OBJ]->isa('Bio::PrimarySeqI') ) {
	    $self->throw("Have a seq_obj which is not Bio::PrimarySeqI compliant");
	} else {
	    return $self->[SEQ_OBJ]->trunc($self->start, $self->end);
	}
    }
    return undef;
}

=head2 entire_seq

 Title   : entire_seq
 Usage   : $whole_seq = $sf->entire_seq()
 Function: gives the entire sequence that this seqfeature is attached to
 Example :
 Returns : a Bio::PrimarySeqI compliant object, or undef if there is no
           sequence attached
 Args    : none


=cut

sub entire_seq {
    my ($self) = shift;
    if(defined $self->[SEQ_OBJ] ) {
	if( ! ref($self->[SEQ_OBJ]) ||
	    ! $self->[SEQ_OBJ]->isa('Bio::PrimarySeqI') ) {
	    $self->throw("Have a seq_obj which is not Bio::PrimarySeqI compliant");
	} else {
	    return $self->[SEQ_OBJ];
	}
    }
    return undef;
}


=head2 seq_id

 Title   : seq_id
 Usage   : $obj->seq_id($newval)
 Function: There are many cases when you make a feature that you
           do know the sequence name, but do not know its actual
           sequence. This is an attribute such that you can store
           the ID (e.g., display_id) of the sequence.

           This attribute should *not* be used in GFF dumping, as
           that should come from the collection in which the seq
           feature was found.
 Returns : value of seq_id
 Args    : newvalue (optional)


=cut

sub seq_id {
    my ($self) = shift;
    if( @_) {
	($self->[SEQ_ID]) = shift @_;
    }
    return $self->[SEQ_ID];
}

=head2 gff_string

 Title   : gff_string
 Usage   : $str = $feat->gff_string;
           $str = $feat->gff_string($gff_formatter);
 Function: Provides the feature information in GFF format.

           The implementation provided here returns GFF2 by default. If you
           want a different version, supply an object implementing a method
           gff_string() accepting a SeqFeatureI object as argument. E.g., to
           obtain GFF1 format, do the following:

                my $gffio = Bio::Tools::GFF->new(-gff_version => 1);
                $gff1str = $feat->gff_string($gff1io);

 Returns : A string
 Args    : Optionally, an object implementing gff_string().


=cut

=head1 Decorating methods

These methods have an implementation provided by Bio::SeqFeatureI,
but can be validly overwritten by subclasses

=head2 spliced_seq

  Title   : spliced_seq

  Usage   : $seq = $feature->spliced_seq()
            $seq = $feature_with_remote_locations->spliced_seq($db_for_seqs)

  Function: Provides a sequence of the feature which is the most
            semantically "relevant" feature for this sequence. A default
            implementation is provided which for simple cases returns just
            the sequence, but for split cases, loops over the split location
            to return the sequence. In the case of split locations with
            remote locations, eg

            join(AB000123:5567-5589,80..1144)

            in the case when a database object is passed in, it will attempt
            to retrieve the sequence from the database object, and "Do the right thing",
            however if no database object is provided, it will generate the correct
            number of N's (DNA) or X's (protein, though this is unlikely).

            This function is deliberately "magical" attempting to second guess
            what a user wants as "the" sequence for this feature.

            Implementing classes are free to override this method with their
            own magic if they have a better idea what the user wants.

  Args    : [optional]
            -db        A L<Bio::DB::RandomAccessI> compliant object if
                       one needs to retrieve remote seqs.
            -nosort    boolean if the locations should not be sorted
                       by start location.  This may occur, for instance,
                       in a circular sequence where a gene span starts
                       before the end of the sequence and ends after the
                       sequence start. Example : join(15685..16260,1..207)
	    -phase     truncates the returned sequence based on the
	               intron phase (0,1,2).
	    
  Returns : A L<Bio::PrimarySeqI> object

=cut

=head2 location

 Title   : location
 Usage   : my $location = $seqfeature->location()
 Function: returns a location object suitable for identifying location
	   of feature on sequence or parent feature
          NOTE: in the implementation location is READ-ONLY!
           and complicated locations can not be represented because
           this is intended to be used with GFF generated locations which will
           always only be start..stop
 Returns : Bio::LocationI object
 Args    : none


=cut

sub location {
   my ($self) = @_;
   if( @_ ) {
       $self->warn("this implementation does not let setting of LOCATION obj\n");
       return undef;
   }
   # somewhat silly - maybe we should cache this?
   Bio::Location::Simple->new(-start => $self->start,
			      -end   => $self->end,
			      -strand=> $self->strand);
}


=head2 primary_id

 Title   : primary_id
 Usage   : $obj->primary_id($newval)
 Function:
 Example :
 Returns : value of primary_id (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

Primary ID is a synonym for the tag 'ID'

=cut

sub primary_id{
    my $self = shift;
    # note from cjm@fruitfly.org:
    # I have commented out the following 2 lines:

    #return $self->{'primary_id'} = shift if @_;
    #return $self->{'primary_id'};

    #... and replaced it with the following; see
    # http://bioperl.org/pipermail/bioperl-l/2003-December/014150.html
    # for the discussion that lead to this change

    if (@_) {
        if ($self->has_tag('ID')) {
            $self->remove_tag('ID');
        }
        $self->add_tag_value('ID', shift);
    }
    my ($id) = $self->get_tagset_values('ID');
    return $id;
}

sub generate_unique_persistent_id {
    # DEPRECATED - us IDHandler
    my $self = shift;
    require "Bio/SeqFeature/Tools/IDHandler.pm";
    Bio::SeqFeature::Tools::IDHandler->new->generate_unique_persistent_id($self);
}

=head1 Bio::RangeI methods

These methods are inherited from RangeI and can be used
directly from a SeqFeatureI interface. Remember that a
SeqFeature is-a RangeI, and so wherever you see RangeI you
can use a feature ($r in the below documentation).

=cut

=head2 start()

 See L<Bio::RangeI>

=cut

sub start {
    my ($self) = shift;
    if( @_) {
	($self->[START]) = shift @_;
    }
    return $self->[START];  
}

=head2 end()

 See L<Bio::RangeI>

=cut

sub end {
    my ($self) = shift;
    if( @_) {
	($self->[STOP]) = shift @_;
    }
    return $self->[STOP];  
}

=head2 strand()

 See L<Bio::RangeI>

=cut

sub strand {
    my ($self) = shift;
    if( @_) {
	($self->[STRAND]) = shift @_;
    }
    return $self->[STRAND];  
}

=head2 length

  Title   : length
  Usage   : $length = $range->length();
  Function: get/set the length of this range
  Returns : the length of this range
  Args    : optionally allows the length to be set
             using $range->length($length)

=cut

sub length {
  my $self = shift;
  if(@_) {
    $self->warn( ref($self). "->length() is read-only");
  }
  return abs($self->end - $self->start) + 1;
}
=head2 overlaps()

 See L<Bio::RangeI>

=head2 contains()

 See L<Bio::RangeI>

=head2 equals()

 See L<Bio::RangeI>

=head2 intersection()

 See L<Bio::RangeI>

=head2 union()

 See L<Bio::RangeI>

=head1 Bio::AnnotatableI methods

=cut

=head2 add_tag_value

 Title   : add_tag_value
 Usage   : $self->add_tag_value('note',"this is a note");
 Returns : TRUE on success
 Args    : tag (string) and one or more values (any scalar(s))


=cut

sub add_tag_value{
    my $self = shift;
    my $tag = shift;
    $self->[TAGS] ||= [];
    push (@{$self->[TAGS]->{$tag}},@_);
}


=head2 create_seqfeature_generic

 Title   : create_seqfeature_generic
 Usage   : my $feat = $slimfeat->create_seqfeature_generic
 Function: Create a Bio::SeqFeature::Generic object from this Slim object
 Returns : L<Bio::SeqFeature::Generic>
 Args    : None


=cut

sub create_seqfeature_generic{
   my ($self) = shift;
   return Bio::SeqFeature::Generic->new(-location   => $self->location,
					-score      => $self->score,
					-source_tag => $self->source_tag,
					-primary_tag=> $self->primary_tag,
					-frame      => $self->frame,
					-tag        => $self->[TAGS],
					-display_name=> $self->display_name,
       );
}

1;
