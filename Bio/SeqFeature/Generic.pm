# $Id$
#
# BioPerl module for Bio::SeqFeature::Generic
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Generic - Generic SeqFeature

=head1 SYNOPSIS

   $feat = new Bio::SeqFeature::Generic ( -start => 10, -end => 100,
				-strand => -1, -primary => 'repeat',
				-source_tag => 'repeatmasker',
				-score  => 1000,
				-tag    => {
				    new => 1,
				    author => 'someone',
				    sillytag => 'this is silly!' } );

   $feat = new Bio::SeqFeature::Generic ( -gff_string => $string );
   # if you want explicitly GFF1
   $feat = new Bio::SeqFeature::Generic ( -gff1_string => $string );

   # add it to an annotated sequence

   $annseq->add_SeqFeature($feat);



=head1 DESCRIPTION

Bio::SeqFeature::Generic is a generic implementation for the
Bio::SeqFeatureI interface, providing a simple object to provide all
the information for a feature on a sequence.

For many Features, this is all you will need to use (for example, this
is fine for Repeats in DNA sequence or Domains in protein
sequence). For other features, which have more structure, this is a
good base class to extend using inheritence to have new things: this
is what is done in the Bio::SeqFeature::Gene,
Bio::SeqFeature::Transcript and Bio::SeqFeature::Exon, which provide
well coordinated classes to represent genes on DNA sequence (for
example, you can get the protein sequence out from a transcript
class).

For many Features, you want to add some piece of information, for
example a common one is that this feature is 'new' whereas other
features are 'old'.  The tag system, which here is implemented using a
hash can be used here.  You can use the tag system to extend the
SeqFeature::Generic programmatically: that is, you know that you have
read in more information into the tag 'mytag' which you can then
retrieve. This means you do not need to know how to write inherieted
Perl to provide more complex information on a feature, and/or, if you
do know but you do not want to write a new class every time you need
some extra piece of information, you can use the tag system to easily
store and then retrieve information.

The tag system can be written in/out of GFF format, and also into EMBL
format via the SeqIO system

=head1 Implemented Interfaces

This class implementes the following interfaces.

=over 4

=item Bio::SeqFeatureI

Note that this includes implementing Bio::RangeI.

=item Bio::AnnotatableI

=item Bio::FeatureHolderI

Features held by a feature are essentially sub-features.

=back

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
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Ewan Birney

Ewan Birney E<lt>birney@sanger.ac.ukE<gt>

=head1 DEVELOPERS

This class has been written with an eye out of inheritence. The fields
the actual object hash are:

   _gsf_tag_hash  = reference to a hash for the tags
   _gsf_sub_array = reference to an array for subfeatures

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Generic;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::SeqFeatureI;
use Bio::AnnotatableI;
use Bio::FeatureHolderI;
use Bio::Annotation::Collection;
use Bio::Location::Simple;
use Bio::Tools::GFF;
use Bio::AlternativeLocationHolderI;
#use Tie::IxHash;

@ISA = qw(Bio::Root::Root Bio::SeqFeatureI 
          Bio::AnnotatableI Bio::FeatureHolderI
	  Bio::AlternativeLocationHolderI);

sub new {
    my ( $caller, @args) = @_;   
    my ($self) = $caller->SUPER::new(@args); 

    $self->{'_parse_h'}       = {};
    $self->{'_gsf_tag_hash'}  = {};
#    tie %{$self->{'_gsf_tag_hash'}}, "Tie::IxHash";

    # bulk-set attributes
    $self->set_attributes(@args);

    # done - we hope
    return $self;
}


=head2 set_attributes

 Title   : set_attributes
 Usage   :
 Function: Sets a whole array of parameters at once.
 Example :
 Returns : none
 Args    : Named parameters, in the form as they would otherwise be passed
           to new(). Currently recognized are:

                  -start          start position
 	 	  -end            end position
 	 	  -strand         strand
 	 	  -primary        primary tag
 	 	  -source         source tag
 	 	  -frame          frame
 	 	  -score          score value
 	 	  -tag            a reference to a tag/value hash
 	 	  -gff_string     GFF v.2 string to initialize from
 	 	  -gff1_string    GFF v.1 string to initialize from
 	 	  -seq_id         the display name of the sequence
 	 	  -annotation     the AnnotationCollectionI object
 	 	  -location       the LocationI object

=cut

sub set_attributes {
    my ($self,@args) = @_;
    my ($start, $end, $strand, $primary_tag, $source_tag, $primary, $source, $frame, 
	$score, $tag, $gff_string, $gff1_string,
	$seqname, $seqid, $annot, $location) =
	    $self->_rearrange([qw(START
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
				  )], @args);
    $location    && $self->location($location);
    $gff_string  && $self->_from_gff_string($gff_string);
    $gff1_string  && do {
	$self->gff_format(Bio::Tools::GFF->new('-gff_version' => 1));
	$self->_from_gff_stream($gff1_string);
    };
    $primary_tag    && $self->primary_tag($primary_tag);
    $source_tag     && $self->source_tag($source_tag);
    $primary        && $self->primary_tag($primary);
    $source         && $self->source_tag($source);
    defined $start  && $self->start($start);
    defined $end    && $self->end($end);
    defined $strand && $self->strand($strand);
    defined $frame  && $self->frame($frame);
    $score          && $self->score($score);
    $annot          && $self->annotation($annot);
    if($seqname) {
	$self->warn("-seqname is deprecated. Please use -seq_id instead.");
	$seqid = $seqname unless $seqid;
    }
    $seqid          && $self->seq_id($seqid);
    $tag            && do {
	foreach my $t ( keys %$tag ) {
	    $self->add_tag_value($t,$tag->{$t});
	}
    };
}


=head2 direct_new

 Title   : direct_new
 Usage   : my $obj = Bio::SeqFeature::Generic->direct_new
 Function: create a blessed hash - for performance improvement in 
           object creation
 Returns : Bio::SeqFeature::Generic object
 Args    : none


=cut

sub direct_new {
    my ( $class) = @_;   
    my ($self) = {};

    bless $self,$class;

    return $self;
}

=head2 location

 Title   : location
 Usage   : my $location = $seqfeature->location()
 Function: returns a location object suitable for identifying location
	   of feature on sequence or parent feature
 Returns : Bio::LocationI object
 Args    : [optional] Bio::LocationI object to set the value to.


=cut

sub location {
    my($self, $value ) = @_;  

    if (defined($value)) {
        unless (ref($value) and $value->isa('Bio::LocationI')) {
	    $self->throw("object $value pretends to be a location but ".
			 "does not implement Bio::LocationI");
        }
        $self->{'_location'} = $value;
    }
    elsif (! $self->{'_location'}) {
        # guarantees a real location object is returned every time
        $self->{'_location'} = Bio::Location::Simple->new();
    }
    return $self->{'_location'};
}

=head2 add_alternative_locations

 Title   : add_alternative_locations
 Usage   : $seqfeature->add_alternative_locations(@locationi)
 Function: adds new alternative location to the object
 Returns : void
 Args    : one or more LocationI-implementing object

 This adds one or more alternative locations to the feature.  These are
 to be viewed as alternative coordinate systems, such as
 assembly-to-assembly alignments, and not as alternative locations in
 the same coordinate space.

=cut

sub add_alternative_locations {
    my $self = shift;
    foreach (@_) {
      $self->throw("object $_ pretends to be a location but ".
			 "does not implement Bio::LocationI")
          unless ref($_) and $_->isa('Bio::LocationI');
      push @{$self->{'_gsf_alternative_locations'}},$_;
    }
}

*add_alternative_location = \&add_alternative_locations;

=head2 alternative_locations

 Title   : alternative_locations
 Usage   : @locations = $seqfeature->alternative_locations([$seq_id])
 Function: returns alternative locations
 Returns : list of alternative locations
 Args    : optionally, a seq_id to filter on

=cut

sub alternative_locations {
    my $self = shift;
    my $seqid_filter = shift;
    return unless $self->{'_gsf_alternative_locations'};
    if ( $seqid_filter ) {
       return grep {$seqid_filter eq $_->seq_id} @{$self->{'_gsf_alternative_locations'}};
    } else {
       return @{$self->{'_gsf_alternative_locations'}};
    }
}

=head2 clear_alternative_locations

 Title   : clear_alternative_locations
 Usage   : $seqfeature->clear_alternative_locations([$seqid])
 Function: clears all alternative locations
 Returns : void
 Args    : optionally, a seq_id to clear locations on

=cut

sub clear_alternative_locations {
    my $self = shift;
    my $seqid_filter = shift;
    return unless $self->{'_gsf_alternative_locations'};
    if ( $seqid_filter ) {
       my @locations = grep {$seqid_filter ne $_->seq_id} @{$self->{'_gsf_alternative_locations'}};
       return $self->{'_gsf_alternative_locations'} = \@locations;
    } else {
       $self->{'_gsf_alternative_locations'} = [];
   }
}

=head2 start

 Title   : start
 Usage   : $start = $feat->start
           $feat->start(20)
 Function: Get/set on the start coordinate of the feature
 Returns : integer
 Args    : none


=cut

sub start {
   my ($self,$value) = @_;
   return $self->location->start($value);
}

=head2 end

 Title   : end
 Usage   : $end = $feat->end
           $feat->end($end)
 Function: get/set on the end coordinate of the feature
 Returns : integer
 Args    : none


=cut

sub end {
   my ($self,$value) = @_;
   return $self->location->end($value);
}

=head2 length

 Title   : length
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


=cut

sub length {
   my ($self) = @_;
   return $self->end - $self->start() + 1;
}

=head2 strand

 Title   : strand
 Usage   : $strand = $feat->strand()
           $feat->strand($strand)
 Function: get/set on strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : none


=cut

sub strand {
   my ($self,$value) = @_;
   return $self->location->strand($value);
}

=head2 seq_id

 Title   : seq_id
 Usage   : $obj->seq_id($newval)
 Function: This method returns the ID of the sequence that the
           start, end and strand are relative to.  It is a simple
           passthrough to $obj->location->seq_id().

           This attribute *should* be used in GFF dumping.
 Returns : value of unique_id
 Args    : newvalue (optional)


=cut

sub seq_id {
    my $obj = shift;
    $obj->location->seq_id(@_);
}

=head2 score

 Title   : score
 Usage   : $score = $feat->score()
           $feat->score($score)
 Function: get/set on score information
 Returns : float
 Args    : none if get, the new value if set


=cut

sub score {
  my ($self,$value) = @_;

  if (defined($value)) {
       if ( $value !~ /^[+-]?\d+\.?\d*(e-\d+)?/ ) {
	   $self->throw("'$value' is not a valid score");
       }
       $self->{'_gsf_score'} = $value;
  }

  return $self->{'_gsf_score'};
}

=head2 frame

 Title   : frame
 Usage   : $frame = $feat->frame()
           $feat->frame($frame)
 Function: get/set on frame information
 Returns : 0,1,2, '.'
 Args    : none if get, the new value if set


=cut

sub frame {
  my ($self,$value) = @_;

  if ( defined $value ) {
       if ( $value !~ /^[0-2.]$/ ) {
	   $self->throw("'$value' is not a valid frame");
       }
       if( $value eq '.' ) { $value = '.'; } 
       $self->{'_gsf_frame'} = $value;
  }
  return $self->{'_gsf_frame'};
}

=head2 primary_tag

 Title   : primary_tag
 Usage   : $tag = $feat->primary_tag()
           $feat->primary_tag('exon')
 Function: get/set on the primary tag for a feature,
           eg 'exon'
 Returns : a string
 Args    : none


=cut

sub primary_tag {
   my ($self,$value) = @_;
   if ( defined $value ) {
       $self->{'_primary_tag'} = $value;
   }
   return $self->{'_primary_tag'};
}

=head2 source_tag

 Title   : source_tag
 Usage   : $tag = $feat->source_tag()
           $feat->source_tag('genscan');
 Function: Returns the source tag for a feature,
           eg, 'genscan'
 Returns : a string
 Args    : none


=cut

sub source_tag {
   my ($self,$value) = @_;

   if( defined $value ) {
       $self->{'_source_tag'} = $value;
   }
   return $self->{'_source_tag'};
}

=head2 has_tag

 Title   : has_tag
 Usage   : $value = $self->has_tag('some_tag')
 Function: Tests wether a feature contaings a tag
 Returns : TRUE if the SeqFeature has the tag,
           and FALSE otherwise.
 Args    : The name of a tag


=cut

sub has_tag {
    my ($self, $tag) = @_;
    return exists $self->{'_gsf_tag_hash'}->{$tag};
}

=head2 add_tag_value

 Title   : add_tag_value
 Usage   : $self->add_tag_value('note',"this is a note");
 Returns : TRUE on success
 Args    : tag (string) and value (a scalar or list)


=cut

sub add_tag_value {
    my ($self, $tag, @value) = @_;
    $self->{'_gsf_tag_hash'}->{$tag} ||= [];
    push(@{$self->{'_gsf_tag_hash'}->{$tag}},@value);
}


=head2 get_tag_values

 Title   : get_tag_values
 Usage   : @values = $gsf->get_tag_values('note');
 Function: Returns a list of all the values stored
           under a particular tag.
 Returns : A list of scalars
 Args    : The name of the tag


=cut

sub get_tag_values {
   my ($self, $tag) = @_;

   if( ! defined $tag ) { return (); }
   if ( ! exists $self->{'_gsf_tag_hash'}->{$tag} ) {
       $self->throw("asking for tag value that does not exist $tag");
   }
   return @{$self->{'_gsf_tag_hash'}->{$tag}};
}


=head2 get_all_tags

 Title   : get_all_tags
 Usage   : @tags = $feat->get_all_tags()
 Function: Get a list of all the tags in a feature
 Returns : An array of tag names
 Args    : none


=cut

sub get_all_tags {
   my ($self, @args) = @_;   
   return keys %{ $self->{'_gsf_tag_hash'}};
}

=head2 remove_tag

 Title   : remove_tag
 Usage   : $feat->remove_tag('some_tag')
 Function: removes a tag from this feature
 Returns : the array of values for this tag before removing it
 Args    : tag (string)


=cut

sub remove_tag {
   my ($self, $tag) = @_;

   if ( ! exists $self->{'_gsf_tag_hash'}->{$tag} ) {
       $self->throw("trying to remove a tag that does not exist: $tag");
   }
   my @vals = @{$self->{'_gsf_tag_hash'}->{$tag}};
   delete $self->{'_gsf_tag_hash'}->{$tag};
   return @vals;
}

=head2 attach_seq

 Title   : attach_seq
 Usage   : $sf->attach_seq($seq)
 Function: Attaches a Bio::Seq object to this feature. This
           Bio::Seq object is for the *entire* sequence: ie
           from 1 to 10000
 Example :
 Returns : TRUE on success
 Args    : a Bio::PrimarySeqI compliant object


=cut

sub attach_seq {
   my ($self, $seq) = @_;

   if ( ! ($seq && ref($seq) && $seq->isa("Bio::PrimarySeqI")) ) {
       $self->throw("Must attach Bio::PrimarySeqI objects to SeqFeatures");
   }

   $self->{'_gsf_seq'} = $seq;

   # attach to sub features if they want it
   foreach ( $self->sub_SeqFeature() ) {
       $_->attach_seq($seq);
   }

   return 1;
}

=head2 seq

 Title   : seq
 Usage   : $tseq = $sf->seq()
 Function: returns the truncated sequence (if there) for this
 Example :
 Returns : sub seq (a Bio::PrimarySeqI compliant object) on attached sequence
           bounded by start & end, or undef if there is no sequence attached
 Args    : none


=cut

sub seq {
   my ($self, $arg) = @_;

   if ( defined $arg ) {
       $self->throw("Calling SeqFeature::Generic->seq with an argument. You probably want attach_seq");
   }

   if ( ! exists $self->{'_gsf_seq'} ) {
       return undef;
   }

   # assumming our seq object is sensible, it should not have to yank
   # the entire sequence out here.

   my $seq = $self->{'_gsf_seq'}->trunc($self->start(), $self->end());


   if ( $self->strand == -1 ) {

       # ok. this does not work well (?)
       #print STDERR "Before revcom", $seq->str, "\n";
       $seq = $seq->revcom;
       #print STDERR "After  revcom", $seq->str, "\n";
   }

   return $seq;
}

=head2 entire_seq

 Title   : entire_seq
 Usage   : $whole_seq = $sf->entire_seq()
 Function: gives the entire sequence that this seqfeature is attached to
 Example :
 Returns : a Bio::PrimarySeqI compliant object, or undef if there is no
           sequence attached
 Args    :


=cut

sub entire_seq {
   my ($self) = @_;

   return $self->{'_gsf_seq'};
}


=head2 unique_id

 Title   : unique_id
 Usage   : $obj->unique_id($newval)
 Function: This is a unique identifier that identifies this object.
           If not set, will return the memory location.

           This attribute should *not* be used in GFF dumping, as
           that should come from the collection in which the seq
           feature was found or from seq_id().
 Returns : value of unique_id
 Args    : newvalue (optional)


=cut

sub unique_id {
    my ($obj,$value) = @_;
    if ( defined $value ) {
	$obj->{'_gsf_unique_id'} = $value;
    }
    return $obj->{'_gsf_unique_id'};
}

=head2 display_name

 Title   : display_name
 Usage   : $featname = $obj->display_name
 Function: Implements the display_name() method, which is a human-readable
           name for the feature. 
 Returns : value of display_name (a string)
 Args    : Optionally, on set the new value or undef 

=cut

sub display_name{
    my $self = shift;

    return $self->{'display_name'} = shift if @_;
    return $self->{'display_name'};
}

=head1 Methods for implementing Bio::AnnotatableI

=cut

=head2 annotation

 Title   : annotation
 Usage   : $obj->annotation($annot_obj)
 Function: Get/set the annotation collection object for annotating this
           feature.

 Example : 
 Returns : A Bio::AnnotationCollectionI object
 Args    : newvalue (optional)


=cut

sub annotation {
    my ($obj,$value) = @_;

    # we are smart if someone references the object and there hasn't been
    # one set yet
    if(defined $value || ! defined $obj->{'annotation'} ) {
        $value = new Bio::Annotation::Collection unless ( defined $value );
        $obj->{'annotation'} = $value;
    }
    return $obj->{'annotation'};
}

=head1 Methods to implement Bio::FeatureHolderI

This includes methods for retrieving, adding, and removing
features. Since this is already a feature, features held by this
feature holder are essentially sub-features.

=cut

=head2 get_SeqFeatures

 Title   : get_SeqFeatures
 Usage   : @feats = $feat->get_SeqFeatures();
 Function: Returns an array of sub Sequence Features
 Returns : An array
 Args    : none


=cut

sub get_SeqFeatures {
    my ($self) = @_;

    if ($self->{'_gsf_sub_array'}) {
        return @{$self->{'_gsf_sub_array'}};
    } else {
        return;
    }
}

=head2 add_SeqFeature

 Title   : add_SeqFeature
 Usage   : $feat->add_SeqFeature($subfeat);
           $feat->add_SeqFeature($subfeat,'EXPAND')
 Function: adds a SeqFeature into the subSeqFeature array.
           with no 'EXPAND' qualifer, subfeat will be tested
           as to whether it lies inside the parent, and throw
           an exception if not.

           If EXPAND is used, the parent's start/end/strand will
           be adjusted so that it grows to accommodate the new
           subFeature
 Returns : nothing
 Args    : An object which has the SeqFeatureI interface


=cut

#'
sub add_SeqFeature{
    my ($self,$feat,$expand) = @_;

    if ( !$feat->isa('Bio::SeqFeatureI') ) {
        $self->warn("$feat does not implement Bio::SeqFeatureI. Will add it anyway, but beware...");
    }

    if($expand && ($expand eq 'EXPAND')) {
        $self->_expand_region($feat);
    } else {
        if ( !$self->contains($feat) ) {
	    $self->throw("$feat is not contained within parent feature, and expansion is not valid");
        }
    }

    $self->{'_gsf_sub_array'} = [] unless exists($self->{'_gsf_sub_array'});
    push(@{$self->{'_gsf_sub_array'}},$feat);

}

=head2 remove_SeqFeatures

 Title   : remove_SeqFeatures
 Usage   : $sf->remove_SeqFeatures
 Function: Removes all sub SeqFeatures

           If you want to remove only a subset, remove that subset from the
           returned array, and add back the rest.

 Example :
 Returns : The array of Bio::SeqFeatureI implementing sub-features that was
           deleted from this feature.
 Args    : none


=cut

sub remove_SeqFeatures {
   my ($self) = @_;

   my @subfeats = @{$self->{'_gsf_sub_array'}};
   $self->{'_gsf_sub_array'} = []; # zap the array implicitly.
   return @subfeats;
}

=head1 GFF-related methods

=cut

=head2 gff_format

 Title   : gff_format
 Usage   : # get:
           $gffio = $feature->gff_format();
           # set (change the default version of GFF2):
           $feature->gff_format(Bio::Tools::GFF->new(-gff_version => 1));
 Function: Get/set the GFF format interpreter. This object is supposed to 
           format and parse GFF. See Bio::Tools::GFF for the interface.

           If this method is called as class method, the default for all
           newly created instances will be changed. Otherwise only this
           instance will be affected.
 Example : 
 Returns : a Bio::Tools::GFF compliant object
 Args    : On set, an instance of Bio::Tools::GFF or a derived object.


=cut

sub gff_format {
    my ($self, $gffio) = @_;

    if(defined($gffio)) {
	if(ref($self)) {
	    $self->{'_gffio'} = $gffio;
	} else {
	    $Bio::SeqFeatureI::static_gff_formatter = $gffio;
	}
    }
    return (ref($self) && exists($self->{'_gffio'}) ?
	    $self->{'_gffio'} : $self->_static_gff_formatter);
}

=head2 gff_string

 Title   : gff_string
 Usage   : $str = $feat->gff_string;
           $str = $feat->gff_string($gff_formatter);
 Function: Provides the feature information in GFF format.

           We override this here from Bio::SeqFeatureI in order to use the
           formatter returned by gff_format().

 Returns : A string
 Args    : Optionally, an object implementing gff_string().


=cut

sub gff_string{
   my ($self,$formatter) = @_;

   $formatter = $self->gff_format() unless $formatter;
   return $formatter->gff_string($self);
}

#  =head2 slurp_gff_file
#
#   Title   : slurp_file
#   Usage   : @features = Bio::SeqFeature::Generic::slurp_gff_file(\*FILE);
#   Function: Sneaky function to load an entire file as in memory objects.
#             Beware of big files.
#
#             This method is deprecated. Use Bio::Tools::GFF instead, which can
#             also handle large files.
#
#   Example :
#   Returns :
#   Args    :
#
#  =cut

sub slurp_gff_file {
   my ($f) = @_;
   my @out;
   if ( !defined $f ) {
       die "Must have a filehandle";
   }

   Bio::Root::Root->warn("deprecated method slurp_gff_file() called in Bio::SeqFeature::Generic. Use Bio::Tools::GFF instead.");
  
   while(<$f>) {

       my $sf = Bio::SeqFeature::Generic->new('-gff_string' => $_);
       push(@out, $sf);
   }

   return @out;

}

=head2 _from_gff_string

 Title   : _from_gff_string
 Usage   :
 Function: Set feature properties from GFF string. 

           This method uses the object returned by gff_format() for the
           actual interpretation of the string. Set a different GFF format
           interpreter first if you need a specific version, like GFF1. (The
           default is GFF2.)
 Example :
 Returns : 
 Args    : a GFF-formatted string


=cut

sub _from_gff_string {
   my ($self, $string) = @_;

   $self->gff_format()->from_gff_string($self, $string);
}


=head2 _expand_region

 Title   : _expand_region
 Usage   : $self->_expand_region($feature);
 Function: Expand the total region covered by this feature to
           accomodate for the given feature.

           May be called whenever any kind of subfeature is added to this
           feature. add_sub_SeqFeature() already does this.
 Returns : 
 Args    : A Bio::SeqFeatureI implementing object.


=cut

sub _expand_region {
    my ($self, $feat) = @_;
    if(! $feat->isa('Bio::SeqFeatureI')) {
	$self->warn("$feat does not implement Bio::SeqFeatureI");
    }
    # if this doesn't have start/end set - forget it!
    if((! defined($self->start())) && (! defined $self->end())) {
	$self->start($feat->start());
	$self->end($feat->end());
	$self->strand($feat->strand) unless defined($self->strand());
    } else {
	my $range = $self->union($feat);
	$self->start($range->start);
	$self->end($range->end);
	$self->strand($range->strand);
    }
}

=head2 _parse

 Title   : _parse
 Usage   :
 Function: Parsing hints
 Example :
 Returns :
 Args    :


=cut

sub _parse {
   my ($self) = @_;

   return $self->{'_parse_h'};
}

=head2 _tag_value

 Title   : _tag_value
 Usage   : 
 Function: For internal use only. Convenience method for those tags that
           may only have a single value.
 Returns : 
 Args    : 


=cut

sub _tag_value {
    my ($self, $tag, $value) = @_;

    if(defined($value) || (! $self->has_tag($tag))) {
	$self->remove_tag($tag) if($self->has_tag($tag));
	$self->add_tag_value($tag, $value);
    }
    return ($self->each_tag_value($tag))[0];
}

#######################################################################
# aliases for methods that changed their names in an attempt to make  #
# bioperl names more consistent                                       #
#######################################################################

sub seqname {
    my $self = shift;
    $self->warn("SeqFeatureI::seqname() is deprecated. Please use seq_id() instead.");
    return $self->seq_id(@_);
}

sub display_id {
    my $self = shift;
    $self->warn("SeqFeatureI::display_id() is deprecated. Please use display_name() instead.");
    return $self->display_name(@_);
}

# this is towards consistent naming
sub each_tag_value { return shift->get_tag_values(@_); }
sub all_tags { return shift->get_all_tags(@_); }

# we revamped the feature containing property to implementing
# Bio::FeatureHolderI
*sub_SeqFeature = \&get_SeqFeatures;
*add_sub_SeqFeature = \&add_SeqFeature;
*flush_sub_SeqFeatures = \&remove_SeqFeatures;
# this one is because of inconsistent naming ...
*flush_sub_SeqFeature = \&remove_SeqFeatures;


1;
