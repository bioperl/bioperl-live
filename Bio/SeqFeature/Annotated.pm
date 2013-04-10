#
# BioPerl module for Bio::SeqFeature::Annotated
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Allen Day <allenday at ucla.edu>
#
# Copyright Allen Day
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Annotated - PLEASE PUT SOMETHING HERE

=head1 SYNOPSIS

    # none yet, complain to authors

=head1 DESCRIPTION

None yet, complain to authors.

=head1 Implemented Interfaces

This class implements the following interfaces.

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
the bugs and their resolution.  Bug reports can be submitted via 
the web:

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR - Allen Day

Allen Day E<lt>allenday at ucla.eduE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


package Bio::SeqFeature::Annotated;

use strict;

use Bio::Annotation::Collection;
use Bio::Annotation::OntologyTerm;
use Bio::Annotation::Target;
use Bio::LocatableSeq;
use Bio::Location::Simple;
use Bio::Ontology::OntologyStore;
use Bio::Tools::GFF;
use Bio::SeqFeature::AnnotationAdaptor;
use Data::Dumper;
use URI::Escape;

use base qw(Bio::Root::Root
    Bio::SeqFeature::TypedSeqFeatureI
    Bio::AnnotatableI
    Bio::FeatureHolderI);

our %tagclass = (
  comment        => 'Bio::Annotation::Comment',
  dblink         => 'Bio::Annotation::DBLink',
  description    => 'Bio::Annotation::SimpleValue',
  gene_name      => 'Bio::Annotation::SimpleValue',
  ontology_term  => 'Bio::Annotation::OntologyTerm',
  reference      => 'Bio::Annotation::Reference',
  __DEFAULT__    => 'Bio::Annotation::SimpleValue',
);

our %tag2text = (
  'Bio::Annotation::Comment'        => 'text',
  'Bio::Annotation::DBLink'         => 'primary_id',
  'Bio::Annotation::SimpleValue'    => 'value',
  'Bio::Annotation::SimpleValue'    => 'value',
  'Bio::Annotation::OntologyTerm'   => 'name',
  'Bio::Annotation::Reference'      => 'title',
  __DEFAULT__                       => 'value',
);

######################################
#get_SeqFeatures
#display_name
#primary_tag
#source_tag                       x with warning
#has_tag
#get_tag_values
#get_tagset_values
#get_all_tags
#attach_seq
#seq                              x
#entire_seq                       x
#seq_id
#gff_string
#_static_gff_handler
#start                            x
#end                              x
#strand                           x
#location
#primary_id

=head1 PREAMBLE

Okay, where to start...

The original idea for this class appears to lump all SeqFeatureI data
(primary_tag, source_tag, etc) into AnnotationI objects into an
Bio::Annotation::Collection. The type is then checked against SOFA.

There have been several requests to have type checking be optionally run. 

Bio::FeatureHolderI::create_hierarchy_from_ParentIDs
Bio::FeatureHolderI::feature_count
Bio::FeatureHolderI::get_all_SeqFeatures
Bio::FeatureHolderI::set_ParentIDs_from_hierarchy
Bio::RangeI::contains
Bio::RangeI::disconnected_ranges
Bio::RangeI::equals
Bio::RangeI::intersection
Bio::RangeI::offsetStranded
Bio::RangeI::overlap_extent
Bio::RangeI::overlaps
Bio::RangeI::subtract
Bio::RangeI::union
Bio::SeqFeature::Annotated::Dumper
Bio::SeqFeature::Annotated::MAX_TYPE_CACHE_MEMBERS
Bio::SeqFeature::Annotated::add_Annotation
Bio::SeqFeature::Annotated::add_SeqFeature
Bio::SeqFeature::Annotated::add_tag_value
Bio::SeqFeature::Annotated::add_target
Bio::SeqFeature::Annotated::annotation
Bio::SeqFeature::Annotated::attach_seq
Bio::SeqFeature::Annotated::display_name
Bio::SeqFeature::Annotated::each_target
Bio::SeqFeature::Annotated::end
Bio::SeqFeature::Annotated::entire_seq
Bio::SeqFeature::Annotated::frame
Bio::SeqFeature::Annotated::from_feature
Bio::SeqFeature::Annotated::get_Annotations
Bio::SeqFeature::Annotated::get_SeqFeatures
Bio::SeqFeature::Annotated::get_all_tags
Bio::SeqFeature::Annotated::get_tag_values
Bio::SeqFeature::Annotated::get_tagset_values
Bio::SeqFeature::Annotated::has_tag
Bio::SeqFeature::Annotated::length
Bio::SeqFeature::Annotated::location
Bio::SeqFeature::Annotated::name
Bio::SeqFeature::Annotated::new
Bio::SeqFeature::Annotated::phase
Bio::SeqFeature::Annotated::primary_tag
Bio::SeqFeature::Annotated::remove_Annotations
Bio::SeqFeature::Annotated::remove_SeqFeatures
Bio::SeqFeature::Annotated::remove_tag
Bio::SeqFeature::Annotated::score
Bio::SeqFeature::Annotated::seq
Bio::SeqFeature::Annotated::seq_id
Bio::SeqFeature::Annotated::source
Bio::SeqFeature::Annotated::source_tag
Bio::SeqFeature::Annotated::start
Bio::SeqFeature::Annotated::strand
Bio::SeqFeature::Annotated::type
Bio::SeqFeature::Annotated::uri_escape
Bio::SeqFeature::Annotated::uri_unescape
Bio::SeqFeature::TypedSeqFeatureI::croak
Bio::SeqFeature::TypedSeqFeatureI::ontology_term
Bio::SeqFeatureI::generate_unique_persistent_id
Bio::SeqFeatureI::gff_string
Bio::SeqFeatureI::primary_id
Bio::SeqFeatureI::spliced_seq

=cut

sub new {
    my ( $caller, @args) = @_;
    my ($self) = $caller->SUPER::new(@args); 

    $self->_initialize(@args);

    return $self;
}

sub _initialize {
  my ($self,@args) = @_;
  my ($start, $end, $strand, $frame, $phase, $score,
      $name, $annot, $location,
      $display_name, # deprecate
      $seq_id, $type,$source,$feature
     ) =
        $self->_rearrange([qw(START
                              END
                              STRAND
                              FRAME
                              PHASE
                              SCORE
                              NAME
                              ANNOTATION
                              LOCATION
                              DISPLAY_NAME
                              SEQ_ID
                              TYPE
                              SOURCE
			      FEATURE
                             )], @args);
  defined $start        && $self->start($start);
  defined $end          && $self->end($end);
  defined $strand       && $self->strand($strand);
  defined $frame        && $self->frame($frame);
  defined $phase        && $self->phase($phase);
  defined $score        && $self->score($score);
  defined $source       && ref($source) ? $self->source($source) : $self->source_tag($source);
  defined $type         && ref($type) ? $self->type($type) : $self->primary_tag($type);
  defined $location     && $self->location($location);
  defined $annot        && $self->annotation($annot);
  defined $feature      && $self->from_feature($feature);

  if( defined($display_name) && defined($name) ){
	  $self->throw('Cannot define (-id and -seq_id) or (-name and -display_name) attributes');
  }
  defined $seq_id                   && $self->seq_id($seq_id);
  defined ($name || $display_name)  && $self->name($name || $display_name);
}

=head1 ATTRIBUTE ACCESSORS FOR Bio::SeqFeature::Annotated

=cut

=head2 from_feature

  Usage: $obj->from_feature($myfeature);
  Desc : initialize this object with the contents of another feature
         object.  Useful for converting objects like
         L<Bio::SeqFeature::Generic> to this class
  Ret  : nothing meaningful
  Args : a single object of some other feature type,
  Side Effects: throws error on failure
  Example:

=cut

sub from_feature {
    my ($self,$feat,%opts) = @_;
  
    # should deal with any SeqFeatureI implementation (i.e. we don't want to
    # automatically force a OO-heavy implementation on all classes)
    ref($feat) && ($feat->isa('Bio::SeqFeatureI')) 
      or $self->throw('invalid arguments to from_feature');
  
    #TODO: add overrides in opts for these values, so people don't have to screw up their feature object
    #if they don't want to
  
    ### set most of the data
    foreach my $fieldname (qw/ start end strand frame score location seq_id source_tag primary_tag/) {
      #no strict 'refs'; #using symbolic refs, yes, but using them for methods is allowed now
      $self->$fieldname( $feat->$fieldname );
    }

    # now pick up the annotations/tags of the other feature
    # We'll use AnnotationAdaptor to convert everything over

    my %no_copy = map {$_ => 1} qw/seq_id source type frame phase score/;
    my $adaptor = Bio::SeqFeature::AnnotationAdaptor->new(-feature => $feat);
    for my $key ( $adaptor->get_all_annotation_keys() ) {
        next if $no_copy{$key};
        my @values = $adaptor->get_Annotations($key);
        @values = _aggregate_scalar_annotations(\%opts,$key,@values);
        foreach my $val (@values) {
            $self->add_Annotation($key,$val)
        }
    }
}
#given a key and its values, make the values into
#Bio::Annotation::\w+ objects

sub _aggregate_scalar_annotations {
  my ($opts,$key,@values) = @_;

  #anything that's not an object, make it a SimpleValue
  @values = map { ref($_) ? $_ : Bio::Annotation::SimpleValue->new(-value => $_) } @values;

  #try to make Target objects
  if($key eq 'Target' && (@values == 3 || @values == 4)
     && @values == grep {$_->isa('Bio::Annotation::SimpleValue')} @values
    ) {
    @values = map {$_->value} @values;
    #make a strand if it doesn't have one, enforcing start <= end
    if(@values == 3) {
      if($values[1] <= $values[2]) {
	$values[3] = '+';
      } else {
	@values[1,2] = @values[2,1];
	$values[3] = '-';
      }
    }
    return ( Bio::Annotation::Target->new( -target_id => $values[0],
					   -start     => $values[1],
					   -end       => $values[2],
					   -strand    => $values[3],
					 )
	   );
  }
  #try to make DBLink objects
  elsif($key eq 'dblink' || $key eq 'Dbxref') {
    return map {
      if( /:/ ) { #convert to a DBLink if it has a colon in it
	my ($db,$id) = split /:/,$_->value;
	Bio::Annotation::DBLink->new( -database   => $db,
				      -primary_id => $id,
				    );
      } else { #otherwise leave as a SimpleValue
	$_
      }
    } @values;
  }
  #make OntologyTerm objects
  elsif($key eq 'Ontology_term') {
    return map { Bio::Annotation::OntologyTerm->new(-identifier => $_->value) } @values
  }
  #make Comment objects
  elsif($key eq 'comment') {
    return map { Bio::Annotation::Comment->new( -text => $_->value ) } @values;
  }

  return @values;
}


=head2 seq_id()

 Usage   : $obj->seq_id($newval)
 Function: holds a string corresponding to the unique
           seq_id of the sequence underlying the feature
           (e.g. database accession or primary key).
 Returns : string representing the seq_id.
 Args    : on set, some string or a Bio::Annotation::SimpleValue object.

=cut

sub seq_id {
  my($self,$val) = @_;
  if (defined($val)) {
      my $term = undef;
      if (!ref($val)) {
	  $term = Bio::Annotation::SimpleValue->new(-value => uri_unescape($val));
      } elsif (ref($val) && $val->isa('Bio::Annotation::SimpleValue')) {
	  $term = $val;
      }
      if (!defined($term) || ($term->value =~ /^>/)) {
	  $self->throw('give seq_id() a scalar or Bio::Annotation::SimpleValue object, not '.$val);
      }
      $self->remove_Annotations('seq_id');
      $self->add_Annotation('seq_id', $term);
  }

  $self->seq_id('.') unless $self->get_Annotations('seq_id'); # make sure we always have something

  return ($self->get_Annotations('seq_id'))[0]->value;
}

=head2 name()

 Usage   : $obj->name($newval)
 Function: human-readable name for the feature.
 Returns : value of name (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub name {
  my($self,$val) = @_;
  $self->{'name'} = $val if defined($val);
  return $self->{'name'};
}

=head2 type()

 Usage   : $obj->type($newval)
 Function: a SOFA type for the feature.
 Returns : Bio::Annotation::OntologyTerm object representing the type.
           NB: to get a string, use primary_tag().
 Args    : on set, Bio::Annotation::OntologyTerm object.
           NB: to set a string (SOFA name or identifier), use primary_tag()

=cut

use constant MAX_TYPE_CACHE_MEMBERS => 20;
sub type {
  my($self,$val) = @_;
  if(defined($val)){
    my $term = undef;

    if(!ref($val)){
      $self->throw("give type() a Bio::Annotation::OntologyTerm object, not a string");
    }
    elsif(ref($val) && $val->isa('Bio::Annotation::OntologyTerm')){
      $term = $val;
    }
    else {
      #we have the wrong type of object
      $self->throw('give type() a SOFA term name, identifier, or Bio::Annotation::OntologyTerm object, not '.$val);
    }
    $self->remove_Annotations('type');
    $self->add_Annotation('type',$term);
  }
  
  return $self->get_Annotations('type');
}

=head2 source()

 Usage   : $obj->source($newval)
 Function: holds the source of the feature.
 Returns : a Bio::Annotation::SimpleValue representing the source.
           NB: to get a string, use source_tag()
 Args    : on set, a Bio::Annotation::SimpleValue object.
           NB: to set a string, use source_tag()

=cut

sub source {
  my($self,$val) = @_;

  if (defined($val)) {
      my $term;
      if (!ref($val)) {
        $self->throw("give source() a Bio::Annotation::SimpleValue object, not a string");
        #$term = Bio::Annotation::SimpleValue->new(-value => uri_unescape($val));
      } elsif (ref($val) && $val->isa('Bio::Annotation::SimpleValue')) {
	  $term = $val;
      } else {
	  $self->throw('give source() a scalar or Bio::Annotation::SimpleValue object, not '.$val);
      }
      $self->remove_Annotations('source');
      $self->add_Annotation('source', $term);
  }
  
  unless ($self->get_Annotations('source')) {
    $self->source(Bio::Annotation::SimpleValue->new(-value => '.'));
  }
  return $self->get_Annotations('source');
}

=head2 score()

 Usage   : $score = $feat->score()
           $feat->score($score)
 Function: holds a value corresponding to the score of the feature.
 Returns : a string representing the score.
 Args    : on set, a scalar or a Bio::Annotation::SimpleValue object.

=cut

sub score {
  my $self = shift;
  my $val = shift;

  if(defined($val)){
      my $term = undef;
      if (!ref($val)) {
	  $term = Bio::Annotation::SimpleValue->new(-value => $val);
      } elsif (ref($val) && $val->isa('Bio::Annotation::SimpleValue')) {
	  $term = $val;
      }

      if ($term->value ne '.' &&
           (!defined($term) || ($term->value !~ /^[+-]?\d+\.?\d*(e-\d+)?/))) {
	  $self->throw("'$val' is not a valid score");
      }
      $self->remove_Annotations('score');
      $self->add_Annotation('score', $term);
  }

  $self->score('.') unless scalar($self->get_Annotations('score')); # make sure we always have something

  return ($self->get_Annotations('score'))[0]->display_text;
}

=head2 phase()

 Usage   : $phase = $feat->phase()
           $feat->phase($phase)
 Function: get/set on phase information
 Returns : a string 0,1,2,'.'
 Args    : on set, one of 0,1,2,'.' or a Bio::Annotation::SimpleValue
           object holding one of 0,1,2,'.' as its value.

=cut

sub phase {
  my $self = shift;
  my $val = shift;

  if(defined($val)){
      my $term = undef;
      if (!ref($val)) {
	  $term = Bio::Annotation::SimpleValue->new(-value => $val);
      } elsif (ref($val) && $val->isa('Bio::Annotation::SimpleValue')) {
	  $term = $val;
      }
      if (!defined($term) || ($term->value !~ /^[0-2.]$/)) {
	  $self->throw("'$val' is not a valid phase");
      }
      $self->remove_Annotations('phase');
      $self->add_Annotation('phase', $term);
  }

  $self->phase('.') unless $self->get_Annotations('phase'); # make sure we always have something
  
  return ($self->get_Annotations('phase'))[0]->value;
}


=head2 frame()

 Usage   : $frame = $feat->frame()
           $feat->frame($phase)
 Function: get/set on phase information
 Returns : a string 0,1,2,'.'
 Args    : on set, one of 0,1,2,'.' or a Bio::Annotation::SimpleValue
           object holding one of 0,1,2,'.' as its value.

=cut

sub frame {
  my $self = shift;
  my $val = shift;

  if(defined($val)){
      my $term = undef;
      if (!ref($val)) {
	  $term = Bio::Annotation::SimpleValue->new(-value => $val);
      } elsif (ref($val) && $val->isa('Bio::Annotation::SimpleValue')) {
	  $term = $val;
      }
      if (!defined($term) || ($term->value !~ /^[0-2.]$/)) {
	  $self->throw("'$val' is not a valid frame");
      }
      $self->remove_Annotations('frame');
      $self->add_Annotation('frame', $term);
  }

  $self->frame('.') unless $self->get_Annotations('frame'); # make sure we always have something
  
  return ($self->get_Annotations('frame'))[0]->value;
}

############################################################

=head1 SHORTCUT METHODS TO ACCESS Bio::AnnotatableI INTERFACE METHODS

=cut

=head2 add_Annotation()

 Usage   :
 Function: $obj->add_Annotation() is a shortcut to $obj->annotation->add_Annotation
 Returns : 
 Args    :

=cut

sub add_Annotation {
  my ($self,@args) = @_;
  return $self->annotation->add_Annotation(@args);
}

=head2 remove_Annotations()

 Usage   :
 Function: $obj->remove_Annotations() is a shortcut to $obj->annotation->remove_Annotations
 Returns : 
 Args    :

=cut

sub remove_Annotations {
  my ($self,@args) = @_;
  return $self->annotation->remove_Annotations(@args);
}

############################################################

=head1 INTERFACE METHODS FOR Bio::SeqFeatureI

Note that no methods are deprecated.  Any SeqFeatureI methods must return
strings (no objects).

=cut

=head2 display_name()

=cut

sub display_name {
  my $self = shift;
  return $self->name(@_);
}

=head2 primary_tag()

=cut

sub primary_tag {
  my $self = shift;
  if (@_) {
    my $val = shift;
    my $term;
    if(!ref($val) && $val){
      #we have a plain text annotation coming in.  try to map it to SOFA.

      our %__type_cache; #a little cache of plaintext types we've already seen

      #clear our cache if it gets too big
      if(scalar(keys %__type_cache) > MAX_TYPE_CACHE_MEMBERS) {
        %__type_cache = ();
      }

      #set $term to either a cached value, or look up a new one, throwing
      #up if not found
      my $anntext = $val;
      if ($__type_cache{$anntext}) {
        $term = $__type_cache{$anntext};
      } else {
        my $sofa = Bio::Ontology::OntologyStore->get_instance->get_ontology('Sequence Ontology OBO');
        my ($soterm) = $anntext =~ /^\D+:\d+$/ #does it look like an ident?
          ? ($sofa->find_terms(-identifier => $anntext))[0] #yes, lookup by ident
          : ($sofa->find_terms(-name => $anntext))[0];      #no, lookup by name
        #throw if it's not in SOFA
        unless($soterm){
          $self->throw("couldn't find a SOFA term matching type '$val'.");
        }
        my $newterm = Bio::Annotation::OntologyTerm->new;
        $newterm->term($soterm);
        $term = $newterm;
      }
      
      $self->type($term);
    }
  }
  
  my $t = $self->type() || return;
  return $t->name;
}

=head2 source_tag()

=cut

sub source_tag {
  my $self = shift;
  if (@_) {
    my $val = shift;
    if(!ref($val) && $val){
      my $term = Bio::Annotation::SimpleValue->new(-value => uri_unescape($val));
      $self->source($term);
    }
  }
  my $t = $self->source() || return;
  return $t->display_text;
}


=head2 attach_seq()

 Usage   : $sf->attach_seq($seq)
 Function: Attaches a Bio::Seq object to this feature. This
           Bio::Seq object is for the *entire* sequence: ie
           from 1 to 10000
 Returns : TRUE on success
 Args    : a Bio::PrimarySeqI compliant object

=cut

sub attach_seq {
   my ($self, $seq) = @_;

   if ( ! ($seq && ref($seq) && $seq->isa("Bio::PrimarySeqI")) ) {
       $self->throw("Must attach Bio::PrimarySeqI objects to SeqFeatures");
   }

   $self->{'seq'} = $seq;

   # attach to sub features if they want it
   foreach ( $self->get_SeqFeatures() ) {
       $_->attach_seq($seq);
   }
   return 1;
}

=head2 seq()

 Usage   : $tseq = $sf->seq()
 Function: returns a truncated version of seq() with bounds matching this feature
 Returns : sub seq (a Bio::PrimarySeqI compliant object) on attached sequence
           bounded by start & end, or undef if there is no sequence attached
 Args    : none

=cut

sub seq {
  my ($self) = @_;

  return unless defined($self->entire_seq());

  my $seq = $self->entire_seq->trunc($self->start(), $self->end());

  if ( defined $self->strand && $self->strand == -1 ) {
    $seq = $seq->revcom;
  }

  return $seq;
}

=head2 entire_seq()

 Usage   : $whole_seq = $sf->entire_seq()
 Function: gives the entire sequence that this seqfeature is attached to
 Returns : a Bio::PrimarySeqI compliant object, or undef if there is no
           sequence attached
 Args    : none

=cut

sub entire_seq {
  return shift->{'seq'};
}

############################################################

=head1 INTERFACE METHODS FOR Bio::RangeI

 as inherited via Bio::SeqFeatureI

=cut

=head2 length()

 Usage   : $feature->length()
 Function: Get the feature length computed as $feat->end - $feat->start + 1
 Returns : integer
 Args    : none

=cut

sub length {
  my $self = shift;
  return $self->end() - $self->start() + 1;
}

=head2 start()

 Usage   : $obj->start($newval)
 Function: Get/set on the start coordinate of the feature
 Returns : integer
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub start {
  my ($self,$value) = @_;
  return $self->location->start($value);
}

=head2 end()

 Usage   : $obj->end($newval)
 Function: Get/set on the end coordinate of the feature
 Returns : integer
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub end {
  my ($self,$value) = @_;
  return $self->location->end($value);
}

=head2 strand()

 Usage   : $strand = $feat->strand($newval)
 Function: get/set on strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub strand {
  my $self = shift;
  return $self->location->strand(@_);
}


############################################################

=head1 INTERFACE METHODS FOR Bio::FeatureHolderI

This includes methods for retrieving, adding, and removing
features. Since this is already a feature, features held by this
feature holder are essentially sub-features.

=cut

=head2 get_SeqFeatures

 Usage   : @feats = $feat->get_SeqFeatures();
 Function: Returns an array of Bio::SeqFeatureI objects
 Returns : An array
 Args    : none

=cut

sub get_SeqFeatures {
  return @{ shift->{'sub_array'} || []};
}

=head2 add_SeqFeature()

 Usage   : $feat->add_SeqFeature($subfeat);
           $feat->add_SeqFeature($subfeat,'EXPAND')
 Function: adds a SeqFeature into the subSeqFeature array.
           with no 'EXPAND' qualifer, subfeat will be tested
           as to whether it lies inside the parent, and throw
           an exception if not.

           If EXPAND is used, the parent''s start/end/strand will
           be adjusted so that it grows to accommodate the new
           subFeature
 Example :
 Returns : nothing
 Args    : a Bio::SeqFeatureI object

=cut

sub add_SeqFeature {
  my ($self,$val, $expand) = @_;

  return unless $val;

  if ((!ref($val)) || !$val->isa('Bio::SeqFeatureI') ) {
      $self->throw((ref($val) ? ref($val) : $val)
                   ." does not implement Bio::SeqFeatureI.");
  }

  if($expand && ($expand eq 'EXPAND')) {
      $self->_expand_region($val);
  } else {
      if ( !$self->contains($val) ) {
	  $self->warn("$val is not contained within parent feature, and expansion is not valid, ignoring.");
	  return;
      }
  }

  push(@{$self->{'sub_array'}},$val);
}

=head2 remove_SeqFeatures()

 Usage   : $obj->remove_SeqFeatures
 Function: Removes all sub SeqFeatures.  If you want to remove only a subset,
           remove that subset from the returned array, and add back the rest.
 Returns : The array of Bio::SeqFeatureI implementing sub-features that was
           deleted from this feature.
 Args    : none

=cut

sub remove_SeqFeatures {
  my ($self) = @_;

  my @subfeats = @{$self->{'sub_array'} || []};
  $self->{'sub_array'} = []; # zap the array.
  return @subfeats;
}

############################################################

=head1 INTERFACE METHODS FOR Bio::AnnotatableI

=cut

=head2 annotation()

 Usage   : $obj->annotation($annot_obj)
 Function: Get/set the annotation collection object for annotating this
           feature.
 Returns : A Bio::AnnotationCollectionI object
 Args    : newvalue (optional)

=cut

sub annotation {
    my ($obj,$value) = @_;

    # we are smart if someone references the object and there hasn't been
    # one set yet
    if(defined $value || ! defined $obj->{'annotation'} ) {
        $value = Bio::Annotation::Collection->new() unless ( defined $value );
        $obj->{'annotation'} = $value;
    }
    return $obj->{'annotation'};
}

############################################################

=head2 location()

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
    $self->{'location'} = $value;
  }
  elsif (! $self->{'location'}) {
    # guarantees a real location object is returned every time
    $self->{'location'} = Bio::Location::Simple->new();
  }
  return $self->{'location'};
}

=head2 add_target()

 Usage   : $seqfeature->add_target(Bio::LocatableSeq->new(...));
 Function: adds a target location on another reference sequence for this feature
 Returns : true on success
 Args    : a Bio::LocatableSeq object

=cut

sub add_target {
  my ($self,$seq) = @_;
  $self->throw("$seq is not a Bio::LocatableSeq, bailing out") unless ref($seq) and seq->isa('Bio::LocatableSeq');
  push @{ $self->{'targets'} }, $seq;
  return $seq;
}

=head2 each_target()

 Usage   : @targets = $seqfeature->each_target();
 Function: Returns a list of Bio::LocatableSeqs which are the locations of this object.
           To obtain the "primary" location, see L</location()>.
 Returns : a list of 0..N Bio::LocatableSeq objects
 Args    : none

=cut

sub each_target {
  my ($self) = @_;
  return $self->{'targets'} ? @{ $self->{'targets'} } : ();
}

=head2 _expand_region

 Title   : _expand_region
 Usage   : $self->_expand_region($feature);
 Function: Expand the total region covered by this feature to
           accommodate for the given feature.

           May be called whenever any kind of subfeature is added to this
           feature. add_SeqFeature() already does this.
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
#        $self->strand($feat->strand) unless $self->strand();
    } else {
        my $range = $self->union($feat);
        $self->start($range->start);
        $self->end($range->end);
        $self->strand($range->strand);
    }
}

=head2 get_Annotations

 Usage   : my $parent   = $obj->get_Annotations('Parent');
           my @parents = $obj->get_Annotations('Parent');
 Function: a wrapper around Bio::Annotation::Collection::get_Annotations().
 Returns : returns annotations as
           Bio::Annotation::Collection::get_Annotations() does, but
           additionally returns a single scalar in scalar context
           instead of list context so that if an annotation tag
           contains only a single value, you can do:

           $parent = $feature->get_Annotations('Parent');

           instead of:

           ($parent) = ($feature->get_Annotations('Parent'))[0];

           if the 'Parent' tag has multiple values and is called in a
           scalar context, the number of annotations is returned.

 Args    : an annotation tag name.

=cut

sub get_Annotations {
    my $self = shift;

    my @annotations = $self->annotation->get_Annotations(@_);

    if(wantarray){
        return @annotations;
    } elsif(scalar(@annotations) == 1){
        return $annotations[0];
    } else {
        return scalar(@annotations);
    }
}

=head1 Bio::SeqFeatureI implemented methods

These are specialized implementations of SeqFeatureI methods which call the
internal Bio::Annotation::AnnotationCollection object. Just prior to the 1.5
release the below methods were moved from Bio::SeqFeatureI to Bio::AnnotatableI,
and having Bio::SeqFeatureI inherit Bio::AnnotatableI. This behavior forced all
Bio::SeqFeatureI-implementing classes to use Bio::AnnotationI objects for any
data. It is the consensus of the core developers that this be rolled back in
favor of a more flexible approach by rolling back the above changes and making
this class Bio::AnnotatableI. The SeqFeatureI tag-related methods are
reimplemented in order to approximate the same behavior as before.

The methods below allow mapping of the "get_tag_values()"-style annotation
access to Bio::AnnotationCollectionI. These need not be implemented in a
Bio::AnnotationCollectionI compliant class, as they are built on top of the
methods.  For usage, see Bio::SeqFeatureI.

=cut

=head2 has_tag

=cut

sub has_tag {
  my ($self,$tag) = @_;
  return scalar($self->annotation->get_Annotations($tag));
}

=head2 add_tag_value

=cut

sub add_tag_value {
  my ($self,$tag,@vals) = @_;

  foreach my $val (@vals){
    my $class = $tagclass{$tag}   || $tagclass{__DEFAULT__};
    my $slot  = $tag2text{$class};

    my $a = $class->new();
    $a->$slot($val);

    $self->annotation->add_Annotation($tag,$a);
  }

  return 1;
}

=head2 get_tag_values

 Usage   : @annotations = $obj->get_tag_values($tag)
 Function: returns annotations corresponding to $tag
 Returns : a list of scalars
 Args    : tag name

=cut

sub get_tag_values {
    my ($self,$tag) = @_;
    if(!$tagclass{$tag} && $self->annotation->get_Annotations($tag)){
        #new tag, haven't seen it yet but it exists.  add to registry
        my($proto) = $self->annotation->get_Annotations($tag);
        # we can only register if there's a method known for obtaining the value
        if (exists($tag2text{ref($proto)})) {
            $tagclass{$tag} = ref($proto);
        }
    }

    my $slot  = $tag2text{ $tagclass{$tag} || $tagclass{__DEFAULT__} };
    
    return map { $_->$slot } $self->annotation->get_Annotations($tag);
}

=head2 get_tagset_values

 Usage   : @annotations = $obj->get_tagset_values($tag1,$tag2)
 Function: returns annotations corresponding to a list of tags.
           this is a convenience method equivalent to multiple calls
           to get_tag_values with each tag in the list.
 Returns : a list of Bio::AnnotationI objects.
 Args    : a list of tag names

=cut

sub get_tagset_values {
  my ($self,@tags) = @_;
  my @r = ();
  foreach my $tag (@tags){
    my $slot  = $tag2text{ $tagclass{$tag} || $tagclass{__DEFAULT__} };
    push @r, map { $_->$slot } $self->annotation->get_Annotations($tag);
  }
  return @r;
}

=head2 get_all_tags

 Usage   : @tags = $obj->get_all_tags()
 Function: returns a list of annotation tag names.
 Returns : a list of tag names
 Args    : none

=cut

sub get_all_tags {
  my ($self,@args) = @_;
  return $self->annotation->get_all_annotation_keys(@args);
}

=head2 remove_tag

 Usage   : See remove_Annotations().
 Function:
 Returns : 
 Args    : 
 Note    : Contrary to what the name suggests, this method removes
           all annotations corresponding to $tag, not just a
           single anntoation.

=cut

sub remove_tag {
  my ($self,@args) = @_;
  return $self->annotation->remove_Annotations(@args);
}

1;
