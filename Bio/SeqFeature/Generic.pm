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

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Generic;
use vars qw( @ISA );
use strict;

use Bio::SeqFeature::SimpleSegment;
use Bio::SeqFeatureI;
use Bio::AnnotatableI;
use Bio::AlternativeLocationHolderI;
@ISA = qw( Bio::SeqFeature::SimpleSegment
           Bio::SeqFeatureI
           Bio::AnnotatableI
           Bio::AlternativeLocationHolderI );

use Bio::DB::GFF::Util::Rearrange; # for 'rearrange'
#use Tie::IxHash;
use Bio::Annotation::Collection;
use Bio::Location::Simple;
use Bio::Tools::GFF;

use overload 
  '==' => 'equals';

sub new {
  my ( $caller, @args) = @_;   
  my ( $self ) = $caller->SUPER::new( @args ); 
  
  $self->{ '_gsf_parse_h' }       = {};
  $self->{ '_gsf_tag_hash' }  = {};
  #tie %{$self->{'_gsf_tag_hash'}}, "Tie::IxHash";
  
  # bulk-set attributes
  $self->set_attributes( @args );
  
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
                  -unique_id      unique id
 	 	  -primary        primary tag
 	 	  -source         source tag
 	 	  -frame          frame
 	 	  -score          score value
 	 	  -tag            a reference to a tag/value hash
                  -type           a TypeI or a string representing the type
 	 	  -gff_string     GFF v.2 string to initialize from
 	 	  -gff1_string    GFF v.1 string to initialize from
 	 	  -seq_id         the display name of the sequence
 	 	  -annotation     the AnnotationCollectionI object
 	 	  -location       the LocationI object

=cut

sub set_attributes {
    my ($self,@args) = @_;
    my ($start, $end, $strand, $unique_id, $display_name, $primary_tag, $source_tag, $primary, $source, $frame, 
	$score, $tag, $type, $gff_string, $gff1_string,
	$seqid, $annot, $location, $orientation_policy) =
	    rearrange( [ [ qw( START BEGIN LOW ) ],
                         [ qw( STOP END HIGH ) ],
                         qw( STRAND ),
                         [ qw( UNIQUE_ID PRIMARY_ID ) ],
                         [ qw( DISPLAY_NAME DISPLAY_ID NAME ) ],
                         qw( PRIMARY_TAG ),
			 qw( SOURCE_TAG ),
                         qw( PRIMARY ),
                         qw( SOURCE ),
                         [ qw( FRAME PHASE ) ],
                         qw( SCORE ),
                         [ qw( TAG ATTRIBUTES ) ],
                         qw( TYPE ),
                         qw( GFF_STRING ),
                         qw( GFF1_STRING ),
                         [ qw( SEQ_ID SEQ SEQID ID SEQNAME SEQ_NAME REF ) ],
                         qw( ANNOTATION ),
                         qw( LOCATION ),
                         [ qw( ORIENTATION_POLICY ORIENTATIONPOLICY
                               ORIENTATION POLICY ) ]
		       ], @args );
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
    $type           && $self->source_tag($type);
    defined $unique_id && $self->unique_id( $unique_id );
    defined $strand && $self->strand($strand);
    defined $start && $self->start( $start );
    defined $end && $self->end( $end );
    defined $orientation_policy &&
      $self->orientation_policy( $orientation_policy );
    defined $display_name && $self->display_name( $display_name );
    defined $frame  && $self->frame($frame);
    $score          && $self->score($score);
    $annot          && $self->annotation($annot);
    $seqid          && $self->seq_id($seqid);
    $tag            && do {
	foreach my $t ( keys %$tag ) {
	    $self->add_tag_value($t,$tag->{$t});
	}
    };

    $self->ensure_orientation();
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

=head2 new_from_feature

 Title   : new_from_feature
 Usage   : my $new_feature =
             Bio::SeqFeature::Generic->new_from_feature( $copy_from );
 Function: Create a new Bio::SeqFeature::Generic object by copying
           values from another Generic object.
 Returns : A new L<Bio::SeqFeature::Generic> object
 Args    : Another L<Bio::SeqFeature::Generic> object
 Status  : Protected

  This is a special copy constructor.  It forces the new feature into
  the L<Bio::SeqFeature::Generic> package, regardless of the
  package that it is called from.  This causes subclass-specfic
  information to be dropped.

  This also does not copy into the new feature the features held in
  the existing feature.  If you would like the new feature to hold the
  same features you must explicitly add them, like so:
    $new_feature->add_features( $copy_from->features() );

  As a special bonus you may also pass an existing hash and it will be the
  blessed an anointed object that is returned, like so:
    $new_feature =
      Bio::SeqFeature::Generic->new_from_feature(
        $copy_from,
        $new_feature
      );

=cut

sub new_from_feature {
  my $pack = shift; # ignored
  my $feature = shift || $pack;
  my $new_feature = shift;
  $new_feature =
    Bio::SeqFeature::SimpleSegment->new_from_segment(
      $feature,
      $new_feature
    );
  foreach my $key ( grep { /^_gsf_/ } keys %$feature ) {
    next if( $key eq '_gsf_unique_id' );
    next if( $key eq '_gsf_tag_hash' );
    $new_feature->{ $key } = $feature->{ $key };
  }
  foreach my $key ( keys %{ $feature->{ '_gsf_tag_hash' } } ) {
    next unless $feature->{ '_gsf_tag_hash' }->{ $key };
    foreach my $value ( @{ $feature->{ '_gsf_tag_hash' }->{ $key } } ) {
      push( @{ $new_feature->{ '_gsf_tag_hash' }->{ $key } }, $value );
    }
  }
  return bless $new_feature, __PACKAGE__;
} # new_from_feature(..)

=head2 location

 Title   : location
 Usage   : my $location = $seqfeature->location()
 Function: returns a location object suitable for identifying location
	   of feature on sequence or parent feature
 Returns : Bio::LocationI object
 Args    : [optional] Bio::LocationI object to set the value to.


=cut

sub location {
  my $self = shift;
  my $new_location = shift;

  # TODO: The problem with this approach is that if the user changes
  # the returned location, the location of this feature doesn't
  # change (or maybe this isn't a problem, I dunno. -Paul
  my $old_location = $self->{ '_gsf_location' } ||
    Bio::Location::Simple->new(
      '-seq_id' => $self->seq_id(),
      '-start' => $self->start(),
      '-end' => $self->end(),
      '-strand' => $self->strand()
    );
  if( defined( $new_location ) ) {
    unless( ref( $new_location ) and $new_location->isa( 'Bio::LocationI' ) ) {
      $self->throw( "object $new_location pretends to be a location but ".
                    "does not implement Bio::LocationI" );
    }
    $self->{ '_gsf_location' } = $new_location;
  }
  return $old_location;
} # location(..)

# Internal overridable getter/setter for the actual stored value of
# seq_id.  Delegates to the location object at $self->{ '_gsf_location' }
# if there is one; otherwise delegates to the superclass _seq_id(..).
sub _seq_id {
  my $self = shift;
  if( $self->{ '_gsf_location' } ) {
    my ( $new_val ) = @_;
    my $old_val = $self->{ '_gsf_location' }->seq_id();
    if( defined $new_val ) {
      $self->{ '_gsf_location' }->seq_id( $new_val );
    }
    return $old_val;
  } else {
    return $self->SUPER::_seq_id( @_ );
  }
} # _seq_id(..)

# Internal overridable getter/setter for the actual stored value of
# start.  Delegates to the location object at $self->{ '_gsf_location' }
# if there is one; otherwise delegates to the superclass _start(..).
sub _start {
  my $self = shift;
  if( $self->{ '_gsf_location' } ) {
    my ( $new_val ) = @_;

    my $old_val = $self->{ '_gsf_location' }->start();
    if( defined $new_val ) {
      $self->{ '_gsf_location' }->start( $new_val );
    }
    return $old_val;
  } else {
    return $self->SUPER::_start( @_ );
  }
} # _start(..)

# Internal overridable getter/setter for the actual stored value of
# end.  Delegates to the location object at $self->{ '_gsf_location' }
# if there is one; otherwise delegates to the superclass _end(..).
sub _end {
  my $self = shift;
  if( $self->{ '_gsf_location' } ) {
    my ( $new_val ) = @_;
    my $old_val = $self->{ '_gsf_location' }->end();
    if( defined $new_val ) {
      $self->{ '_gsf_location' }->end( $new_val );
    }
    return $old_val;
  } else {
    return $self->SUPER::_end( @_ );
  }
} # _end(..)

# Internal overridable getter/setter for the actual stored value of
# strand.  Delegates to the location object at $self->{ '_gsf_location' }
# if there is one; otherwise delegates to the superclass _strand(..).
sub _strand {
  my $self = shift;
  if( $self->{ '_gsf_location' } ) {
    my ( $new_val ) = @_;
    my $old_val = $self->{ '_gsf_location' }->strand();
    if( defined $new_val ) {
      $self->{ '_gsf_location' }->strand( $new_val );
    }
    return $old_val;
  } else {
    return $self->SUPER::_strand( @_ );
  }
} # _strand(..)

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

## Just use the one inherited from SeqFeatureI.
#=head2 primary_tag
#
# Title   : primary_tag
# Usage   : $tag = $feat->primary_tag()
#           $feat->primary_tag('exon')
# Function: get/set on the primary tag for a feature,
#           eg 'exon'
# Returns : a string
# Args    : none if get, the new value if set
#
#=cut
#
#sub primary_tag {
#   my ($self,$value) = @_;
#   my $type = $self->type;
#   if ( defined $value ) {
#       if (!$type) {
#           require "Bio/Ontology/TermFactory.pm";
#           my $factory = Bio::Ontology::TermFactory->new(-type => 'Bio::Ontology::Term');
#           $type = 
#             $factory->create_object(-name => $value,
#                                     -category => 'adhoc_sequence');
#           $self->type($type);
#       }
#       else {
#           $type->name($value);
#       }
#   }
#   if( defined $type ) {
#     if( ref( $type ) && $type->isa( 'Bio::Ontology::Term' ) ) {
#       return $type->name();
#     } else {
#       return $type;
#     }
#   } else {
#     return;
#   }
#}
#
=head2 type_string

 Title   : type_string
 Usage   : $tag = $feat->type_string()
           $feat->type_string('exon')
 Function: synonym for 'primary_tag'
           get/set on the primary tag for a feature,
           eg 'exon'
 Returns : a string
 Args    : none if get, the new value if set

=cut

sub type_string {
  shift->primary_tag(@_);
}

=head2 type

 Title   : type
 Usage   : $seq_feature->type( [$new_type] )
 Function: Getter/setter for the type of this feature.
 Returns : the current (or former, if used as a set method) type (either a
           string or a Bio::SeqFeature::TypeI)
 Args    : (optional) A new type (either a string or, preferably, a
           Bio::SeqFeature::TypeI)
 Status  : Public

=cut

sub type {
  my $self = shift;
  my ( $new_value ) = @_;

  my $old_value = $self->{ '_gsf_type' };
  if( defined( $new_value ) ) {
    if( $new_value eq 'undef' ) {
      undef $self->{ '_gsf_type' };
    } else {
      $self->{ '_gsf_type' } = $new_value;
    }
  }
  return $old_value;
} # type(..)

=head2 resolve_type

 Title   : resolve_type
 Usage   :
 Function: resolves a $feat->type() (a Bio::OntologyTermI object)
           using the current type_string() [primary_tag()]
 Example :
 Returns : Bio::OntologyTermI on success, undef on failure
 Args    : Graph (currently must be a GO::Model::Graph object)
           Hash mapping [OPTIONAL] - a hashref og mapping of ad-hoc termnames to 
                                     ontology termnames

=cut

sub resolve_type {
  my ( $self, $graph, $mapping ) = @_;

  my $tname = $self->type_string;
  $tname = $mapping->{ $tname } if $mapping->{ $tname };
  my $term = $graph->get_term_by_name( $tname );
  if( $term ) {
    # turn GO object into Bio object
    $self->type_string( $term->name() );
    if( ref( $self->type() ) &&
        $self->type()->isa( 'Bio::Ontology::TermI' ) ) {
      $self->type()->identifier( $term->acc );
      $self->type()->category( 'sequence' );
      $term = $self->type();
    }
  }
  return $term;
} # resolve_type(..)


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
  my ( $self, $value ) = @_;

  if( defined $value ) {
    $self->{'_gsf_source_tag'} = $value;
  }
  return $self->{'_gsf_source_tag'};
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

   unless( ( defined $self->seq_id() ) && ref( $self->seq_id() ) ) {
     $self->seq_id( $seq );
   }

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

  if ( !exists $self->{'_gsf_seq'} ) {
    return undef;
  }

  # assumming our seq object is sensible, it should not have to yank
  # the entire sequence out here.
  my $seq = $self->{'_gsf_seq'}->trunc( $self->abs_start( 'plus' ), $self->abs_end( 'plus' ) );

  if ( $self->abs_strand() < 0 ) {
    # ok. this does not work well (?)
    ## TODO: REMOVE
    #print STDERR "Before revcom: ".$seq->seq()."\n";
    $seq = $seq->revcom();
    ## TODO: REMOVE
    #print STDERR "After  revcom: ".$seq->seq()."\n";
  }

  ## TODO: REMOVE
  #print STDERR "For $self: ".$self->stack_trace_dump();

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
 Usage   : my $unique_id = $feature->unique_id( [$newval] )
 Function: This is a unique identifier that identifies this object.
           If not set, will return undef per L<Bio::LocallyIdentifiableI.pm>
           If a value is given, the unique_id will be set to it, unless that
           value is the string 'undef', in which case the unique_id will
           become undefined.

           This attribute should *not* be used in GFF dumping, as
           that should come from the collection in which the seq
           feature was found or from seq_id().
 Returns : The current (or former, if used as a set method) value of unique_id
 Args    : newvalue (optional)

=cut

# Note that this method notifies observers.
sub unique_id {
  my ($self,$value) = @_;
  my $current_value = $self->{'_gsf_unique_id'};
  if ( defined $value ) {
    if( !$value || ( $value eq 'undef' ) ) {
      $self->{'_gsf_unique_id'} = undef;
    } else {
      $self->{'_gsf_unique_id'} = $value;
    }
    $self->notify_observers( 'unique_id', '-old'=>$current_value, '-new'=>$value );
  }
  return $current_value;
} # unique_id()

=head2 display_name

 Title   : display_name
 Usage   : $featname = $feature->display_name
 Function: Implements the display_name() method, which is a human-readable
           name for the feature. 
 Returns : value of display_name (a string)
 Args    : Optionally, on set the new value or undef 

=cut

sub display_name {
    my $self = shift;

    return $self->{'_gsf_display_name'} = shift if @_;
    return $self->{'_gsf_display_name'};
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
        $value = Bio::Annotation::Collection->new() unless ( defined $value );
        $obj->{'annotation'} = $value;
    }
    return $obj->{'annotation'};
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
	    $self->{'_gsf_gffio'} = $gffio;
	} else {
	    $Bio::SeqFeatureI::static_gff_formatter = $gffio;
	}
    }
    return (ref($self) && exists($self->{'_gsf_gffio'}) ?
	    $self->{'_gsf_gffio'} : $self->_static_gff_formatter);
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

sub gff_string {
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

   return $self->{'_gsf_parse_h'};
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

sub _expand_region {
  shift->adjust_bounds( 0, @_ );
} # _expand_region(..)

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


# Aliases because we can never quite settle on the right method names.
sub sub_SeqFeature {
  shift->Bio::SeqFeature::SimpleSegment::get_SeqFeatures( @_ );
}
sub all_sub_SeqFeature_types {
  shift->Bio::SeqFeature::SimpleSegment::types( @_ );
}
sub add_SeqFeature {
  shift->add_SeqFeatures( @_ );
}
sub add_sub_SeqFeature {
  shift->add_SeqFeatures( @_ );
}
sub flush_sub_SeqFeatures {
  my $self = shift;
  my @features = @_ || $self->features();
  $self->remove_SeqFeatures( @features );
}
sub flush_sub_SeqFeature {
  shift->flush_sub_SeqFeatures( @_ );
}

=head2 equals

  Title   : equals
  Usage   : if( $seq_feature->equals( $another_seq_feature ) )
  Function: Return true iff the other feature is 'equal' to this one
            (see below); if the argument is a RangeI but not a
            SeqFeatureI, then return true iff it has the same start,
            end, length as this SeqFeatureI does; if it is any other
            object, return false.
  Args    : a Bio::SeqFeatureI
  Returns : true iff they are describing the same feature

   This method is implemented in the interface to use the following
   equality test for other SeqFeatureIs:

   * If either SeqFeatureI returns something other than undef from its
     unique_id() method, then return false unless they are the same
     (by eq).
   * Return false unless they have the same range (using the RelRangeI
     superclass equals method).
   * If they have the same range and either returns something other
     than undef from its display_name() method, then return false
     unless they are the same (by eq).
   * If all other tests do not fail (ie. they have the same range and
     no unique_id or display_name), return true.

  Implementing classes may override or augment this behavior.

=cut

sub equals {
  my ( $self, $other, @args ) = @_;

  return unless( defined $other );

  if( ref( $other ) && $other->isa( 'Bio::SeqFeatureI' ) ) {
    my $my_unique_id = $self->unique_id();
    my $other_unique_id = $other->unique_id();
    if( defined( $my_unique_id ) ) {
      if( defined( $other_unique_id ) ) {
        ## TODO: REMOVE
        #if( $my_unique_id eq $other_unique_id ) {
        #  print STDERR "$self and $other ARE equal because their unique_ids are the same ($my_unique_id)\n";
        #} else {
        #  print STDERR "$self and $other are not equal because their unique_ids are different.\n";
        #}
        return ( $my_unique_id eq $other_unique_id );
      }
      ## TODO: REMOVE
      #print STDERR "$self and $other are not equal because $self has a unique_id but $other does not.\n";
      return 0;
    }
    if( defined( $other_unique_id ) ) {
      ## TODO: REMOVE
      #print STDERR "$self and $other are not equal because $other has a unique_id but $self does not.\n";
      return 0;
    }
    ## Okay, neither has a unique_id.
    unless( $self->Bio::RelRangeI::equals( $other ) ) {
      ## TODO: REMOVE
      #print STDERR "$self and $other are not equal because RelRangeI says so.\n";
      return 0;
    }
    my $my_display_name = $self->display_name();
    my $other_display_name = $other->display_name();
    if( defined( $my_display_name ) ) {
      if( defined( $other_display_name ) ) {
        ## TODO: REMOVE
        #if( $my_display_name eq $other_display_name ) {
        #  print STDERR "$self and $other ARE equal because their display_names are the same ($my_display_name)\n";
        #} else {
        #  print STDERR "$self and $other are not equal because their display_names are different.\n";
        #}
        return ( $my_display_name eq $other_display_name );
      }
      ## TODO: REMOVE
      #print STDERR "$self and $other are not equal because $self has a display_name but $other does not.\n";
      return 0;
    }
    if( defined( $other_display_name ) ) {
      ## TODO: REMOVE
      #print STDERR "$self and $other are not equal because $other has a display_name but $self does not.\n";
      return 0;
    }
    ## Okay, neither has a display_name
    ## TODO: REMOVE
    #print STDERR "$self and $other ARE equal because neither has a display_name.\n";
    return 1;
  } else {
    ## TODO: REMOVE
    #print STDERR "$self and $other are not equal because $other isn't even a feature.\n";
    return 0;
  }
} # equals(..)

1;

__END__
