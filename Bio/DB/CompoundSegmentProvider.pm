package Bio::DB::CompoundSegmentProvider;

# $Id$
# An implementation of Bio::DB::SegmentProviderI that provides the
# data (the Segments) from another implementation in addition to its
# own.

=head1 NAME

Bio::DB::CompoundSegmentProvider -- An implementation of
Bio::DB::SegmentProviderI that provides the data (the Segments) from
another implementation in addition to its own.

=head1 SYNOPSIS

TODO

=head1 DESCRIPTION

TODO

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
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR

Paul Edlefsen E<lt>paul@systemsbiology.orgE<gt>.
Copyright (c) 2003 Institute for Systems Biology

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=head1 CONTRIBUTORS

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...
use strict;
use vars qw( @ISA );

use Bio::DB::SimpleSegmentProvider;
@ISA = qw( Bio::DB::SimpleSegmentProvider );

use vars '$VERSION';
$VERSION = '1.00';

use Bio::DB::GFF::Util::Rearrange; # for 'rearrange'
use Bio::SeqFeature::CompoundSegment;
use Bio::SeqFeature::CompoundIterator;

sub new {
  my( $class, @args ) = @_;

  my $self = $class->SUPER::new( @args );
  $self->_initialize_compound_segment_provider( @args );
  return $self;
} # new(..)

sub _initialize_compound_segment_provider {
  my $self = shift;
  my @args = @_;

  return if( $self->{ '_compound_segment_provider_initialized' } );

  $self->_initialize_simple_segment_provider( @args );

  $self->{ '_next_providers' } = [];

  $self->{ '_compound_segment_provider_initialized' }++;
  return $self;
} # _initialize_compound_segment_provider(..)

## TODO: Document this.
sub add_next_provider {
  my $self = shift;

  foreach my $next_provider ( @_ ) {
    push( @{ $self->{ '_next_providers' } }, $next_provider );
  }
  return;
} # add_next_provider(..)

## TODO: Document this.
sub get_next_providers {
  my $self = shift;

  return @{ $self->{ '_next_providers' } };
} # get_next_providers()

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 get_collection

 Title    : get_collection
 Usage    : my $segment = $segmentprovider->get_collection( %args );
            OR
            my $segment = $segmentprovider->get_collection( @types );
            OR
            my @segments = $segmentprovider->get_collection( %args );
            OR
            my @segments = $segmentprovider->get_collection( @types );
 Returns  : A L<Bio::SeqFeature::SegmentI> object or a list thereof.
 Args     : see below
 Status   : Public
 Exception: "Those features do not share a common sequence" if this
            method is called in scalar context and the features that
            would otherwise be included in the resulting segment do
            not all fall on the same sequence.

NOTE: This method is (almost) identical to the get_collection() method
from L<Bio::DB::FeatureProviderI> that it overrides.  The entire
documentation follows, but first a brief summary of the changes:
  * This method returns L<Bio::SeqFeature::SegmentI> objects instead
    of mere CollectionI objects.  SegmentI objects are CollectionI
    objects, so this is an additional constraint on the interface.
    The returned SegmentI objects will have as their range the range
    searched, if any, or the smallest range that encloses the returned
    features.
  * This method will return a list of objects if called in list
    context; one L<Bio::SeqFeature::SegmentI> object per root sequence
    of the requested features.  Each returned SegmentI will have as
    its seq_id the common sequences' unique_id() or primary_id().
  * This method will throw an exception if called in scalar context
    and the features that would be included in the resulting SegmentI
    do not all share a common sequence.

This routine will retrieve one or more L<Bio::SeqFeature::SegmentI>
objects based on feature type, location or attributes.  The
SeqFeatureI objects in the returned SegmentIs may or may not be newly
instantiated by this request.  They will have as their range the range
searched, if any, or the smallest range that encloses the returned
features.  They will have as their seq_id() the unique_id() or
primary_id() of the returned features' common sequence.  If this
method is called in list context then one SegmentI object will be
returned per root sequence.  If this method is called in scalar
context and the returned features do not share a common sequence then
an exception will be thrown.

If you make a modification to a feature you must call
update_collection with a collection that contains that feature to
ensure that the data provider is in sync with your change.  You may
not, however, assume that modifications to the feature do not
auto-sync (they might!).

If a range is specified using the -range argument then this range will
 be used to narrow the results, according to the specified -rangetype
 and -strandtype arguments.

-rangetype is one of:
   "overlaps"      return all features that overlap the range (default)
   "contains"      return features completely contained within the range
   "contained_in"  return features that completely contain the range

-strandmatch is one of:
   "strong"        ranges must have the same strand
   "weak"          ranges must have the same strand or no strand (default)
   "ignore"        ignore strand information

Two types of argument lists are accepted.  In the positional argument
form, the arguments are treated as a list of feature types (as if they
were given as -types => \@_).  In the named parameter form, the
arguments are a series of -name=E<gt>value pairs.  Note that the table
below is not exhaustive; implementations must support these but may
support other arguments as well (and are responsible for documenting the
difference).

  Argument       Description
  --------       ------------

  -type          A type name or an object of type L<Bio::SeqFeature::TypeI>
  -types         An array reference to multiple type names or TypeI objects

  -unique_id     A (string) unique_id.  See also -namespace.
  -unique_ids    An array reference to multiple unique_id values.

  -name          A (string) display_name or unique_id.  See also -namespace.
  -names         An array reference to multiple display_name/unique_id values.

  -namespace     A (string) namespace qualifier to help resolve the name/id(s)
  -class         same as -namespace

  -attributes    A hashref containing a set of attributes to match.  See
                 below.

  -baserange     A L<Bio::RangeI> object defining the range to which
                 the -range argument is relative.  There may be a
                 default -baserange.  If this SegmentProviderI is also a
                 L<Bio::RangeI>, then the default -baserange should be
                 itself.

  -range         A L<Bio::RangeI> object defining the range to search.
                 See also -strandmatch, -rangetype, and -baserange.
  -ranges        An array reference to multiple ranges.

  -rangetype     One of "overlaps", "contains", or "contained_in".

  -strandmatch   One of "strong", "weak", or "ignore".  Note that the strand
                 attribute of a given -range must be non-zero for this to work
                 (a 0/undef strand forces a 'weak' strandmatch to become
                 'ignore' and cripples the 'strong' strandmatch).

All plural arguments are interchangeable with their singular counterparts.

The -attributes argument is a hashref containing one or more
attributes to match against:

  -attributes => { Gene => 'abc-1',
                   Note => 'confirmed' }

Attribute matching is simple string matching, and multiple attributes
are ANDed together.

The -unique_ids argument is a reference to a list of strings.  Every
returned feature must have its unique_id value in this list or, if a
feature has no defined unique_id, then its display_name value in the
list if the list is provided.  A -unique_id argument is treated as a
single-element list of unique_ids.

The -names argument is a reference to a list of strings.  Every
returned feature must have its display_name or its unique_id value in this
list if the list is provided.  A -name argument is treated as a
single-element list of names.

If a -namespace is provided then names and ids (both queries and
targets) will be prepended with "$namespace:" as a bonus.  So
if you do features( -names => [ 'foo', 'bar' ], -namespace => 'ns' )
then any feature with the display_name or unique_id 'foo', 'ns:foo',
'bar', or 'ns:bar' will be returned.

=cut

sub get_collection {
  my $self = shift;

  ## TODO: REMOVE
  #warn $self.'::get_collection( '.join( ', ', @_ )." )";

  ## NOTE/TODO: Sorry about this, but I can't figure out how to design
  ## the ISA hierarchy so that this is automatic.  When I'm a
  ## CompoundSegment, my first ISA is CompoundSegmentProvider, and its
  ## first SUPER is SimpleSegmentProvider, skipping the get_collection
  ## stuff from SimpleSegment.
  my @segments;
  if( $self->isa( 'Bio::SeqFeature::SimpleSegment' ) ) {
    ## TODO: There's a serious problem here when there's subclasses of SimpleSegment that are cut out of the loop.
    @segments = $self->Bio::SeqFeature::SimpleSegment::get_collection( @_ );
  } else {
    @segments = $self->SUPER::get_collection( @_ );
  }
  my %seq_ids;
  foreach my $segment ( @segments ) {
    ## If we've managed to get a compound segment, then it already
    ## includes everything in itself..
    if( $segment->isa( 'Bio::SeqFeature::CompoundSegment' ) ) {
      return $segment;
    }
    ## TODO: REMOVE
    #warn "From self, got $segment, on ".$segment->abs_seq_id().".";
    ## TODO: REMOVE.  Testing.
    #unless( $segment->seq_id() ) {
    #  warn "Oh, mang, shoot.  the args were ( ".join( ', ', @_ )." ).  Types were [ ".join( ', ', @{ $_[ 1 ] } )." ].  self is $self, a ".ref( $self ).".";
    #}
    #next unless( $segment->feature_count() );
    unless( $segment->feature_count() ) {
      ## TODO: REMOVE
      #warn "..but there's no features on this one.";
      next;
    }
    $seq_ids{ $segment->abs_seq_id() } = $segment;
  }
  foreach my $next_provider ( $self->get_next_providers() ) {
    @segments = $next_provider->get_collection( @_ );
    ## TODO: REMOVE
    #warn "from $next_provider, got ( " . join( ', ', @segments ) . " )";
    foreach my $segment ( @segments ) {
      ## TODO: REMOVE
      #warn "From $next_provider, got $segment, on ".$segment->abs_seq_id().".";
      if( exists $seq_ids{ $segment->abs_seq_id() } ) {
        ## TODO: REMOVE
        #warn "Yes, we already know about ".$segment->abs_seq_id().".";
        #next unless( $segment->feature_count() );
        unless( $segment->feature_count() ) {
          ## TODO: REMOVE
          #warn "..but there's no features on this one.";
          next;
        }
        # Merge them.
        ## TODO: REMOVE
        #warn "Merging the two segments now...";
        my $existing_segment = $seq_ids{ $segment->abs_seq_id() };
        my $abs_range = $existing_segment->abs_range();
        if( !ref( $abs_range ) ||
            !$abs_range->isa( 'Bio::PrimarySeqI' ) ||
            !$abs_range->length()
          ) {
          $abs_range = $segment->abs_range();
          ## TODO: REMOVE
          #warn "Okay, so the new abs_range is $abs_range, a ".ref( $abs_range )."; we got it from $next_provider.";
        } else {
          ## TODO: REMOVE
          #warn "Okay, so the new abs_range is $abs_range, a ".ref( $abs_range )."; we got it from the existing segment.";
        }
        my $union_strand = $existing_segment->strand();
        if( defined( $union_strand ) ) {
          if( $union_strand != $segment->strand() ) {
            $union_strand = 0;
          }
        } else {
          $union_strand = $segment->strand();
        }
        my $low = $existing_segment->abs_low();
        if( !defined( $low ) or ( $low > $segment->abs_low() ) ) {
          $low = $segment->abs_low();
        }
        my $high = $existing_segment->abs_high();
        if( !defined( $high ) or ( $high < $segment->abs_high() ) ) {
          $high = $segment->abs_high();
        }
        if( $existing_segment->isa( 'Bio::SeqFeature::CompoundSegment' ) ) {
          $existing_segment->seq_id( $abs_range );
          $existing_segment->strand( $union_strand );
          $existing_segment->start( $low );
          $existing_segment->end( $high );
          $existing_segment->add_next_segment( $segment );
          ## TODO: REMOVE
          warn "Using existing compound segment $existing_segment";
          warn $self->stack_trace_dump();
        } else {
          my $compound_segment =
            Bio::SeqFeature::CompoundSegment->new(
              '-seq_id' => $abs_range,
              '-start'  => $low,
              '-end'    => $high,
              '-strand' => $union_strand
            );
          $compound_segment->add_next_segment( $existing_segment );
          $compound_segment->add_next_segment( $segment );
          ## TODO: REMOVE
          #warn "Created compound segment $compound_segment";
          $seq_ids{ $segment->abs_seq_id() } = $compound_segment;
        }
      } else {
        ## TODO: REMOVE
        #warn "got $segment from $next_provider, on ".$segment->abs_seq_id().", using these args: ( ".join( ', ', @_ )." )";
        $seq_ids{ $segment->abs_seq_id() } = $segment;
      }
    }
  }
  @segments = values %seq_ids;
  ## TODO: REMOVE
  #warn "CompoundSegmentProvider::get_collection(..): done.";
  if( wantarray ) {
    return @segments;
  } elsif( scalar( @segments ) > 1 ) {
    ## TODO: REMOVE
    warn "Uh-oh.  Attempting to return multiple segments, but we were called in a scalar context: segments are ( ".join( ', ', @segments )." )";
    $self->throw( "Those features do not share a common sequence" );
  } else { # 1 segment. 
    return $segments[ 0 ];
  }
} # get_collection(..)

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 insert_or_update_collection

 Title   : insert_or_update_collection
 Usage   : $collectionprovider->insert_or_update($collection);
 Function: Attempts to update all the features of a collection.  If
           a feature doesn\'t exist it inserts it automatically.
 Returns : None
 Args    : L<Bio::SeqFeature::CollectionI> object

=cut

#sub insert_or_update_collection {
#  ## TODO
#  shift->throw_not_implemented( @_ );
#}

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 insert_collection

 Title   : insert_collection
 Usage   : $collectionprovider->insert_collection($collection);
 Function: Insert all the features of a collection.  If any features
           already exist throw an exception. 
 Returns : None
 Args    : L<Bio::SeqFeature::CollectionI> object

=cut

#sub insert_collection {
#  ## TODO
#  shift->throw_not_implemented( @_ );
#}

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 update_collection

 Title   : update_collection
 Usage   : $collectionprovider->update_collection($collection);
 Function: Updates all the features of a collection.  If any do not
           already exist throw an exception.
 Returns : Return the updated collection upon success or undef
           upon failure.
 Args    : L<Bio::SeqFeature::CollectionI> object

  If you make a modification to a feature you must call
  update_collection with a collection that contains that feature to
  ensure that the data provider is in sync with your change.  You may
  not, however, assume that modifications to the feature do not
  auto-sync (they might!).

=cut

#sub update_collection {
#  ## TODO
#  shift->throw_not_implemented( @_ );
#}

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 remove_collection

 Title   : remove_collection
 Usage   : $collectionprovider->remove_collection($collection);
 Function: Removes all the features in a collection.  If any features 
           do not exists throw an exception.
 Returns : None
 Args    : L<Bio::SeqFeature::CollectionI> object

=cut

#sub remove_collection {
#  ## TODO
#  shift->throw_not_implemented( @_ );
#}

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 types

 Title   : types
 Usage   : my @types = $collectionprovider->types();
           OR
           my %types_and_counts = $collectionprovider->types( -count => 1 );
 Function: Enumerate the feature types provided by this provider, and possibly
           count the features in each type.
 Returns : a list of L<Bio::SeqFeature::TypeI> objects
           OR
           a hash mapping type id strings to integer counts
 Args    : see below

This routine returns a list of feature types known to the provider.
If the -count argument is given, it returns a hash of known types
mapped to their occurrence counts in this provider.  Note that the
returned list (or the keys of the returned hash) may include types for
which the count is 0.  Also note that the hierarchy of TypeI objects
is ignored, so if there are 4 features of type 'foo' which is a child
of type 'bar', and only 1 feature (explicitly) of type 'bar', then the
count for 'bar' will be 1, not 5.

Arguments are -option=E<gt>value pairs as follows:

  -count aka -enumerate  if true, count the features

The returned value will be a list of L<Bio::SeqFeature::TypeI> objects
or a hash with the string values of these objects as keys.

=cut

sub types {
  my $self = shift;

  my $count;
  if( scalar( @_ ) && $_[ 0 ] =~ /^-/ ) {
    ( $count ) =
      rearrange( [ qw( COUNT COUNTS ENUM ENUMERATE ) ], @_ );
  }
  my @args = @_;
  unless( $count ) {
    push( @args, ( '-count' => 1 ) );
  }
  my %types = $self->SUPER::types( @args );
  foreach my $next_provider ( $self->get_next_providers() ) {
    my %next_types = $next_provider->types( @_ );
    foreach my $type ( keys %next_types ) {
      $types{ $type } += $next_types{ $type };
    }
  }
  if( $count ) {
    return %types;
  } else {
    return keys %types;
  }
} # types(..)

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 seq_ids

 Title   : seq_ids
 Usage   : my @seq_ids = $segmentprovider->seq_ids();
           OR
           my %seq_ids_and_counts =
               $segmentprovider->seq_ids( -count => 1 );
 Function: Enumerate all root seq_ids of features provided by this
           provider, and all seq_ids of sequences provided by this
           provider, and possibly count the features with each seq_id.
 Returns : a list of strings
           OR
           a hash mapping seq_id strings to integer counts
 Args    : see below

This routine returns a list of feature root seq_ids known to the
provider.  If the -count argument is given, it returns a hash of known
seq_ids mapped to their occurrence counts in this provider.  Note that
the returned list (or the keys of the returned hash) may include
seq_ids for which the count is 0, which indicates that the sequence is
provided but there are no features on it.

Arguments are -option=E<gt>value pairs as follows:

  -count aka -enumerate  if true, count the features

=cut

sub seq_ids {
  my $self = shift;

  my $count;
  if( scalar( @_ ) && $_[ 0 ] =~ /^-/ ) {
    ( $count ) =
      rearrange( [ qw( COUNT COUNTS ENUM ENUMERATE ) ], @_ );
  }
  my @args = @_;
  unless( $count ) {
    push( @args, ( '-count' => 1 ) );
  }
  my %seq_ids = $self->SUPER::seq_ids( @args );
  foreach my $next_provider ( $self->get_next_providers() ) {
    my %next_seq_ids = $next_provider->seq_ids( @_ );
    foreach my $type ( keys %next_seq_ids ) {
      $seq_ids{ $type } += $next_seq_ids{ $type };
    }
  }
  if( $count ) {
    return %seq_ids;
  } else {
    return keys %seq_ids;
  }
} # seq_ids(..)

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 get_all_primary_ids

 Title   : get_all_primary_ids
 Usage   : my @primary_ids = $provider->get_all_primary_ids()
 Function: Returns an array of all the primary_ids of the sequence
           objects in this data store. These maybe ids (display style)
           or accession numbers or something else completely different
           - they may be anything so long as each sequence has a
           different one.  Note that although some sequences may have
           undefined primary_ids (bad!), the returned list will not
           include undef.
 Returns : an array of strings
 Args    : none
 Status  : Public

=cut

#sub get_all_primary_ids {
#  ## TODO
#  shift->throw_not_implemented();
#} # get_all_primary_ids(..)

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 unique_ids

 Title     : unique_ids
 Usage     : my @unique_ids = $sequenceprovider->unique_ids();
 Function  : Return a list of the unique_ids of all sequences provided
             by this SequenceProvider.  Note that although some
             sequences may have undefined unique_ids, the returned
             list will not include undef.
 Returns   : an array of strings.
 Args      : none
 Status    : Public

=cut

#sub unique_ids {
#  ## TODO
#  shift->throw_not_implemented();
#} # unique_ids(..)

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 sequence_count

 Title   : sequence_count
 Usage   : $collection->sequence_count()
 Function: Return the number of L<Bio::PrimarySeqI> objects that would
           be returned by a call to sequences() with no arguments.
 Returns : integer representing the number of sequence objects
 Args    : None

  This method is implemented in the interface to return
    scalar( $self->sequences() )
  Because this is not particularly efficient, implementers are
  encouraged to override it, but the result should of course be the
  same.

=cut

#sub sequence_count {
#  ## TODO
#  shift->throw_not_implemented();
#} # sequence_count()

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 add_sequences

 Title   : add_sequences
 Usage   : my @added = $collection->add_sequences( @sequence_list );
 Function: Adds the given sequences to this provider.
 Returns : The sequences added (or their count, in scalar context).
 Args    : An array of L<Bio::PrimarySeqI>s
 Status  : Public

=cut

#sub add_sequences {
#  ## TODO
#  shift->throw_not_implemented();
#} # add_sequences(..)

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 remove_sequences

 Title   : remove_sequences
 Usage   : my @removed = $collection->remove_sequences( @sequence_list )
 Function: Removes the requested sequences from this provider.
 Returns : The removed sequences (or their count, in scalar context)
 Args    : An array of L<Bio::PrimarySeqI>s or their ids (see below)
 Status  : Public

  If any argument is a string, it will be taken as either the
  unique_id or the primary_id or the accession of a sequence to be
  removed.  The return list will contain the actual removed sequence
  if there was a match.

=cut

#sub remove_sequences {
#  ## TODO
#  shift->throw_not_implemented();
#} # remove_sequences()

=head2 sequences

 Title   : sequences
 Usage   : my @seqs = $provider->sequences( @names );
           OR
           my @seqs = $provider->sequences( %args );
 Function: Retrieves a list of L<Bio::PrimarySeqI> objects.
 Returns : a list of L<Bio::PrimarySeqI> objects
           OR
           (when the -iterator option is true)
             a L<Bio::Seq::IteratorI> object
           OR
           (when the -callback option is true) true iff the callbacks
             completed.
 Args    : see below
 Status  : Public
 Throws  : "$some_text does not exist" ($some_text might be anything)
            if a version is given and a matching sequence exists,
            but not of that version.

Two types of argument lists are accepted.  In the positional argument
form, the arguments are treated as a list of sequence names (as if they
were given as -names => \@_).  In the named parameter form, the
arguments are a series of -name=E<gt>value pairs.  Note that the table
below is not exhaustive; implementations must support these but may
support next arguments as well (and are responsible for documenting the
difference).

  Argument       Description
  --------       ------------

  -unique_id     A (string) unique_id.  See also -namespace and -id.
  -unique_ids    An array reference to multiple unique_id values.

  -primary_id    A (string) primary_id.  See also -namespace and -id.
  -primary_ids   An array reference to multiple primary_id values.

  -display_id    A (string) display_id.  See also -namespace and -id.
  -display_ids   An array reference to multiple display_id values.

  -id            A (string) unique_id or primary_id or display_id.
                 See also -namespace.
  -ids           An array reference to multiple id values.

  -accession     A (string) accession.  See also -namespace and -name.
  -accessions    An array reference to multiple accession values.

  -name          A (string) accession or id.  See also -namespace and -id.
  -names         An array reference to multiple accession/id values.

  -namespace     A (string) namespace qualifier to help resolve the names.
  -class         same as -namespace

  -iterator      Return a L<Bio::Seq::IteratorI>
  -stream        same as -iterator

  -callback      A callback to invoke on each sequence

All plural arguments are interchangeable with their singular counterparts.

The -unique_ids argument is a reference to a list of strings.  Every
returned sequence must have its unique_id value in this list if the
list is provided*.  Iff a sequence has no unique_id value then its
primary_id will be used instead.  A -unique_id argument is treated as
a single-element list of unique_ids.

The -ids argument is a reference to a list of strings.  Every returned
sequence must have its unique_id, primary_id, or display_id value in
this list if the list is provided*.  An -id argument is treated as a
single-element list of ids.

The -accessions argument is a reference to a list of strings.  Every
returned sequence must have its accession value in this list if the
list is provided*.  An -accession argument is treated as a
single-element list of names.  If the accession value contains a dot
('.') then the suffix after the final dot will be interpreted as a
version number.  If the accession is available but not in the given
version then an exception (ending in ' does not exist') will be
thrown.  An empty version value is acceptable and means 'latest version',
which is also the default, so -accession => 'foo' and -accession =>
'foo.' mean the same thing.

NOTE: If your accession value contains a dot (unrelated to version
number) then you should postpend the value with a dot to disambiguate:
'accession.that.contains.dots' will be interpreted as the accession
'accession.that.contains' with version 'dots' unless you postpend it
with a dot, like so: -accession => 'accession.that.contains.dots.'

The -names argument is a reference to a list of strings.  Every
returned sequence must have its accession, unique_id, primary_id,
or display_id value in this list if the list is provided*.  A -name
argument is treated as a single-element list of names.

If a -namespace is provided then names and ids (both queries and
targets) will be prepended with "$namespace:" as a bonus.  So if you
do sequences( -names => [ 'foo', 'bar' ], -namespace => 'ns' ) then any
sequence with the accession, unique_id, primary_id, or display_id
'foo', 'ns:foo', 'bar', or 'ns:bar' will be returned.

If -iterator is true, then the method returns an object of type
L<Bio::Seq::IteratorI>.  Each call to next_seq() on this
object returns a L<Bio::PrimarySeqI> object from this provider.

If -callback is passed a code reference, the code reference will be
invoked on each sequence returned.  The code will be passed two
arguments consisting of the current sequence and this SequenceProviderI
object, and must return a true value.  If the code returns a false
value, sequence retrieval will be aborted.

-callback and -iterator are mutually exclusive options.  If -iterator
is defined, then -callback is ignored.

Footnote (*): -unique_ids, -primary_ids, -display_ids, -ids,
-accessions, and -names will be ORed together if they appear together,
so if you do $provider->sequences( -unique_id => 'foo', -accession =>
'bar' ) then you will get all sequences with unique_id 'foo' and
(also) all sequences with accession 'bar'.

NOTE: the following methods all build on top of sequences(), and do not
need to be explicitly implemented.

    get_Seq_by_id()
    get_Seq_by_primary_id()
    get_Seq_by_accession()
    get_Seq_by_version()
    get_Seq_stream()
    get_PrimarySeq_stream()

=cut

sub sequences {
  ## TODO: ERE I AM
  my $self = shift;

  my ( $iterator, $callback );
  if( scalar( @_ ) && $_[ 0 ] =~ /^-/ ) {
    ( $iterator, $callback ) =
      rearrange(
        [ [ qw( ITERATOR STREAM ) ],
          [ qw( CALLBACK CALL_BACK ) ]
        ],
        @_
      );
  }

  ## -iterator and -callback are mutually exclusive.
  if( $iterator && $callback ) {
    $self->throw( "The -iterator and -callback options are mutually exclusive, and yet both have been specified.  Perhaps you could apply your callback method to each element returned by the iterator?" );
  }

  ## TODO: We should make sure that the same sequence is not returned twice.

  if( $callback ) {
    my $callbacks_completed = $self->SUPER::sequences( @_ );
    unless( $callbacks_completed ) {
      return $callbacks_completed;
    }
    foreach my $next_provider ( $self->get_next_providers() ) {
      $next_provider->sequences( @_ );
      unless( $callbacks_completed ) {
        return $callbacks_completed;
      }
    }
  } elsif( $iterator ) {
    my $compound_iterator =
      Bio::SeqFeature::CompoundIterator->new();
    $compound_iterator->add_next_iterator( $self->SUPER::sequences( @_ ) );
    foreach my $next_provider ( $self->get_next_providers() ) {
      $compound_iterator->add_next_iterator( $next_provider->sequences( @_ ) );
    }
    return $compound_iterator;
  } else {
    my @sequences = $self->SUPER::sequences( @_ );
    foreach my $next_provider ( $self->get_next_providers() ) {
      push( @sequences, $next_provider->sequences( @_ ) );
    }
    return @sequences;
  }
} # sequences(..)

=head2 sequence_count

 Title   : sequence_count
 Usage   : $collection->sequence_count()
 Function: Return the number of L<Bio::PrimarySeqI> objects that would
           be returned by a call to sequences() with no arguments.
 Returns : integer representing the number of sequence objects
 Args    : None

=cut

sub sequence_count {
  my $self = shift;

  ## There shouldn't be any arguments, but jic we'll pass them along.
  my $count = $self->SUPER::sequence_count( @_ );
  ## TODO: We should make sure that the same sequence isn't counted twice.
  foreach my $next_provider ( $self->get_next_providers() ) {
    $count += $next_provider->sequence_count( @_ );
  }
  return $count;
} # sequence_count()

1;

__END__
