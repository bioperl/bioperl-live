# 

=head1 NAME

Bio::DB::GFF::Segment -- Simple DNA segment object

=head1 SYNOPSIS

See L<Bio::DB::GFF>.

=head1 DESCRIPTION

Bio::DB::GFF::Segment provides the basic representation of a range of
DNA contained in a GFF database.  It is the base class from which the
Bio::DB::GFF::Segment and Bio::DB::GFF::Feature classes are
derived.

Generally, you will not create or manipulate Bio::DB::GFF::Segment
objects directly, but use those that are returned by the Bio::DB::GFF
module.

=cut

package Bio::DB::GFF::Segment;

use strict;
use Bio::SeqFeature::SimpleSegment;
use Bio::SeqI;

use vars qw($VERSION @ISA);
@ISA = qw( Bio::SeqFeature::SimpleSegment Bio::SeqI );
$VERSION = '0.31';

## Other imports
use Bio::Annotation::Collection;
use Bio::DB::GFF::Util::Rearrange; # for &rearrange
use Bio::RelRange qw( &absRange &absSeqId &absStart &absEnd &absStrand );

use overload 
  '""'     => 'asString',
  eq       => 'equals',
  fallback => 1;

=head1 API

The remainder of this document describes the API for
Bio::DB::GFF::Segment.

=cut

=head2 new

 Title   : new
 Usage   : $s = Bio::DB::GFF::Segment->new(@args)
 Function: create a new segment
 Returns : a new Bio::DB::GFF::Segment object
 Args    : see below
 Status  : Public

This method creates a new Bio::DB::GFF::Segment object.  Generally
this is called automatically by the Bio::DB::GFF module and
derivatives.

Arguments may either be in named argument form or positional.
  There are five positional arguments:
  
   $factory      a Bio::DB::GFF::Adaptor to use for database access
   $sourceclass  class of the source sequence
   $sourceseq    ID of the source sequence
   $start        start of the desired segment relative to source sequence
   $stop         stop of the desired segment relative to source sequence

  Named arguments are:
   -factory    => factory and DBI interface
   -seq        => seq_id, the sequence or segment that this one is relative to
   -class      => The namespace of the -seq sequence
   -start      => start relative to seq_id
   -stop       => stop relative to seq_id
   -offset     => 0-based offset relative to seq_id (instead of -start, -stop)
   -length     => length of segment (instead of -stop; required with -offset)
   -ref        => new_seq_id, range to re-relativize to after creation
   -refclass   => The namespace of new_seq_id
   -absolute   => use absolute coordinate addressing
   -nocheck    => turn off checking, force segment to be constructed

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
                 the -range argument is relative.  The default
                 baserange depends on the value of the absolute()
                 flag.  If absolute() is true then the default is the
                 value of the abs_seq_id() method.  If absolute() is
                 false then the default is $self.  Note that the
                 baserange affects the sort order.  See also
                 -absolute.

  -absolute      If -absolute is given and true then all behavior will be as
                 if this SegmentI's absolute() flag was set to true,
                 even if it isn't.  If -absolute is given and false
                 then all behavior will be as if this SegmentI's
                 absolute() flag was set to false, even if it isn't.
                 Note that -baserange can still be given and can force
                 relativeness, and that takes precedence over -absolute.

  -range         A L<Bio::RangeI> object defining the range to search.
                 See also -strandmatch, -rangetype, and -baserange.
  -ranges        An array reference to multiple ranges.

  -rangetype     One of "overlaps", "contains", or "contained_in".

  -strandmatch   One of "strong", "weak", or "ignore".  Note that the strand
                 attribute of a given -range must be non-zero for this to work
                 (a 0/undef strand forces a 'weak' strandmatch to become
                 'ignore' and cripples the 'strong' strandmatch).

=cut
 
sub new {
  my $pack = shift;
  $pack = ref $pack if ref $pack;

  my ($parent,$seq_id,$start,$stop,$new_seq_id,$class,$new_class,$offset,$length,$force_absolute,$nocheck,$closure_options);
  if( scalar( @_ ) && $_[ 0 ] =~ /^-/ ) {
    ($parent,$seq_id,$start,$stop,$new_seq_id,$class,$new_class,$offset,$length,$force_absolute,$nocheck,$closure_options) =
      rearrange([
                  [qw(FACTORY PARENT)],
                  [qw(NAME SEQ SEQUENCE SOURCESEQ SEQ_ID SEQID)],
                  [qw(START BEGIN)],
                  [qw(STOP END)],
                  [qw(REFSEQ REF REFNAME NEW_SEQ_ID NEWSEQID)],
                  [qw(CLASS SEQCLASS SEQ_ID_CLASS SEQIDCLASS)],
                  [qw(REFCLASS NEW_SEQ_ID_CLASS NEWSEQIDCLASS)],
                  [qw(OFFSET OFF)],
                  [qw(LENGTH LEN)],
                  [qw(ABSOLUTE)],
                  [qw(NOCHECK FORCE)],
                  [qw(CLOSURE OPTIONS CLOSURE_OPTIONS CLOSURE_OPTIONS)]
                ], @_ );
  } else {
    ( $parent, $class, $seq_id, $start, $stop ) = @_;
    $force_absolute = 1;
  }

  unless( defined( $parent ) &&
          ref( $parent ) &&
          $parent->isa( 'Bio::DB::GFF' )
        ) {
    if( defined( $parent ) &&
        ref( $parent ) &&
        $parent->isa( 'Bio::DB::GFF::Segment' ) ) {
      $parent = $parent->factory();
    } else {
      $pack->throw( "new(): you must provide a factory argument that isa Bio::DB::GFF object.".( defined( $parent ) ? ( ref( $parent ) ? '  The object given is a '.ref( $parent ). '.' : "  The value given, '$parent', is not an object." ) : '  No factory was given.' ) );
    }
  }

  # support for Featname objects
  if( ref( $seq_id ) && $seq_id->isa( 'Bio::DB::GFF::Featname' ) ) {
    $class = $seq_id->class();
    $seq_id = $seq_id->name();
  }
  # if the class of the landmark is not specified then default to 'Sequence'
  $class ||= 'Sequence';

  # We allow people to use segments as sources.  This is the same as
  # if the user had done $seq_id->subseq( $start, $stop ) instead,
  # so we just do that for her.
  if( ref( $seq_id ) && $seq_id->isa( 'Bio::DB::GFF::Segment' ) ) {
    unless( defined $start ) {
      $start = 1;
    }
    unless( defined $stop ) {
      $stop = $seq_id->length();
    }
    #return $seq_id->subseq( $start, $stop );
    my $result = $seq_id->subseq( $start, $stop );
    unless( $result->isa( 'Bio::DB::GFF::Segment' ) ) {
      $pack->throw( 'Internal error: The segment created by \$seq_id->subseq( $start, $stop ) is not a Bio::DB::GFF::Segment!  It is a '.ref( $result ).'.  \$seq_id is $seq_id, a '.ref( $seq_id ).'.' );
    }
  }

  if( $nocheck ) {
    $force_absolute++;
    $start = 1;
  }

  # Clean up range info
  my $strand = 1;
  if( defined $offset ) {
    if( defined $start ) {
      warn "new(): bad idea to call new() with both a start and an offset";
    }
    $start = $offset + 1;
  }
  # An explicit length overrides start and stop
  if( defined $length ) {
    if( defined $stop ) {
      warn "new(): bad idea to call new() with both a stop and a length";
    }
    $stop = ( $start + $length - 1 );
  }
  if( $stop < $start ) {
    # Flip and note it in the strand.
    my $temp = $start;
    $start = $stop;
    $stop = $temp;
    $strand = -1;
  }

  # It is possible that $seq_id is a range, in which
  # case all of our start, end, etc. are relative to that.  This next
  # section uses the database to get the absolute coords of the seq_id
  # if it is not a range; in
  # the case that seq_id refers to a sequence, this will just be
  # 1..sequence_length.  It is even possible that the database will
  # return multiple possible ranges/sequences that the given
  # seq_id/class
  # combination could refer to.  In this case we'll return a list of
  # multiple new Segment objects, one for each possibility.
  my @seq_ids;
  if( ref( $seq_id ) && $seq_id->isa( 'Bio::RangeI' ) ) {
    push( @seq_ids, $seq_id );
  } else {
    # @abscoords is an array ref, each element of which is
    # ($abs_ref,$abs_class,$abs_start,$abs_stop,$abs_strand,$sname).
    # That's the format returned by GFF's abscoords(..) method.
    my @abscoords;
  
    if( $force_absolute && defined( $start ) ) {
      # absolute position is given to us
      @abscoords = ( [ $seq_id, $class, $start, $stop, '+' ] );
    } else {
      my $result =
        $parent->abscoords(
          $seq_id,
          $class,
          ( $force_absolute ? $seq_id : () )
        );
      return unless $result; ## TODO: Complain?  This seems like a big deal 2 me.
      @abscoords = @$result;
    }
  
    foreach ( @abscoords ) {
      # This is a possible range for $seq_id, so eg. $seq_id_seq_id is
      # one possibility for $seq_id's $seq_id.
      my ( $seq_id_seq_id, $seq_id_class, $seq_id_start,
           $seq_id_stop, $seq_id_strand, $seq_id_sname ) =
        @$_;
      ## TODO: REMOVE?
      unless( defined( $seq_id_seq_id ) ) {
        warn "Broken result from the abscoords(..) method: ( ".
          join( ", ", @$_ )." )";
        next;
      }

      unless( defined( $seq_id_sname ) ) {
        $seq_id_sname = $seq_id; # Yes, $seq_id, not $seq_id_seq_id.
      }
      $seq_id_strand ||= '+';

      ## Segment_NamedRelRange is an inner class, defined below.
      my $a_possible_seq_id = Bio::DB::GFF::Segment_NamedRelRange->new(
        '-name' => $seq_id_sname,
        '-seq_id' => $seq_id_seq_id,
        '-start' => $seq_id_start,
        '-end' => $seq_id_stop,
        '-strand' => $seq_id_strand
      );

      # We're gonna store some extra info in there.  This is kinda
      # hacky, but hey, this is Perl.  We'll precede all additions
      # with _gff_, to avoid name clashes.

      # This is how we'll hackily keep track of the class of it.
      $a_possible_seq_id->{ '_gff_class' } = $seq_id_class;
      
      # Apparently "this allows an SQL optimization way down deep" --
      # when we do end up making a $self, we'll store this temporary
      # '_gff_whole' value in $self as $self->{ 'whole' }.
      if( ( $seq_id_seq_id eq $seq_id_sname ) and
          !defined( $start ) and
          !defined( $stop )
        ) {
        $a_possible_seq_id->{ '_gff_whole' }++;
      }
      push( @seq_ids, $a_possible_seq_id );
    } # End foreach @abscoords
  } # End building up the list of possible @seq_ids. (if..else..)

  # So we get to return one $self for each possible $seq_id.
  my @object_results;
  foreach my $a_possible_seq_id ( @seq_ids ) {
    # Localize, man.
    my $a_possible_start = $start;
    my $a_possible_stop = $stop;

    # If we haven't been given start/stop info we make it the
    # entirety of the seq.
    unless( defined $a_possible_start ) {
      $a_possible_start = 1;
    }
    unless( defined $a_possible_stop ) {
      $a_possible_stop = $a_possible_seq_id->length();
    }

    # When $force_absolute is given, the coords provided are
    # interpreted as absolute coords. (I think).
    # TODO: rerelativize the coords.

    # also TODO: This is what it used to be, and it just doesn't make
    # any sense to me.
#    if( $force_absolute ) {
#      ( $a_possible_start, $a_possible_stop ) =
#        ( $seq_id->abs_start(), $seq_id->abs_end() );
#      $a_possible_self->absolute( 1 );
#    }

    # Make sure that our possible seq_id is rooted at a real sequence.
    my $sequence_range = $a_possible_seq_id;
    if( ref( $a_possible_seq_id ) &&
        $a_possible_seq_id->isa( 'Bio::DB::GFF::Segment_NamedRelRange' ) &&
        ( $a_possible_seq_id->name() ne $a_possible_seq_id->seq_id() )
      ) {
      my $a =
        $parent->abscoords( $a_possible_seq_id->seq_id(), 'Sequence' );
      my $seq_start = $a->[ 0 ][ 2 ];
      my $seq_stop  = $a->[ 0 ][ 3 ];
      ## Paul's question: why would $seq_start ever not be 1?  How
      ## should that affect absolute positions?
      $a_possible_seq_id->seq_id(
        Bio::RelRange->new(
          '-seq_id' => $a->[ 0 ][ 0 ],
          '-start' => $a->[ 0 ][ 2 ],
          '-end' => $a->[ 0 ][ 3 ],
          '-strand' => 1,
          '-orientation_policy' => 'dependent'
        )
      );
      $sequence_range = $a_possible_seq_id->seq_id();
    } # End if we need to get the root sequence's bounds

    my $a_possible_self = $pack->SUPER::new(
      '-seq_id' => $a_possible_seq_id,
      '-start' => $a_possible_start,
      '-end' => $a_possible_stop,
      '-strand' => $strand,
      '-parent' => $parent,
      '-orientation_policy' => 'dependent'
    );

    $a_possible_self->{ 'class' } = $class;

    # Handle truncation in either direction.  This only happens if the
    # segment runs off the end of the reference sequence.
    if( $parent->strict_bounds_checking() ) {
      my %truncated;
      # Skip this possible seq_id if we are completely off the end of
      # the sequence, 'cause that's, like, impossible.
      if( ( $a_possible_self->abs_start() > $sequence_range->end() ) ||
          ( $a_possible_self->abs_end() < $sequence_range->start() ) ) {
        next;
      }
      if( $a_possible_self->abs_start() < $sequence_range->start() ) {
        $a_possible_self->start(
          $a_possible_self->abs2rel( $sequence_range->start() )
        );
        $truncated{ 'start' }++;
      }
      if( $a_possible_self->abs_end() > $sequence_range->end() ) {
        $a_possible_self->end(
          $a_possible_self->abs2rel( $sequence_range->end() )
        );
        $truncated{ 'stop' }++;
      }
      $a_possible_self->{ 'truncated' } = %truncated;
    }

    # We (hackily) stored the 'whole' value in the possible seq_id.
    $a_possible_self->{ 'whole' } = $a_possible_seq_id->{ '_gff_whole' };
    undef $a_possible_seq_id->{ '_gff_whole' };

    # If they want to rerelativize, do it.
    if( defined( $new_class ) ) {
      if( defined( $new_seq_id ) &&
          ref( $new_seq_id ) &&
          $new_seq_id->isa( 'Bio::DB::GFF:Featname' )
        ) {
        $new_seq_id->class( $new_class );
      } elsif( defined( $new_seq_id ) ) {
        $new_seq_id =
          Bio::DB::GFF::Featname->new( $new_class, $new_seq_id );
      } else {
        $a_possible_seq_id->{ '_gff_class' } = $new_class;
      }
    }
    if( defined( $new_seq_id ) ) {
      $a_possible_self->seq_id( $new_seq_id );
    }

    if( $force_absolute ) {
      $a_possible_self->absolute( 1 );
    }

    ## TODO: REMOVE.  Testing.
    if( $closure_options ) {
      my @co = %{ $closure_options };
      ## TODO: REMOVE
      #warn "Segment::new(): Using these closure options: { ".join( ', ', @co )." }";
      my ( $co_types, $co_unique_ids, $co_namespace, $co_names, $co_attributes, $co_baserange, $co_ranges, $co_strandmatch, $co_rangetype, $co_sort, $co_sparse, $co_merge, $more_closure_options ) =
        rearrange(
          [ [ qw( TYPE TYPES ) ],
            [ qw( UNIQUE_IDS UNIQUEIDS IDS UNIQUE_ID UNIQUEID ID ) ],
            [ qw( NAMESPACE NAME_SPACE CLASS ) ],
            [ qw( NAME NAMES DISPLAY_NAME DISPLAY_NAMES DISPLAYNAME DISPLAYNAMES ) ],
            [ qw( ATTRIBUTE ATTRIBUTES ) ],
            [ qw( BASERANGE BASE_RANGE BASELOCATION BASE_LOCATION BASE ) ],
            [ qw( RANGE RANGES LOCATION LOCATIONS LOC ) ],
            [ qw( STRANDMATCH STRAND_MATCH ) ],
            [ qw( RANGETYPE RANGE_TYPE ) ],
            [ qw( SORT SORTED ) ],
            [ qw( RARE SPARSE ) ],
            [ qw( MERGE AUTOMERGE ) ]
          ],
          @co
        );
  
      $closure_options = $more_closure_options;

      ## Fix up types.
      if( $co_types ) {
        unless( ref $co_types eq 'ARRAY' ) {
          ## The incoming value might be a simple scalar instead of a list.
          $co_types = [ $co_types ];
        }
        ## Just in case it's an array ref to an empty string:
        if(
           !scalar( @$co_types ) ||
           ( ( scalar( @$co_types ) == 1 ) && !( $co_types->[ 0 ] ) )
          ) {
          undef $co_types;
        }
      }
      if( $co_types ) {
        $closure_options->{ 'types' } = $co_types;
      }
      
      ## Fix up unique_ids.
      if( $co_unique_ids ) {
        unless( ref $co_unique_ids eq 'ARRAY' ) {
          ## The incoming value might be a simple scalar instead of a list.
          $co_unique_ids = [ $co_unique_ids ];
        }
        ## Just in case it's an array ref to an empty string:
        if(
           !scalar( @$co_unique_ids ) ||
           ( ( scalar( @$co_unique_ids ) == 1 ) && !( $co_unique_ids->[ 0 ] ) )
          ) {
          undef $co_unique_ids;
        }
      }
      if( $co_unique_ids ) {
        $closure_options->{ 'unique_ids' } = $co_unique_ids;
      }
      
      ## Fix up names.
      if( $co_names ) {
        unless( ref $co_names eq 'ARRAY' ) {
          ## The incoming value might be a simple scalar instead of a list.
          $co_names = [ $co_names ];
        }
        ## Just in case it's an array ref to an empty string:
        if(
           !scalar( @$co_names ) ||
           ( ( scalar( @$co_names ) == 1 ) && !( $co_names->[ 0 ] ) )
          ) {
          undef $co_names;
        }
      }
      if( $co_names ) {
        $closure_options->{ 'names' } = $co_names;
      }
      
      ## Attributes better be a hash ref if it is anything.
      if( $co_attributes ) {
        unless( ref( $co_attributes ) eq 'HASH' ) {
          $a_possible_self->throw(
            "The -attributes argument must be a HASH REF."
          );
        }
        $closure_options->{ 'attributes' } = $co_attributes;
      }
      
      ## Fix up ranges.
      if( $co_ranges ) {
        unless( ref $co_ranges eq 'ARRAY' ) {
          ## The incoming value might be a simple scalar instead of a list.
          $co_ranges = [ $co_ranges ];
        }
        ## Just in case it's an array ref to an empty string:
        if(
           !scalar( @$co_ranges ) ||
           ( ( scalar( @$co_ranges ) == 1 ) && !( $co_ranges->[ 0 ] ) )
          ) {
          undef $co_ranges;
        }
      }
      
      ## Derelativize, man.
      if( $co_ranges && @$co_ranges && $co_baserange ) {
        my @new_ranges;
        foreach my $co_range ( @$co_ranges ) {
          unless( ref( $co_range ) && $co_range->isa( 'Bio::RangeI' ) ) {
            $a_possible_self->throw( "Expected the -ranges argument to be a reference to a list of Bio::RangeI objects, but it contains something incompatible: " . overload::StrVal( $co_range ) );
          }
          unless( defined $co_range->seq_id() ) {
            $co_range = Bio::RelRange->new(
                       -seq_id => $co_baserange,
                       -start =>  $co_range->start(),
                       -end =>    $co_range->end(),
                       -strand => $co_range->strand(),
                       -orientation_policy => 'dependent'
                     );
          } elsif( !$co_range->isa( 'Bio::RelRangeI' ) ) {
            $co_range = Bio::RelRange->new(
                       -seq_id => $co_range->seq_id(),
                       -start =>  $co_range->start(),
                       -end =>    $co_range->end(),
                       -strand => $co_range->strand(),
                       -orientation_policy => 'dependent'
                     );
          }
          $co_range->absolute( 1 );
          push( @new_ranges, $co_range );
        }
        $co_ranges = \@new_ranges;
      } # End derelativizing ranges.  Now they're absolute.
      if( $co_ranges ) {
        $closure_options->{ 'ranges' } = $co_ranges;
      }

      if( $co_namespace ) {
        $closure_options->{ 'namespace' } = $co_namespace;
      }
      if( $co_baserange ) {
        $closure_options->{ 'baserange' } = $co_baserange;
      }
      if( $co_strandmatch ) {
        $closure_options->{ 'strandmatch' } = $co_strandmatch;
      }
      if( $co_rangetype ) {
        $closure_options->{ 'rangetype' } = $co_rangetype;
      }
      if( $co_sort ) {
        $closure_options->{ 'sort' } = $co_sort;
      }
      if( $co_sparse ) {
        $closure_options->{ 'sparse' } = $co_sparse;
      }
      if( $co_merge ) {
        $closure_options->{ 'merge' } = $co_merge;
      }

      #warn "!+++++++ Adding closure options.";
      $a_possible_self->{ '_closure_options' } = $closure_options;
    }
    push( @object_results, $a_possible_self );
  } # End creating all of the possible $selfs.

  # If they only want one, we just give 'em one.
  return ( wantarray ? @object_results : $object_results[ 0 ] );
} # new(..)

=head2 new_from_segment

 Title   : new_from_segment
 Usage   : $s = Bio::DB::GFF::Segment->new_from_segment( $copy_from )
 Function: create a new L<Bio::DB::GFF::Segment>
 Returns : A new L<Bio::DB::GFF::Segment> object
 Args    : Another L<Bio::DB::GFF::Segment> object
 Status  : Protected

  This constructor is used internally by the subseq() method.  It forces
  the new segment into the L<Bio::DB::GFF::Segment> package, regardless
  of the package that it is called from.  This causes subclass-specific
  information, such as feature types, to be dropped when a subsequence
  is created.

  This also does not copy into the new segment the features held in
  the existing segment.  If you would like the new segment to hold the
  same features you must explicitly add them, like so:
    $new_segment->add_features( $copy_from->features() );

  As a special bonus you may also pass an existing hash and it will be the
  blessed an anointed object that is returned, like so:
    $new_segment =
      Bio::DB::GFF::Segment->new_from_segment(
        $copy_from,
        $new_segment
      );

=cut

sub new_from_segment {
  my $pack = shift; # ignored
  my $copy_from = shift || $pack;
  my $new_segment = shift;
  $new_segment =
    Bio::SeqFeature::SimpleSegment->new_from_segment(
      $copy_from,
      $new_segment
    );
  @{ $new_segment }{ qw( _class ) } = @{ $copy_from }{ qw( _class ) };
  bless $new_segment, __PACKAGE__;
  $new_segment->orientation_policy( 'dependent' );
  $new_segment->ensure_orientation();
  return $new_segment;
} # new_from_segment

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 get_collection

 Title   : get_collection
 Usage   : my $segment = $segment->get_collection( %args );
           OR
           my $segment = $segment->get_collection( @types );
 Returns : A L<Bio::SeqFeature::SegmentI> object
 Args    : see below
 Status  : Public

This routine will retrieve a L<Bio::SeqFeature::SegmentI> object based
on feature type, location or attributes.  The SeqFeatureI objects in
the returned SegmentI may or may not be newly instantiated by this
request.  They will have as their range the range searched, if any, or
the smallest range that encloses the returned features.  They will
have as their seq_id() the -baserange used here (if the baserange is
absolute by any means then their seq_id() will be this SegmentI's
abs_seq_id()).

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
                   (default ONLY when -strand is specified and non-zero)
   "weak"          ranges must have the same strand or no strand
   "ignore"        ignore strand information
                   (default unless -strand is specified and non-zero)

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
                 the -range argument is relative.  The default
                 baserange depends on the value of the absolute()
                 flag.  If absolute() is true then the default is the
                 value of the abs_seq_id() method.  If absolute() is
                 false then the default is $self.  Note that the
                 baserange affects the sort order.  See also
                 -absolute.

  -absolute      If -absolute is given and true then all behavior will be as
                 if this SegmentI's absolute() flag was set to true,
                 even if it isn't.  If -absolute is given and false
                 then all behavior will be as if this SegmentI's
                 absolute() flag was set to false, even if it isn't.
                 Note that -baserange can still be given and can force
                 relativeness, and that takes precedence over -absolute.

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

  ## TODO: REMOVE.  Testing.
  #return $self->factory()->segment( @_ );

  ## NEW: Testing..

  ## TODO: Add offset stuff...  See GFF.pm's segment() method.
  my ( $types, $unique_ids, $namespace, $names, $attributes, $baserange, $ranges, $strandmatch, $rangetype, $sort, $sparse, $merge, $seq_id, $start, $end, $refclass, $absolute, $force );
  if( scalar( @_ ) && $_[ 0 ] =~ /^-/ ) {
    ( $types, $unique_ids, $namespace, $names, $attributes, $baserange, $ranges, $strandmatch, $rangetype, $sort, $sparse, $merge, $seq_id, $start, $end, $refclass, $absolute, $force ) =
      rearrange(
        [ [ qw( TYPE TYPES ) ],
          [ qw( UNIQUE_IDS UNIQUEIDS IDS UNIQUE_ID UNIQUEID ID ) ],
          [ qw( NAMESPACE NAME_SPACE CLASS SEQCLASS SEQ_CLASS ) ],
          [ qw( NAME NAMES DISPLAY_NAME DISPLAY_NAMES DISPLAYNAME DISPLAYNAMES ) ],
          [ qw( ATTRIBUTE ATTRIBUTES ) ],
          [ qw( BASERANGE BASE_RANGE BASELOCATION BASE_LOCATION BASE ) ],
          [ qw( RANGE RANGES LOCATION LOCATIONS LOC ) ],
          [ qw( STRANDMATCH STRAND_MATCH ) ],
          [ qw( RANGETYPE RANGE_TYPE ) ],
          [ qw( SORT SORTED ) ],
          [ qw( RARE SPARSE ) ],
          [ qw( MERGE AUTOMERGE ) ],
          [ qw( SEQ_ID SEQID SEQ ID REF REF_SEQ REFSEQ ) ],
          [ qw( START BEGIN LOW ) ],
          [ qw( STOP END HIGH ) ],
          [ qw( REFCLASS REF_CLASS ) ],
          [ qw( ABSOLUTE ABS ) ],
          [ qw( FORCE NOCHECK NO_CHECK ) ]
        ],
        @_
      );
  } else {
    ## Types.
    $types = \@_;
  }

  my $closure_options = $self->{ '_closure_options' } || {};

  unless( defined $namespace ) {
    $namespace = $closure_options->{ 'namespace' };
  }

  ## Fix up types.
  if( $types ) {
    unless( ref $types eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $types = [ $types ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$types ) ||
       ( ( scalar( @$types ) == 1 ) && !( $types->[ 0 ] ) )
      ) {
      undef $types;
    }
  }

  if( $closure_options->{ 'types' } ) {
    if( $types ) {
      unshift( @$types, @$closure_options->{ 'types' } );
    } else {
      $types = $closure_options->{ 'types' };
    }
  }

  ## Fix up unique_ids.
  if( $unique_ids ) {
    unless( ref $unique_ids eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $unique_ids = [ $unique_ids ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$unique_ids ) ||
       ( ( scalar( @$unique_ids ) == 1 ) && !( $unique_ids->[ 0 ] ) )
      ) {
      undef $unique_ids;
    }
  }

  if( $closure_options->{ 'unique_ids' } ) {
    if( $unique_ids ) {
      unshift( @$unique_ids, @$closure_options->{ 'unique_ids' } );
    } else {
      $unique_ids = $closure_options->{ 'unique_ids' };
    }
  }

  ## Now we need to remove from $unique_ids the id of this Segment, if
  ## it is in there.
  if( $self->unique_id() && $unique_ids ) {
    my $unique_id = $self->unique_id();

    # Remove any entries that are the same as the unique_id, perhaps
    # taking into account the namespace..
    @$unique_ids =
      grep { ( $_ ne $unique_id ) && ( $namespace ? ( ( ( $namespace.$_ ) ne $unique_id ) && ( $_ ne ( $namespace.$unique_id ) ) ) : 1 ) } @$unique_ids;
  }

  ## Fix up names.
  if( $names ) {
    unless( ref $names eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $names = [ $names ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$names ) ||
       ( ( scalar( @$names ) == 1 ) && !( $names->[ 0 ] ) )
      ) {
      undef $names;
    }
  }

  if( $closure_options->{ 'names' } ) {
    if( $names ) {
      unshift( @$names, @$closure_options->{ 'names' } );
    } else {
      $names = $closure_options->{ 'names' };
    }
  }

  ## Alas, now we need to remove from $names the name of this Segment,
  ## if that name is in there.
  my $unique_id = $self->unique_id();
  my $display_name = eval{ $self->display_name() };
  my $display_id = $self->display_id();
  my $primary_id = $self->primary_id();
  my $accession_number = $self->accession_number();
  if( $names && ( $unique_id || $display_name || $display_id || $primary_id || $accession_number ) ) {

    # Remove any entries that are the same as the name of this segment,
    # perhaps taking into account the namespace..
    @$names =
      grep { ( $_ ne $unique_id ) && ( $namespace ? ( ( ( $namespace.$_ ) ne $unique_id ) && ( $_ ne ( $namespace.$unique_id ) ) ) : 1 ) && ( $_ ne $display_name ) && ( $namespace ? ( ( ( $namespace.$_ ) ne $display_name ) && ( $_ ne ( $namespace.$display_name ) ) ) : 1 ) && ( $_ ne $display_id ) && ( $namespace ? ( ( ( $namespace.$_ ) ne $display_id ) && ( $_ ne ( $namespace.$display_id ) ) ) : 1 ) && ( $_ ne $primary_id ) && ( $namespace ? ( ( ( $namespace.$_ ) ne $primary_id ) && ( $_ ne ( $namespace.$primary_id ) ) ) : 1 ) && ( $_ ne $accession_number ) && ( $namespace ? ( ( ( $namespace.$_ ) ne $accession_number ) && ( $_ ne ( $namespace.$accession_number ) ) ) : 1 ) } @$names;
  }

  ## Attributes better be a hash ref if it is anything.
  if( $attributes ) {
    unless( ref( $attributes ) eq 'HASH' ) {
      $self->throw( "The -attributes argument must be a HASH REF." );
    }
  }

  if( $closure_options->{ 'attributes' } ) {
    if( $attributes ) {
      unshift( @$attributes, @$closure_options->{ 'attributes' } );
    } else {
      $attributes = $closure_options->{ 'attributes' };
    }
  }

  ## Fix up ranges.
  if( $ranges ) {
    unless( ref $ranges eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $ranges = [ $ranges ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$ranges ) ||
       ( ( scalar( @$ranges ) == 1 ) && !( $ranges->[ 0 ] ) )
      ) {
      undef $ranges;
    }
  } elsif( defined( $seq_id ) || defined( $start ) || defined( $end ) ) {
    $ranges = [ Bio::RelRange->new(
                  '-seq_id' => $seq_id,
                  '-start' => $start,
                  '-end' => $end,
                   -orientation_policy => 'dependent'
                ) ];
  }

  if( $closure_options->{ 'ranges' } ) {
    if( $ranges ) {
      unshift( @$ranges, @$closure_options->{ 'ranges' } );
    } else {
      $ranges = $closure_options->{ 'ranges' };
    }
  }

  if( !$baserange && $closure_options->{ 'baserange' } ) {
    $baserange = $closure_options->{ 'baserange' };
  }

  ## We can use ourselves as the baserange, if we are a Bio::RangeI.
  if( $self->isa( 'Bio::RangeI' ) && !$baserange ) {
    $baserange = $self;
  }

  ## Derelativize, man.
  if( $ranges && @$ranges && $baserange ) {
    my @new_ranges;
    foreach my $range ( @$ranges ) {
      unless( ref( $range ) && $range->isa( 'Bio::RangeI' ) ) {
        $self->throw( "Expected the -ranges argument to be a reference to a list of Bio::RangeI objects, but it contains something incompatible: " . overload::StrVal( $range ) );
      }
      unless( defined $range->seq_id() ) {
        $range = Bio::RelRange->new(
                   -seq_id => $baserange,
                   -start =>  $range->start(),
                   -end =>    $range->end(),
                   -strand => $range->strand(),
                   -orientation_policy => 'dependent'
                 );
      } elsif( !$range->isa( 'Bio::RelRangeI' ) ) {
        $range = Bio::RelRange->new(
                   -seq_id => $range->seq_id(),
                   -start =>  $range->start(),
                   -end =>    $range->end(),
                   -strand => $range->strand(),
                   -orientation_policy => 'dependent'
                 );
      }
      $range->absolute( 1 );
      push( @new_ranges, $range );
    }
    $ranges = \@new_ranges;
  } # End derelativizing ranges.  Now they're absolute.

  ## -strandmatch defaults to 'weak'
  unless( defined $strandmatch ) {
    $strandmatch = $closure_options->{ 'strandmatch'} || 'weak';
  }

  ## -rangetype defaults to 'overlaps'
  unless( defined $rangetype ) {
    $rangetype = $closure_options->{ 'rangetype'} || 'overlaps';
  }

  ## -sort is redundant if sorted() is true.
  unless( defined( $sort ) ) {
    $sort = $closure_options->{ 'sort' };
  }
  if( $sort && $self->sorted() ) {
    undef $sort;
  }

  unless( defined $sparse ) {
    $sparse = $closure_options->{ 'sparse' };
  }
  unless( defined $merge ) {
    $merge = $closure_options->{ 'merge' };
  }
  unless( defined $absolute ) {
    $absolute = $closure_options->{ 'absolute' };
  }
  unless( defined $force ) {
    $force = $closure_options->{ 'force' };
  }

  ## TODO: REMOVE.  Testing.
  ## TODO: Fix this up a bit.
  my %args = (
    '-seq' => ( $ranges && @$ranges ? $ranges->[ 0 ]->seq_id() : undef ) ||
      $self->abs_seq_id(),
    '-refseq' => ( $ranges && @$ranges ? $ranges->[ 0 ]->seq_id() : undef ) ||
      $self->abs_seq_id(),
    '-class' => $namespace || $self->class(),
    '-start' => ( $ranges && @$ranges ? $ranges->[ 0 ]->start() : undef ) ||
      ( $self->{ 'whole' } ? undef : $self->abs_start() ),
    '-end' => ( $ranges && @$ranges ? $ranges->[ 0 ]->end() : undef ) ||
      ( $self->{ 'whole' } ? undef : $self->abs_end() ),
    '-closure_options' =>
      {
       '-types' => $types,
       '-unique_ids' => $unique_ids,
       '-names' => $names,
       '-namespace' => $namespace,
       '-attributes' => $attributes,
       '-baserange' => $baserange,
       '-absolute'  => $absolute,
       '-ranges'    => $ranges,
       '-range_type' => $rangetype,
       '-strand_match' => $strandmatch,
       '-sparse' => $sparse,
       '-merge' => $merge,
       '-force' => $force
      },
    '-sparse' => $sparse,
    '-merge'  => $merge,
    '-parent' => $self,
    '-absolute' => $absolute,
    '-force' => $force
  );
  ## TODO: REMOVE
  #warn "Args to segment() of ".$self->factory(). " are ( ".join( ', ', ( my @args = %args ) )." ).  Closure options are { ".join( ', ', ( my @co = %{ $args{ '-closure_options' } } ) )." }.  names are { ".($names?join( ', ', @$names ):'')." }.  types are { ".($types?join( ', ', @$types ):''). " }.";
  #eval { $self->throw( "Args to segment() of ".$self->factory(). " are ( ".join( ', ', ( my @args = %args ) )." )"."-names are ".($names&&@$names?'[ '.join( ', ', @$names ).' ]':'undef')."." ) };
  #warn $@ if $@;

  ## TODO: REMOVE.  Testing.
  #return $self->factory()->segment( %args );
  if( wantarray ) {
    my @r = $self->factory()->segment( %args );
    ## TODO: REMOVE
    #warn "Got ( ".join( ', ', @r )." )\n";
    return @r;
  } else {
    my $r = $self->factory()->segment( %args );
    ## TODO: REMOVE
    #warn "Got $r\n";
    return $r;
  }

  # END NEW. Testing...

  # OLD:

  my ( $types, $absolute, $baserange, $ranges );
  my %args;
  if( scalar( @_ ) && $_[ 0 ] =~ /^-/ ) {
    ( $types, $absolute, $baserange, $ranges ) =
      rearrange(
        [ [ qw( TYPE TYPES ) ],
          [ qw( ABSOLUTE ABS ) ],
          [ qw( BASERANGE BASE_RANGE BASELOCATION BASE_LOCATION BASE ) ],
          [ qw( RANGE RANGES LOCATION LOCATIONS LOC ) ]
        ],
        @_
      );
    %args = @_;
  } elsif( scalar( @_ ) ) {
    ## Types.
    $types = \@_;
    %args = ( '-types' => $types );
  }

  # If -baserange is given but is not a RangeI then the
  # -baserange argument to the superclass should be
  # $self->abs_seq_id().  If -baserange is not given then the argument
  # to the superclass method should be $self unless -absolute is given
  # and true or absolute() is true.
  if( not defined( $baserange ) ) {
    if( $absolute || $self->absolute() ) {
      $baserange = $self->abs_range();
      $absolute = 1;
    } else {
      $baserange = $self;
    }
  } elsif( defined( $baserange ) &&
           ( not ref( $baserange ) || not $baserange->isa( 'Bio::RangeI' ) )
         ) {
    $baserange = $self->abs_range();
  }
  $args{ '-baserange' } = $baserange;
  # Oh yeah and make sure it's the only baserange argument..
  ## TODO: What if they're uppercase?
  delete @args{ qw( -base_range -baselocation -base_location -base
                    baserange base_range baselocation base_location base ) };

  ## TODO: Dehackify
  @_ = %args;

  my $seq_id = $self->seq_id();

  # If the user named the sequence itself in their request, give it to 'em.
  # This entails removing from the request any args that name
  # the sequence.
  ## TODO: Dehackify
  my ( $unique_ids, $namespace, $names );
  if( scalar( @_ ) && $_[ 0 ] =~ /^-/ ) {
    ( $unique_ids, $namespace, $names ) =
      rearrange(
        [ [ qw( UNIQUE_IDS UNIQUEIDS IDS UNIQUE_ID UNIQUEID ID ) ],
          [ qw( NAMESPACE NAME_SPACE CLASS ) ],
          [ qw( NAME NAMES DISPLAY_NAME DISPLAY_NAMES DISPLAYNAME DISPLAYNAMES ) ],
        ],
        @_
      );
  }
  ## Fix up unique_ids.
  if( $unique_ids ) {
    unless( ref $unique_ids eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $unique_ids = [ $unique_ids ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$unique_ids ) ||
       ( ( scalar( @$unique_ids ) == 1 ) && !( $unique_ids->[ 0 ] ) )
      ) {
      undef $unique_ids;
    }
  }
  if( $unique_ids ) {
    my @to_be_removed;
    for( my $i = 0; $i < scalar( @$unique_ids ); $i++ ) {
      if( ( $unique_ids->[ $i ] eq $seq_id ) ||
          ( $namespace.':'.$unique_ids->[ $i ] eq $seq_id ) ||
          ( $unique_ids->[ $i ] eq $namespace.':'.$seq_id ) ) {
        ## TODO: REMOVE
        #warn "\$unique_ids->[ $i ] is $seq_id" if Bio::Graphics::Browser::DEBUG;
        push( @to_be_removed, $i );
      }
    }
    my $offset = 0;
    foreach my $i ( @to_be_removed ) {
      splice( @$unique_ids, $i + $offset++, 1 );
    }
    # If we removed any unique_ids, modify @_.
    if( $offset ) {
      if( scalar( @$unique_ids ) ) {
        $args{ '-unique_ids' } = $unique_ids;
      } else {
        delete $args{ '-unique_ids' };
      }
      # Oh yeah and make sure it's the only unique_ids argument..
      ## TODO: What if they're uppercase?
      [ qw( UNIQUE_IDS UNIQUEIDS IDS UNIQUE_ID UNIQUEID ID ) ],
      delete @args{ qw( -unique_id -ids -id -uniqueids -uniqueid
                        uniqueids unique_id ids id uniqueids uniqueid ) };
      ## TODO: REMOVE
      #warn "\$unique_ids is now [ ".join( ', ', @$unique_ids )." ]" if Bio::Graphics::Browser::DEBUG;
    }
  }
  ## Fix up names.
  if( $names ) {
    unless( ref $names eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $names = [ $names ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$names ) ||
       ( ( scalar( @$names ) == 1 ) && !( $names->[ 0 ] ) )
      ) {
      undef $names;
    }
  }
  if( $names ) {
    my @to_be_removed;
    for( my $i = 0; $i < scalar( @$names ); $i++ ) {
      if( ( $names->[ $i ] eq $seq_id ) ||
          ( $namespace.':'.$names->[ $i ] eq $seq_id ) ||
          ( $names->[ $i ] eq $namespace.':'.$seq_id ) ) {
        ## TODO: REMOVE
        #warn "\$names->[ $i ] is $seq_id" if Bio::Graphics::Browser::DEBUG;
        push( @to_be_removed, $i );
      }
    }
    my $offset = 0;
    foreach my $i ( @to_be_removed ) {
      splice( @$names, $i + $offset++, 1 );
    }
    # If we removed any names, modify @_.
    if( $offset ) {
      if( scalar( @$names ) ) {
        $args{ '-names' } = $names;
      } else {
        delete $args{ '-names' };
      }
      # Oh yeah and make sure it's the only names argument..
      ## TODO: What if they're uppercase?
      delete @args{ qw( -name -display_name -display_names
                        -displaynames -displaynames
                        name names display_name display_names
                        displaynames displaynames ) };
      ## TODO: REMOVE
      #warn "\$names is now [ ".join( ', ', @$names )." ]" if Bio::Graphics::Browser::DEBUG;
    }
  }
  
  unless( %args ) {
    ## TODO: REMOVE
    warn "Returning \$self.";
    return $self;
  }

  ## TODO: Dehackify
  @_ = %args;
  
  ## Use our own features() method to do the hard work..
  my @features = $self->features( @_ );
  return unless( @features ); ## NOTE: We return if there's no features...
  ## TODO: REMOVE
  #warn "Bio::DB::GFF::Segment [$self]: got these features: [ ".join( ', ', @features )." ]\n";
  my $segment = $self->_create_collection( \@_, @features );
  ## TODO: REMOVE?  Testing.
  ## TODO: REMOVE
  #warn "bouts to adjust bounds of $segment, which has these features: [ ".join( ', ', $segment->features() )." ]\n";
  $segment->adjust_bounds();
  return $segment;
} # get_collection(..)

=head2 features

 Title   : features
 Usage   : @features = $segment->features( %args );
           OR
           @features = $segment->features( @types );
 Returns : a list of L<Bio::SeqFeatureI> objects,
           OR
           (when the -iterator option is true) an L<Bio::SeqFeature::IteratorI>
           OR
           (when the -callback argument is given) true iff the callbacks
             completed.
 Args    : see below
 Status  : Public

This routine will retrieve features associated with this segment
object.  It can be used to return all features, or a subset based on
their type, location, or attributes.  Features that are returned in
relative mode (relative either to this SegmentI or to a given RangeI)
will be returned with coordinates that are relative.  Features that
are returned in absolute mode will be returned with absolute
coordinates.  The mode is determined by the -baserange and -absolute
arguments and by the absolute() flag, in that precedence order.

If ranges are specified using the -ranges argument, then these ranges
will be used to narrow the results, according to the specified
-rangetype and -strandtype arguments.

If no ranges are specified but the -rangetype argument is given then a
special and strange thing happens: the method call is delegated to the
parent_segment_provider.  If it is a SegmentI then its features()
method will be called with all the same arguments but with *this*
segment as the -range argument.  If the parent_segment_provider is a
L<Bio::DB::SegmentProviderI> (but not a SegmentI) then the same thing
will happen, but to the SegmentI returned by its get_collection()
method with no arguments.  If the parent_segment_provider is null then
no features will be returned.

If a -baserange is specified then unqualified ranges given with the
-ranges argument will be interpreted as relative to that baserange,
and qualified ranges will be re-relativized to the baserange.  If no
-baserange is given then a default will be provided that will depend
on the value of the -absolute argument or the absolute() flag.  If
-absolute is given and true or if absolute() is true then the default
baserange is the value returned by the abs_seq_id() method; if
( -absolute || absolute() ) is false then the default is this SegmentI
object ($self).  You may force absolute range interpretations by
giving a -baserange that is not a L<Bio::RangeI> (such as the string
'absolute', though any string will do the trick), by providing a true
value to the -absolute argument, or by setting the absolute() flag to
true.

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
                 the -range argument is relative.  The default
                 baserange depends on the value of the absolute()
                 flag.  If absolute() is true then the default is the
                 value of the abs_seq_id() method.  If absolute() is
                 false then the default is $self.  Note that the
                 baserange affects the sort order.  See also
                 -absolute.

  -absolute      If -absolute is given and true then all behavior will be as
                 if this SegmentI's absolute() flag was set to true,
                 even if it isn't.  If -absolute is given and false
                 then all behavior will be as if this SegmentI's
                 absolute() flag was set to false, even if it isn't.
                 Note that -baserange can still be given and can force
                 relativeness, and that takes precedence over -absolute.

  -range         A L<Bio::RangeI> object defining the range to search.
                 See also -strandmatch, -rangetype, and -baserange.
  -ranges        An array reference to multiple ranges.

  -rangetype     One of "overlaps", "contains", or "contained_in".  If no
                 range is given then a strange thing happens (it is
                 described above).

  -strandmatch   One of "strong", "weak", or "ignore".  Note that the strand
                 attribute of a given -range must be non-zero for this to work
                 (a 0/undef strand forces a 'weak' strandmatch to become
                 'ignore' and cripples the 'strong' strandmatch).

  -iterator      Return a L<Bio::SeqFeature::IteratorI>

  -callback      A callback to invoke on each feature

  -sort          Return the features in order (of their start positions).
                 Note that if sorted() is true, then this argument is
                 redundant (the features will be returned in order
                 regardless).  If the baserange (see -baserange) has a
                 negative strand then the sort order will be reversed.

All plural arguments are interchangeable with their singular counterparts.

The -attributes argument is a hashref containing one or more
attributes to match against:

  -attributes => { Gene => 'abc-1',
                   Note => 'confirmed' }

Attribute matching is simple string matching, and multiple attributes
are ANDed together.  More complex filtering can be performed using the
-callback option (see below).

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

If -iterator is true, then the method returns an object of type
Bio::SeqFeature::IteratorI.  Each call to next_seq() on this
object returns a Bio::SeqFeatureI object from this collection.

If -callback is passed a code reference, the code reference will be
invoked on each feature returned.  The code will be passed two
arguments consisting of the current feature and this SegmentI
object, and must return a true value. If the code returns a false
value, feature retrieval will be aborted.

-callback and -iterator are mutually exclusive options.  If -iterator
is defined, then -callback is ignored.

-callback and -sort are mutually exclusive options.  If -sort is
defined, then -callback is ignored.  If you want to do a sorted
callback, set the sorted() flag of this collection to true.

If -sort or sorted() is true then the features will be returned in
order of the features' start positions.  This order will be reversed
if the baserange has a negative strand (remember that the default
baserange depends upon the value of the absolute() flag, but this may
be overridden by the -baserange argument).

Note that no guarantees are made by the SegmentI interface about
the order of the features, except when the sorted() flag is true or
when the -sort option is given to the features method.  Therefore
the implementation may choose to reorder the underlying data structure
to better accomodate -sorted feature requests as a result of a
features() call.  When this happens the SegmentI's sorted() flag
should be set to true, so that the client can detect that the -sorted
argument to features() is now irrelevant.

NOTE: the following methods all build on top of features(), and do not
need to be explicitly implemented.

    features_in_range()
    overlapping_features()
    contained_features()
    contained_in()
    get_feature_stream()
    get_feature_by_name()
    get_feature_by_id()
    get_feature_by_attribute()

=cut

sub features {
  my $self = shift;

  my ( $types, $unique_ids, $namespace, $names, $attributes, $baserange, $ranges, $strandmatch, $rangetype, $iterator, $callback, $sort, $sparse, $merge, $seq_id, $start, $end );
  if( scalar( @_ ) && $_[ 0 ] =~ /^-/ ) {
    ( $types, $unique_ids, $namespace, $names, $attributes, $baserange, $ranges, $strandmatch, $rangetype, $iterator, $callback, $sort, $sparse, $merge, $seq_id, $start, $end ) =
      rearrange(
        [ [ qw( TYPE TYPES ) ],
          [ qw( UNIQUE_IDS UNIQUEIDS IDS UNIQUE_ID UNIQUEID ID ) ],
          [ qw( NAMESPACE NAME_SPACE CLASS ) ],
          [ qw( NAME NAMES DISPLAY_NAME DISPLAY_NAMES DISPLAYNAME DISPLAYNAMES ) ],
          [ qw( ATTRIBUTE ATTRIBUTES ) ],
          [ qw( BASERANGE BASE_RANGE BASELOCATION BASE_LOCATION BASE ) ],
          [ qw( RANGE RANGES LOCATION LOCATIONS LOC ) ],
          [ qw( STRANDMATCH STRAND_MATCH ) ],
          [ qw( RANGETYPE RANGE_TYPE ) ],
          [ qw( ITERATOR STREAM ) ],
          [ qw( CALLBACK CALL_BACK ) ],
          [ qw( SORT SORTED ) ],
          [ qw( RARE SPARSE ) ],
          [ qw( MERGE AUTOMERGE ) ],
          [ qw( SEQ_ID SEQID SEQ ID REF REF_SEQ REFSEQ ) ],
          [ qw( START BEGIN LOW ) ],
          [ qw( STOP END HIGH ) ]
        ],
        @_
      );
  } else {
    ## Types.
    $types = \@_;
  }

  my $closure_options = $self->{ '_closure_options' } || {};

  ## Fix up types.
  if( $types ) {
    unless( ref $types eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $types = [ $types ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$types ) ||
       ( ( scalar( @$types ) == 1 ) && !( $types->[ 0 ] ) )
      ) {
      undef $types;
    }
  }

  if( $closure_options->{ 'types' } ) {
    if( $types ) {
      unshift( @$types, @$closure_options->{ 'types' } );
    } else {
      $types = $closure_options->{ 'types' };
    }
  }

  ## Fix up unique_ids.
  if( $unique_ids ) {
    unless( ref $unique_ids eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $unique_ids = [ $unique_ids ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$unique_ids ) ||
       ( ( scalar( @$unique_ids ) == 1 ) && !( $unique_ids->[ 0 ] ) )
      ) {
      undef $unique_ids;
    }
  }

  if( $closure_options->{ 'unique_ids' } ) {
    if( $unique_ids ) {
      unshift( @$unique_ids, @$closure_options->{ 'unique_ids' } );
    } else {
      $unique_ids = $closure_options->{ 'unique_ids' };
    }
  }

  ## Now we need to remove from $unique_ids the id of this Segment, if
  ## it is in there.
  if( $self->unique_id() && $unique_ids ) {
    my $unique_id = $self->unique_id();

    # Remove any entries that are the same as the unique_id, perhaps
    # taking into account the namespace..
    @$unique_ids =
      grep { ( $_ ne $unique_id ) && ( $namespace ? ( ( ( $namespace.$_ ) ne $unique_id ) && ( $_ ne ( $namespace.$unique_id ) ) ) : 1 ) } @$unique_ids;
  }

  ## Fix up names.
  if( $names ) {
    unless( ref $names eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $names = [ $names ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$names ) ||
       ( ( scalar( @$names ) == 1 ) && !( $names->[ 0 ] ) )
      ) {
      undef $names;
    }
  }

  if( $closure_options->{ 'names' } ) {
    if( $names ) {
      unshift( @$names, @$closure_options->{ 'names' } );
    } else {
      $names = $closure_options->{ 'names' };
    }
  }

  ## Alas, now we need to remove from $names the name of this Segment,
  ## if that name is in there.
  my $unique_id = $self->unique_id();
  my $display_name = eval{ $self->display_name() };
  my $display_id = $self->display_id();
  my $primary_id = $self->primary_id();
  my $accession_number = $self->accession_number();
  if( $names && ( $unique_id || $display_name || $display_id || $primary_id || $accession_number ) ) {

    # Remove any entries that are the same as the name of this segment,
    # perhaps taking into account the namespace..
    @$names =
      grep { ( $_ ne $unique_id ) && ( $namespace ? ( ( ( $namespace.$_ ) ne $unique_id ) && ( $_ ne ( $namespace.$unique_id ) ) ) : 1 ) && ( $_ ne $display_name ) && ( $namespace ? ( ( ( $namespace.$_ ) ne $display_name ) && ( $_ ne ( $namespace.$display_name ) ) ) : 1 ) && ( $_ ne $display_id ) && ( $namespace ? ( ( ( $namespace.$_ ) ne $display_id ) && ( $_ ne ( $namespace.$display_id ) ) ) : 1 ) && ( $_ ne $primary_id ) && ( $namespace ? ( ( ( $namespace.$_ ) ne $primary_id ) && ( $_ ne ( $namespace.$primary_id ) ) ) : 1 ) && ( $_ ne $accession_number ) && ( $namespace ? ( ( ( $namespace.$_ ) ne $accession_number ) && ( $_ ne ( $namespace.$accession_number ) ) ) : 1 ) } @$names;
  }

  ## Attributes better be a hash ref if it is anything.
  if( $attributes ) {
    unless( ref( $attributes ) eq 'HASH' ) {
      $self->throw( "The -attributes argument must be a HASH REF." );
    }
  }

  if( $closure_options->{ 'attributes' } ) {
    if( $attributes ) {
      unshift( @$attributes, @$closure_options->{ 'attributes' } );
    } else {
      $attributes = $closure_options->{ 'attributes' };
    }
  }

  ## Fix up ranges.
  if( $ranges ) {
    unless( ref $ranges eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $ranges = [ $ranges ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$ranges ) ||
       ( ( scalar( @$ranges ) == 1 ) && !( $ranges->[ 0 ] ) )
      ) {
      undef $ranges;
    }
  } elsif( defined( $seq_id ) || defined( $start ) || defined( $end ) ) {
    $ranges = [ Bio::RelRange->new(
                  '-seq_id' => $seq_id,
                  '-start' => $start,
                  '-end' => $end,
                   -orientation_policy => 'dependent'
                ) ];
  }

  if( $closure_options->{ 'ranges' } ) {
    if( $ranges ) {
      unshift( @$ranges, @$closure_options->{ 'ranges' } );
    } else {
      $ranges = $closure_options->{ 'ranges' };
    }
  }

  if( !$baserange && $closure_options->{ 'baserange' } ) {
    $baserange = $closure_options->{ 'baserange' };
  }

  ## We can use ourselves as the baserange, if we are a Bio::RangeI.
  if( $self->isa( 'Bio::RangeI' ) && !$baserange ) {
    $baserange = $self;
  }

  ## Derelativize, man.
  if( $ranges && @$ranges && $baserange ) {
    my @new_ranges;
    foreach my $range ( @$ranges ) {
      unless( ref( $range ) && $range->isa( 'Bio::RangeI' ) ) {
        $self->throw( "Expected the -ranges argument to be a reference to a list of Bio::RangeI objects, but it contains something incompatible: " . overload::StrVal( $range ) );
      }
      unless( defined $range->seq_id() ) {
        $range = Bio::RelRange->new(
                   -seq_id => $baserange,
                   -start =>  $range->start(),
                   -end =>    $range->end(),
                   -strand => $range->strand(),
                   -orientation_policy => 'dependent'
                 );
      } elsif( !$range->isa( 'Bio::RelRangeI' ) ) {
        $range = Bio::RelRange->new(
                   -seq_id => $range->seq_id(),
                   -start =>  $range->start(),
                   -end =>    $range->end(),
                   -strand => $range->strand(),
                   -orientation_policy => 'dependent'
                 );
      }
      $range->absolute( 1 );
      push( @new_ranges, $range );
    }
    $ranges = \@new_ranges;
  } # End derelativizing ranges.  Now they're absolute.

  ## -strandmatch defaults to 'weak'
  unless( defined $strandmatch ) {
    $strandmatch = $closure_options->{ 'strandmatch'} || 'weak';
  }

  ## -rangetype defaults to 'overlaps'
  unless( defined $rangetype ) {
    $rangetype = $closure_options->{ 'rangetype'} || 'overlaps';
  }

  ## -iterator and -callback are mutually exclusive.
  if( $iterator && $callback ) {
    $self->throw( "The -iterator and -callback options are mutually exclusive, and yet both have been specified.  Perhaps you could apply your callback method to each element returned by the iterator?" );
  }

  ## -sort is redundant if sorted() is true.
  unless( defined( $sort ) || $callback ) {
    $sort = $closure_options->{ 'sort' };
  }
  if( $sort && $self->sorted() ) {
    undef $sort;
  }

  ## -sort and -callback are mutually exclusive.
  if( $sort && $callback ) {
    $self->throw( "The -sort and -callback options are mutually exclusive, and yet both have been specified.  Perhaps you could first set the sorted() property of this collection to true, and then try again." );
  }

  unless( defined $namespace ) {
    $namespace = $closure_options->{ 'namespace' };
  }

  unless( defined $sparse ) {
    $sparse = $closure_options->{ 'sparse' };
  }

  unless( defined $merge ) {
    $merge = $closure_options->{ 'merge' };
  }

  ## We actually implement the sorting on-the-fly anyway, so the above
  ## redundancy check is silly (but we keep it for the mutex check)
  $sort ||= $self->sorted();
  ## Do a reverse sort iff the baserange has a negative strand.
  my $reverse_sort = ( $baserange && ( $baserange->strand() < 0 ) );

  my %args = (
    '-ref' => ( $ranges && @$ranges ? $ranges->[ 0 ]->seq_id() : undef ) ||
      $self->abs_seq_id(),
    '-class' => $namespace || $self->class(),
    '-start' => ( $ranges && @$ranges ? $ranges->[ 0 ]->low( 'plus' ) : undef ) ||
      ( $self->{ 'whole' } ? undef : $self->abs_start( 'plus' ) ),
    '-end' => ( $ranges && @$ranges ? $ranges->[ 0 ]->high( 'plus' ) : undef ) ||
      ( $self->{ 'whole' } ? undef : $self->abs_end( 'plus' ) ),
    '-range_type' => $rangetype,
    '-types' => $types,
    '-unique_ids' => $unique_ids,
    '-names' => $names,
    '-attributes' => $attributes,
    '-sparse' => $sparse,
    '-merge' => $merge,
    '-iterator' => $iterator,
    '-callback' => $callback,
    '-parent' => $self
  );
  ## TODO: REMOVE
  #warn "Args to features() of ".$self->factory(). " are ( ".join( ', ', ( my @args = %args ) )." ).  names are { ".($names?join( ', ', @$names ):'')." }.  types are { ".($types?join( ', ', @$types ):''). " }.";
  #eval { $self->throw( "Args to features() of ".$self->factory(). " are ( ".join( ', ', ( my @args = %args ) )." )"."-names are ".($names&&@$names?'[ '.join( ', ', @$names ).' ]':'undef')."." ) };
  #warn $@ if $@;
  if( $iterator ) {
    ## TODO: Do what we do for the list, but in an iterator.
    ## TODO: REMOVE
    warn "May not be filtering appropriately\n";
    return $self->factory()->features_in_range( %args );
  } elsif( $callback ) {
    ## TODO: Do what we do for the list, but in a callback.
    ## TODO: REMOVE
    warn "May not be filtering appropriately\n";
    return $self->factory()->features_in_range( %args );
  }
  ## TODO: This doesn't account for the whole picture.  Update.
  my @features = $self->factory()->features_in_range( %args );

  ## The implementation may not do the filtering correctly, so we have
  ## to do it again, alas.
  my @features_to_return;
  ## This next line is fancy speak for "for each feature we've got,
  ## perhaps sorted, perhaps reverse-sorted.."
  foreach my $feature ( @features ) {
    ## TODO: REMOVE
    unless( ref( $feature ) && $feature->isa( 'Bio::SeqFeature::SegmentI' ) ) {
      print STDERR "--------------Hey, man.  $feature ain't a feature!=======\n";
    }

    ## Filter them:
    my $passes_filter = 1;

    ## Filter on types
    if( $passes_filter && $types ) {
      # Failure until proven success:
      $passes_filter = 0;
      my $feature_type = $feature->type();
      foreach my $type ( @$types ) {
        if( ref( $feature_type ) &&
            $feature_type->isa( 'Bio::DB::GFF::Typename' ) &&
            $feature_type->match( $type ) ) {
          # success
          $passes_filter = 1;
          last;
        } elsif( ref( $type ) && $type->isa( 'Bio::SeqFeature::TypeI' ) ) {
          if( ( $type eq $feature_type ) ||
              $type->is_descendent( $feature_type ) ) {
            # success
            $passes_filter = 1;
            last;
          }
        } elsif( ref( $type ) && $type->isa( 'Bio::DB::GFF::Typename' ) ) {
          ## I've added this to support the Typename stuff used by GFF.
          if( $type->match( $feature_type ) ) {
            # success
            $passes_filter = 1;
            last;
          }
        } else {
          ## $type is a string.
          if( ref( $feature_type ) &&
              $feature_type->isa( 'Bio::SeqFeature::TypeI' ) ) {
            if( ( $feature_type eq $type ) ||
                $feature_type->is_ancestor( $type ) ) {
              # success
              $passes_filter = 1;
              last;
            }
          } else {
            ## $feature_type is also a string.
            if( $type =~ /$feature_type/i ) { # ignore case (should we? ok)
              # success
              $passes_filter = 1;
              last;
            }
          }
        } # End if $type is a TypeI .. else ..
      } # End for each $type
    } # End if we're filtering on types

    ## Filter on unique_ids
    ## If they've provided a $namespace, also try "$namespace:$feature_id"
    ## and "$namespace:$id"
    if( $passes_filter && $unique_ids ) {
      # Failure until proven success:
      $passes_filter = 0;
      my $feature_id = $feature->unique_id() || $feature->display_name();
      foreach my $unique_id ( @$unique_ids ) {
        if( ( $feature_id eq $unique_id ) ||
            ( $namespace && ( ( $namespace.':'.$feature_id ) eq $unique_id ) ) ||
            ( $namespace && ( $feature_id eq ( $namespace.':'.$unique_id ) ) )
          ) {
          # success
          $passes_filter = 1;
          last;
        }
      } # End for each $unique_id
    } # End if we're filtering on unique_ids

    ## Filter on names
    ## If they've provided a $namespace, also try "$namespace:$feature_name"
    ## and "$namespace:$name" and "$namespace:$feature_id"
    if( $passes_filter && $names ) {
      # Failure until proven success:
      $passes_filter = 0;
      my $feature_name = $feature->display_name();
      my $feature_id = $feature->unique_id();
      foreach my $name ( @$names ) {
        if( ( $feature_name eq $name ) ||
            ( $namespace && ( ( $namespace.':'.$feature_name ) eq $name ) ) ||
            ( $namespace && ( $feature_name eq ( $namespace.':'.$name ) ) ) ||
            ( $feature_id eq $name ) ||
            ( $namespace && ( ( $namespace.':'.$feature_id ) eq $name ) ) ||
            ( $namespace && ( $feature_id eq ( $namespace.':'.$name ) ) )
          ) {
          # success
          $passes_filter = 1;
          last;
        }
      } # End for each $name
    } # End if we're filtering on names

    ## Filter on attributes
    if( $passes_filter && $attributes ) {
      foreach my $tag ( keys %$attributes ) {
        if( !$feature->has_tag( $tag ) ) {
          $passes_filter = 0;
          last;
        }
        $passes_filter = 0;
        foreach my $value ( $feature->get_tag_values( $tag ) ) {
          if( $value eq $attributes->{ $tag } ) {
            $passes_filter = 1;
            last;
          }
        } # End for each $value
        last unless( $passes_filter );
      } # End for each $tag
    } # End if we're filtering on attributes

    ## Filter on range
    if( $passes_filter && $ranges ) {
      foreach my $range ( @$ranges ) {
        if( $range->seq_id() &&
            $feature->seq_id() &&
            !( $feature->seq_id() eq $range->seq_id() )
          ) {
          ## If they both define seq_id(), and they don't match, then
          ## the ranges don't match.
          $passes_filter = 0;
        } elsif(
            ( $rangetype eq 'overlaps' ) &&
            !$range->overlaps( $feature, $strandmatch )
          ) {
          $passes_filter = 0;
        } elsif(
            ( $rangetype eq 'contains' ) &&
            !$range->contains( $feature, $strandmatch )
          ) {
          $passes_filter = 0;
        } elsif(
            ( $rangetype eq 'contained_in' ) &&
            !$feature->contains( $range, $strandmatch )
          ) {
          $passes_filter = 0;
        }
        last unless $passes_filter;
      } # End foreach $range in @$ranges
    } # End if we're filtering on range

    next unless( $passes_filter );

    if( defined $callback ) {
      unless( $callback->( $feature, $self ) ) {
        return 0; # return false when the callback fails.
      }
    } else {
      push( @features_to_return, $feature );
    }
  } # End for each $feature

  return @features_to_return;
} # features(..)

=head2 types

 Title   : types
 Usage   : @types = $s->types( [ '-enumerate' => 1 ] )
 Function: list feature types that overlap this segment
 Returns : a list of Bio::DB::GFF::Typename objects or a hash
 Args    : see below
 Status  : Public

The types() method will return a list of Bio::DB::GFF::Typename
objects, each corresponding to a feature that overlaps the segment.
If the optional -enumerate parameter is set to a true value, then the
method will return a hash in which the keys are the type names and the 
values are the number of times a feature of that type is present on
the segment.  For example:

  %count = $s->types(-enumerate=>1);

=cut 

# wrapper for lower-level types() call.
sub types {
  my $self = shift;

  my @args;
  if( @_ && $_[ 0 ] !~ /^-/ ) {
    @args = ( '-type' => \@_ );
  } else {
    @args = @_;
  }
  $self->factory->types(
    '-ref'   => $self->abs_seq_id(),
    '-class' => $self->class(),
    '-start' => $self->abs_start( 'plus' ),
    '-stop'  => $self->abs_end( 'plus' ),
    @args
  );
} # types(..)

sub feature_count {
  my $self = shift;

  # count all feature types in the segment
  my %type_counts = $self->types( '-enumerate' => 1 );

  my $count = 0;
  map { $count += $type_counts{ $_ } } keys %type_counts;
  return $count;
} # feature_count(..)

=head2 stop

 Title   : stop
 Usage   : $s->stop
 Function: stop of segment
 Returns : integer
 Args    : none
 Status  : Public

This is an alias for end(), provided for AcePerl compatibility.

=cut

sub stop {
  shift->end( @_ );
} # stop(..)

=head2 sourceseq

 Title   : sourceseq
 Usage   : $s->sourceseq
 Function: get the segment source
 Returns : a string
 Args    : none
 Status  : Public

An alias for seq_id.

=cut

sub sourceseq {
  shift->seq_id( @_ );
}

=head2 class

 Title   : class
 Usage   : $s->class([$newclass])
 Function: get the source sequence class
 Returns : a string
 Args    : new class (optional)
 Status  : Public

Gets or sets the class for the source sequence for this segment.

=cut

sub class     { 
  my $self = shift;
  my $d = $self->{class};
  $self->{class} = shift if @_;
  $d;
}

=head2 subseq

 Title   : subseq
 Usage   : $s->subseq( $new_start, $new_end )
 Function: Generate a sub-segment, relative to the same refseq as this segment.
 Returns : A L<Bio::DB::GFF::Segment> object
 Args    : Start and end of subsegment, relative to this segment.
 Status  : Public

This method generates a new segment from the start and end positions
given in the arguments.  If end E<lt> start, then the strand is
reversed.  Note that although the coordinates are given relative to
this segment, the returned segment will have start and end positions
relative to this segment's current seq_id (and it will have the same seq_id).

=cut

sub subseq {
  my $self = shift;
  my ( $new_start, $new_end ) = @_;

  # We deliberately force subseq to return objects of type Segment
  # Otherwise, when we get a subsequence from a Feature object,
  # its method and source go along for the ride, which is incorrect.
  my $new = $self->new_from_segment( $self );
  if( $new_start <= $new_end ) {
    my $temp = $new_start;
    $new_start = $new_end;
    $new_end = $temp;
    $new->strand( ( $self->strand() || 1 ) * -1 );
  }

  $new->start( 'plus', ( $self->start( 'plus' ) + $new_start - 1 ) );
  $new->end( 'plus', ( $self->start( 'plus' ) + $new_end - 1 ) );

  return $new;
} # subseq(..)

=head2 seq

 Title   : seq
 Usage   : $s->seq
 Function: get the sequence string for this segment
 Returns : a string
 Args    : none
 Status  : Public

Returns the sequence for this segment as a simple string.  (-) strand
segments are automatically reverse complemented

This method is also called dna() and protein() for backward
compatibility with AceDB.

=cut

sub seq {
  my $self = shift;
  return $self->factory()->dna(
    '-id'    => $self->abs_seq_id(),
    '-start' => $self->abs_start( 'plus' ),
    '-end'   => $self->abs_end( 'plus' ),
    '-class' => $self->class()
  );
} # seq(..)

*protein = *dna = \&seq;

=head2 equals

  Title   : equals
  Usage   : $s->equals( $d )
  Function: segment equality
  Returns : true, if two segments are equal
  Args    : another segment
  Status  : Public

  Returns true if the two segments have the same abs_seq_id,
  abs_start, and abs_end.

=cut

sub equals {
  my $self = shift;
  my $peer = shift;
  return unless defined $peer;
  if( ref( $peer ) && $peer->isa( 'Bio::DB::GFF::Segment' ) ) {
    return ( ( $self->abs_start() eq $peer->abs_start() ) &&
             ( $self->abs_end() eq $peer->abs_end() ) &&
             #( $self->seq_id() eq $peer->seq_id() ) &&
             ( $self->abs_seq_id() eq $peer->abs_seq_id() ) );
  } else {
    return ( $self->asString() eq $peer );
  }
} # equals(..)

=head2 asString

 Title   : asString
 Usage   : $s->asString
 Function: human-readable string for segment
 Returns : a string
 Args    : none
 Status  : Public

Returns a human-readable string representing this sequence.  Format
is:

   sourceseq:start,stop

=cut

sub asString {
  ## TODO: REMOVE
  #return shift->toString();

  my $self = shift;
  my $label = $self->refseq();
  my $start = $self->start( 'plus' );
  my $stop  = $self->stop( 'plus' );
  return "$label:$start,$stop";
}

=head2 clone

 Title   : clone
 Usage   : $copy = $s->clone
 Function: make a copy of this segment
 Returns : a L<Bio::DB::GFF::Segment> object
 Args    : none
 Status  : Public

This method creates a copy of the segment and returns it.

=cut

# deep copy of the thing
sub clone {
  my $self = shift;
  my %h = %$self;
  return bless \%h,ref($self);
}

=head2 error

 Title   : error
 Usage   : $error = $s->error([$new_error])
 Function: get or set the last error
 Returns : a string
 Args    : an error message (optional)
 Status  : Public

In case of a fault, this method can be used to obtain the last error
message.  Internally it is called to set the error message.

=cut

sub error {
  my $self = shift;
  my $g = $self->{error};
  $self->{error} = shift if @_;
  $g;
}

=head2 abs_ref

 Title   : abs_ref
 Usage   : $s->abs_ref
 Function: the reference sequence for this segment
 Returns : a string
 Args    : none
 Status  : Public

This is an alias to abs_seq_id().

=cut

sub abs_ref {
  shift->abs_seq_id( @_ );
}

# Internal overridable getter/setter for the actual stored value of seq_id.
# This one can take an optional class argument as well, so the usage becomes:
# Usage   : my $seq_id = $segment->_seq_id( [$new_seq_id] );
#           OR
#           my $old_seq_id = $segment->_seq_id( $new_seq_class, $new_seq_id );
#
# This method will get or set the reference sequence.  Called with no
# arguments, it returns the current reference sequence.  Called with
# either a class and sequence ID, a Bio::RelRangeI object (or
# subclass), or a Bio::DB::GFF::Featname object, it will set the current
# reference sequence and return the previous one.
# 
# The method will generate an exception if you attempt to set the
# reference sequence to a sequence that isn't contained in the database,
# or one that has a different source sequence from the segment.
sub _seq_id {
  my $self = shift;
  my $old_val = $self->{ '_seq_id' };
  if( @_ ) {
    my ( $newref, $newclass );
    if( @_ == 2 ) {
      $newclass = shift;
      $newref   = shift;
    } else {
      $newref   = shift;
      $newclass = 'Sequence';
    }

    # support for Featname objects
    if( ref( $newref ) && $newref->isa( 'Bio::DB::GFF::Featname' ) ) {
      $newclass = $newref->class();
    }

    my $abs_seq_id = $self->abs_seq_id();
    my $new_seq_id;
    if( ref( $newref ) && $newref->isa( 'Bio::RangeI' ) ) {
      if( $self->{ '_seq_id' } &&
          ref( $self->{ '_seq_id' } ) &&
          $self->{ '_seq_id' }->isa( 'Bio::RangeI' ) &&
          $self->{ '_seq_id' } eq $newref ) {

        # Attempt to set the seq_id to what it already is.  Who am I
        # to complain?  Go ahead.
        # (We've bothered to test for this because the elsif clause
        # can grab this case sometimes and it'd be rude to throw an
        # exception).
      } elsif( $self eq $newref ) {
        $self->throw( "Unable to set the reference sequence to \$self: \$self is $self, attempted new seq_id is $newref." );
      }
      # Avoid circularity by forcing new_seq_id to be defined relative
      # to the absolute range.
      my $name =
        ( $newref->isa( 'Bio::DB::GFF::Segment_NamedRelRange' ) ?
          $newref->name() :
          ( $newref->isa( 'Bio::RangeI' ) ?
            $newref->seq_id() :
            undef ) );
      ## TODO: Perhaps we don't *always* have to do this.. but this
      ## works for now.
      $new_seq_id = Bio::DB::GFF::Segment_NamedRelRange->new(
        '-name' => $name,
        '-seq_id' => absSeqId( $newref ),#absRange( $newref ),
        '-start' => absStart( $newref, 'plus' ),
        '-end' => absEnd( $newref, 'plus' ),
        '-strand' => absStrand( $newref )
      );
      ## TODO: REMOVE
      #warn "CREATED new_seq_id $new_seq_id (RelRangeString is ".$new_seq_id->toRelRangeString( 'both', 'plus' ).") from newref $newref (RelRangeString is ".$newref->toRelRangeString( 'both', 'plus' ).")";
    } else {
      my ( $newref_abs_seq_id, $refclass, $refstart, $refstop, $refstrand, $sname );
      my $coords = $self->factory()->abscoords( $newref, $newclass );
      foreach ( @$coords ) { # find the appropriate one
	( $newref_abs_seq_id, $refclass, $refstart, $refstop, $refstrand, $sname ) =
          @$_;
	last if( !defined( $abs_seq_id ) ||
                 ( $newref_abs_seq_id ne $abs_seq_id ) );
      }
      $new_seq_id = Bio::DB::GFF::Segment_NamedRelRange->new(
        '-name' => $sname,
        '-seq_id' => $newref_abs_seq_id,
        '-start' => $refstart,
        '-end' => $refstop,
        '-strand' => $refstrand
      );

      # Make sure that our new seq_id is rooted at a real sequence.
      if( $new_seq_id->isa( 'Bio::DB::GFF::Segment_NamedRelRange' ) &&
          ( $new_seq_id->name() ne $new_seq_id->seq_id() ) 
        ) {
        my $a =
          $self->factory()->abscoords( $new_seq_id->seq_id(), 'Sequence' );
        my $seq_start = $a->[ 0 ][ 2 ];
        my $seq_stop  = $a->[ 0 ][ 3 ];
        ## Paul's question: why would $seq_start ever not be 1?  How
        ## should that affect absolute positions?
        $new_seq_id->seq_id(
          Bio::RelRange->new(
            '-seq_id' => $a->[ 0 ][ 0 ],
            '-start' => $a->[ 0 ][ 2 ],
            '-end' => $a->[ 0 ][ 3 ],
            '-strand' => 1,
            '-orientation_policy' => 'dependent'
          )
        );
      } # End if we need to get the root sequence's bounds
    }

    if( defined( $abs_seq_id ) &&
        ( $new_seq_id->abs_seq_id() ne $abs_seq_id )
      ) {
      $self->throw( "Can't set reference sequence: $newref and $self are on different sequences ($newref is on '".$new_seq_id->abs_seq_id()."', $self is on '".$abs_seq_id."')." );
    }

    # Our strategy is to get absolute coords for $self and then
    # rerelativize them to the new_seq_id.
    # We can use another range that is relative to $new_seq_id to
    # convert the coords.  Anything will do, so long as it starts at
    # 1 and is on the + strand.
    my $temp_rel_range =
      Bio::RelRange->new(
        '-seq_id' => $new_seq_id,
        '-start' => 1,
        '-end' => 1, # we could use $new_seq_id->length, but it's unnecessary.
        '-strand' => 1,
        '-orientation_policy' => 'dependent'
      );
    my $abs_self = new Bio::RelRange( $self );
    $abs_self->absolute( 1 );
    $self->{ '_start' } =
      $temp_rel_range->abs2rel( $abs_self->start( 'plus' ) );
    $self->{ '_end' } =
      $temp_rel_range->abs2rel( $abs_self->end( 'plus' ) );
    $self->{ '_strand' } =
      $temp_rel_range->abs2rel_strand( $abs_self->strand() );

    if( $self->{ '_end' } < $self->{ '_start' } ) {
      my $tmp = $self->{ '_start' };
      $self->{ '_start' } = $self->{ '_end' };
      $self->{ '_end' } = $tmp;
    }

    $self->{ '_seq_id' } = $new_seq_id;
    ## TODO: ERE I AM.  If we're going to make it so that segments with string seq_ids have 'plus' position_policies (we've added code to new(..) for this, BTW), then we're going to want to switch now if the new seq_id is not a string and the old one was a string.  Actually we're going to want to make GFF work in a plus mode at all times, I think.  We're going to want to be able to know that from outside I think, so that the constructors are compatible... or something.  Actually the constructor should only be called from _create_segment(..) in other non-GFF classes, so then it's just a matter of making GFF's methods obey the contract, which should explicitly state a stranded assumption.  Anyways we're going to have to alter RelRange and RelRangeI so that you can be sure that the stored values are always stranded (or better, always 'plus' so that we don't need to know the length of the parent sequence), and also to do the right thing on setting the values as well as on retrieving them.
  }
  ## TODO: This is what this was from RelSegment::refseq(..):  why sourceseq?
  ## TODO: Note that presently we won't allow setting in absolute mode anyways.
  # return $self->absolute ? $self->sourceseq : $old_value;
  return $old_val;
} # _seq_id(..)

=head2 refseq

 Title   : refseq
 Usage   : $s->refseq
 Function: get or set the reference sequence
 Returns : a string
 Args    : none
 Status  : Public

Examine or change the reference sequence. This is an alias to
seq_id().

=cut

sub refseq {
  shift->sourceseq( @_ );
}

=head2 ref

 Title   : ref
 Usage   : $s->refseq
 Function: get or set the reference sequence
 Returns : a string
 Args    : none
 Status  : Public

An alias for refseq()

=cut

sub ref { shift->refseq(@_) }

=head2 truncated

 Title   : truncated
 Usage   : $truncated = $s->truncated
 Function: Flag indicating that the segment was truncated during creation
 Returns : A boolean flag
 Args    : none
 Status  : Public

This indicates that the sequence was truncated during creation.  The
returned flag is undef if no truncation occured.  If truncation did
occur, the flag is actually an array ref in which the first element is
true if truncation occurred on the left, and the second element
occurred if truncation occurred on the right.

=cut

sub truncated {
  my $self = shift;
  my $hash = $self->{truncated} or return;
  CORE::ref($hash) eq 'HASH' or return [1,1];  # paranoia -- not that this would ever happen ;-)
  return [$hash->{start},$hash->{stop}];
}

=head2 Bio::SeqI implementation

=cut

=head2 primary_id

 Title   : primary_id
 Usage   : $unique_implementation_key = $obj->primary_id;
 Function: Returns the unique id for this object in this
           implementation. This allows implementations to manage their
           own object ids in a way the implementaiton can control
           clients can expect one id to map to one object.

           For sequences with no accession number, this method should
           return a stringified memory location.

 Returns : A string
 Args    : None
 Status  : Virtual


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

sub display_id { shift->seq_id }

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
    return 'unknown';
}

=head2 alphabet

 Title   : alphabet
 Usage   : if( $obj->alphabet eq 'dna' ) { /Do Something/ }
 Function: Returns the type of sequence being one of
           'dna', 'rna' or 'protein'. This is case sensitive.

           This is not called <type> because this would cause
           upgrade problems from the 0.5 and earlier Seq objects.

 Returns : a string either 'dna','rna','protein'. NB - the object must
           make a call of the type - if there is no type specified it
           has to guess.
 Args    : none
 Status  : Virtual


=cut

sub alphabet{
    return 'dna'; # no way this will be anything other than dna!
}

=head2 desc

 Title   : desc
 Usage   : $seqobj->desc($string) or $seqobj->desc()
 Function: Sets or gets the description of the sequence
 Example :
 Returns : The description
 Args    : The description or none


=cut

sub desc { shift->asString }

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
   if( defined $value || ! defined $obj->{'annotation'} ) {
       $value = new Bio::Annotation::Collection() unless defined $value;
      $obj->{'annotation'} = $value;
    }
    return $obj->{'annotation'};

}

## Inner class ##############################################################
#============================================================================
# Bio::DB::GFF::Segment_NamedRelRange: A Bio::RelRange that has a
# name and uses it for overloading stringification and equality.
#============================================================================
package Bio::DB::GFF::Segment_NamedRelRange;
use Bio::RelRange;
use vars qw( @ISA );

@ISA = qw( Bio::RelRange );

use overload 
#  '""'     => 'name',
# For debugging, try using:
  '""'     => 'toString',
  eq       => 'equals',
  fallback => 1;

=head2 new

  Title   : new
  Usage   : $range = Bio::DB::GFF::Segment_NamedRelRange->new(
                       $another_range_to_copy_from,
                       $name
                     );
            OR
            $range = Bio::DB::GFF::Segment_NamedRelRange->new(
                       -name => $name,
                       -seq_id => $another_range_or_a_sequence_id,
                       -start => 100,
                       -end => 200,
                       -strand => +1,
                       -absolute => 1
                     );
  Function: generates a new Bio::DB::GFF:Segment_NamedRelRange object
  Returns : a new Bio::DB::GFF::Segment_NamedRelRange object
  Args    : a L<Bio::RangeI> to copy from and a name
            OR
              -name (required)
              two of (-start, -end, -length) - the third is calculated
              -strand (defaults to 0)
              -seq_id (not required but highly recommended)
              -absolute (defaults to 0)

    Note that if you pass a RelRangeI that is in absolute() mode, the
    copied values will be absolute, and the relative values will not
    be preserved.  To work around this, set absolute() to false before
    passing it in.

    If the -absolute argument is true, the other values will be
    interpreted relative to the given -seq_id, but then absolute()
    will be set to true.

=cut

sub new {
  my $caller = shift;
  my $self = $caller->SUPER::new( @_ );

  my $name;
  if( scalar( @_ ) && ( $_[ 0 ] =~ /^-/ ) ) {
    ( $name ) =
      $self->_rearrange( [ qw( NAME ) ], @_ );
  } else {
    ( undef, $name ) = @_;
  }
  unless( defined( $name ) ) {
    $self->throw( "You must provide a -name argument to the Segment_NamedRelRange constructor.\n" );
  }
  $self->name( $name );

  ## GFF uses the dependent orientation policy.
  $self->orientation_policy( 'dependent' );
  $self->ensure_orientation();

  return $self;
} # new(..)

=head2 name

 Title   : name
 Usage   : my $name = $range->name()
 Function: Get/set the name of this range.
 Returns : The current (or former, if used as a set method) name.
 Args    : [optional] a new name
 Status  : Public

=cut

sub name {
  my $self = shift;
  my $new_value = shift;
  my $old_value = $self->{ '_name' };

  if( defined( $new_value ) ) {
    $self->{ '_name' } = $new_value;
  }
  return $old_value;
} # name(..)

=head2 equals

  Title   : equals
  Usage   : if( $r1->equals( $r2 ) ) { do something }
  Function: Test whether $r1 has the same abs_start, abs_end, length,
            and abs_seq_id as $r2.  If $r2 is a string, test whether it
            is equal to name.
  Args    : arg #1 = a L<Bio::RangeI> or string to compare this to (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') (optional)
  Returns : true iff they are describing the same range

  If either range has no defined (abs) seq_id then (abs) seq_id will
  be ignored in the test.

  The behavior of this method differs from its behavior in
  L<Bio::RelRangeI> when the argument is a string.  When the argument
  is a string it will be tested for equality (using eq) with this
  range's name.

=cut
# '

sub equals {
  my $self = shift;
  my ( $other, $strand_option ) = @_;
  if( defined( $other ) && !ref( $other ) && ( ref( \$other ) eq 'SCALAR' ) ) {
    return ( $self->name() eq $other );
  } else {
    return $self->SUPER::equals( @_ );
  }
} # equals(..)

#============================================================================
## This is the end of Segment_NamedRelRange, an inner class of Segment.
#============================================================================
## End Inner class ##########################################################

1;

__END__

=head1 BUGS

Report them please.

=head1 SEE ALSO

L<bioperl>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.  

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 CONTRIBUTORS

Jason Stajich E<lt>jason@bioperl.orgE<gt>.

=cut

