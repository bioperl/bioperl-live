# $Id$
#
# BioPerl module for Bio::SeqFeatureI
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeatureI - Abstract interface of a Sequence Feature

=head1 SYNOPSIS

    # get a seqfeature somehow, eg,

    foreach $feat ( $seq->top_SeqFeatures() ) {
            print "Feature from ", $feat->start, "to ", 
	          $feat->end, " Primary tag  ", $feat->primary_tag, 
	          ", produced by ", $feat->source_tag(), "\n";

            if( $feat->strand == 0 ) {
		print "Feature applicable to either strand\n";
            } else {
                print "Feature on strand ", $feat->strand,"\n"; # -1,1
            }

            foreach $tag ( $feat->all_tags() ) {
		print "Feature has tag ", $tag, "with values, ",
		      join(' ',$feat->each_tag_value($tag)), "\n";
            }
	    print "new feature\n" if $feat->has_tag('new');
	    # features can have sub features
	    my @subfeat = $feat->get_SeqFeatures();
	}

=head1 DESCRIPTION

This interface is the functions one can expect for any Sequence
Feature, whatever its implementation or whether it is a more complex
type (eg, a Gene). This object doesn't actually provide any
implemention, it just provides the definitions of what methods one can
call. See Bio::SeqFeature::Generic for a good standard implementation
of this object

=head1 FEEDBACK

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

=head1 CONTRIBUTORS

Lincoln Stein E<lt>lstein@cshl.orgE<gt>
Paul Edlefsen E<lt>paul@systemsbiology.orgE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::SeqFeatureI;
use vars qw( @ISA );
use overload 
  '""' => 'toString',
  cmp   => '_cmp',
  '==' => 'equals';
use strict;

use Bio::RangeI;
use Bio::SeqFeature::CollectionI;
use Bio::LocallyIdentifiableI;
@ISA = qw( Bio::RangeI Bio::SeqFeature::CollectionI Bio::LocallyIdentifiableI );

use Bio::Seq;
use Carp;

=head1 Identification methods

=head2 unique_id (from Bio::LocallyIdentifiableI)

 Title   : unique_id
 Usage   : $id = $feat->unique_id( [$new_id] )
 Function: Getter/setter for the unique id for this feature
 Returns : the current (or former, if used as a set method) unique identifier (or undef)
 Args    : (optional) A new unique identifier, or "undef" if there is none
 Status  : Public

  This method will return undef if a unique identifier has not been
  set for this feature.  If the argument is the string "undef" then
  the unique_id will become undefined.  Note that the unique_id may
  not be changed (if, for instance, the implementing class does not
  allow unique_id changes).

=cut

sub unique_id {
  shift->throw_not_implemented();
}

=head2 display_name

 Title   : display_name
 Usage   : $name = $feat->display_name( [$new_name] )
 Function: Getter/setter for the human-readable name of the feature
           for displays.
 Returns : the current (or former, if used as a set method) display_name
 Args    : (optional) A new display name
 Status  : Public

=cut

sub display_name { 
  shift->throw_not_implemented();
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
  shift->throw_not_implemented();
}

=head1 Tag methods

 SeqFeatures may have arbitrary name/value pairs called tags.
 Convenience methods are provided for the most common of these
 ('primary' and 'source').  Additionally the 'display_name' and
 'identifier' attributes are available as tags.  All available tags
 may be retrieved using get_all_tags(..), and the presence of a
 particular tag may be tested using has_tag(..).

=head2 primary_tag

 Title   : primary_tag
 Usage   : $tag = $feat->primary_tag()
 Function: Returns the string value of this feature's type().
 Returns : a string 
 Args    : none
 Status  : Public

 This method is implemented in the interface in terms of type().  You
 do not have to override it, but you may.

=cut

sub primary_tag {
  my $self = shift;
  my $type = $self->type();
  if( ref( $type ) && $type->isa( 'Bio::SeqFeature::TypeI' ) ) {
    return $type->toString();
  } else {
    return "$type";
  }
} # primary_tag()

=head2 type_string

 Title   : type_string
 Usage   : $tag = $feat->type_string()
 Function: Returns the feature type as a string.
           This is a synonym for 'primary_tag'
 Returns : a string 
 Args    : none

=cut

sub type_string {
  shift->primary_tag( @_ );
}

=head2 source_tag

 Title   : source_tag
 Usage   : $tag = $feat->source_tag()
 Function: Returns the source tag for a feature,
           eg, 'genscan' 
 Returns : a string 
 Args    : none
 Status  : Public

=cut

sub source_tag {
  shift->throw_not_implemented();
}

=head2 has_tag

 Title   : has_tag
 Usage   : $tag_exists = $self->has_tag('some_tag')
 Function: 
 Returns : TRUE if the specified tag exists, and FALSE otherwise
 Args    : a tag name
 Status  : Public

=cut

sub has_tag {
  shift->throw_not_implemented();
}

=head2 get_tag_values

 Title   : get_tag_values
 Usage   : @values = $self->get_tag_values( 'some_tag' );
 Function: 
 Returns : An array comprising the values of the specified tag.
 Args    :


=cut

sub get_tag_values {
  shift->throw_not_implemented();
}

=head2 get_all_tags

 Title   : get_all_tags
 Usage   : @tags = $feat->get_all_tags()
 Function: gives all tags for this feature
 Returns : an array of strings
 Args    : none
 Status  : Public

=cut

sub get_all_tags {
  shift->throw_not_implemented();
}

=head2 get_SeqFeatures

 Title   : get_SeqFeatures
 Usage   : @feats = $feat->get_SeqFeatures();
 Function: Returns an array of sub Sequence Features
 Returns : An array
 Args    : none


=cut

sub get_SeqFeatures {
  shift->throw_not_implemented();
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
  shift->throw_not_implemented();
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
  shift->throw_not_implemented();
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
  shift->throw_not_implemented();
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

sub gff_string {
   my ($self,$formatter) = @_;

   $formatter = $self->_static_gff_formatter unless $formatter;
   return $formatter->gff_string($self);
}

my $static_gff_formatter = undef;
sub _static_gff_formatter {
   my ($self,@args) = @_;

   if( !defined $static_gff_formatter ) {
       $static_gff_formatter = Bio::Tools::GFF->new('-gff_version' => 2);
   }
   return $static_gff_formatter;
}

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
            what a user wants as "the" sequence for this feature

            Implementing classes are free to override this method with their
            own magic if they have a better idea what the user wants

  Args    : [optional] A Bio::DB::RandomAccessI compliant object
  Returns : A Bio::Seq

   This method is implemented in the interface.

=cut

sub spliced_seq {
    my ($self,$db) = shift;

    if( ! $self->location->isa("Bio::Location::SplitLocationI") ) {
	return $self->seq(); # nice and easy!
    }

    # redundant test, but the above ISA is probably not ideal.
    if( ! $self->location->isa("Bio::Location::SplitLocationI") ) {
	$self->throw("not atomic, not split, yikes, in trouble!");
    }

    my $seqstr;
    my $seqid = $self->entire_seq->display_id;
    # This is to deal with reverse strand features
    # so we are really sorting features 5' -> 3' on their strand
    # i.e. rev strand features will be sorted largest to smallest
    # as this how revcom CDSes seem to be annotated in genbank.
    # Might need to eventually allow this to be programable?    
    # (can I mention how much fun this is NOT! --jason)
    
    my ($mixed,$fstrand) = (0);
    if( $self->isa('Bio::Das::SegmentI') &&
	! $self->absolute ) { 
	$self->warn("Calling spliced_seq with a Bio::Das::SegmentI which does have absolute set to 1 -- be warned you may not be getting things on the correct strand");
    }
    
    my @locs = map { $_->[0] }
    # sort so that most negative is first basically to order
    # the features on the opposite strand 5'->3' on their strand
    # rather than they way most are input which is on the fwd strand

    sort { $a->[1] <=> $b->[1] } # Yes Tim, Schwartzian transformation
    map { 
	$fstrand = $_->strand unless defined $fstrand;
	$mixed = 1 if defined $_->strand && $fstrand != $_->strand;
	[ $_, $_->start* ($_->strand || 1)];	    
    } $self->location->each_Location; 
    
    if ( $mixed ) { 
	$self->warn("Mixed strand locations, spliced seq using the input order rather than trying to sort");    
	@locs = $self->location->each_Location; 
    }

    foreach my $loc ( @locs  ) {
	if( ! $loc->isa("Bio::Location::Atomic") ) {
	    $self->throw("Can only deal with one level deep locations");
	}
	my $called_seq;
	if( $fstrand != $loc->strand ) {
	    $self->warn("feature strand is different from location strand!");
	}
	# deal with remote sequences

	if( $loc->seq_id ne $seqid ) {
	    if( defined $db ) {
		my $sid = $loc->seq_id;
		$sid =~ s/\.\d+//g;
		eval {
		    $called_seq = $db->get_Seq_by_acc($sid);
		};
		if( $@ ) {
		    $self->warn("In attempting to join a remote location, sequence $sid was not in database. Will provide padding N's. Full exception \n\n$@");
		    $called_seq = undef;
		}
	    } else {
		$called_seq = undef;
	    }
	    if( !defined $called_seq ) {
		$seqstr .= 'N' x $self->length;
		next;
	    }
	} else {
	    $called_seq = $self->entire_seq;
	}
	
	if( $self->isa('Bio::Das::SegmentI') ) {
	    my ($s,$e) = ($loc->start,$loc->end);	    
	    $seqstr .= $called_seq->subseq($s,$e)->seq();
	} else { 
	    # This is dumb subseq should work on locations...
	    if( $loc->strand == 1 ) {
		$seqstr .= $called_seq->subseq($loc->start,$loc->end);
	    } else {
		$seqstr .= $called_seq->trunc($loc->start,$loc->end)->revcom->seq();
	    }
	}
    }
    my $out = Bio::Seq->new( -id => $self->entire_seq->display_id . "_spliced_feat",
				      -seq => $seqstr);
    
    return $out;
}

=head1 Bio::SeqFeature::CollectionI methods

List of methods inherited from Bio::SeqFeature::CollectionI (see
L<Bio::SeqFeature::CollectionI> for details).  To be sure, do go check
the pod for these.  There's unfortunately a shortcoming in pod that
prevents the inheritance of documentation for inherited methods.

=cut

# NOTE: We've copied the pod to here because the pod code is not clever
# enough to do it for us.  (Please kind developers, give us a super-pod
# capable of following inheritance relationships!)
#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 types (from Bio::SeqFeature::CollectionProviderI)

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

=head2 add_features

 Title   : add_features
 Usage   : $seq_feature->add_features( @feature_list );
 Function: Adds the given features to this Collection.
 Returns : The features added (or their count, in scalar context).
 Args    : An array of SeqFeatures or their ids
 Status  : Public

=head2 remove_features

 Title   : remove_features
 Usage   : $seq_feature->remove_features( @feature_list )
 Function: Removes the requested sequence features
 Returns : The removed features (or their count, in scalar context)
 Args    : An array of SeqFeatures or their ids
 Status  : Public

=head2 features

 Title   : features
 Usage   : @features = $collection->features( %args );
           OR
           @features = $collection->features( @types );
 Returns : a list of L<Bio::SeqFeatureI> objects,
           OR
           (when the -iterator option is true) an L<Bio::SeqFeature::IteratorI>
           OR
           (when the -callback argument is given) true iff the callbacks
             completed.
 Args    : see below
 Status  : Public

This routine will retrieve features associated with this collection
object.  It can be used to return all features, or a subset based on
their type, location, or attributes.

If ranges are specified using the -ranges argument, then these ranges
will be used to narrow the results, according to the specified
-rangetype and -strandtype arguments.

If a -baserange is specified or is provided by default* then
unqualified ranges given with the -ranges argument will be interpreted
as relative to that baserange.  Note that this only applies to
unqualified ranges, ie. ranges that have no defined seq_id.  You may
force absolute range interpretations by giving a -baserange that is
not a L<Bio::RangeI> (such as the string 'absolute') or by qualifying all
given ranges.

Footnote (*): All implementing classes that also implement L<Bio::RangeI>
              B<must> provide $self as the default baserange!

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
                 default -baserange.  If this CollectionI is also a
                 L<Bio::RangeI>, then the default -baserange should be
                 itself.  Note that the baserange affects the sort order.

  -range         A L<Bio::RangeI> object defining the range to search.
                 See also -strandmatch, -rangetype, and -baserange.
  -ranges        An array reference to multiple ranges.

  -rangetype     One of "overlaps", "contains", or "contained_in".

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
arguments consisting of the current feature and this CollectionI
object, and must return a true value. If the code returns a false
value, feature retrieval will be aborted.

-callback and -iterator are mutually exclusive options.  If -iterator
is defined, then -callback is ignored.

-callback and -sort are mutually exclusive options.  If -sort is
defined, then -callback is ignored.  If you want to do a sorted
callback, set the sorted() flag of this collection to true.

If -sort or sorted() is true then the features will be returned in
order of the features' start positions.  This order will be reversed
if the baserange has a negative strand (remember that a CollectionI
implementation that is also a L<Bio::RangeI> must provide itself as
the default baserange, but this may be overridden by the -baserange
argument).

Note that no guarantees are made by the CollectionI interface about
the order of the features, except when the sorted() flag is true or
when the -sort option is given to the features method.  Therefore
the implementation may choose to reorder the underlying data structure
to better accomodate -sorted feature requests as a result of a
features() call.  When this happens the CollectionI's sorted() flag
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

=head2 overlapping_features

 Title   : overlapping_features
 Usage   : @features = $collection->overlapping_features( %args )
 Function: get features that overlap the range of this collection
 Returns : a list of L<Bio::SeqFeatureI> objects,
           OR
           an iterator (when the -iterator option is true)
           OR
           (when the -callback option is given) true iff the callbacks
             completed.
 Args    : same as features()
 Status  : Public

This method is identical to features() except that it defaults to
finding overlapping features.

NOTE: This is defined in the interface in terms of features().  You do not
have to implement it.

=head2 contained_features

 Title   : contained_features
 Usage   : @features = $collection->contained_features( %args )
 Function: get features that are contained in the range of this collection
 Returns : a list of L<Bio::SeqFeatureI> objects,
           OR
           an iterator (when the -iterator option is true)
           OR
           (when the -callback option is given) true iff the callbacks
             completed.
 Args    : same as features()
 Status  : Public

This method is identical to features() except that it defaults to
a range type of 'contains'.

NOTE: This is defined in the interface in terms of features().  You do not
have to implement it.

=head2 contained_in

 Title   : contained_in
 Usage   : @features = $collection->contained_in( %args )
 Function: get features that contain the range of this collection
 Returns : a list of L<Bio::SeqFeatureI> objects,
           OR
           an iterator (when the -iterator option is true)
           OR
           (when the -callback option is given) true iff the callbacks
             completed.
 Args    : same as features()
 Status  : Public

This method is identical to features() except that it defaults to
a range type of 'contained_in'.

NOTE: This is defined in the interface in terms of features().  You do not
have to implement it.

=head2 get_feature_stream

 Title   : get_feature_stream
 Usage   : $iterator = $collection->get_feature_stream( %args )
 Function: get an iterator over the features in this collection
 Returns : a Bio::SeqFeature::IteratorI
 Args    : same as features()
 Status  : Public

This method is identical to features() except that it always generates
an iterator.

NOTE: This is defined in the interface in terms of features().  You do not
have to implement it.

=head2 features_in_range

 Title   : features_in_range
 Usage   : @features = $collection->features_in_range( $range );
             OR
           @features = $collection->features_in_range( %args );
 Function: Retrieves a list of features which were contained or overlap the
           the requested range
 Returns : a list of L<Bio::SeqFeatureI> objects,
           OR
           an iterator (when the -iterator option is true)
           OR
           (when the -callback option is given) true iff the callbacks
             completed.
 Args    : same as features(), or a single L<Bio::RangeI> argument
 Status  : Public

This method is identical to features() except that its first argument, if it is a RangeI object, will be used as the -range argument.

NOTE: This is defined in the interface in terms of features().  You do not
have to implement it.

=head2 get_feature_by_name

 Title   : get_feature_by_name
 Usage   : my @features = $collection->get_feature_by_name( $name )
           OR
           my @features = $collection->get_feature_by_name( $namespace, $name )
           OR
           my @features = $collection->get_feature_by_name( %args )
 Function: fetch features by their name
 Returns : a list of L<Bio::SeqFeatureI> objects,
           OR
           an iterator (when the -iterator option is true)
           OR
           (when the -callback option is given) true iff the callbacks
             completed.
 Args    : the string name or array ref of names
           OR
           the string namespace and the string name or array ref of names
           OR
           a hash, same as features()
 Status  : Public

This method is identical to features() except that it can take as
unnamed arguments the -name/-names value (one argument) OR the
-namespace value and the -name/-names value (two arguments).

Again, here's the deal:
  1) one argument: the argument is treated as the -name/-names value
  2) two arguments: the arguments are treated as the -namespace value and
       the -name/-names value, in that order.
     (note: this uses _rearrange() so the first argument must not
     begin with a hyphen or it will be interpreted as a named
     argument).
  3) an args hash as with features()

NOTE: This is defined in the interface in terms of features().  You do not
have to implement it.

=head2 get_feature_by_id

 Title   : get_feature_by_id
 Usage   : my @features = $collection->get_feature_by_id( $unique_id )
           OR
           my @features = $collection->get_feature_by_id( %args )
 Function: fetch features by their unique_ids
 Returns : a list of L<Bio::SeqFeatureI> objects,
           OR
           an iterator (when the -iterator option is true)
           OR
           (when the -callback option is given) true iff the callbacks
             completed.
 Args    : the string unique_id or array ref of unique_ids
           OR
           a list of string unique_ids
           OR
           a hash, same as features()
 Status  : Public

This method is identical to features() except that it can take as
an unnamed argument(s) the -unique_id/-unique_ids value(s).

NOTE: This is defined in the interface in terms of features().  You do not
have to implement it.

=head2 get_feature_by_attribute

 Title   : get_feature_by_attribute
 Usage   : my @features = $collection->get_feature_by_attribute( %attrs )
           OR
           my @features = $collection->get_feature_by_attribute( $attrs_ref, %args )
 Function: fetch features by their attributes
 Returns : a list of L<Bio::SeqFeatureI> objects,
           OR
           an iterator (when the -iterator option is true)
           OR
           (when the -callback option is given) true iff the callbacks
             completed.
 Args    : a hash, as would be passed to features() as -attributes => \%attrs
           OR
           a hash ref, as would be passed to features() as -attributes => $ref,
             then some other args to the features() method.
 Status  : Public

This method is identical to features() except that it assumes that the
given hash is the value meant for the -attributes argument, or (if the
first argument does not begin with '-') that the first argument is
meant as the -attributes value and the rest of them are usual
features() arguments.

NOTE: This is defined in the interface in terms of features().  You do not
have to implement it.

=head1 Bio::RangeI methods

List of methods inherited from Bio::RangeI (see L<Bio::RangeI> for
details).  To be sure, do go check the pod for these.  There's
unfortunately a shortcoming in pod that prevents the inheritance of
documentation for inherited methods.

=cut

# NOTE: We've copied the pod to here because the pod code is not clever
# enough to do it for us.  (Please kind developers, give us a super-pod
# capable of following inheritance relationships!)
#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 seq_id

  Title   : seq_id
  Usage   : my $seq_id = $range->seq_id( [new_seq_id] );
  Function: Get/Set a unique_id or primary_id of a L<Bio::PrimarySeqI>
            or another L<Bio::RangeI> that this RangeI is defined
            over or relative to.
  Returns : The current (or former, if used as a set method) value of
            the seq_id.
  Args    : [optional] A new (string or L<Bio::RangeI> seq_id value

  Ranges may have no defined seq_id, but this should be considered
  deprecated.  The concept of a 'range' requires that it is a range
  over some sequence; this method returns (and optionally sets) that
  sequence.  It is also possible to specify another range, to support
  relative ranges.  If the value of seq_id is another L<Bio::RangeI>,
  then this RangeI's positions are relative to that RangeI's
  positions.  If seq_id is the id of a sequence then it should provide
  enough information for a user of a RangeI to retrieve that sequence;
  ideally it should be a L<Bio::GloballyIdentifiableI> unique_id.

=head2 start

  Title   : start
  Usage   : $start = $range->start();
  Function: get/set the start of this range
  Returns : the start of this range
  Args    : optionaly allows the start to be set
            using $range->start($start)

=head2 end

  Title   : end
  Usage   : $end = $range->end();
  Function: get/set the end of this range
  Returns : the end of this range
  Args    : optionaly allows the end to be set
            using $range->end($end)

=head2 length

  Title   : length
  Usage   : $length = $range->length();
  Function: get/set the length of this range
  Returns : the length of this range
  Args    : optionaly allows the length to be set
             using $range->length($length)

=head2 strand

  Title   : strand
  Usage   : $strand = $range->strand();
  Function: get/set the strand of this range
  Returns : the strandidness (-1, 0, +1)
  Args    : optionaly allows the strand to be set
            using $range->strand($strand)

=head2 overlaps

  Title   : overlaps
  Usage   : if($r1->overlaps($r2)) { do stuff }
  Function: tests if $r2 overlaps $r1
  Args    : arg #1 = a range to compare this one to (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') (optional)
  Returns : true if the ranges overlap, false otherwise

=head2 contains

  Title   : contains
  Usage   : if($r1->contains($r2) { do stuff }
  Function: tests whether $r1 totally contains $r2 
  Args    : arg #1 = a range to compare this one to (mandatory)
	             alternatively, integer scalar to test
            arg #2 = strand option ('strong', 'weak', 'ignore') (optional)
  Returns : true if the argument is totaly contained within this range

=head2 intersection

  Title   : intersection
  Usage   : my $intersection_range = $r1->intersection( $r2 ) (scalar context)
            OR
            my ( $start, $end, $strand ) = $r1->intersection( $r2 )
             (list context)
  Function: gives the range that is contained by both ranges
  Args    : arg #1 = a range to compare this one to (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') (optional)
  Returns : undef if they do not overlap,
            or new range object containing the overlap
            or (in list context) the start, end, and strand of that range.

=head2 union

  Title   : union
  Usage   : my $union_range = $r1->union( @other_ranges ); (scalar context)
            OR
            my ( $start, $end, $strand ) = $r1->union( @other_ranges );
              (list context)
            OR
            my $union_range = Bio::RangeI->union( @ranges );
              (scalar context)
            OR
            my ( $start, $end, $strand ) = Bio::RangeI->union( @ranges );
              (list context)
  Function: finds the minimal range that contains all of the ranges
  Args    : a range or list of ranges to find the union of
  Returns : a new range object that contains all of the given ranges, or
            (in list context) the start, end, and strand of that range object.

=head2 overlap_extent

 Title   : overlap_extent
 Usage   : ($a_unique,$common,$b_unique) = $a->overlap_extent($b)
 Function: Provides actual amount of overlap between two different
           ranges.
 Example :
 Returns : array of values for 
           - the amount unique to a
           - the amount common to both
           - the amount unique to b
 Args    : a range


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
   * Return false unless they have the same range (using the RangeI
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

  if( $other->isa( 'Bio::SeqFeatureI' ) ) {
    my $my_unique_id = $self->unique_id();
    my $other_unique_id = $self->unique_id();
    if( defined( $my_unique_id ) ) {
      if( defined( $other_unique_id ) ) {
        return ( $my_unique_id eq $other_unique_id );
      }
      return 0;
    }
    if( defined( $other_unique_id ) ) {
      return 0;
    }
    ## Okay, neither has a unique_id.
    unless( $self->start() == $other->start() and
	    $self->end()   == $other->end()       ) {
      return 0;
    }
    my $my_display_name = $self->display_name();
    my $other_display_name = $self->display_name();
    if( defined( $my_display_name ) ) {
      if( defined( $other_display_name ) ) {
        return ( $my_display_name eq $other_display_name );
      }
      return 0;
    }
    if( defined( $other_display_name ) ) {
      return 0;
    }
    ## Okay, neither has a display_name
    return 1;
  } else {
    return 0;
  }
} # equals(..)

=head2 toString

 Title   : toString
 Usage   : $str_val = $feature->toString()
 Function: returns $self->unique_id() || $overload->StrVal( $self )
 Returns : a String
 Args    : None
 Status  : Public

  This method is a hack.

=cut

sub toString {
  my $self = shift;

  return $self->unique_id() || overload::StrVal( $self );
} # toString()

## method for overload for comparing two SeqFeature objects.  Uses toString().
sub _cmp {
  my $self = shift;
  my ( $b, $reversed ) = @_;
  my $a = $self->toString();
  ( $a, $b ) = ( $b, $a ) if $reversed;
  return ( $a cmp $b );
}

1;

__END__
