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

use Bio::SeqFeature::SegmentI;
@ISA = qw( Bio::SeqFeature::SegmentI );

use Bio::Seq;
use Carp;

=head1 Identification methods

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
  shift->throw_not_implemented( @_ );
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
  shift->throw_not_implemented( @_ );
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
  my $current_string_value;
  if( ref( $type ) && $type->isa( 'Bio::SeqFeature::TypeI' ) ) {
    $current_string_value = $type->toString();
  } else {
    $current_string_value = "$type";
  }
  if( @_ ) {
    $self->type( @_ );
  }
  return $current_string_value;
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
  shift->throw_not_implemented( @_ );
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
  shift->throw_not_implemented( @_ );
}

=head2 get_tag_values

 Title   : get_tag_values
 Usage   : @values = $self->get_tag_values( 'some_tag' );
 Function: 
 Returns : An array comprising the values of the specified tag.
 Args    :


=cut

sub get_tag_values {
  shift->throw_not_implemented( @_ );
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
  shift->throw_not_implemented( @_ );
}

=head2 get_SeqFeatures

 Title   : get_SeqFeatures
 Usage   : @feats = $feat->get_SeqFeatures();
 Function: Returns an array of sub Sequence Features
 Returns : An array
 Args    : none


=cut

sub get_SeqFeatures {
  shift->throw_not_implemented( @_ );
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
  shift->throw_not_implemented( @_ );
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
  shift->throw_not_implemented( @_ );
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
  shift->throw_not_implemented( @_ );
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

=head1 SeqFeatureI-unique variations on Bio::SeqFeature::CollectionI methods

=head2 add_SeqFeatures

 Title   : add_SeqFeatures
 Usage   : my @added_features = $feature->add_SeqFeatures( $sub_feature );
           OR
           $feature->add_SeqFeatures( $sub_feature, $type );
           OR
           $feature->add_SeqFeatures( $sub_feature, $expand, $type );
 Function: Adds subfeatures to this feature, optionally with the given
           type, optionally expanding this feature to accommodate them.
 Returns : nothing
 Args    : 1 or 2 or 3 arguments: the first is a L<Bio::SeqFeatureI>
           object or a reference to a list of them; the second if
           there's 3 arguments is a boolean EXPAND value; the third if
           there's 3 or second if there's 2 arguments is a
           L<Bio::SeqFeature::TypeI> object or a type id string.

  Note that the $sub_feature argument may be a reference to a list of them.
  If a type is given then before the feature is added its type(..)
  method will be called with the given value as an argument.  If a
  true value for $expand is given (when there's 3 arguments) then this
  feature will grow to accomodate the new sub_feature.  Note that if
  only 2 arguments are given but the second argument is the string
  'EXPAND', the value will be interpreted as a true $expand argument,
  and the type of the feature will not be changed.

=cut

sub add_SeqFeatures {
  my $self = shift;
  my @sub_features = ( shift );

  if( ref( $sub_features[ 0 ] ) eq 'ARRAY' ) {
    @sub_features = @{ $sub_features[ 0 ] };
  }
  my $type = shift;
  my $expand;
  if( $_[ 0 ] ) {
    $expand = $type;
    $type = shift;
  } elsif( $type eq 'EXPAND' ) {
    $expand = $type;
    undef $type;
  }

  # Add them, and replace the list with those actually added.
  @sub_features = $self->add_features( @sub_features );

  ## TODO: REMOVE
  #print STDERR "add_SeqFeatures(..): \@sub_features are ( ", join( ", ", @sub_features ), " ), \$expand is $expand, \$type is $type\n";

  if( defined( $type ) ) {
    foreach my $sub_feature ( @sub_features ) {
      $sub_feature->type( $type );
    } # End foreach $sub_feature, 
  } # End if a $type was given.

  # Adjust for the features if requested
  if( $expand ) {
    ## TODO: REMOVE
    #warn "Before expanding, I am ".$self.", starting at ".$self->start();
    # Adjust the bounds, don't shrink.
    $self->adjust_bounds( 0, @sub_features );
    ## TODO: REMOVE
    #warn "After expanding, I am ".$self.", starting at ".$self->start();
  }
  return ( wantarray ? @sub_features : $sub_features[ 0 ] );
} # add_SeqFeatures(..)

=head2 remove_SeqFeatures

 Title   : remove_SeqFeatures
 Usage   : $feature->remove_SeqFeatures();
           OR
           $feature->remove_SeqFeatures( @types );
 Function: Removes all sub features or all sub features of the given types
 Returns : The array of L<Bio::SeqFeatureI> objects removed.
 Args    : none

=cut

sub remove_SeqFeatures {
  my $self = shift;
  return $self->remove_features( $self->features( @_ ) );
} # remove_SeqFeatures(..)

=head1 Bio::SeqFeature::SegmentI methods

List of methods inherited from Bio::SeqFeature::SegmentI (see
L<Bio::SeqFeature::SegmentI> for details).  To be sure, do go check
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
 Usage   : my @types = $seq_feature->types();
           OR
           my %types_and_counts = $seq_feature->types( -count => 1 );
 Function: Enumerate the feature types of all contained features, and possibly
           count the features in each type.
 Returns : a list of L<Bio::SeqFeature::TypeI> objects
           OR
           a hash mapping type id strings to integer counts
 Args    : see below

This routine returns a list of the types of sub_features contained in
this feature.  If the -count argument is given, it returns a hash of
types mapped to their occurrence counts.  Note that the hierarchy of
TypeI objects is ignored, so if there are 4 features of type 'foo'
which is a child of type 'bar', and only 1 feature (explicitly) of
type 'bar', then the count for 'bar' will be 1, not 5.

Arguments are -option=E<gt>value pairs as follows:

  -count aka -enumerate  if true, count the features for each type

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
 Usage   : @features = $seq_feature->features( %args );
           OR
           @features = $seq_feature->features( @types );
 Returns : a list of L<Bio::SeqFeatureI> objects,
           OR
           (when the -iterator option is true) an L<Bio::SeqFeature::IteratorI>
           OR
           (when the -callback argument is given) true iff the callbacks
             completed.
 Args    : see below
 Status  : Public

This routine will retrieve features contained within this feature
object.  It can be used to return all features, or a subset based on
their type, location, or attributes.  Features that are returned in
relative mode (relative either to this SeqFeatureI or to a given RangeI)
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
( -absolute || absolute() ) is false then the default is this SeqFeatureI
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
                 if this SeqFeatureI's absolute() flag was set to true,
                 even if it isn't.  If -absolute is given and false
                 then all behavior will be as if this SeqFeatureI's
                 absolute() flag was set to false, even if it isn't.
                 Note that -baserange can still be given and can force
                 relativeness, and that takes precedence over -absolute.

  -range         A L<Bio::RangeI> object defining the range to search.
                 See also -strandmatch, -rangetype, and -baserange.
  -ranges        An array reference to multiple ranges.

  -rangetype     One of "overlaps", "contains", or "contained_in".  If no
                 range is given then a strange thing happens (it is
                 described above).

  -strandmatch   One of "strong", "weak", or "ignore".  Note that the
                 strand attribute of a given -range must be non-zero
                 for this to work (a 0/undef strand forces a 'weak'
                 strandmatch to become 'ignore' and cripples the
                 'strong' strandmatch).

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
arguments consisting of the current feature and this SeqFeatureI
object, and must return a true value. If the code returns a false
value, feature retrieval will be aborted.

-callback and -iterator are mutually exclusive options.  If -iterator
is defined, then -callback is ignored.

-callback and -sort are mutually exclusive options.  If -sort is
defined, then -callback is ignored.  If you want to do a sorted
callback, set the sorted() flag of this feature to true.

If -sort or sorted() is true then the features will be returned in
order of the features' start positions.  This order will be reversed
if the baserange has a negative strand (remember that the default
baserange depends upon the value of the absolute() flag, but this may
be overridden by the -baserange argument).

Note that no guarantees are made by the SeqFeatureI interface about
the order of the features, except when the sorted() flag is true or
when the -sort option is given to the features method.  Therefore
the implementation may choose to reorder the underlying data structure
to better accomodate -sorted feature requests as a result of a
features() call.  When this happens the SeqFeatureI's sorted() flag
should be set to true, so that the client can detect that the -sorted
argument to features() is now irrelevant.

=head2 overlapping_features

 Title   : overlapping_features
 Usage   : @features = $seq_feature->overlapping_features( %args )
 Function: get features that overlap the range of this feature
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
 Usage   : @features = $seq_feature->contained_features( %args )
 Function: get features that are contained in the range of this feature
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
 Usage   : @features = $seq_feature->contained_in( %args )
 Function: get features that contain the range of this feature
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
 Usage   : $iterator = $seq_feature->get_feature_stream( %args )
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
 Usage   : @features = $seq_feature->features_in_range( $range );
             OR
           @features = $seq_feature->features_in_range( %args );
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
 Usage   : my @features = $seq_feature->get_feature_by_name( $name )
           OR
           my @features = $seq_feature->get_feature_by_name( $namespace, $name )
           OR
           my @features = $seq_feature->get_feature_by_name( %args )
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
 Usage   : my @features = $seq_feature->get_feature_by_id( $unique_id )
           OR
           my @features = $seq_feature->get_feature_by_id( %args )
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
 Usage   : my @features = $seq_feature->get_feature_by_attribute( %attrs )
           OR
           my @features = $seq_feature->get_feature_by_attribute( $attrs_ref, %args )
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

=head2 get_collection

 Title   : get_collection
 Usage   : my $segment = $seq_feature->get_collection( %args );
           OR
           my $segment = $seq_feature->get_collection( @types );
 Returns : A L<Bio::SeqFeature::SegmentI> object
 Args    : see below
 Status  : Public

This routine will retrieve a L<Bio::SeqFeature::SegmentI> object based
on feature type, location or attributes.  The SeqFeatureI objects in
the returned SegmentI may or may not be newly instantiated by this
request.  They will have as their range the range searched, if any, or
the smallest range that encloses the returned features.  They will
have as their seq_id() the -baserange used here (if the baserange is
absolute by any means then their seq_id() will be this SeqFeatureI's
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
                 if this SeqFeatureI's absolute() flag was set to true,
                 even if it isn't.  If -absolute is given and false
                 then all behavior will be as if this SeqFeatureI's
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

=head2 parent_segment_provider

 Title   : parent_segment_provider
 Usage   : my $parent = $seq_feature->parent_segment_provider();
 Function: Return the SegmentProviderI that is the parent of this feature.
 Returns : a L<Bio::DB::SegmentProviderI> or undef if there is none
 Args    : none

=head2 absolute

  Title   : absolute
  Usage   : my $absolute_flag = $seq_feature->absolute( [$new_absolute_flag] );
  Function: Get/set the absolute flag.
  Returns : The current (or former, if used as a set method) value of the
            absolute flag.
  Args    : [optional] a new value for the absolute flag.

  If the absolute() flag is set then the start(), end(), and strand()
  methods will behave like the abs_start(), abs_end(), and abs_strand()
  methods, meaning that they will return values relative to abs_seq_id()
  rather than to seq_id().

=head2 abs_seq_id

  Title   : abs_seq_id
  Usage   : my $abs_seq_id = $seq_feature->abs_seq_id();
  Function: Get the unique_id or primary_id of the L<Bio::PrimarySeqI>
            that this SeqFeatureI is defined over.
  Returns : The root seq_id, or undef if there is none.
  Args    : none

  Features may have no defined abs_seq_id, but this should be considered
  deprecated.  The concept of a feature requires that it is a range
  over some sequence; this method returns that sequence.  If the value
  of seq_id() is a string (the unique_id or primary_id of a
  L<Bio::PrimarySeqI>) then this method will be identical to seq_id().
  If the value of seq_id() is another L<Bio::RangeI>, then this method
  will return its seq_id() if that is a string, or keep searching up the
  tree until a string (or undef) is reached.

=head2 abs_start

  Title   : abs_start
  Usage   : my $abs_start = $seq_feature->abs_start();
  Function: Get the absolute start position of this feature.
  Returns : The current start position of this feature, relative to the
            abs_seq_id.
  Args    : none

  Note the interdependence of abs_start() and start().  Changing start() will
  change abs_start().

  Note the interdependence of abs_start() and length().  Changing length() will
  change abs_start().

=head2 abs_end

  Title   : abs_end
  Usage   : my $abs_end = $seq_feature->abs_end();
  Function: Get the absolute end position of this feature.
  Returns : The current absolute end position of this feature, relative
            to the abs_seq_id.
  Args    : none

  Note the interdependence of abs_end() and end().  Changing end() will
  change abs_end().

  Note the interdependence of abs_end() and length().  Changing length() will
  change abs_end().

=head2 abs_strand

  Title   : abs_strand
  Usage   : my $abs_strand = $seq_id->abs_strand();
  Function: Get the absolute strandedness (-1, 0, or 1) of this feature.
  Returns : The current absolute strand value of this feature.
  Args    : none

=head2 abs_low

  Title   : abs_low
  Usage   : my $abs_low = $seq_id->abs_low();
  Function: Get the least-valued absolute position of this feature.
  Returns : The current lowest position of this feature, relative to the
            abs_seq_id.
  Args    : none

  This will return either abs_start() or abs_end(), depending on which
  is lower.

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

=head2 abs_high

  Title   : abs_high
  Usage   : my $abs_high = $seq_id->abs_high();
  Function: Get the greatest-valued absolute position of this feature.
  Returns : The current highest position of this feature, relative to the
            abs_seq_id.
  Args    : none

  This will return either abs_start() or abs_end(), depending on which
  is higher.

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

=head2 rel2abs

  Title   : rel2abs
  Usage   : my @abs_coords = $seq_id->rel2abs( @rel_coords );
  Function: Convert relative coordinates into absolute coordinates
  Returns : a list of absolute coordinates
  Args    : a list of relative coordinates

  This function takes a list of positions in relative coordinates
  (relative to seq_id()), and converts them into absolute coordinates.

  Note that if absolute() is true this method still interprets
  incoming coordinates as if they were relative to what seq_id() would
  be if absolute() were false.

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.  Note that this implementation
  uses abs_start() and abs_strand(), so these methods should not be
  defined in terms of rel2abs(), lest a vicious cycle occur.

=head2 abs2rel

  Title   : abs2rel
  Usage   : my @rel_coords = $seq_id->abs2rel( @abs_coords )
  Function: Convert absolute coordinates into relative coordinates
  Returns : a list of relative coordinates
  Args    : a list of absolute coordinates

  This function takes a list of positions in absolute coordinates
  and converts them into relative coordinates (relative to seq_id()).

  Note that if absolute() is true this method still produces
  coordinates relative to what seq_id() would be if absolute() were
  false.

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.  Note that this implementation
  uses abs_start() and abs_strand(), so these methods should not be
  defined in terms of abs2rel(), lest a vicious cycle occur.

=head2 rel2abs_strand

  Title   : rel2abs_strand
  Usage   : my $abs_strand = $seq_id->rel2abs_strand( $rel_strand );
  Function: Convert a strand that is relative to seq_id() into one that
            is relative to abs_seq_id().
  Returns : a strand value (-1, 0, or 1).
  Args    : a strand value (-1, 0, or 1).

  This function takes a strand value that is relative to seq_id()
  and converts it so that it is absolute (ie. relative to abs_seq_id()).

  Note that if absolute() is true this method still interprets
  the argument strand as it were relative to what seq_id() would
  be if absolute() were false.

=head2 abs2rel_strand

  Title   : abs2rel_strand
  Usage   : my $rel_strand = $seq_id->abs2rel_strand( $abs_strand )
  Function: Convert a strand that is relative to abs_seq_id() into one that
            is relative to seq_id().
  Returns : a strand value (-1, 0, or 1).
  Args    : a strand value (-1, 0, or 1).

  This function takes a strand value that is absolute (ie. relative to
  abs_seq_id()) and converts it so that it is relative to seq_id().

  Note that if absolute() is true this method still returns the strand
  relative to what seq_id() would be if absolute() were false.

  This method turns out to be identical to rel2abs_strand, so it is
  implemented in the interface as a (glob ref) alias for
  rel2abs_strand.

=head2 seq_id

  Title   : seq_id
  Usage   : my $seq_id = $seq_id->seq_id( [new_seq_id] );
  Function: Get/Set a unique_id or primary_id of a L<Bio::PrimarySeqI>
            or another L<Bio::RangeI> that this SeqFeatureI is defined
            over or relative to.  If absolute() is true, this will be
            identical to abs_seq_id().
  Returns : The current (or former, if used as a set method) value of
            the seq_id.
  Args    : [optional] A new (string or L<Bio::RangeI> seq_id value

  Features may have no defined seq_id, but this should be considered
  deprecated.  The concept of a 'feature' requires that it is a range
  over some sequence; this method returns (and optionally sets) that
  sequence.  It is also possible to specify another range, to support
  relative ranges.  If the value of seq_id is another L<Bio::RangeI>,
  then this SeqFeatureI's positions are relative to that RangeI's
  positions (unless absolute() is true, in which case they are
  relative to the root seq_id).  If seq_id is the id of a sequence then
  it should provide enough information for a user of a SeqFeatureI to
  retrieve that sequence; ideally it should be a
  L<Bio::GloballyIdentifiableI> unique_id.

  You may not set the seq_id when absolute() is true.

=head2 start

  Title   : start
  Usage   : my $start = $seq_id->start( [$new_start] );
  Function: Get/set the start of this feature.
  Returns : The current (or former, if used as a set method) start position
            of this feature.  If absolute() is true then this value will
            be relative to the abs_seq_id; otherwise it will be
            relative to the seq_id.
  Args    : [optional] a new start position

  Note the interdependence of start() and abs_start().  Changing start() will
  change abs_start().

  You may not set start() when absolute() is true.

  Note the interdependence of start() and length().  Changing start() will
  change length().

=head2 end

  Title   : end
  Usage   : my $end = $seq_id->end( [$new_end] );
  Function: Get/set the end of this feature.
  Returns : The current (or former, if used as a set method) end position
            of this feature.  If absolute() is true then this value will
            be relative to the abs_seq_id; otherwise it will be
            relative to the seq_id.
  Args    : [optional] a new end position

  Note the interdependence of end() and abs_end().  Changing end() will
  change abs_end().

  You may not set end() when absolute() is true.

  Note the interdependence of end() and length().  Changing one will
  change the other.

=head2 strand

  Title   : strand
  Usage   : my $strand = $seq_id->strand( [$new_strand] );
  Function: Get/set the strandedness (-1, 0, or 1) of this feature.
  Returns : The current (or former, if used as a set method) strand value
            of this feature.  If absolute() is true then this value will
            be absolute.  Otherwise it will be relative to the
            strandedness (if any) of seq_id.
  Args    : [optional] a new strand value.

  You may not set strand() when absolute() is true.

=head2 length

  Title   : length
  Usage   : my $length = $seq_id->length( [$new_length] );
  Function: Get/set the length of this feature.
  Returns : The current (or former, if used as a set method) length
            of this feature.
  Args    : [optional] a new length

  length = ( ( end - start ) + 1 ) = ( ( abs_high - abs_low ) + 1 ).

  Note the interdependence of start()|end()|abs_start()|abs_end() and
  length().  Changing start() or end() will change the length.
  Changing the length will change the end() (and consequently abs_end()).

  You may not set the length when absolute() is true.

=head2 overlaps

  Title   : overlaps
  Usage   : if( $r1->overlaps( $r2 ) ) { do stuff }
  Function: tests if $r2 overlaps $r1
  Args    : arg #1 = a L<Bio::RangeI> to compare this one to (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') (optional)
  Returns : true if the ranges overlap, false otherwise

  The second argument's values may be:
   "strong"        ranges must have the same strand
   "weak"          ranges must have the same strand or no strand
   "ignore"        ignore strand information (default)

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

=head2 contains

  Title   : contains
  Usage   : if( $r1->contains( $r2 ) ) { do stuff }
  Function: tests if $r2 is totally contained within $r1
  Args    : arg #1 = a L<Bio::RangeI> to compare this one to,
                     or an integer position (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') (optional)
  Returns : true iff this range wholly contains the given range

  The second argument's values may be:
   "strong"        ranges must have the same strand
   "weak"          ranges must have the same strand or no strand
   "ignore"        ignore strand information (default)

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

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

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

=head2 union

  Title   : union
  Usage   : my $union_range = $r1->union( @other_ranges ); (scalar context)
            OR
            my ( $start, $end, $strand ) = $r1->union( @other_ranges );
              (list context)
            OR
            my $union_range = Bio::RelRangeI->union( @ranges );
              (scalar context)
            OR
            my ( $start, $end, $strand ) = Bio::RelRangeI->union( @ranges );
              (list context)
  Function: finds the minimal range that contains all of the ranges
  Args    : a range or list of ranges to find the union of
  Returns : a new range object that contains all of the given ranges, or
            (in list context) the start, end, and strand of that range object.

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

=head2 overlap_extent

 Title   : overlap_extent
 Usage   : my ( $a_unique, $common, $b_unique ) = $a->overlap_extent( $b );
 Function: Provides actual amount of overlap between two different ranges.
 Returns : 3-tuple consisting of:
           - the number of positions unique to a
           - the number of positions common to both
           - the number of positions unique to b
 Args    : a L<Bio::RangeI> object

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

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
    unless( $self->Bio::RelRangeI::equals( $other ) ) {
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
 Function: returns $self->unique_id() || $self->display_name() ||
           overload::StrVal( $self )
 Returns : a String
 Args    : None
 Status  : Public

  This method is a hack.

=cut

sub toString {
  my $self = shift;

  return $self->unique_id() || $self->display_name() ||
         overload::StrVal( $self );
} # toString()

## method for overload for comparing two SeqFeature objects.
sub _cmp {
  my $self = shift;
  my ( $b, $reversed ) = @_;
  my $a = ( $self->unique_id() ||
            $self->display_name() ||
          overload::StrVal( $self ) );
  ( $a, $b ) = ( $b, $a ) if $reversed;
  return ( $a cmp $b );
}

1;

__END__
