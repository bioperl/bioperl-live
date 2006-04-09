package Bio::DB::SeqFeature::Store;

# $Id$

=head1 NAME

Bio::DB::SeqFeature::Store -- Storage and retrieval of sequence annotation data

=head1 SYNOPSIS

  use Bio::DB::SeqFeature::Store

  # Open the sequence database
  my $db      = Bio::DB::SeqFeature::Store->new( -adaptor => 'DBI::mysq',
                                                 -dsn     => 'dbi:mysql:elegans');

  # get a feature from somewhere
  my $feature = Bio::SeqFeature::Generic->new(...);

  # store it
  $db->store($feature) or die "Couldn't store!";

  # primary ID of the feature is changed to indicate its primary ID
  # in the database...
  my $id = $feature->primary_id;

  # get the feature back out
  my $f  = $db->fetch($id);

  # change the feature and update it
  $f->start(100);
  $db->update($f) or die "Couldn't update!";

  # searching...
  # ...by id
  my @features = $db->fetch_many(@list_of_ids);

  # ...by name
  @features = $db->get_features_by_name('ZK909');

  # ...by alias
  @features = $db->get_features_by_alias('sma-3');

  # ...by type
  @features = $db->get_features_by_name('gene');

  # ...by location
  @features = $db->get_features_by_location(-seq_id=>'Chr1',-start=>4000,-end=>600000);

  # ...by attribute
  @features = $db->get_features_by_attribute({description => 'protein kinase'})

  # ...by the GFF "note" field
  @features = $db->search_notes('kinase');

  # ...by arbitrary combinations of selectors
  @features = $db->features(-name => $name,
                            -type => $types,
                            -seq_id => $seqid,
                            -start  => $start,
                            -end    => $end,
                            -attributes => $attributes);

  # ...using an iterator
  my $iterator = $db->get_seq_stream(-name => $name,
                                     -type => $types,
                                     -seq_id => $seqid,
                                     -start  => $start,
                                     -end    => $end,
                                     -attributes => $attributes);

  while (my $feature = $iterator->next_seq) {
    # do something with the feature
  }

  # ...limiting the search to a particular region
  my $segment  = $db->segment('Chr1',5000=>6000);
  my @features = $db->features(-type=>['mRNA','match']);

  # getting & storing sequence information
  # warning -- this is a string not a PrimarySeq object
  $db->insert_sequence('Chr1','GATCCCCCGGGATTCCAAAA...');
  my $sequence = $db->fetch_sequence('Chr1',5000=>6000);

=head1 DESCRIPTION

Bio::DB::SeqFeature::Store implements the Bio::SeqFeature::CollectionI
interface to allow you to persistently store Bio::SeqFeatureI objects
in a database and to later to retrieve them by a variety of
searches. This module is similar to the older Bio::DB::GFF module,
with the following differences:

=over 4

=item 1. No limitation on Bio::SeqFeatureI implementations

Unlike Bio::DB::GFF, Bio::DB::SeqFeature::Store works with
any Bio::SeqFeatureI object.

=item 2. No limitation on nesting of features & subfeatures

Bio::DB::GFF is limited to features that have at most one
level of subfeature. Bio::DB::SeqFeature::Store can work with features
that have unlimited levels of nesting.

=item 3. No aggregators

The aggregator architecture, which was necessary to impose order on
the GFF2 files that Bio::DB::GFF works with, does not apply to
Bio::DB::SeqFeature::Store. It is intended to store features that obey
well-defined ontologies, such as the Sequence Ontology
(http://song.sourceforge.net).

=item 4. No relative locations

All locations defined by this module are relative to an absolute
sequence ID, unlike Bio::DB::GFF which allows you to define the
location of one feature relative to another.

=back

We'll discuss major concepts in Bio::DB::SeqFeature::Store and then
describe how to use the module.

=head2 Adaptors

Bio::DB::SeqFeature::Store is designed to work with a variety of
storage back ends called "adaptors." Adaptors are specialized
subclasses of Bio::DB::SeqFeature::Store and provide the interface
between the store() and fetch() methods and the physical
database. Currently the number of adaptors is quite limited, but the
number will grow soon.

=over 4

=item DBI::mysql

A full-featured implementation on top of the MySQL relational database
system.

=item bdb

A partial implementation that runs on top of the BerkeleyDB
database. The fetch() and store() methods are implemented, but the
various search functions (e.g. get_features_by_name()) are not.

=back

If you do not explicitly specify the adaptor, then DBI::mysql will be
used by default.

=head2 Serializers

When Bio::DB::SeqFeature::Store stores a Bio::SeqFeatureI object into
the database, it serializes it into binary or text form. When it later
fetches the feature from the database, it unserializes it. Two
serializers are available: Recent versions of 

=over 4

=item Storable

This is a fast binary serializer. It is available in Perl versions
5.8.7 and higher and is used when available.

=item Data::Dumper

This is a slow text serializer that is available in Perl 5.8.0 and
higher. It is used when Storable is unavailable.

=back

If you do not specify the serializer, then Storable will be used if
available; otherwise Data::Dumper.

=head2 Loaders and Lazy Features

The Bio::DB::SeqFeature::GFF3Loader parses a GFF3-format file and
loads the annotations and sequence data into the database of your
choice. The script bp_seqfeature_load.pl (found in the
scripts/Bio-SeqFeature-Store/ subdirectory) is a thin front end to the
GFF3Loader. Other loaders may be written later.

Although Bio::DB::SeqFeature::Store should work with any
Bio::SeqFeatureI object, there are some disadvantages to using
Bio::SeqFeature::Generic and other vanilla implementations. The major
issue is that if two vanilla features share the same subfeature
(e.g. two transcripts sharing an exon), the shared subfeature will be
cloned when stored into the database.

The special-purpose L<Bio::DB::SeqFeature::LazyTableFeature> class is
able to normalize its subfeatures in the database, so that shared
subfeatures are stored only once. This minimizes wasted storage
space. In addition, when combined with the
L<Bio::DB::SeqFeature::Store::Cacher> module, each shared subfeature
will usually occupy only a single memory location upon restoration.

=cut


use strict;

use base 'Bio::SeqFeature::CollectionI';
use Carp 'croak';
use Bio::DB::GFF::Util::Rearrange;
use Bio::DB::SeqFeature::Segment;
use Scalar::Util 'blessed';

*dna = *get_dna = *get_sequence = \&fetch_sequence;
*get_SeqFeatures = \&fetch_SeqFeatures;

=head1 Methods for Connecting and Initializating a Database

=head2 new

 Title   : new
 Usage   : $db = $db->new(@options)
 Function: connect to a database
 Returns : A descendent of Bio::DB::Seqfeature::Store
 Args    : several - see below
 Status  : public

This class method creates a new database connection. The following
-name=>$value arguments are accepted:http://iowg.brcdevel.org/gff3.html#a_fasta

 Name               Value
 ----               -----

 -adaptor           The name of the Adaptor class (default DBI::mysql)

 -serializer        The name of the serializer class (default Storable)

 -index_subfeatures Whether or not to make subfeatures searchable
                    (default true)
      
The B<-index_subfeatures> argument, if true, tells the module to
create indexes for a feature and all its subfeatures (and its
subfeatues' subfeatures). Indexing subfeatures means that you will be
able to search for the gene, its mRNA subfeatures and the exons inside
each mRNA. It also means when you search the database for all features
contained within a particular location, you will get the gene, the
mRNAs and all the exons as individual objects as well as subfeatures
of each other.

The new() method of individual adaptors recognize additional
arguments. The default DBI::mysql adaptor recognizes the following
ones:

 Name               Value
 ----               -----

 -dsn               DBI data source (default dbi:mysql:test)

 -autoindex         A flag that controls whether or not to update
                    all search indexes whenever a feature is stored
                    or updated (default true).

 -namespace         A string that will be used to qualify each table,
                    thereby allowing you to store several independent
                    sequence feature databases in a single Mysql
                    database.

 -dumpdir           The path to a temporary directory that will be
                    used during "fast" loading. See
		    L<Bio::DB::SeqFeature::GFF3Loader> for a
		    description of this. Default is the current
                    directory.

=cut

### 
# object constructor
#
sub new {
  my $self      = shift;
  my ($adaptor,$serializer,$index_subfeatures,$args);
  if (@_ == 1) {
    $args = {DSN => shift}
  }
  else {
    ($adaptor,$serializer,$index_subfeatures,$args) =
      rearrange(['ADAPTOR',
		 'SERIALIZER',
		 'INDEX_SUBFEATURES'
		],@_);
  }
  $adaptor ||= 'DBI::mysql';

  my $class = "Bio::DB::SeqFeature::Store::$adaptor";
  eval "require $class " or croak $@;
  my $obj = $class->new_instance();
  $obj->init($args);
  $obj->serializer($serializer)               if defined $serializer;
  $obj->index_subfeatures($index_subfeatures) if defined $index_subfeatures;
  $obj;
}

=head2 init_database

 Title   : init_database
 Usage   : $db->init_database([$erase_flag])
 Function: initialize a database
 Returns : true
 Args    : (optional) flag to erase current data
 Status  : public

Call this after Bio::DB::SeqFeature::Store->new() to initialize a new
database. In the case of a DBI database, this method installs the
schema but does B<not> create the database. You have to do this
offline using the appropriate command-line tool. In the case of the
"bdb" BerkeleyDB adaptor, this creates an empty BTREE database.

If there is any data already in the database, init_database() called
with no arguments will have no effect. To permanently erase the data
already there and prepare to receive a fresh set of data, pass a true
argument.

=cut

###
# wipe database clean and reinstall schema
#
sub init_database {
  my $self = shift;
  $self->_init_database(@_);
}


=head2 store

 Title   : store
 Usage   : $success = $db->store(@features)
 Function: store one or more features into the database
 Returns : true if successful
 Args    : list of Bio::SeqFeatureI objects
 Status  : public

This method stores a list of features into the database. Each feature
is updated so that its primary_id becomes the primary ID of the
serialized feature stored in the database. If all features were
successfully stored, the method returns true. In the DBI
implementation, the store is performed as a single transaction and the
transaction is rolled back if one or more store operations failed.

You can find out what the primary ID of the feature has become by
calling the feature's primary_id() method:

  $db->store($my_feature) or die "Oh darn";
  my $id = $my_feature->primary_id;

If the feature contains subfeatures, they will all be stored
recursively. In the case of Bio::DB::SeqFeature::LazyFeature and
Bio::DB::SeqFeature::LazyTableFeature, the subfeatures will be stored
in a normalized way so that each subfeature appears just once in the
database.

Subfeatures will be indexed for separate retrieval based on the
current value of index_subfeatures().

=cut

###
# store one or more Bio::SeqFeatureI objects
#      if they already have a primary_id will replace into the database
#      otherwise will insert and primary_id will be added
#

# this version stores the object and flags it to be indexed
# for search via attributes, name, type or location

sub store {
  my $self = shift;
  $self->_store(1,@_);
}

=head2 store_noindex

 Title   : store_noindex
 Usage   : $success = $db->store_noindex(@features)
 Function: store one or more features into the database without indexing
 Returns : true if successful
 Args    : list of Bio::SeqFeatureI objects
 Status  : public

This method stores a list of features into the database but does not
make them searchable. The only way to access the features is via their
primary IDs. This method is ordinarily only used internally to store
subfeatures that are not indexed.

=cut

# this version stores the object and flags it so that it is
# not searchable via attributes, name, type or location
# (typically used only for subfeatures)
sub store_noindex {
  my $self = shift;
  $self->_store(0,@_);
}

=head2 delete

 Title   : delete
 Usage   : $success = $db->delete(@features)
 Function: delete a list of feature from the database
 Returns : true if successful
 Args    : list of features
 Status  : public

This method looks up the primary IDs from a list of features and
deletes them from the database, returning true if all deletions are
successful.

WARNING: The current DBI::mysql implementation has some issues that
need to be resolved, namely (1) normalized subfeatures are NOT
recursively deleted; and (2) the deletions are not performed in a
transaction.

=cut

sub delete {
  my $self   = shift;
  my $success = 1;
  for my $object (@_) {
    my $id = $object->primary_id;
    $success &&= $self->_deleteid($id);
  }
  $success;
}

=head2 update

 Title   : update
 Usage   : $success = $db->update($feature)
 Function: update a single feature
 Returns : true if successful
 Args    : feature to be updated
 Status  : public

This method writes a modified feature back into the database. The
feature must contain a valid primary_id.

When using normalized feature classes such as
Bio::DB::SeqFeature::LazyTableFeature the subfeatures are not
recursively updated when you update the parent feature. You must
manually update each of the subfeatures you have changed.

=cut

sub update {
  my $self = shift;
  my $object = shift;
  defined (my $primary_id = eval { $object->primary_id})
    or $self->throw("$object has no primary ID: $@");
  $self->_update($object,$primary_id);
}



=head2 fetch

 Title   : fetch
 Usage   : $feature = $db->fetch($primary_id)
 Function: fetch a feature from the database using its primary ID
 Returns : a feature
 Args    : primary ID of desired feature
 Status  : public

This method returns a previously-stored feature from the database
using its primary ID. If the primary ID is invalid, it returns undef.

=cut

###
# Fetch a Bio::SeqFeatureI from database using its primary_id
#
sub fetch {
  my $self       = shift;
  @_ or croak "usage: fetch(\$primary_id)";
  my $primary_id = shift;
  $self->_fetch($primary_id);
}

=head2 fetch_many

 Title   : fetch_many
 Usage   : @features = $db->fetch_many($primary_id,$primary_id,$primary_id...)
 Function: fetch many features from the database using their primary ID
 Returns : list of features
 Args    : a list of primary IDs or an array ref of primary IDs
 Status  : public

Same as fetch() except that you can pass a list of primary IDs or a
ref to an array of IDs.

=cut

###
# Efficiently fetch a series of IDs from the database
# Can pass an array or an array ref
#
sub fetch_many {
  my $self       = shift;
  @_ or croak 'usage: fetch_many($id1,$id2,$id3...)';
  my @ids = map {ref($_) ? @$_ : $_} @_ or return;
  $self->_fetch_many(@ids);
}

=head2 get_seq_stream

 Title   : get_seq_stream
 Usage   : $iterator = $db->get_seq_stream(@args)
 Function: return an iterator across all features in the database
 Returns : a Bio::DB::SeqFeature::Store::Iterator object
 Args    : (optional) the feature() method
 Status  : public

When called without any arguments this method will return an iterator
object that will traverse all indexed features in the database. Call
the iterator's next_seq() method to step through them (in no
particular order):

  my $iterator = $db->get_seq_stream;
  while (my $feature = $iterator->next_seq) {
    print $feature->primary_tag,' ',$feature->display_name,"\n";
  }

You can select a subset of features by passing a series of filter
arguments. The arguments are identical to those accepted by
$db->features().

=cut

###
# Return an iterator across all features that are indexable
#
sub get_seq_stream {
  my $self = shift;
  $self->_features(-iterator=>1,@_);
}

=head2 get_features_by_name

 Title   : get_features_by_name
 Usage   : @features = $db->get_features_by_name($name)
 Function: looks up features by their display_name
 Returns : a list of matching features
 Args    : the desired name
 Status  : public

This method searches the display_name of all features for matches
against the provided name. GLOB style wildcares ("*", "?") are
accepted, but may be slow.

The method returns the list of matches, which may be zero, 1 or more
than one features. Be prepared to receive more than one result, as
display names are not guaranteed to be unique.

For backward compatibility with gbrowse, this method is also known as
get_feature_by_name().

=cut

# backward compatibility for gbrowse
sub get_feature_by_name { shift->get_features_by_name(@_) }

###
# get_feature_by_name() return 0 or more features using a name lookup
# uses the Bio::DB::GFF API
#
sub get_features_by_name {
  my $self   = shift;
  my ($class,$name,$types,$allow_alias);

  if (@_ == 1) {  # get_features_by_name('name');
    $name = shift;
  } else {        # get_features_by_name('class'=>'name'), get_feature_by_name(-name=>'name')
    ($class,$name,$allow_alias,$types) = rearrange([qw(CLASS NAME ALIASES),[qw(TYPE TYPES)]],@_);
  }

  # hacky workaround for assumption in Bio::DB::GFF that unclassed reference points were of type "Sequence"
  undef $class if $class && $class eq 'Sequence';

  $self->_features(-name=>$name,-class=>$class,-aliases=>$allow_alias,-type=>$types);
}

=head2 get_features_by_alias

 Title   : get_features_by_alias
 Usage   : @features = $db->get_features_by_alias($name)
 Function: looks up features by their display_name or alias
 Returns : a list of matching features
 Args    : the desired name
 Status  : public

This method is similar to get_features_by_name() except that it will
also search through the feature aliases.  Aliases can be created by
storing features that contain one or more Alias tags. Wildards are
accepted.

=cut

sub get_features_by_alias {
  my $self = shift;
  my @args = @_;
  if (@_ == 1) {
    @args  = (-name=>shift);
  }
  push @args,(-aliases=>1);
  $self->get_features_by_name(@args);
}

=head2 get_features_by_type

 Title   : get_features_by_type
 Usage   : @features = $db->get_features_by_type(@types)
 Function: looks up features by their primary_tag
 Returns : a list of matching features
 Args    : list of primary tags
 Status  : public

This method will return a list of features that have any of the
primary tags given in the argument list. For compatibility with
gbrowse and Bio::DB::GFF, types can be qualified using a colon:

  primary_tag:source_tag

in which case only features that match both the primary_tag B<and> the
indicated source_tag will be returned. If the database was loaded from
a GFF3 file, this corresponds to the third and second columns of the
row, in that order.

For example, given the GFF3 lines:

  ctg123 geneFinder exon 1300 1500 . + . ID=exon001
  ctg123 fgenesH    exon 1300 1520 . + . ID=exon002

exon001 and exon002 will be returned by searching for type "exon", but
only exon001 will be returned by searching for type "exon:fgenesH".

=cut

sub get_features_by_type {
  my $self = shift;
  my @types = @_;
  $self->_features(-type=>\@types);
}

=head2 get_features_by_location

 Title   : get_features_by_location
 Usage   : @features = $db->get_features_by_location(@args)
 Function: looks up features by their location
 Returns : a list of matching features
 Args    : see below
 Status  : public

This method fetches features based on a location range lookup. You
call it using a positional list of arguments, or a list of
(-argument=>$value) pairs.

The positional form is as follows:

 $db->get_features_by_location($seqid [[,$start,]$end])

The $seqid is the name of the sequence on which the feature resides,
and start and end are optional endpoints for the match. If the
endpoints are missing then any feature on the indicated seqid is
returned.

Examples:

 get_features_by_location('chr1');      # all features on chromosome 1
 get_features_by_location('chr1',5000); # features between 5000 and the end
 get_features_by_location('chr1',5000,8000); # features between 5000 and 8000

Location lookups are overlapping. A feature will be returned if it
partially or completely overlaps the indicated range.

The named argument form gives you more control:

  Argument       Value
  --------       -----

  -seq_id        The name of the sequence on which the feature resides
  -start         Start of the range
  -end           End of the range
  -strand        Strand of the feature
  -range_type    Type of range to search over

The B<-strand> argument, if present, can be one of "0" to find
features that are on both strands, "+1" to find only plus strand
features, and "-1" to find only minus strand features. Specifying a
strand of undef is the same as not specifying this argument at all,
and retrieves all features regardless of their strandedness.

The B<-range_type> argument, if present, can be one of "overlaps" (the
default), to find features whose positions overlap the indicated
range, "contains," to find features whose endpoints are completely
contained within the indicated range, and "contained_in" to find
features whose endpoints are both outside the indicated range.

=cut

sub get_features_by_location {
  my $self = shift;
  my ($seqid,$start,$end,$strand,$rangetype) = 
    rearrange([['SEQ_ID','SEQID','REF'],'START',['STOP','END'],'STRAND','RANGE_TYPE'],@_);
  $self->_features(-seqid=>$seqid,
		   -start=>$start||undef,
		   -end=>$end||undef,
		   -strand=>$strand||undef,
		   -range_type=>$rangetype);
}

=head2 get_features_by_attribute

 Title   : get_features_by_attribute
 Usage   : @features = $db->get_features_by_attribute(@args)
 Function: looks up features by their attributes/tags
 Returns : a list of matching features
 Args    : see below
 Status  : public

MORE DOCUMENTATION PENDING

=cut

sub get_features_by_attribute {
  my $self       = shift;
  my $attributes = shift;
  $attributes  or croak "Usage: get_feature_by_attribute({attribute_hash})";
  $self->_features(-attributes=>$attributes);
}
###
# features() call -- main query interface
#

# documentation of args
#   my ($seq_id,$start,$end,$strand,
#       $name,$class,$allow_aliases,
#       $types,
#       $attributes,
#       $range_type,
#       $iterator,
#      ) = rearrange([['SEQID','SEQ_ID','REF'],'START',['STOP','END'],'STRAND',
# 		    'NAME','CLASS','ALIASES',
# 		    ['TYPES','TYPE','PRIMARY_TAG'],
# 		    ['ATTRIBUTES','ATTRIBUTE'],
# 		    'RANGE_TYPE',
# 		   ],@_);
#   $range_type ||= 'overlaps';
sub features {
  my $self = shift;
  my @args;
  if (@_ == 1) {
    @args = (-type=>shift);
  } else {
    @args = @_;
  }
  $self->_features(@args);
}

###
# search_notes()
#
sub search_notes {
  my $self = shift;
  my ($search_string,$limit) = @_;
  return $self->_search_notes($search_string,$limit);
}

###
# insert_sequence()
#
# insert a bit of primary sequence into the database
#
sub insert_sequence {
  my $self = shift;
  my ($seqid,$seq,$offset) = @_;
  $offset ||= 0;
  $self->_insert_sequence($seqid,$offset,$seq);
}

###
# get_sequence()
#
# equivalent to old Bio::DB::GFF->dna() method
#
sub fetch_sequence {
  my $self = shift;
  my ($seqid,$start,$end,$class) = rearrange([['NAME','SEQID','SEQ_ID'],'START',['END','STOP'],'CLASS'],@_);
  $seqid = "$seqid:$class" if defined $class;
  $self->_fetch_sequence($seqid,$start,$end);
}

=head2 segment

 Title   : segment
 Usage   : $segment = $db->segment($seq_id [,$start] [,$end])
 Function: restrict the database to a sequence range
 Returns : a Bio::DB::SeqFeature::Segment object
 Args    : sequence id, start and end ranges (optional)
 Status  : public

This is a convenience method that can be used when you are interested
in the contents of a particular sequence landmark, such as a
contig. Specify the ID of a sequence or other landmark in the database
and optionally a start and endpoint relative to that landmark. The
method will look up the region and return a
Bio::DB::SeqFeature::Segment object that spans it. You can then use
this segment object to make location-restricted queries on the database.

Example:

 $segment  = $db->segment('contig23',1,1000);  # first 1000 bp of contig23
 my @mRNAs = $segment->features('mRNA');       # all mRNAs that overlap segment

Although you will usually want to fetch segments that correspond to
physical sequences in the database, you can actually use any feature
in the database as the sequence ID. The segment() method will perform
a get_features_by_name() internally and then transform the feature
into the appropriate coordinates.

=cut

###
# Replacement for Bio::DB::GFF->segment() method
#
sub segment {
  my $self = shift;
  my (@features,@args);

  if (@_ == 1 && blessed($_[0])) {
    @features = @_;
    @args = ();
  }
  else {
    @args     = $self->setup_segment_args(@_);
    @features = $self->get_features_by_name(@args);
  }
  if (!wantarray && @features > 1) {
    $self->throw(<<END);
segment() called in a scalar context but multiple features match.
Either call in a list context or narrow your search using the -types or -class arguments
END
  }
  my ($rel_start,$rel_end) = rearrange(['START',['STOP','END']],@args);
  $rel_start = 1 unless defined $rel_start;

  my @segments;
  for my $f (@features) {
    my $seqid  = $f->seq_id;
    my $strand = $f->strand;
    my ($start,$end);
    $rel_end = $f->end - $f->start + 1 unless defined $rel_end;

    if ($strand >= 0) {
      $start = $f->start + $rel_start - 1;
      $end   = $f->start + $rel_end   - 1;
    }
    else {
      $start = $f->end - $rel_end   + 1;
      $end   = $f->end - $rel_start + 1;
    }
    push @segments,Bio::DB::SeqFeature::Segment->new($self,$seqid,$start,$end,$strand);
  }
  return wantarray ? @segments : $segments[0];
}

sub setup_segment_args {
  my $self = shift;
  return @_ if defined $_[0] && $_[0] =~ /^-/;
  return (-name=>$_[0],-start=>$_[1],-end=>$_[2]) if @_ == 3;
  return (-class=>$_[0],-name=>$_[1])              if @_ == 2;
  return (-name=>$_[0])                            if @_ == 1;
  return;
}

###
# force reindexing
#
sub reindex {
  my $self = shift;

  my $count = 0;
  my $now;
  my $last_time = time();

  $self->_start_reindexing;

  my $iterator = $self->get_seq_stream;
  while (my $f = $iterator->next_seq) {
    if (++$count %1000 == 0) {
      $now = time();
      my $elapsed = sprintf(" in %5.2fs",$now - $last_time);
      $last_time = $now;
      print STDERR "$count features indexed$elapsed...",' 'x60;
      print STDERR -t STDOUT && !$ENV{EMACS} ? "\r" : "\n";
    }
    $self->_update_indexes($f);
  }

  $self->_end_reindexing;
}

sub _load_class {
  my $self = shift;
  my $obj  = shift;
  return if $self->{class_loaded}{ref $obj}++;
  unless ($obj && $obj->can('primary_id')) {
    my $class = ref $obj;
    eval "require $class";
  }
}

sub start_bulk_update  { shift->_start_bulk_update(@_) }
sub finish_bulk_update { shift->_finish_bulk_update(@_) }

=head1 Changing the Behavior of the Database

These methods allow you to modify the behavior of the database.

=head2 debug

 Title   : debug
 Usage   : $debug_flag = $db->debug([$new_flag])
 Function: set the debug flag
 Returns : current debug flag
 Args    : new debug flag
 Status  : public

This method gets/sets a flag that turns on verbose progress
messages. Currently this will not do very much.

=cut

sub debug {
  my $self = shift;
  my $d = $self->{debug};
  $self->{debug} = shift if @_;
  $d;
}

=head2 serializer

 Title   : serializer
 Usage   : $serializer = $db->serializer([$new_serializer])
 Function: get/set the name of the serializer
 Returns : the name of the current serializer class
 Args    : (optional) the name of a new serializer
 Status  : public

You can use this method to set the serializer, but do not attempt to
change the serializer once the database is initialized and populated.

=cut

###
# serializer
#
sub serializer {
  my $self = shift;
  my $d    = $self->setting('serializer');
  if (@_) {
    my $serializer = shift;
    eval "require $serializer; 1" or croak $@;
    $self->setting(serializer=>$serializer);
    $Storable::forgive_me =1 if $serializer eq 'Storable';
  }
  $d;
}

=head2 index_subfeatures

 Title   : index_subfeatures
 Usage   : $flag = $db->index_subfeatures([$new_value])
 Function: flag whether to index subfeatures
 Returns : current value of the flag
 Args    : (optional) new value of the flag
 Status  : public

If true, the store() method will add a searchable index to both the
top-level feature and all its subfeatures, allowing the search
functions to return features at any level of the conainment
hierarchy. If false, only the top level feature will be indexed,
meaning that you will only be able to get at subfeatures by fetching
the top-level feature and then traversing downward using
get_SeqFeatures().

You are free to change this setting at any point during the creation
and population of a database. One database can contain both indexed
and unindexed subfeatures.

=cut

###
# whether to index subfeatures by default
#
sub index_subfeatures {
  my $self = shift;
  my $d    = $self->setting('index_subfeatures');
  $self->setting('index_subfeatures'=>shift) if @_;
  $d;
}

################################# TIE interface ####################

sub TIEHASH {
  my $class = shift;
  return $class->new(@_);
}

sub STORE {
  my $self = shift;
  my ($key,$feature) = @_;
  $key =~ /^\d+$/ && $key > 0 or croak "keys must be positive integers";
  $self->_load_class($feature);
  $feature->primary_id($key);
  $self->store($feature);
}

sub FETCH {
  my $self = shift;
  $self->fetch(@_);
}

sub FIRSTKEY {
  my $self = shift;
  $self->_firstid;
}

sub NEXTKEY {
  my $self    = shift;
  my $lastkey = shift;
  $self->_nextid($lastkey);
}

sub EXISTS {
  my $self = shift;
  my $key  = shift;
  $self->existsid($key);
}

sub DELETE {
  my $self = shift;
  my $key  = shift;
  $self->_deleteid($key);
}

sub CLEAR {
  my $self = shift;
  $self->_clearall;
}

sub SCALAR {
  my $self = shift;
  $self->_featurecount;
}

=head1 Private Methods

These methods are private to Bio::DB::SeqFeature::Store and its
descendents (the adaptors). Some methods in this section B<must> be
implemented by adaptors in order to have full functionality. They are
preceded by an understore (_).

=head2 init

 Title   : init
 Usage   : $db->init(@args)
 Function: initialize object
 Returns : none
 Args    : Arguments passed to new()
 Status  : private

This method is called internally by new() to initialize a
newly-created object using the arguments passed to new(). It is to be
overridden by Bio::DB::SeqFeature::Store adaptors.

=cut

sub init {
  my $self = shift;
  $self->default_settings();
}

=head2 default_settings

 Title   : default_settings
 Usage   : $db->default_settings()
 Function: set up default settings for the adaptor
 Returns : none
 Args    : none
 Status  : private

This method is may be overridden by adaptors. It is responsible for
setting up object default settings.

=cut

###
# default settings -- set up whatever are the proper default settings
#
sub default_settings {
  my $self = shift;
  $self->serializer($self->default_serializer);
  $self->index_subfeatures(1);
}

=head2 default_serializer

 Title   : default_serializer
 Usage   : $serializer = $db->default_serializer
 Function: finds an available serializer
 Returns : the name of an available serializer
 Args    : none
 Status  : private

This method returns the name of an available serializer module.

=cut

###
# choose a serializer
#
sub default_serializer {
  my $self = shift;
  # try Storable
  eval "require Storable; 1"     and return 'Storable';
  eval "require Data::Dumper; 1" and return 'Data::Dumper';
  croak "Unable to load either Storable or Data::Dumper. Please provide a serializer using -serializer";
}

=head2 setting

 Title   : setting
 Usage   : $value = $db->setting('setting_name' [=> $new_value])
 Function: get/set the value of a setting
 Returns : the value of the current setting
 Args    : the name of the setting and optionally a new value for the setting
 Status  : private

This is a low-level procedure for persistently storing database
settings. It can be overridden by adaptors.

=cut

# persistent settings
# by default we store in the object
sub setting {
  my $self  = shift;
  my $variable_name = shift;
  my $d    = $self->{setting}{$variable_name};
  $self->{setting}{$variable_name} = shift if @_;
  $d;
}

=head2 subfeatures_are_indexed

 Title   : subfeatures_are_indexed
 Usage   : $flag = $db->subfeatures_are_indexed([$new_value])
 Function: flag whether subfeatures are indexed
 Returns : a flag indicating that all subfeatures are indexed
 Args    : (optional) new value of the flag
 Status  : private

This method is used internally by the
Bio::DB::SeqFeature::LazyTableFeature class to optimize some of its
operations. It returns true if all of the subfeatures in the database
are indexed; it returns false if at least one of the subfeatures is
not indexed. Do not attempt to change the value of this setting unless
you are writing an adaptor.

=cut

###
# whether subfeatures are all indexed
#
sub subfeatures_are_indexed {
  my $self = shift;
  my $d    = $self->setting('subfeatures_are_indexed');
  $self->setting(subfeatures_are_indexed => shift) if @_;
  $d;
}


=head2 add_SeqFeature

 Title   : add_SeqFeature
 Usage   : $db->add_SeqFeature($parent,@children)
 Function: add a series of subfeatures to a parent feature
 Returns : true if successful
 Args    : parent feature and list of subfeatures
 Status  : private

This method is used internally by the
Bio::DB::SeqFeature::LazyTableFeature class to store its
parent/children relationships in the database.

=cut

###
# Add a subparts to a feature. Both feature and all subparts must already be in database.
#
sub add_SeqFeature {
  my $self     = shift;
  my $parent   = shift;
  my @children = @_;

  $self->_add_SeqFeature($parent,@children);
}

=head2 fetch_SeqFeatures

 Title   : fetch_SeqFeatures
 Usage   : @children = $db->fetch_SeqFeatures($parent)
 Function: retrieve the subfeatures of a parent feature
 Returns : list of Bio::SeqFeatureI objects
 Args    : parent feature
 Status  : private

This method is used internally by the
Bio::DB::SeqFeature::LazyTableFeature class to retreive its
parent/children relationships from the database.

=cut

sub fetch_SeqFeatures {
  my $self   = shift;
  my $parent = shift;
  my @types  = @_;
  $self->_add_SeqFeatures($parent,@types);
}


###################### TO BE IMPLEMENTED BY ADAPTOR ##########

# DOC THIS!!!!!
sub _new { shift->throw_not_implemented}

sub _init_database { shift->throw_not_implemented }

# _store($indexed,@objs)
sub _store {
  my $self    = shift;
  my $indexed = shift;
  my @objs    = @_;
  $self->throw_not_implemented;
}

# _store($indexed,@objs)
sub _update {
  my $self       = shift;
  my $object     = shift;
  my $primary_id = shift;
  $self->throw_not_implemented;
}

# this is called to index a feature
sub _update_indexes { shift->throw_not_implemented }

# these do not necessary have to be overridden
# they are called at beginning and end of reindexing process
sub _start_reindexing {}
sub _end_reindexing   {}

# _fetch($id)
sub _fetch { shift->throw_not_implemented }

# _fetch_many(@ids)
# this one will fall back to many calls on fetch() if you don't
# override it
sub _fetch_many {
  my $self = shift;
  return map {$self->_fetch($_)} @_;
}

# bottleneck query generator
sub _features { shift->throw_not_implemented }

sub _search_notes { shift->throw_not_implemented }

# return true here if the storage engine is prepared to store parent/child
# relationships using _add_SeqFeature and return them using _fetch_SeqFeatures
sub can_store_parentage { return; }

# these two are called only if _can_store_subFeatures() returns true
# _add_SeqFeature ($parent,@children)
sub _add_SeqFeature { shift->throw_not_implemented }

# _get_SeqFeatures($parent,@list_of_child_types)
sub _fetch_SeqFeatures {shift->throw_not_implemented }

# _fetch_sequence() is similar to old dna() method
sub _insert_sequence   { shift->throw_not_implemented }
sub _fetch_sequence    { shift->throw_not_implemented }

# for full TIE() interface  - not necessary to implement in most cases
sub _firstid  { shift->throw_not_implemented }
sub _nextid   { shift->throw_not_implemented }
sub _existsid { shift->throw_not_implemented }
sub _deleteid { shift->throw_not_implemented }
sub _clearall { shift->throw_not_implemented }
sub _featurecount { shift->throw_not_implemented }

# Optional flags to change behavior to optimize bulk updating.
sub _start_bulk_update { }
sub _finish_bulk_update { }

=head1 Internal Methods

These methods are used internally and will not ordinarily be of use to
either application-level scripts or those implementing adaptors.

=head2 new_instance

 Title   : new_instance
 Usage   : $db = $db->new_instance()
 Function: class constructor
 Returns : A descendent of Bio::DB::SeqFeature::Store
 Args    : none
 Status  : internal

This method is called internally by new() to create a new
uninitialized instance of Bio::DB::SeqFeature::Store. It is used
internally and should not be called by application software.

=cut

sub new_instance {
  my $class = shift;
  return bless {},ref($class) || $class;
}


#################################### Internal methods ####################

sub _freeze {
  my $self = shift;
  my $obj  = shift;

  # Bio::SeqFeature::Generic contains cleanup methods, so we need to
  # localize the methods to undef temporarily so that we can serialize
  local $obj->{'_root_cleanup_methods'} if exists $obj->{'_root_cleanup_methods'};

  my ($id,$store);
  $id    = $obj->primary_id();
  $obj->primary_id(undef);     # don't want primary ID to be stored in object
  eval {
    $store = $obj->object_store;
    $obj->object_store(undef);   # don't want a copy of the store in the object
  };
  my $serializer = $self->serializer;
  my $data;
  if ($serializer eq 'Data::Dumper') {
    my $d    = Data::Dumper->new([$obj]);
    $d->Terse(1);
    $d->Deepcopy(1);
    $data = $d->Dump;
  } elsif ($serializer eq 'Storable') {
    $data = Storable::freeze($obj);
  }

  $obj->primary_id($id);       # restore to original state
  eval {
    $obj->object_store($store);
  };

  return $data;
}

sub _thaw {
  my $self               = shift;
  my ($obj,$primary_id)  = @_;
  my $serializer = $self->serializer;
  my $object;
  if ($serializer eq 'Data::Dumper') {
    $object = $obj;
  } elsif ($serializer eq 'Storable') {
    $object = Storable::thaw($obj);
  }

  # remember the primary ID of this object as well as the
  # identity of the store, so that we can do lazy loading;
  # both of these are wrapped in an eval because not all
  # bioseqfeatures support them (or want to)
  $self->_load_class($object);
  eval {
    $object->primary_id($primary_id);
    $object->object_store($self);
  };
  $object;
}

1;
