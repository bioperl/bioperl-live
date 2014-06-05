package Bio::DB::SeqFeature::Store;


=head1 NAME

Bio::DB::SeqFeature::Store -- Storage and retrieval of sequence annotation data

=head1 SYNOPSIS

  use Bio::DB::SeqFeature::Store;

  # Open the feature database
  my $db = Bio::DB::SeqFeature::Store->new( -adaptor => 'DBI::mysql',
                                            -dsn     => 'dbi:mysql:test',
                                            -create  => 1 );

  # Get a feature from somewhere
  my $feature = Bio::SeqFeature::Generic->new(...);

  # Store it
  $db->store($feature) or die "Couldn't store!";

  # If absent, a primary ID is added to the feature when it is stored in the
  # database. Retrieve the primary ID
  my $id = $feature->primary_id;

  # Get the feature back out
  my $feature = $db->fetch($id);

  # .... which is identical to
  my $feature = $db->get_feature_by_primary_id($id);

  # Change the feature and update it
  $f->start(100);
  $db->store($f) or die "Couldn't update!";

  # Get all features at once
  my @features = $db->features( );

  # Retrieve multiple features by primary id
  my @features = $db->fetch_many(@list_of_ids);

  # ...by name
  @features = $db->get_features_by_name('ZK909');

  # ...by alias
  @features = $db->get_features_by_alias('sma-3');

  # ...by type
  @features = $db->get_features_by_type('gene');

  # ...by location
  @features = $db->get_features_by_location(-seq_id=>'Chr1',-start=>4000,-end=>600000);

  # ...by attribute
  @features = $db->get_features_by_attribute({description => 'protein kinase'})

  # ...by the GFF "Note" field
  @result_list = $db->search_notes('kinase');

  # ...by arbitrary combinations of selectors
  @features = $db->features(-name => $name,
                            -type => $types,
                            -seq_id => $seqid,
                            -start  => $start,
                            -end    => $end,
                            -attributes => $attributes);

  # Loop through the features using an iterator
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
  my @features = $segment->features(-type=>['mRNA','match']);

  # Getting coverage statistics across a region
  my $summary = $db->feature_summary('Chr1',10_000=>1_110_000);
  my ($bins)  = $summary->get_tag_values('coverage');
  my $first_bin = $bins->[0];

  # Getting & storing sequence information
  # Warning: this returns a string, and not a PrimarySeq object
  $db->insert_sequence('Chr1','GATCCCCCGGGATTCCAAAA...');
  my $sequence = $db->fetch_sequence('Chr1',5000=>6000);

  # What feature types are defined in the database?
  my @types    = $db->types;

  # Create a new feature in the database
  my $feature = $db->new_feature(-primary_tag => 'mRNA',
                                 -seq_id      => 'chr3',
                                 -start      => 10000,
                                 -end        => 11000);

  # Load an entire GFF3 file, using the GFF3 loader...
  my $loader = Bio::DB::SeqFeature::Store::GFF3Loader->new(-store    => $db,
							   -verbose  => 1,
							   -fast     => 1);

  $loader->load('./my_genome.gff3');



=head1 DESCRIPTION

Bio::DB::SeqFeature::Store implements the Bio::SeqFeature::CollectionI
interface to allow you to persistently store Bio::SeqFeatureI objects
in a database and to later to retrieve them by a variety of
searches. This module is similar to the older Bio::DB::GFF module,
with the following differences:

=over 4

=item 1.

No limitation on Bio::SeqFeatureI implementations

Unlike Bio::DB::GFF, Bio::DB::SeqFeature::Store works with
any Bio::SeqFeatureI object.

=item 2.

No limitation on nesting of features & subfeatures

Bio::DB::GFF is limited to features that have at most one
level of subfeature. Bio::DB::SeqFeature::Store can work with features
that have unlimited levels of nesting.

=item 3.

No aggregators

The aggregator architecture, which was necessary to impose order on
the GFF2 files that Bio::DB::GFF works with, does not apply to
Bio::DB::SeqFeature::Store. It is intended to store features that obey
well-defined ontologies, such as the Sequence Ontology
(http://song.sourceforge.net).

=item 4.

No relative locations

All locations defined by this module are relative to an absolute
sequence ID, unlike Bio::DB::GFF which allows you to define the
location of one feature relative to another.

=back

We'll discuss major concepts in Bio::DB::SeqFeature::Store and then
describe how to use the module.

=head2 Adaptors

Bio::DB::SeqFeature::Store is designed to work with a variety of
storage back ends called "adaptors." Adaptors are subclasses of
Bio::DB::SeqFeature::Store and provide the interface between the
store() and fetch() methods and the physical database. Currently the
number of adaptors is quite limited, but the number will grow soon.

=over 4

=item memory

An implementation that stores all data in memory. This is useful for
small data sets of no more than 10,000 features (more or less,
depending on system memory).

=item DBI::mysql

A full-featured implementation on top of the MySQL relational database
system.

=item berkeleydb

A full-feature implementation that runs on top of the BerkeleyDB
database. See L<Bio::DB::SeqFeature::Store::berkeleydb>.


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

=head2 Loaders and Normalized Features

The Bio::DB::SeqFeature::Store::GFF3Loader parses a GFF3-format file
and loads the annotations and sequence data into the database of your
choice. The script bp_seqfeature_load.pl (found in the
scripts/Bio-SeqFeature-Store/ subdirectory) is a thin front end to the
GFF3Loader. Other loaders may be written later.

Although Bio::DB::SeqFeature::Store should work with any
Bio::SeqFeatureI object, there are some disadvantages to using
Bio::SeqFeature::Generic and other vanilla implementations. The major
issue is that if two vanilla features share the same subfeature
(e.g. two transcripts sharing an exon), the shared subfeature will be
cloned when stored into the database.

The special-purpose L<Bio::DB::SeqFeature> class is able to normalize
its subfeatures in the database, so that shared subfeatures are stored
only once. This minimizes wasted storage space. In addition, when
in-memory caching is turned on, each shared subfeature will usually
occupy only a single memory location upon restoration.

=cut


use strict;
use warnings;

use base 'Bio::SeqFeature::CollectionI';
use Carp 'croak';
use Bio::DB::GFF::Util::Rearrange;
use Bio::DB::SeqFeature::Segment;
use Scalar::Util 'blessed';

# this probably shouldn't be here
use Bio::DB::SeqFeature;

*dna = *get_dna = *get_sequence = \&fetch_sequence;
*get_SeqFeatures = \&fetch_SeqFeatures;

# local version
sub api_version { 1.2 }

=head1 Methods for Connecting and Initializating a Database

## TODO: http://iowg.brcdevel.org/gff3.html#a_fasta is a dead link

=head2 new

 Title   : new
 Usage   : $db = Bio::DB::SeqFeature::Store->new(@options)
 Function: connect to a database
 Returns : A descendent of Bio::DB::Seqfeature::Store
 Args    : several - see below
 Status  : public

This class method creates a new database connection. The following
-name=E<gt>$value arguments are accepted:

 Name               Value
 ----               -----

 -adaptor           The name of the Adaptor class (default DBI::mysql)

 -serializer        The name of the serializer class (default Storable)

 -index_subfeatures Whether or not to make subfeatures searchable
                    (default false)

 -cache             Activate LRU caching feature -- size of cache

 -compress          Compresses features before storing them in database
                    using Compress::Zlib

 -create            (Re)initialize the database.

The B<-index_subfeatures> argument, if true, tells the module to
create indexes for a feature and all its subfeatures (and its
subfeatures' subfeatures). Indexing subfeatures means that you will be
able to search for the gene, its mRNA subfeatures and the exons inside
each mRNA. It also means when you search the database for all features
contained within a particular location, you will get the gene, the
mRNAs and all the exons as individual objects as well as subfeatures
of each other. NOTE: this option is only honored when working with a
normalized feature class such as Bio::DB::SeqFeature.

The B<-cache> argument, if true, tells the module to try to create a
LRU (least-recently-used) object cache using the Tie::Cacher
module. Caching will cause two objects that share the same primary_id
to (often, but not always) share the same memory location, and may
improve performance modestly. The argument is taken as the desired
size for the cache. If you pass "1" as the cache value, a reasonable
default cache size will be chosen. Caching requires the Tie::Cacher
module to be installed. If the module is not installed, then caching
will silently be disabled.

The B<-compress> argument, if true, will cause the feature data to be
compressed before storing it. This will make the database somewhat
smaller at the cost of decreasing performance.

The B<-create> argument, if true, will either initialize or
reinitialize the database. It is needed the first time a database is
used.

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
		    L<Bio::DB::SeqFeature::Store::GFF3Loader> for a
		    description of this. Default is the current
                    directory.

 -write             Make the database writable (implied by -create)

 -fasta             Provide an alternative DNA accessor object or path.

By default the database will store DNA sequences internally. However,
you may override this behavior by passing either a path to a FASTA
file, or any Perl object that recognizes the seq($seqid,$start,$end)
method. In the former case, the FASTA path will be passed to
Bio::DB::Fasta, possibly causing an index to be constructed. Suitable
examples of the latter type of object include the Bio::DB::Sam and
Bio::DB::Sam::Fai classes.

=cut

###
# object constructor
#
sub new {
  my $self      = shift;
  my ($adaptor,$serializer,$index_subfeatures,$cache,$compress,$debug,$create,$fasta,$args);
  if (@_ == 1) {
    $args = {DSN => shift}
  }
  else {
    ($adaptor,$serializer,$index_subfeatures,$cache,$compress,$debug,$create,$fasta,$args) =
      rearrange(['ADAPTOR',
		 'SERIALIZER',
		 'INDEX_SUBFEATURES',
		 'CACHE',
		 'COMPRESS',
		 'DEBUG',
		 'CREATE',
		 'FASTA',
		],@_);
  }
  $adaptor ||= 'DBI::mysql';
  $args->{WRITE}++  if $create;
  $args->{CREATE}++ if $create;

  my $class = "Bio::DB::SeqFeature::Store::$adaptor";
  eval "require $class " or croak $@;
  $cache &&= eval "require Tie::Cacher; 1";
  my $obj = $class->new_instance();
  $obj->debug($debug) if defined $debug;
  $obj->init($args);
  $obj->init_cache($cache) if $cache;
  $obj->do_compress($compress);
  $obj->serializer($serializer)               if defined $serializer;
  $obj->index_subfeatures($index_subfeatures) if defined $index_subfeatures;
  $obj->seqfeature_class('Bio::DB::SeqFeature');
  $obj->set_dna_accessor($fasta)              if defined $fasta;
  $obj->post_init($args);
  $obj;
}

=head2 init_database

 Title   : init_database
 Usage   : $db->init_database([$erase_flag])
 Function: initialize a database
 Returns : true
 Args    : (optional) flag to erase current data
 Status  : public

Call this after Bio::DB::SeqFeature::Store-E<gt>new() to initialize a
new database. In the case of a DBI database, this method installs the
schema but does B<not> create the database. You have to do this
offline using the appropriate command-line tool. In the case of the
"berkeleydb" adaptor, this creates an empty BTREE database.

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

=head2 post_init

This method is invoked after init_database for use by certain adaptors
(currently only the memory adaptor) to do automatic data loading after
initialization. It is passed a copy of the init_database() args.

=cut

sub post_init { }


=head2 add_features

 Title   : add_features
 Usage   : $success = $db->add_features(\@features)
 Function: store one or more features into the database
 Returns : true if successful
 Args    : array reference of Bio::SeqFeatureI objects
 Status  : public

=cut

sub add_features {
  my ($self, $feats) = @_;
  my $result = $self->store_and_cache(1, @$feats);
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

In most cases, you should let the database assign the primary id. If
the object you store already has a primary_id, then the ID must adhere
to the datatype expected by the adaptor: an integer in the
case of the various DB adaptors, and a string in the case of the
memory and berkeley adaptors.

You can find out what the primary ID of the feature has become by
calling the feature's primary_id() method:

  $db->store($my_feature) or die "Oh darn";
  my $id = $my_feature->primary_id;

If the feature contains subfeatures, they will all be stored
recursively. In the case of Bio::DB::SeqFeature and
Bio::DB::SeqFeature::Store::NormalizedFeature, the subfeatures will be
stored in a normalized way so that each subfeature appears just once
in the database.

Subfeatures will be indexed for separate retrieval based on the
current value of index_subfeatures().

If you call store() with one or more features that already have valid
primary_ids, then any existing objects will be B<replaced>. Note that
when using normalized features such as Bio::DB::SeqFeature, the
subfeatures are not recursively updated when you update the parent
feature. You must manually update each subfeatures that has changed.

=cut

###
# store one or more Bio::SeqFeatureI objects
#      if they already have a primary_id will replace into the database
#      otherwise will insert and primary_id will be added
#

# this version stores the object and flags it to be indexed
# for search via attributes, name, type or location

sub store {
  my ($self, @feats) = @_;
  for my $feat (@feats) {
    if ( (not ref $feat) || (not $feat->isa('Bio::SeqFeatureI')) ) {
      die "Cannot store non-Bio::SeqFeatureI object '$feat'\n";
    }
  }
  my $result = $self->store_and_cache(1,@feats);
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
  $self->store_and_cache(0,@_);
}

=head2 no_blobs

 Title   : no_blobs
 Usage   : $db->no_blobs(1);
 Function: decide if objects should be stored in the database as blobs.
 Returns : boolean (default false)
 Args    : boolean (true to no longer store objects; when the corresponding
           feature is retrieved it will instead be a minimal representation of
           the object that was stored, as some simple Bio::SeqFeatureI object)
 Status  : dubious (new)

This method saves lots of space in the database, which may in turn lead to large
performance increases in extreme cases (over 7 million features in the db).

=cut

sub no_blobs {
    my $self = shift;
    if (@_) { $self->{no_blobs} = shift }
    return $self->{no_blobs} || 0;
}

=head2 new_feature

 Title   : new_feature
 Usage   : $feature = $db->new_feature(@args)
 Function: create a new Bio::DB::SeqFeature object in the database
 Returns : the new seqfeature
 Args    : see below
 Status  : public

This method creates and stores a new Bio::SeqFeatureI object using the
specialized Bio::DB::SeqFeature class. This class is able to store its
subfeatures in a normalized fashion, allowing subfeatures to be shared
among multiple parents (e.g. multiple exons shared among several
mRNAs).

The arguments are the same as for Bio::DB::SeqFeature-E<gt>new(), which in
turn are similar to Bio::SeqFeature::Generic-E<gt>new() and
Bio::Graphics::Feature-E<gt>new(). The most important difference is the
B<-index> option, which controls whether the feature will be indexed
for retrieval (default is true). Ordinarily, you would only want to
turn indexing off when creating subfeatures, because features stored
without indexes will only be reachable via their primary IDs or their
parents.

Arguments are as follows:

  -seq_id       the reference sequence
  -start        the start position of the feature
  -end          the stop position of the feature
  -display_name the feature name (returned by seqname)
  -primary_tag  the feature type (returned by primary_tag)
  -source       the source tag
  -score        the feature score (for GFF compatibility)
  -desc         a description of the feature
  -segments     a list of subfeatures (see Bio::Graphics::Feature)
  -subtype      the type to use when creating subfeatures
  -strand       the strand of the feature (one of -1, 0 or +1)
  -phase        the phase of the feature (0..2)
  -url          a URL to link to when rendered with Bio::Graphics
  -attributes   a hashref of tag value attributes, in which the key is the tag
                  and the value is an array reference of values
  -index        index this feature if true

Aliases:

  -id           an alias for -display_name
  -seqname      an alias for -display_name
  -display_id   an alias for -display_name
  -name         an alias for -display_name
  -stop         an alias for end
  -type         an alias for primary_tag

You can change the seqfeature implementation generated by new() by
passing the name of the desired seqfeature class to
$db-E<gt>seqfeature_class().

=cut

sub new_feature {
  my $self = shift;
  return $self->seqfeature_class->new(-store=>$self,@_);
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
    if ( not defined $id ) {
      warn "Could not delete feature without primary_id: $object";
      $success = 0;
      next;
    }
    my $result = $self->_deleteid($id);
    warn "Could not delete feature with id=$id" unless $result;
    $success &&= $result;
  }
  $success;
}

=head2 fetch / get_feature_by_id / get_feature_by_primary_id

 Title   : fetch
           get_feature_by_id
           get_feature_by_primary_id
 Usage   : $feature = $db->fetch($primary_id)
 Function: fetch a feature from the database using its primary ID
 Returns : a feature
 Args    : primary ID of desired feature
 Status  : public

This method returns a previously-stored feature from the database
using its primary ID. If the primary ID is invalid, it returns undef.
Use fetch_many() to rapidly retrieve multiple features.

=cut

###
# Fetch a Bio::SeqFeatureI from database using its primary_id
#
sub fetch {
  my $self       = shift;
  @_ or croak "usage: fetch(\$primary_id)";
  my $primary_id = shift;
  if (my $cache = $self->cache()) {
    return $cache->fetch($primary_id) if $cache->exists($primary_id);
    my $object = $self->_fetch($primary_id);
    $cache->store($primary_id,$object);
    return $object;
  }
  else {
    return $self->_fetch($primary_id);
  }
}
*get_feature_by_id = *get_feature_by_primary_id = \&fetch;

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
 Args    : feature filters (optional)
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
$db-E<gt>features().

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

=head2 get_feature_by_name

 Title   : get_feature_by_name
 Usage   : @features = $db->get_feature_by_name($name)
 Function: looks up features by their display_name
 Returns : a list of matching features
 Args    : the desired name
 Status  : Use get_features_by_name instead.

This method is provided for backward compatibility with gbrowse.

=cut

sub get_feature_by_name { shift->get_features_by_name(@_) }

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
(-argument=E<gt>$value) pairs.

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

This implements a simple tag filter. Pass a list of tag names and
their values. The module will return a list of features whose tag
names and values match. Tag names are case insensitive. If multiple
tag name/value pairs are present, they will be ANDed together. To
match any of a list of values, use an array reference for the value.

Examples:

 # return all features whose "function" tag is "GO:0000123"
 @features = $db->get_features_by_attribute(function => 'GO:0000123');

 # return all features whose "function" tag is "GO:0000123" or "GO:0000555"
 @features = $db->get_features_by_attribute(function => ['GO:0000123','GO:0000555']);

 # return all features whose "function" tag is "GO:0000123" or "GO:0000555"
 # and whose "confirmed" tag is 1
 @features = $db->get_features_by_attribute(function  => ['GO:0000123','GO:0000555'],
                                            confirmed => 1);

=cut

sub get_features_by_attribute {
  my $self       = shift;
  my %attributes = ref($_[0]) ? %{$_[0]} : @_;
  %attributes  or $self->throw("Usage: get_feature_by_attribute(attribute_name=>\$attribute_value...)");
  $self->_features(-attributes=>\%attributes);
}
###
# features() call -- main query interface
#

=head2 features

 Title   : features
 Usage   : @features = $db->features(@args)
 Function: generalized query & retrieval interface
 Returns : list of features
 Args    : see below
 Status  : Public

This is the workhorse for feature query and retrieval. It takes a
series of -name=E<gt>$value arguments filter arguments. Features that
match all the filters are returned.

  Argument       Value
  --------       -----

 Location filters:
  -seq_id        Chromosome, contig or other DNA segment
  -seqid         Synonym for -seq_id
  -ref           Synonym for -seqid
  -start         Start of range
  -end           End of range
  -stop          Synonym for -end
  -strand        Strand
  -range_type    Type of range match ('overlaps','contains','contained_in')

 Name filters:
  -name          Name of feature (may be a glob expression)
  -aliases       If true, match aliases as well as display names
  -class         Archaic argument for backward compatibility.
                  (-class=>'Clone',-name=>'ABC123') is equivalent
                  to (-name=>'Clone:ABC123')

 Type filters:
  -types         List of feature types (array reference) or one type (scalar)
  -type          Synonym for the above
  -primary_tag   Synonym for the above

  -attributes    Hashref of attribute=>value pairs as per
                    get_features_by_attribute(). Multiple alternative values
                    can be matched by providing an array reference.
  -attribute     synonym for -attributes

You may also provide features() with a list of scalar values (the
first element of which must B<not> begin with a dash), in which case
it will treat the list as a feature type filter.

Examples:

All features:
 @features = $db->features( );

All features on chromosome 1:

 @features = $db->features(-seqid=>'Chr1');

All features on chromosome 1 between 5000 and 6000:

 @features = $db->features(-seqid=>'Chr1',-start=>5000,-end=>6000);

All mRNAs on chromosome 1 between 5000 and 6000:

 @features = $db->features(-seqid=>'Chr1',-start=>5000,-end=>6000,-types=>'mRNA');

All confirmed mRNAs and repeats on chromosome 1 that overlap the range 5000..6000:

 @features = $db->features(-seqid     => 'Chr1',-start=>5000,-end=>6000,
                           -types     => ['mRNA','repeat'],
                           -attributes=> {confirmed=>1}
                          );

All confirmed mRNAs and repeats on chromosome 1 strictly contained within the range 5000..6000:

 @features = $db->features(-seqid     => 'Chr1',-start=>5000,-end=>6000,
                           -types     => ['mRNA','repeat'],
                           -attributes=> {confirmed=>1}
                           -range_type => 'contained_in',
                          );

All genes and repeats:

 @features = $db->features('gene','repeat_region');

=cut

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
  if (@_ == 0) {
    @args = ();
  }
  elsif ($_[0] !~/^-/) {
    my @types = @_;
    @args = (-type=>\@types);
  } else {
    @args = @_;
  }
  $self->_features(@args);
}


=head2 get_all_features

 Title   : get_all_features
 Usage   : @features = $db->get_all_features()
 Function: get all feature in the database
 Returns : list of features
 Args    : none
 Status  : Public

=cut

# for compatibility with Bio::SeqFeature::Collection
sub get_all_features {
  shift->features();
}

=head2 seq_ids

 Title   : seq_ids
 Usage   : @ids = $db->seq_ids()
 Function: Return all sequence IDs contained in database
 Returns : list of sequence Ids
 Args    : none
 Status  : public

=cut

sub seq_ids {
  my $self = shift;
  return $self->_seq_ids();
}

=head2 search_attributes

 Title   : search_attributes
 Usage   : @result_list = $db->search_attributes("text search string",[$tag1,$tag2...],$limit)
 Function: Search attributes for keywords occurring in a text string
 Returns : array of results
 Args    : full text search string, array ref of attribute names, and an optional feature limit
 Status  : public

Given a search string, this method performs a full-text search of the
specified attributes and returns an array of results.  You may pass a
scalar attribute name to search the values of one attribute
(e.g. "Note") or you may pass an array reference to search inside
multiple attributes (['Note','Alias','Parent']).Each row of the
returned array is a arrayref containing the following fields:

  column 1     The display name of the feature
  column 2     The text of the note
  column 3     A relevance score.
  column 4     The feature type
  column 5     The unique ID of the feature

NOTE: This search will fail to find features that do not have a display name!

You can use fetch() or fetch_many() with the returned IDs to get to
the features themselves.

=cut

sub search_attributes {
  my $self = shift;
  my ($search_string,$attribute_names,$limit) = @_;
  my $attribute_array   = ref $attribute_names
                      && ref $attribute_names eq 'ARRAY' ? $attribute_names : [$attribute_names];
  return $self->_search_attributes($search_string,$attribute_array,$limit);
}

=head2 search_notes

 Title   : search_notes
 Usage   : @result_list = $db->search_notes("full text search string",$limit)
 Function: Search the notes for a text string
 Returns : array of results
 Args    : full text search string, and an optional feature limit
 Status  : public

Given a search string, this method performs a full-text search of the
"Notes" attribute and returns an array of results.  Each row of the
returned array is a arrayref containing the following fields:

  column 1     The display_name of the feature, suitable for passing to get_feature_by_name()
  column 2     The text of the note
  column 3     A relevance score.
  column 4     The type

NOTE: This is equivalent to $db-E<gt>search_attributes('full text search
string','Note',$limit). This search will fail to find features that do
not have a display name!

=cut

###
# search_notes()
#
sub search_notes {
  my $self = shift;
  my ($search_string,$limit) = @_;
  return $self->_search_attributes($search_string,['Note'],$limit);
}

=head2 types

 Title   : types
 Usage   : @type_list = $db->types
 Function: Get all the types in the database
 Returns : array of Bio::DB::GFF::Typename objects
 Args    : none
 Status  : public

=cut

sub types {
    shift->throw_not_implemented;
}

=head2 insert_sequence

 Title   : insert_sequence
 Usage   : $success = $db->insert_sequence($seqid,$sequence_string,$offset)
 Function: Inserts sequence data into the database at the indicated offset
 Returns : true if successful
 Args    : see below
 Status  : public

This method inserts the DNA or protein sequence fragment
$sequence_string, identified by the ID $seq_id, into the database at
the indicated offset $offset. It is used internally by the GFF3Loader
to load sequence data from the files.

=cut

###
# insert_sequence()
#
# insert a bit of primary sequence into the database
#
sub insert_sequence {
  my $self = shift;
  my ($seqid,$seq,$offset) = @_;
  $offset ||= 0;
  $self->_insert_sequence($seqid,$seq,$offset);
}


=head2 fetch_sequence

 Title   : fetch_sequence
 Usage   : $sequence = $db->fetch_sequence(-seq_id=>$seqid,-start=>$start,-end=>$end)
 Function: Fetch the indicated subsequene from the database
 Returns : The sequence string (not a Bio::PrimarySeq object!)
 Args    : see below
 Status  : public

This method retrieves a portion of the indicated sequence. The arguments are:

  Argument       Value
  --------       -----
  -seq_id        Chromosome, contig or other DNA segment
  -seqid         Synonym for -seq_id
  -name          Synonym for -seq_id
  -start         Start of range
  -end           End of range
  -class         Obsolete argument used for Bio::DB::GFF compatibility. If
                  specified will qualify the seq_id as "$class:$seq_id".
  -bioseq        Boolean flag; if true, returns a Bio::PrimarySeq object instead
                  of a sequence string.

You can call fetch_sequence using the following shortcuts:

 $seq = $db->fetch_sequence('chr3');  # entire chromosome
 $seq = $db->fetch_sequence('chr3',1000);        # position 1000 to end of chromosome
 $seq = $db->fetch_sequence('chr3',undef,5000);  # position 1 to 5000
 $seq = $db->fetch_sequence('chr3',1000,5000);   # positions 1000 to 5000

=cut

###
# fetch_sequence()
#
# equivalent to old Bio::DB::GFF->dna() method
#
sub fetch_sequence {
  my $self = shift;
  my ($seqid,$start,$end,$class,$bioseq) = rearrange([['NAME','SEQID','SEQ_ID'],
						      'START',['END','STOP'],'CLASS','BIOSEQ'],@_);
  $seqid = "$seqid:$class" if defined $class;
  my $seq = $self->seq($seqid,$start,$end);
  return $seq unless $bioseq;

  require Bio::Seq unless Bio::Seq->can('new');
  my $display_id = defined $start ? "$seqid:$start..$end" : $seqid;
  return Bio::Seq->new(-display_id=>$display_id,-seq=>$seq);
}

=head2 segment

 Title   : segment
 Usage   : $segment = $db->segment($seq_id [,$start] [,$end] [,$absolute])
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

The named feature should exist once and only once in the database. If
it exists multiple times in the database and you attempt to call
segment() in a scalar context, you will get an exception. A workaround
is to call the method in a list context, as in:

  my ($segment) = $db->segment('contig23',1,1000);

or

  my @segments  = $db->segment('contig23',1,1000);

However, having multiple same-named features in the database is often
an indication of underlying data problems.

If the optional $absolute argument is a true value, then the specified
coordinates are relative to the reference (absolute) coordinates.

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
  my ($rel_start,$rel_end,$abs) = rearrange(['START',['STOP','END'],'ABSOLUTE'],@args);
  $rel_start = 1 unless defined $rel_start;

  my @segments;
  for my $f (@features) {
    my $seqid  = $f->seq_id;
    my $strand = $f->strand;
    my ($start,$end);
    if ($abs) {
      $start = $rel_start;
      $end   = defined $rel_end ? $rel_end : $start + $f->length - 1;
    }
    else {
      my $re = defined $rel_end ? $rel_end : $f->end - $f->start + 1;

      if ($strand >= 0) {
	$start = $f->start + $rel_start - 1;
	$end   = $f->start + $re   - 1;
      }
      else {
	$start = $f->end - $re   + 1;
	$end   = $f->end - $rel_start + 1;
      }
    }
    my $id = eval{$f->primary_id};
    push @segments,Bio::DB::SeqFeature::Segment->new($self,$seqid,$start,$end,$strand,$id);
  }
  return wantarray ? @segments : $segments[0];
}

=head2 seqfeature_class

 Title   : seqfeature_class
 Usage   : $classname = $db->seqfeature_class([$new_classname])
 Function: get or set the name of the Bio::SeqFeatureI class generated by new_feature()
 Returns : name of class
 Args    : new classname (optional)
 Status  : public

=cut

sub seqfeature_class {
  my $self = shift;
  my $d = $self->{seqfeatureclass};
  if (@_) {
    my $class = shift;
    eval "require $class";
    $self->throw("$class does not implement the Bio::SeqFeatureI interface")
      unless $class->isa('Bio::SeqFeatureI');
    $self->{seqfeatureclass} = $class;
  }
  $d;
}

=head2 reindex

 Title   : reindex
 Usage   : $db->reindex
 Function: reindex the database
 Returns : nothing
 Args    : nothing
 Status  : public

This method will force the secondary indexes (name, location,
attributes, feature types) to be recalculated. It may be useful to
rebuild a corrupted database.

=cut

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

=head2 attributes

 Title   : attributes
 Usage   : @a = $db->attributes
 Function: Returns list of all known attributes
 Returns : Returns list of all known attributes
 Args    : nothing
 Status  : public

=cut

sub attributes {
    my $self = shift;
    shift->throw_not_implemented;
}


=head2 start_bulk_update,finish_bulk_update

 Title   : start_bulk_update,finish_bulk_update
 Usage   : $db->start_bulk_update
           $db->finish_bulk_update
 Function: Activate optimizations for large number of insertions/updates
 Returns : nothing
 Args    : nothing
 Status  : public

With some adaptors (currently only the DBI::mysql adaptor), these
methods signal the adaptor that a large number of insertions or
updates are to be performed, and activate certain optimizations. These
methods are called automatically by the
Bio::DB::SeqFeature::Store::GFF3Loader module.

Example:

  $db->start_bulk_update;
  for my $f (@features) {
    $db->store($f);
  }
  $db->finish_bulk_update;

=cut

sub start_bulk_update  { shift->_start_bulk_update(@_) }
sub finish_bulk_update { shift->_finish_bulk_update(@_) }

=head2 add_SeqFeature

 Title   : add_SeqFeature
 Usage   : $count = $db->add_SeqFeature($parent,@children)
 Function: store a parent/child relationship between a $parent and @children
           features that are already stored in the database
 Returns : number of children successfully stored
 Args    : parent feature or primary ID and children features or primary IDs
 Status  : OPTIONAL; MAY BE IMPLEMENTED BY ADAPTORS

If can_store_parentage() returns true, then some store-aware features
(e.g. Bio::DB::SeqFeature) will invoke this method to store
feature/subfeature relationships in a normalized table.

=cut

# these two are called only if _can_store_subFeatures() returns true
# _add_SeqFeature ($parent,@children)
sub add_SeqFeature  { shift->_add_SeqFeature(@_)   }

=head2 fetch_SeqFeatures

 Title   : fetch_SeqFeatures
 Usage   : @children = $db->fetch_SeqFeatures($parent_feature)
 Function: return the immediate subfeatures of the indicated feature
 Returns : list of subfeatures
 Args    : the parent feature and an optional list of children types
 Status  : OPTIONAL; MAY BE IMPLEMENTED BY ADAPTORS

If can_store_parentage() returns true, then some store-aware features
(e.g. Bio::DB::SeqFeature) will invoke this method to retrieve
feature/subfeature relationships from the database.

=cut

# _get_SeqFeatures($parent,@child_types)
sub fetch_SeqFeatures {
  my ($self, $parent, @child_types) = @_;
  return unless defined $parent->primary_id;
  $self->_fetch_SeqFeatures($parent,@child_types);
}



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
  }
  $d;
}

=head2 dna_accessor

 Title   : dna_accessor
 Usage   : $dna_accessor = $db->dna_accessor([$new_dna_accessor])
 Function: get/set the name of the dna_accessor
 Returns : the current dna_accessor object, if any
 Args    : (optional) the dna_accessor object
 Status  : public

You can use this method to request or set the DNA accessor.

=cut

###
# dna_accessor
#
sub dna_accessor {
  my $self = shift;
  my $d    = $self->{dna_accessor};
  $self->{dna_accessor} = shift if @_;
  $d;
}

sub can_do_seq {
    my $self = shift;
    my $obj  = shift;
    return
	UNIVERSAL::can($obj,'seq') ||
	UNIVERSAL::can($obj,'fetch_sequence');
}

sub set_dna_accessor {
    my $self = shift;
    my $accessor = shift;
    if (-e $accessor) {  # a file, assume it is a fasta file
	eval "require Bio::DB::Fasta" unless Bio::DB::Fasta->can('new');
	my $a = Bio::DB::Fasta->new($accessor)
	    or croak "Can't open FASTA file $accessor: $!";
	$self->dna_accessor($a);
    }

    if (ref $accessor && $self->can_do_seq($accessor)) {
	$self->dna_accessor($accessor);  # already built
    }

    return;
}

sub do_compress {
  my $self = shift;
  if (@_) {
    my $do_compress = shift;
    $self->setting(compress => $do_compress);
  }
  my $d    = $self->setting('compress');
  if ($d) {
    eval "use Compress::Zlib; 1" or croak $@ unless Compress::Zlib->can('compress');
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
functions to return features at any level of the containment
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

=head2 clone

The clone() method should be used when you want to pass the
Bio::DB::SeqFeature::Store object to a child process across a
fork(). The child must call clone() before making any queries.

The default behavior is to do nothing, but adaptors that use the DBI
interface may need to implement this in order to avoid database handle
errors. See the dbi adaptor for an example.

=cut

sub clone { }

################################# TIE interface ####################

=head1 TIE Interface

This module implements a full TIEHASH interface. The keys are the
primary IDs of the features in the database. Example:

 tie %h,'Bio::DB::SeqFeature::Store',-adaptor=>'DBI::mysql',-dsn=>'dbi:mysql:elegans';
 $h{123} = $feature1;
 $h{124} = $feature2;
 print $h{123}->display_name;

=cut

sub TIEHASH {
  my $class = shift;
  return $class->new(@_);
}

sub STORE {
  my $self = shift;
  my ($key,$feature) = @_;
  $key =~ /^\d+$/ && $key > 0 or croak "keys must be positive integers";
  $self->load_class($feature);
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


###################### TO BE IMPLEMENTED BY ADAPTOR ##########

=head2 _init_database

 Title   : _init_database
 Usage   : $success = $db->_init_database([$erase])
 Function: initialize an empty database
 Returns : true on success
 Args    : optional boolean flag to erase contents of an existing database
 Status  : ABSTRACT METHOD; MUST BE IMPLEMENTED BY AN ADAPTOR

This method is the back end for init_database(). It must be
implemented by an adaptor that inherits from
Bio::DB::SeqFeature::Store. It returns true on success. @features = $db->features(-seqid=>'Chr1');

=cut

sub _init_database { shift->throw_not_implemented }

=head2 _store

 Title   : _store
 Usage   : $success = $db->_store($indexed,@objects)
 Function: store seqfeature objects into database
 Returns : true on success
 Args    : a boolean flag indicating whether objects are to be indexed,
           and one or more objects
 Status  : ABSTRACT METHOD; MUST BE IMPLEMENTED BY AN ADAPTOR

This method is the back end for store() and store_noindex(). It should
write the seqfeature objects into the database. If indexing is
requested, the features should be indexed for query and
retrieval. Otherwise the features should be stored without indexing
(it is not required that adaptors respect this).

If the object has no primary_id (undef), then the object is written
into the database and assigned a new primary_id. If the object already
has a primary_id, then the system will perform an update, replacing
whatever was there before.

In practice, the implementation will serialize each object using the
freeze() method and then store it in the database under the
corresponding primary_id. The object is then updated with the
primary_id.

=cut

# _store($indexed,@objs)
sub _store {
  my $self    = shift;
  my $indexed = shift;
  my @objs    = @_;
  $self->throw_not_implemented;
}

=head2 _fetch

 Title   : _fetch
 Usage   : $feature = $db->_fetch($primary_id)
 Function: fetch feature from database
 Returns : feature
 Args    : primary id
 Status  : ABSTRACT METHOD; MUST BE IMPLEMENTED BY AN ADAPTOR

This method is the back end for fetch(). It accepts a primary_id and
returns a feature object. It must be implemented by the adaptor.

In practice, the implementation will retrieve the serialized
Bio::SeqfeatureI object from the database and pass it to the thaw()
method to unserialize it and synchronize the primary_id.

=cut

# _fetch($id)
sub _fetch { shift->throw_not_implemented }

=head2 _fetch_many

 Title   : _fetch_many
 Usage   : $feature = $db->_fetch_many(@primary_ids)
 Function: fetch many features from database
 Returns : feature
 Args    : primary id
 Status  : private -- does not need to be implemented

This method fetches many features specified by a list of IDs. The
default implementation simply calls _fetch() once for each
primary_id. Implementors can override it if needed for efficiency.

=cut

# _fetch_many(@ids)
# this one will fall back to many calls on fetch() if you don't
# override it
sub _fetch_many {
  my $self = shift;
  return map {$self->_fetch($_)} @_;
}

=head2 _update_indexes

 Title   : _update_indexes
 Usage   : $success = $db->_update_indexes($feature)
 Function: update the indexes for a feature
 Returns : true on success
 Args    : A seqfeature object
 Status  : ABSTRACT METHOD; MUST BE IMPLEMENTED BY AN ADAPTOR

This method is called by reindex() to update the searchable indexes
for a feature object that has changed.

=cut

# this is called to index a feature
sub _update_indexes { shift->throw_not_implemented }

=head2 _start_reindexing, _end_reindexing

 Title   : _start_reindexing, _end_reindexing
 Usage   : $db->_start_reindexing()
           $db->_end_reindexing
 Function: flag that a series of reindexing operations is beginning/ending
 Returns : true on success
 Args    : none
 Status  : MAY BE IMPLEMENTED BY AN ADAPTOR (optional)

These methods are called by reindex() before and immediately after a
series of reindexing operations. The default behavior is to do
nothing, but these methods can be overridden by an adaptor in order to
perform optimizations, turn off autocommits, etc.

=cut

# these do not necessary have to be overridden
# they are called at beginning and end of reindexing process
sub _start_reindexing {}
sub _end_reindexing   {}

=head2 _features

 Title   : _features
 Usage   : @features = $db->_features(@args)
 Function: back end for all get_feature_by_*() queries
 Returns : list of features
 Args    : see below
 Status  : ABSTRACT METHOD; MUST BE IMPLEMENTED BY ADAPTOR

This is the backend for features(), get_features_by_name(),
get_features_by_location(), etc. Arguments are as described for the
features() method, except that only the named-argument form is
recognized.

=cut

# bottleneck query generator
sub _features { shift->throw_not_implemented }

=head2 _search_attributes

 Title   : _search_attributes
 Usage   : @result_list = $db->_search_attributes("text search string",[$tag1,$tag2...],$limit)
 Function: back end for the search_attributes() method
 Returns : results list
 Args    : as per search_attributes()
 Status  : ABSTRACT METHOD; MUST BE IMPLEMENTED BY ADAPTOR

See search_attributes() for the format of the results list. The only
difference between this and the public method is that the tag list is
guaranteed to be an array reference.

=cut

sub _search_attributes { shift->throw_not_implemented }

=head2 can_store_parentage

 Title   : can_store_parentage
 Usage   : $flag = $db->can_store_parentage
 Function: return true if this adaptor can store parent/child relationships
 Returns : boolean
 Args    : none
 Status  : OPTIONAL; MAY BE IMPLEMENTED BY ADAPTORS

Override this method and return true if this adaptor supports the
_add_SeqFeature() and _get_SeqFeatures() methods, which are used for
storing feature parent/child relationships in a normalized
fashion. Default is false (parent/child relationships are stored in
denormalized form in each feature).

=cut

# return true here if the storage engine is prepared to store parent/child
# relationships using _add_SeqFeature and return them using _fetch_SeqFeatures
sub can_store_parentage { return; }

=head2 _add_SeqFeature

 Title   : _add_SeqFeature
 Usage   : $count = $db->_add_SeqFeature($parent,@children)
 Function: store a parent/child relationship between $parent and @children
 Returns : number of children successfully stored
 Args    : parent feature and one or more children
 Status  : OPTIONAL; MAY BE IMPLEMENTED BY ADAPTORS

If can_store_parentage() returns true, then some store-aware features
(e.g. Bio::DB::SeqFeature) will invoke this method to store
feature/subfeature relationships in a normalized table.

=cut

sub _add_SeqFeature { shift->throw_not_implemented }

=head2 _fetch_SeqFeatures

 Title   : _fetch_SeqFeatures
 Usage   : @children = $db->_fetch_SeqFeatures($parent_feature)
 Function: return the immediate subfeatures of the indicated feature
 Returns : list of subfeatures
 Args    : the parent feature
 Status  : OPTIONAL; MAY BE IMPLEMENTED BY ADAPTORS

If can_store_parentage() returns true, then some store-aware features
(e.g. Bio::DB::SeqFeature) will invoke this method to retrieve
feature/subfeature relationships from the database.

=cut

# _get_SeqFeatures($parent,@list_of_child_types)
sub _fetch_SeqFeatures {shift->throw_not_implemented }

=head2 _insert_sequence

 Title   : _insert_sequence
 Usage   : $success = $db->_insert_sequence($seqid,$sequence_string,$offset)
 Function: Inserts sequence data into the database at the indicated offset
 Returns : true if successful
 Args    : see below
 Status  : ABSTRACT METHOD; MUST BE IMPLEMENTED BY ADAPTOR

This is the back end for insert_sequence(). Adaptors must implement
this method in order to store and retrieve nucleotide or protein
sequence.

=cut

sub _insert_sequence   { shift->throw_not_implemented }

# _fetch_sequence() is similar to old dna() method

=head2 _fetch_sequence

 Title   : _fetch_sequence
 Usage   : $sequence = $db->_fetch_sequence(-seq_id=>$seqid,-start=>$start,-end=>$end)
 Function: Fetch the indicated subsequence from the database
 Returns : The sequence string (not a Bio::PrimarySeq object!)
 Args    : see below
 Status  : ABSTRACT METHOD; MUST BE IMPLEMENTED BY ADAPTOR

This is the back end for fetch_sequence(). Adaptors must implement
this method in order to store and retrieve nucleotide or protein
sequence.

=cut

sub _fetch_sequence    { shift->throw_not_implemented }

sub seq {
    my $self     = shift;
    my ($seq_id,$start,$end) = @_;
    if (my $a = $self->dna_accessor) {
	return $a->can('seq')           ? $a->seq($seq_id,$start,$end)
	      :$a->can('fetch_sequence')? $a->fetch_sequence($seq_id,$start,$end)
          : undef;
    }
    else {
	return $self->_fetch_sequence($seq_id,$start,$end);
    }
}

=head2 _seq_ids

 Title   : _seq_ids
 Usage   : @ids = $db->_seq_ids()
 Function: Return all sequence IDs contained in database
 Returns : list of sequence Ids
 Args    : none
 Status  : TO BE IMPLEMENTED BY ADAPTOR

This method is invoked by seq_ids() to return all sequence IDs
(coordinate systems) known to the database.

=cut

sub _seq_ids { shift->throw_not_implemented }

=head2 _start_bulk_update,_finish_bulk_update

 Title   : _start_bulk_update, _finish_bulk_update
 Usage   : $db->_start_bulk_update
           $db->_finish_bulk_update
 Function: Activate optimizations for large number of insertions/updates
 Returns : nothing
 Args    : nothing
 Status  : OPTIONAL; MAY BE IMPLEMENTED BY ADAPTOR

These are the backends for start_bulk_update() and
finish_bulk_update(). The default behavior of both methods is to do
nothing.

=cut

# Optional flags to change behavior to optimize bulk updating.
sub _start_bulk_update { }
sub _finish_bulk_update { }


# for full TIE() interface  - not necessary to implement in most cases

=head2 Optional methods needed to implement full TIEHASH interface

The core TIEHASH interface will work if just the _store() and _fetch()
methods are implemented. To support the full TIEHASH interface,
including support for keys(), each(), and exists(), the following
methods should be implemented:

=over 4

=item $id = $db-E<gt>_firstid()

Return the first primary ID in the database. Needed for the each()
function.

=item $next_id = $db-E<gt>_nextid($id)

Given a primary ID, return the next primary ID in the series. Needed
for the each() function.

=item $boolean = $db-E<gt>_existsid($id)

Returns true if the indicated primary ID is in the database. Needed
for the exists() function.

=item $db-E<gt>_deleteid($id)

Delete the feature corresponding to the given primary ID. Needed for
delete().

=item $db-E<gt>_clearall()

Empty the database. Needed for %tied_hash = ().

=item $count = $db-E<gt>_featurecount()

Return the number of features in the database. Needed for scalar
%tied_hash.

=back

=cut

sub _firstid  { shift->throw_not_implemented }
sub _nextid   { shift->throw_not_implemented }
sub _existsid { shift->throw_not_implemented }
sub _deleteid { shift->throw_not_implemented }
sub _clearall { shift->throw_not_implemented }
sub _featurecount { shift->throw_not_implemented }


=head1 Internal Methods

These methods are internal to Bio::DB::SeqFeature::Store and adaptors.

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

This method is used internally by the Bio::DB::SeqFeature class to
optimize some of its operations. It returns true if all of the
subfeatures in the database are indexed; it returns false if at least
one of the subfeatures is not indexed. Do not attempt to change the
value of this setting unless you are writing an adaptor.

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

=head2 subfeature_types_are_indexed

 Title   : subfeature_types_are_indexed
 Usage   : $flag = $db->subfeature_types_are_indexed
 Function: whether subfeatures are indexed by type
 Returns : a flag indicating that all subfeatures are indexed
 Args    : none
 Status  : private

This method returns true if subfeature types are indexed. Default is
to return the value of subfeatures_are_indexed().

=cut

sub subfeature_types_are_indexed {
  my $self = shift;
  return $self->subfeatures_are_indexed;
}

=head2 subfeature_locations_are_indexed

 Title   : subfeature_locations_are_indexed
 Usage   : $flag = $db->subfeature_locations_are_indexed
 Function: whether subfeatures are indexed by type
 Returns : a flag indicating that all subfeatures are indexed
 Args    : none
 Status  : private

This method returns true if subfeature locations are indexed. Default is
to return the value of subfeatures_are_indexed().

=cut

sub subfeature_locations_are_indexed {
  my $self = shift;
  return $self->subfeatures_are_indexed;
}

=head2 setup_segment_args

 Title   : setup_segment_args
 Usage   : @args = $db->setup_segment_args(@args)
 Function: munge the arguments to the segment() call
 Returns : munged arguments
 Args    : see below
 Status  : private

This method is used internally by segment() to translate positional
arguments into named argument=E<gt>value pairs.

=cut

sub setup_segment_args {
  my $self = shift;
  return @_ if defined $_[0] && $_[0] =~ /^-/;
  return (-name=>$_[0],-start=>$_[1],-end=>$_[2]) if @_ == 3;
  return (-class=>$_[0],-name=>$_[1])              if @_ == 2;
  return (-name=>$_[0])                            if @_ == 1;
  return;
}

=head2 store_and_cache

 Title   : store_and_cache
 Usage   : $success = $db->store_and_cache(@features)
 Function: store features into database and update cache
 Returns : number of features stored
 Args    : index the features? (0 or 1) and  list of features
 Status  : private

This private method stores the list of Bio::SeqFeatureI objects into
the database and caches them in memory for retrieval.

=cut

sub store_and_cache {
  my $self = shift;
  my $indexit = shift;
  my $result = $self->_store($indexit,@_);
  if (my $cache = $self->cache) {
    for my $obj (@_) {
      defined (my $id     = eval {$obj->primary_id}) or next;
      $cache->store($id,$obj);
    }
  }
  $result;
}

=head2 init_cache

 Title   : init_cache
 Usage   : $db->init_cache($size)
 Function: initialize the in-memory feature cache
 Returns : the Tie::Cacher object
 Args    : desired size of the cache
 Status  : private

This method is used internally by new() to create the Tie::Cacher
instance used for the in-memory feature cache.

=cut

sub init_cache {
  my $self       = shift;
  my $cache_size = shift;
  $cache_size    = 5000 if $cache_size == 1;   # in case somebody treats it as a flag
  $self->{cache} = Tie::Cacher->new($cache_size) or $self->throw("Couldn't tie cache: $!");
}

=head2 cache

 Title   : cache
 Usage   : $cache = $db->cache
 Function: return the cache object
 Returns : the Tie::Cacher object
 Args    : none
 Status  : private

This method returns the Tie::Cacher object used for the in-memory
feature cache.

=cut

sub cache { shift->{cache} }

=head2 load_class

 Title   : load_class
 Usage   : $db->load_class($blessed_object)
 Function: loads the module corresponding to a blessed object
 Returns : empty
 Args    : a blessed object
 Status  : private

This method is used by thaw() to load the code for a blessed
object. This ensures that all the object's methods are available.

=cut

sub load_class {
  my $self = shift;
  my $obj  = shift;
  return unless defined $obj;
  return if $self->{class_loaded}{ref $obj}++;
  unless ($obj && $obj->can('primary_id')) {
    my $class = ref $obj;
    eval "require $class";
  }
}


#################################### Internal methods ####################

=head2 freeze

 Title   : freeze
 Usage   : $serialized_object = $db->freeze($feature)
 Function: serialize a feature object into a string
 Returns : serialized feature object
 Args    : a seqfeature object
 Status  : private

This method converts a Bio::SeqFeatureI object into a serialized form
suitable for storage into a database. The feature's primary ID is set
to undef before it is serialized. This avoids any potential mismatch
between the primary ID used as the database key and the primary ID
stored in the serialized object.

=cut

sub freeze {
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
    $d->Deparse(1);
    $data = $d->Dump;
  } elsif ($serializer eq 'Storable') {
    local $Storable::forgive_me = 1;
    local $Storable::Deparse = 1;
    $data = Storable::nfreeze($obj);
  }

  $obj->primary_id($id);       # restore to original state
  eval {
    $obj->object_store($store);
  };

  $data = compress($data) if $self->do_compress;
  return $data;
}

=head2 thaw

 Title   : thaw
 Usage   : $feature = $db->thaw($serialized_object,$primary_id)
 Function: unserialize a string into a feature object
 Returns : Bio::SeqFeatureI object
 Args    : serialized form of object from freeze() and primary_id of object
 Status  : private

This method is the reverse of the freeze(). The supplied primary_id
becomes the primary_id() of the returned Bio::SeqFeatureI object. This
implementation checks for a deserialized object in the cache before it
calls thaw_object() to do the actual deserialization.

=cut

sub thaw {
  my $self               = shift;
  my ($obj,$primary_id)  = @_;

  if (my $cache = $self->cache) {
    return $cache->fetch($primary_id) if $cache->exists($primary_id);
    my $object = $self->thaw_object($obj,$primary_id) or return;
    $cache->store($primary_id,$object);
    return $object;
  } else {
    return $self->thaw_object($obj,$primary_id);
  }

}

=head2 thaw_object

 Title   : thaw_object
 Usage   : $feature = $db->thaw_object($serialized_object,$primary_id)
 Function: unserialize a string into a feature object
 Returns : Bio::SeqFeatureI object
 Args    : serialized form of object from freeze() and primary_id of object
 Status  : private

After thaw() checks the cache and comes up empty, this method is
invoked to thaw the object.

=cut

sub thaw_object {
  my $self               = shift;
  my ($obj,$primary_id)  = @_;

  my $serializer = $self->serializer;
  my $object;

  $obj = uncompress($obj) if $self->do_compress;

  if ($serializer eq 'Data::Dumper') {
    $object = eval $obj;
  } elsif ($serializer eq 'Storable') {
    local $Storable::forgive_me = 1;
    local $Storable::Eval = 1;
    $object = Storable::thaw($obj);
  }

  # remember the primary ID of this object as well as the
  # identity of the store, so that we can do lazy loading;
  # both of these are wrapped in an eval because not all
  # bioseqfeatures support them (or want to)
  $self->load_class($object);
  eval {
    $object->primary_id($primary_id);
    $object->object_store($self);
  };
  $object;
}

=head2 feature_names

 Title   : feature_names
 Usage   : ($names,$aliases) = $db->feature_names($feature)
 Function: get names and aliases for a feature
 Returns : an array of names and an array of aliases
 Args    : a Bio::SeqFeatureI object
 Status  : private

This is an internal utility function which, given a Bio::SeqFeatureI
object, returns two array refs. The first is a list of official names
for the feature, and the second is a list of aliases. This is slightly
skewed towards GFF3 usage, so the official names are the
display_name(), plus all tag values named 'Name', plus all tag values
named 'ID'. The aliases are all tag values named 'Alias'.

=cut

sub feature_names {
  my $self = shift;
  my $obj  = shift;

  my $primary_id = $obj->primary_id;
  my @names;
  push @names,$obj->display_name           if defined $obj->display_name;
  push @names,$obj->get_tag_values('Name') if $obj->has_tag('Name');
  push @names,$obj->get_tag_values('ID')   if $obj->has_tag('ID');

  # don't think this is desired behavior
  # @names = grep {defined $_ && $_ ne $primary_id} @names;

  my @aliases = grep {defined} $obj->get_tag_values('Alias') if $obj->has_tag('Alias');

  return (\@names,\@aliases);
}

=head2 feature_summary

 Title   : feature_summary
 Usage   : $summary = $db->feature_summary(@args)
 Function: returns a coverage summary across indicated region/type
 Returns : a Bio::SeqFeatureI object containing the "coverage" tag
 Args    : see below
 Status  : public

This method is used to get coverage density information across a
region of interest. You provide it with a region of interest, optional
a list of feature types, and a count of the number of bins over which
you want to calculate the coverage density. An object is returned
corresponding to the requested region. It contains a tag called
"coverage" that will return an array ref of "bins" length. Each
element of the array describes the number of features that overlap the
bin at this postion.

Arguments:

  Argument       Description
  --------       -----------

  -seq_id        Sequence ID for the region
  -start         Start of region
  -end           End of region
  -type/-types   Feature type of interest or array ref of types
  -bins          Number of bins across region. Defaults to 1000.
  -iterator      Return an iterator across the region

Note that this method uses an approximate algorithm that is only
accurate to 500 bp, so when dealing with bins that are smaller than
1000 bp, you may see some shifting of counts between adjacent bins.

Although an -iterator option is provided, the method only ever returns
a single feature, so this is fairly useless.

=cut


sub feature_summary {
    my $self = shift;
    my ($seq_name,$start,$end,$types,$bins,$iterator) =
	rearrange([['SEQID','SEQ_ID','REF'],'START',['STOP','END'],
		   ['TYPES','TYPE','PRIMARY_TAG'],
		   'BINS',
		   'ITERATOR',
		  ],@_);
    my ($coverage,$tag) = $self->coverage_array(-seqid=> $seq_name,
						-start=> $start,
						-end  => $end,
						-type => $types,
						-bins => $bins) or return;
    my $score = 0;
    for (@$coverage) { $score += $_ }
    $score /= @$coverage;

    my $feature = Bio::SeqFeature::Lite->new(-seq_id => $seq_name,
					     -start  => $start,
					     -end    => $end,
					     -type   => $tag,
					     -score  => $score,
					     -attributes =>
					     { coverage => [$coverage] });
    return $iterator
	   ? Bio::DB::SeqFeature::Store::FeatureIterator->new($feature)
	   : $feature;
}

=head2 coverage_array

 Title   : coverage_array
 Usage   : $arrayref = $db->coverage_array(@args)
 Function: returns a coverage summary across indicated region/type
 Returns : an array reference
 Args    : see below
 Status  : public

This method is used to get coverage density information across a
region of interest. The arguments are identical to feature_summary,
except that instead of returning a Bio::SeqFeatureI object, it returns
an array reference of the desired number of bins. The value of each
element corresponds to the number of features in the bin.

Arguments:

  Argument       Description
  --------       -----------

  -seq_id        Sequence ID for the region
  -start         Start of region
  -end           End of region
  -type/-types   Feature type of interest or array ref of types
  -bins          Number of bins across region. Defaults to 1000.

Note that this method uses an approximate algorithm that is only
accurate to 500 bp, so when dealing with bins that are smaller than
1000 bp, you may see some shifting of counts between adjacent bins.

=cut

sub coverage_array {
    shift->throw_not_implemented;
}


package Bio::DB::SeqFeature::Store::FeatureIterator;

sub new {
    my $self     = shift;
    my @features = @_;
    return bless \@features,ref $self || $self;
}
sub next_seq {
  my $self  = shift;
  return unless @$self;
  return shift @$self;
}

sub begin_work { }# noop
sub commit     { }# noop
sub rollback   { }# noop

1;

__END__

=head1 BUGS

This is an early version, so there are certainly some bugs. Please
use the BioPerl bug tracking system to report bugs.

=head1 SEE ALSO

L<Bio::DB::SeqFeature>,
L<Bio::DB::SeqFeature::Store::GFF3Loader>,
L<Bio::DB::SeqFeature::Segment>,
L<Bio::DB::SeqFeature::Store::DBI::mysql>,
L<Bio::DB::SeqFeature::Store::berkeleydb>
L<Bio::DB::SeqFeature::Store::memory>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2006 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut
