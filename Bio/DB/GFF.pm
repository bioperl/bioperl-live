=head1 NAME

Bio::DB::GFF -- Storage and retrieval of sequence annotation data

=head1 SYNOPSIS

  use Bio::DB::GFF;

  # Open the sequence database
  my $db      = Bio::DB::GFF->new( -adaptor => 'dbi:mysql',
                                   -dsn     => 'dbi:mysql:elegans42');

  # fetch a 1 megabase segment of sequence starting at landmark "ZK909"
  my $segment = $db->segment('ZK909', 1 => 1000000);

  # pull out all transcript features
  my @transcripts = $segment->features('transcript');

  # for each transcript, total the length of the introns
  my %totals;
  for my $t (@transcripts) {
    my @introns = $t->Intron;
    $totals{$t->name} += $_->length foreach @introns;
  }

  # Sort the exons of the first transcript by position
  my @exons = sort {$a->start <=> $b->start} $transcripts[0]->Exon;

  # Get a region 1000 bp upstream of first exon
  my $upstream = $exons[0]->segment(-1000,0);

  # get its DNA
  my $dna = $upstream->dna;

  # and get all curated polymorphisms inside it
  @polymorphisms = $upstream->contained_features('polymorphism:curated');

  # get all feature types in the database
  my @types = $db->types;

  # count all feature types in the segment
  my %type_counts = $segment->types(-enumerate=>1);

  # get an iterator on all curated features of type 'exon' or 'intron'
  my $iterator = $db->features(-type     => ['exon:curated','intron:curated'],
                               -iterator => 1);

  while ($_ = $iterator->next_feature) {
      print $_,"\n";
  }

=head1 DESCRIPTION

Bio::DB::GFF provides fast indexed access to a sequence annotation
database.  It supports multiple database types (ACeDB, relational),
and multiple schemas through a system of adaptors and aggregators.

The following operations are supported by this module:

  - retrieving a segment of sequence based on the ID of a landmark
  - retrieving the DNA from that segment
  - finding all annotations that overlap with the segment
  - finding all annotations that are completely contained within the
    segment
  - retrieving all annotations of a particular type, either within a
    segment, or globally
  - conversion from absolute to relative coordinates and back again,
    using any arbitrary landmark for the relative coordinates
  - using a sequence segment to creatie new segments based on relative 
    offsets

The data model used by Bio::DB::GFF is compatible with the GFF flat
file format (http://www.sanger.ac.uk/software/GFF).  The module can
load a set of GFF files into the database, and serves objects that
have methods corresponding to GFF fields.

The objects returned by Bio::DB::GFF are compatible with the
SeqFeatureI interface, allowing their use by the Bio::Graphics and
Bio::DAS modules.

=head2 GFF Fundamentals

The GFF format is a flat tab-delimited file, each line of which
corresponds to an annotation, or feature.  Each annotation has the
following attributes:

=over 4

=item reference sequence

This is the ID of the sequence that is used to establish the
coordinate system of the annotation.

=item start position

The start of the annotation relative to the reference sequence.  

=item stop position

The stop of the annotation relative to the reference sequence.  Start
is always less than or equal to stop.

=item method

The annotation method.  This field describes the general nature of the
annotation, such as "tRNA".

=item source

The source of the annotation.  This field describes how the annotation
was derived, such as "tRNAScanSE-1.2".  Together the method and source
describe the annotation type.

=item strand

For those annotations which are strand-specific, this field is the
strand on which the annotation resides.

=item phase

For annotations that are linked to proteins, this field describes the
phase of the annotation on the codons.

=item score

For annotations that are associated with a numeric score (for example,
a sequence similarity), this field describes the score.

=item group

GFF provides a simple way of generating annotation hierarchies ("is
composed of" relationships) by providing a group field.  The group
field contains the class and ID of an annotation which is the logical
parent of the current one.

The group field is also used to store information about the target of
sequence similarity hits, and miscellaneous notes.

=back

The sequences used to establish the coordinate system for annotations
can correspond to sequenced clones, clone fragments, contigs or
super-contigs.  Thus, this module can be used throughout the lifecycle
of a sequencing project.

In addition to a group ID, the GFF format allows annotations to have a
group class.  For example, in the ACeDB representation, RNA
interference experiments have a class of "RNAi" and an ID that is
unique among the RNAi experiments.  Since not all databases support
this notion, the class is optional in all calls to this module, and
defaults to "Sequence" when not provided.

This module has no relationship to the GFF.pm module.

=head2 Adaptors and Aggregators

This module uses a system of adaptors and aggregators in order to make
it adaptable to use with a variety of databases.

=over 4

=item Adaptors

The core of the module handles the user API, annotation coordinate
arithmetic, and other common issues.  The details of fetching
information from databases is handled by an adaptor, which is
specified during Bio::DB::GFF construction.  The adaptor encapsulates
database-specific information such as the schema, user authentication
and access methods.

Currently there are two adaptors: 'dbi:mysql' and 'dbi:mysqlopt'.  The 
former is an interface to a simple Mysql schema.  The latter is an
optimized version of dbi:mysql which uses a binning scheme to
accelerate range queries and the Bio::DB::Fasta module for rapid
retrieval of sequences.

=item Aggregators

The GFF format uses a "group" field to indicate aggregation properties
of individual features.  For example, a set of exons and introns may
share a common transcript group, and multiple transcripts may share
the same gene group.  

Aggregators are small modules that use the group information to
rebuild the hierarchy.  When a Bio::DB::GFF object is created, you
indicate that it use a set of one or more aggregators.  Each
aggregator provides a new composite annotation type.  Before the
database query is generated each aggregator is called to
"disaggregate" its annotation type into list of component types
contained in the database.  After the query is generated, each
aggregator is called again in order to build composite annotations
from the returned components.

For example, during disaggregation, the standard "transcript"
aggregator generates a list of component feature types including
"intron", "exon", "CDS", and "3'UTR".  Later, it aggregates these
features into a set of annotations of type "transcript".

During aggregation, the list of aggregators is called in reverse
order.  This allows aggregators to collaborate to create multi-level
structures: the transcript aggregator assembles transcripts from
introns and exons; the gene aggregator then assembles genes from sets
of transcripts.

Three aggregators are currently provided:

      - transcript   assembles transcripts
      - clone        assembles clones from Clone_end features
      - alignment    assembles gapped alignments from similarity
		     features

The existing aggregators are easily customized.

=back

=head1 API

The following is the API for Bio::DB::GFF.

=cut

package Bio::DB::GFF;

use strict;

use Bio::DB::GFF::Util::Rearrange;
use Bio::DB::GFF::RelSegment;
use Bio::DB::GFF::Feature;
use Bio::Root::RootI;

use vars qw($VERSION @ISA);
@ISA = qw(Bio::Root::RootI);

$VERSION = '0.30';

=head2 new

 Title   : new
 Usage   : my $db = new Bio::DB::GFF(@args);
 Function: create a new Bio::DB::GFF object
 Returns : new Bio::DB::GFF object
 Args    : lists of adaptors and aggregators
 Status  : Public

These are the arguments:

 -adaptor      Name of the adaptor module to use.  If none
               provided, defaults to "dbi:mysqlopt".

 -aggregator   Array reference to a list of aggregators
               to apply to the database.  If none provided,
	       defaults to ['transcript','clone','alignment'].

  <other>      Any other named argument pairs are passed to
               the adaptor for processing.

The adaptor argument must correspond to a module contained within the
Bio::DB::GFF::Adaptor namespace.  For example, the
Bio::DB::GFF::Adaptor::dbi::mysql adaptor is loaded by specifying
'dbi:mysql'.  By Perl convention, the adaptors names are lower case
because they are loaded at run time.

The aggregator array may contain a list of aggregator names, or a list 
of initialized aggregator objects.  For example, if you wish to change 
the components aggregated by the transcript aggregator, you could
pass it to the GFF constructor this way:

  my $transcript = 
     Bio::DB::Aggregator::transcript->new(-parts=>[qw(exon intron utr
                                                      polyA spliced_leader)]);
  my $db = Bio::DB::GFF->new(-aggregator=>[$transcript,'clone','alignment],
                             -adaptor   => 'dbi:mysql',
                              -dsn      => 'dbi:mysql:elegans42');

=cut

sub new {
  my $package   = shift;
  my ($adaptor,$aggregators,$args) = rearrange([
						[qw(ADAPTOR FACTORY)],
						[qw(AGGREGATOR AGGREGATORS)]
						],@_);

  $adaptor    ||= 'dbi::mysqlopt';
  my $class = "Bio::DB::GFF::Adaptor::\L${adaptor}\E";
  eval "require $class";
  $package->throw("Unable to load $adaptor adaptor: $@") if $@;

  my $self = $class->new($args);

  # handle the aggregators.
  # aggregators are responsible for creating complex multi-part features
  # from the GFF "group" field.  If none are provided, then we provide a
  # list of the two used in WormBase.
  # Each aggregator can be a scalar or a ref.  In the former case
  # it is treated as a class name to call new() on.  In the latter
  # the aggreator is treated as a ready made object.
  $aggregators = $self->default_aggregators unless defined $aggregators;
  my @a = ref($aggregators) ? @$aggregators : $aggregators;
  my @aggregators;
  for my $a (@a) {
    $self->add_aggregator($a);
  }
  $self;
}

=head2 load

 Title   : load
 Usage   : $db->load($file|$directory|$filehandle);
 Function: load GFF data into database
 Returns : count of records loaded
 Args    : a directory, a file, a list of files, 
           or a filehandle
 Status  : Public

This method takes a single overloaded argument, which can be any of:

=over 4

=item 1. a scalar corresponding to a GFF file on the system

A pathname to a local GFF file.  Any files ending with the .gz, .Z, or
.bz2 suffixes will be transparently decompressed with the appropriate
command-line utility.

=item 2. an array reference containing a list of GFF files on the
system

For example ['/home/gff/gff1.gz','/home/gff/gff2.gz']

=item 3. path to a directory

The indicated directory will be searched for all files ending in the
suffixes .gz, .Z or .bz2.

=item 4. a filehandle

An open filehandle from which to read the GFF data.

=item 5. a pipe expression

A pipe expression will also work. For example, a GFF file on a remote
web server can be loaded with an expression like this:

  $db->load("lynx -dump -source http://stein.cshl.org/gff_test |");

=back

If successful, the method will return the number of GFF lines
successfully loaded.

=cut

sub load {
  my $self      = shift;
  my $file_or_directory = shift || '.';

  local @ARGV;  # to play tricks with reader

  if (-d $file_or_directory) {
    @ARGV = glob("$file_or_directory/*.{gff,gff.gz,gff.Z,gff,bz2}");
  } elsif (my $fd = fileno($file_or_directory)) {
    open SAVEIN,"<&STDIN";
    open STDIN,"<&=$fd" or $self->throw("Can't dup STDIN");
    @ARGV = '-';
  } elsif (ref $file_or_directory) {
    @ARGV = @$file_or_directory;
  } else {
    @ARGV = $file_or_directory;
  }

  return unless @ARGV;
  foreach (@ARGV) {
    if (/\.gz$/) {
      $_ = "gunzip -c $_ |";
    } elsif (/\.Z$/) {
      $_ = "uncompress -c $_ |";
    } elsif (/\.bz2$/) {
      $_ = "bunzip2 -c $_ |";
    }
  }

  my $result = $self->load_gff;

  open STDIN,"<&SAVEIN";  # restore STDIN
  return $result;
}

=head2 initialize

 Title   : initialize
 Usage   : $db->initialize($erase);
 Function: initialize a GFF database
 Returns : true if initialization successful
 Args    : an optional flag indicating that existing
           contents should be wiped clean
 Status  : Public

This method can be used to initialize an empty database.  It will not
overwrite existing data unless a true $erase flag is present.

=cut

sub initialize {
    shift->do_initialize(@_);
}

=head2 error

 Title   : error
 Usage   : $db->error( [$new error] );
 Function: read or set error message
 Returns : error message
 Args    : an optional argument to set the error message
 Status  : Public

This method can be used to retrieve the last error message.  Errors
are not reset to empty by successful calls, so contents are only valid
immediately after an error condition has been detected.

=cut

sub error {
  my $self = shift;
  my $g = $self->{error};
  $self->{error} = shift if @_;
  $g;
}

=head2 debug

 Title   : debug
 Usage   : $db->debug( [$flag] );
 Function: read or set debug flag
 Returns : current value of debug flag
 Args    : new debug flag (optional)
 Status  : Public

This method can be used to turn on debug messages.  The exact nature
of those messages depends on the adaptor in use.

=cut

sub debug {
  my $self = shift;
  my $g = $self->{debug};
  $self->{debug} = shift if @_;
  $g;
}

=head2 segment

 Title   : segment
 Usage   : $db->segment(@args);
 Function: create a segment object
 Returns : a segment object
 Args    : numerous, see below
 Status  : public

This method generates a segment object, which is a Perl object
subclassed from Bio::DB::GFF::Segment.  The segment can be used to
find overlapping features and the raw DNA.  

When making the segment() call, you specify the ID of a sequence
landmark (e.g. an accession number, a clone or contig), and a
positional range relative to the landmark.  If no range is specified,
then the entire extent of the landmark is used to generate the
segment.

You may also provide the ID of a "reference" sequence, which will set
the coordinate system and orientation used for all features contained
within the segment.  The reference sequence can be changed later.  If
no reference sequence is provided, then the coordinate system is based
on the landmark.

Arguments:

 -seq          ID of the landmark sequence.

 -class        Database object class for the landmark sequence.
               "Sequence" assumed if not specified.  This is
               irrelevant for databases which do not recognize
               object classes.

 -start        Start of the segment relative to landmark.  Positions
               follow standard 1-based sequence rules.  If not specified,
               defaults to the beginning of the landmark.

 -stop         Stop of the segment relative to the landmark.  If not specified,
               defaults to the end of the landmark.

 -offset       For those who prefer 0-based indexing, the offset specifies the
               position of the new segment relative to the start of the landmark.

 -length       For those who prefer 0-based indexing, the length specifies the
               length of the new segment.

 -refseq       Specifies the ID of the reference landmark used to establish the
               coordinate system for the newly-created segment.

 -refclass     Specifies the class of the reference landmark, for those databases
               that distinguish different object classes.  Defaults to "Sequence".

 -name,-sequence,-sourceseq   Aliases for -seq.

 -begin,-end   Aliases for -start and -stop

 -off,-len     Aliases for -offset and -length

 -seqclass     Alias for -class

Here's an example to explain how this works:

  my $db = Bio::DB::GFF->new(-dsn => 'dbi:mysql:human',-adaptor=>'dbi:mysql');

If successful, $db will now hold the database accessor object.  We now
try to fetch the fragment of sequence whose ID is A0000182 and class
is "Accession."

  my $segment = $db->segment(-name=>'A0000182',-class=>'Accession');

If successful, $segment now holds the entire segment corresponding to
this accession number.  By default, the sequence is used as its own
reference sequence, so its first base will be 1 and its last base will
be the length of the accession.

Assuming that this sequence belongs to a longer stretch of DNA, say a
contig, we can fetch this information like so:

  my $sourceseq = $segment->sourceseq;

and find the start and stop on the source like this:

  my $start = $segment->abs_start;
  my $stop = $segment->abs_stop;

If we had another segment, say $s2, which is on the same contiguous
piece of DNA, we can pass that to the refseq() method in order to
establish it as the coordinat reference point:

  $segment->refseq($s2);

Now calling start() will return the start of the segment relative to
the beginning of $s2, accounting for differences in strandedness:

  my $rel_start = $segment->start;

=cut

sub segment {
  my $self = shift;
  # (see Ace::Sequence::DBI::Segment for all the arguments)
  return $_[0] =~ /^-/ ? Bio::DB::GFF::RelSegment->new(-factory => $self,@_)
                       : Bio::DB::GFF::RelSegment->new($self,@_);
}

=head2 types

 Title   : types
 Usage   : $db->types(@args)
 Function: return list of feature types in range or database
 Returns : a list of Bio::DB::GFF::Typename objects
 Args    : see below
 Status  : public

This routine returns a list of feature types known to the database.
The list can be database-wide or restricted to a region.  It is also
possible to find out how many times each feature occurs.

For range queries, it is usually more convenient to create a
Bio::DB::GFF::Segment object, and then invoke it's types() method.

Arguments are as follows:

  -ref        ID of reference sequence
  -class      class of reference sequence
  -start      start of segment
  -stop       stop of segment
  -enumerate  if true, count the features

The returned value will be a list of Bio::DB::GFF::Typename objects,
which if evaluated in a string context will return the feature type in 
"method:source" format.  This object class also has method() and
source() methods for retrieving the like-named fields.

If -enumerate is true, then the function returns a hash (not a hash
reference) in which the keys are type names in "method:source" format
and the values are the number of times each feature appears in the
database or segment.

The argument -end is a synonum for -stop, and -count is a synonym for
-enumerate.

=cut

sub types {
  my $self = shift;
  my ($refseq,$start,$stop,$enumerate) = rearrange ([
						     [qw(REF REFSEQ)],
						     qw(START),
						     [qw(STOP END)],
						     [qw(ENUMERATE COUNT)],
						     ],@_);
  $self->get_types($refseq,$start,$stop,$enumerate);
}

=head2 dna

 Title   : dna
 Usage   : $db->dna($id,$class,$start,$stop)
 Function: return the raw DNA string for a segment
 Returns : a raw DNA string
 Args    : id of the sequence, its class, start and stop positions
 Status  : public

This method is invoked by Bio::DB::GFF::Segment to fetch the raw DNA
sequence.

NOTE: you will probably prefer to create a Segment and then invoke its
dna() method.

=cut

# call to return the DNA string for the indicated region
# real work is done by get_dna()
sub dna {
  my $self = shift;
  my ($id,$class,$start,$stop) = rearrange([
					    [qw(NAME ID REF REFSEQ)],
					    'CLASS',
					    qw(START),
					    [qw(STOP END)],
					   ],@_);
  return unless defined $start && defined $stop;
  $self->get_dna($id,$class,$start,$stop);
}

=head2 overlapping_features

 Title   : overlapping_features
 Usage   : $db->overlapping_features(@args)
 Function: get features that overlap the indicated range
 Returns : a list of Bio::DB::GFF::Feature objects
 Args    : see below
 Status  : public

This method is invoked by Bio::DB::GFF::Segment->features() to find
the list of features that overlap a given range.  It is generally
preferable to create the Segment first, and then fetch the features.

This method takes set of named arguments:

  -refseq    ID of the reference sequence
  -class     Class of the reference sequence
  -start     Start of the desired range in refseq coordinates
  -stop      Stop of the desired range in refseq coordinates
  -types     List of feature types to return.  Argument is an array
	     reference containing strings of the format "method:source"
  -parent    A parent Bio::DB::GFF::Segment object, used to create
	     relative coordinates in the generated features.
  -merge     Whether to apply aggregators to the generated features.
  -iterator  Whether to return an iterator across the features.

If -iterator is true, then the method returns a single scalar value
consisting of a Bio::SeqIO object.  You can call next_seq() repeatedly
on this object to fetch each of the features in turn.  If iterator is
false or absent, then all the features are returned as a list.

Currently aggregation is disabled when iterating over a series of
features.

Types are indicated using the nomenclature "method:source".  Either of
these fields can be omitted, in which case a wildcard is used for the
missing field.  Type names without the colon (e.g. "exon") are
interpreted as the method name and a source wild card.  Regular
expressions are allowed in either field, as in: "similarity:BLAST.*".

=cut

# call to return the features that overlap the named region
# real work is done by get_features
sub overlapping_features {
  my $self = shift;
  my ($refseq,$class,$start,$stop,$types,$parent,$automerge,$iterator) =
    rearrange([
	       [qw(REF REFSEQ)],
	       qw(CLASS),
	       qw(START),
	       [qw(STOP END)],
	       [qw(TYPE TYPES)],
	       qw(PARENT),
	       [qw(MERGE AUTOMERGE)],
	       'ITERATOR'
	      ],@_);

  # return unless defined $start && defined $stop;
  $automerge = 1 unless defined $automerge;
  $self->_features(0,$refseq,$class,$start,$stop,$types,$parent,$automerge,$iterator);
}


=head2 contained_features

 Title   : contained_features
 Usage   : $db->contained_features(@args)
 Function: get features that are contained within the indicated range
 Returns : a list of Bio::DB::GFF::Feature objects
 Args    : see overlapping_features()
 Status  : public

This call is similar to overlapping_features(), except that it only
retrieves features whose end points are completely contained within
the specified range.

Generally you will want to fetch a Bio::DB::GFF::Segment object and
call its contained_features() method rather than call this directly.

=cut

# The same, except that it only returns features that are completely contained within the
# range (much faster usually)
sub contained_features {
  my $self = shift;
  my ($refseq,$class,$start,$stop,$types,$parent,$automerge,$iterator) = 
    rearrange([
	       [qw(REF REFSEQ)],
	       qw(CLASS),
	       qw(START),
	       [qw(STOP END)],
	       [qw(TYPE TYPES)],
	       qw(PARENT),
	       [qw(MERGE AUTOMERGE)],
	       'ITERATOR'
	      ],@_);

  # return unless defined $start && defined $stop;
  $automerge = 1 unless defined $automerge;
  $self->_features(1,$refseq,$class,$start,$stop,$types,$parent,$automerge,$iterator);
}

=head2 features

 Title   : features
 Usage   : $db->features(@args)
 Function: get all features, possibly filtered by type
 Returns : a list of Bio::DB::GFF::Feature objects
 Args    : see below
 Status  : public

This routine will retrieve features in the database regardless of
position.  It can be used to return all features, or a subset based on
their method and source.

Arguments are as follows:

  -types     List of feature types to return.  Argument is an array
	     reference containing strings of the format "method:source"
  -merge     Whether to apply aggregators to the generated features.
  -iterator  Whether to return an iterator across the features.

If -iterator is true, then the method returns a single scalar value
consisting of a Bio::SeqIO object.  You can call next_seq() repeatedly
on this object to fetch each of the features in turn.  If iterator is
false or absent, then all the features are returned as a list.

Currently aggregation is disabled when iterating over a series of
features.

Types are indicated using the nomenclature "method:source".  Either of
these fields can be omitted, in which case a wildcard is used for the
missing field.  Type names without the colon (e.g. "exon") are
interpreted as the method name and a source wild card.  Regular
expressions are allowed in either field, as in: "similarity:BLAST.*".

=cut

# The same, except that it fetches all features of a particular type regardless
# of position
sub features {
  my $self = shift;
  my ($types,$automerge,$iterator);
  if ($_[0] =~ /^-/) {
    ($types,$automerge,$iterator) = rearrange([
					       [qw(TYPE TYPES)],
					       [qw(MERGE AUTOMERGE)],
					       'ITERATOR'
					      ],@_);
  } else {
    $types = \@_;
  }

  $automerge = 1 unless defined $automerge;
  $self->_features(1,undef,undef,undef,undef,$types,undef,$automerge,$iterator);
}

=head2 add_aggregator

 Title   : add_aggregator
 Usage   : $db->add_aggregator($aggregator)
 Function: add an aggregator to the list
 Returns : nothing
 Args    : an aggregator
 Status  : public

This method will append an aggregator to the end of the list of
registered aggregators.

=cut

sub add_aggregator {
  my $self       = shift;
  my $aggregator = shift;
  my $list = $self->{aggregators} ||= [];
  if (ref $aggregator) { # an object
    push @$list,$a;
  } else {
    my $class = "Bio::DB::GFF::Aggregator::\L${aggregator}\E";
    eval "require $class";
    $self->throw("Unable to load $aggregator aggregator: $@") if $@;
    push @$list,$class->new();
  }
}

=head2 aggregators

 Title   : aggregators
 Usage   : $db->aggregators;
 Function: retrieve list of aggregators
 Returns : list of aggregators
 Args    : none
 Status  : public

This method will return a list of aggregators currently assigned to
the object.

=cut

sub aggregators {
  my $self = shift;
  return unless $self->{aggregators};
  return @{$self->{aggregators}};
}

=head2 abscoords

 Title   : abscoords
 Usage   : $db->abscoords($name,$class)
 Function: finds position of a landmark in reference coordinates
 Returns : ($ref,$class,$start,$stop,$strand)
 Args    : name and class of landmark
 Status  : public

This method is called by Bio::DB::GFF::RelSegment to obtain the
absolute coordinates of a sequence landmark.  The arguments are the
name and class of the landmark.  If successful, abscoords() returns
the ID of the reference sequence, its class, its start and stop
positions, and the orientation of the reference sequence's coordinate
system ("+" for forward strand, "-" for reverse strand).

=cut

# given a sequence class and name, return its coordinates in format (reference,start,stop,strand)
sub abscoords {
  my $self = shift;
  my ($name,$class) = @_;
  $class ||= 'Sequence';
  $self->get_abscoords($name,$class);
}

=head1 Protected API

The following methods are not intended for public consumption, but are
intended to be overridden/implemented by adaptors.

=head2 default_aggregators

 Title   : default_aggregators
 Usage   : $db->default_aggregators;
 Function: retrieve list of aggregators
 Returns : array reference containing list of aggregator names
 Args    : none
 Status  : protected

This method (which is intended to be overridden by adaptors) returns a
list of standard aggregators to be applied when no aggregators are
specified in the constructor.

=cut

sub default_aggregators {
  my $self = shift;
  return ['transcript','clone','alignment'];
}

=head2 load_gff

 Title   : load_gff
 Usage   : $db->load_gff
 Function: load a GFF input stream
 Returns : number of features loaded
 Args    : none
 Status  : protected

This method is called to load a GFF data stream.  The method will read
GFF features from <> and load them into the database.  On exit the
method must return the number of features loaded.

Note that the method is responsible for parsing the GFF lines.  This
is to allow for differences in the interpretation of the "group"
field, which are legion.

=cut

# load from <>
sub load_gff {
  my $self = shift;
  $self->setup_load();

  while (<>) {
    my ($ref,$source,$method,$start,$stop,$score,$strand,$phase,$group) = split "\t";
    next if /^\#/;

    # handle group parsing
    $group =~ s/(\"[^\"]*);([^\"]*\")/$1$;$2/g;  # protect embedded semicolons in the group
    my @groups = split(/\s*;\s*/,$group);
    foreach (@groups) { s/$;/;/g }

    my ($gclass,$gname,$tstart,$tstop,$notes) = $self->_split_group(@groups);

    # call subclass to do the dirty work
    $self->load_gff_line($ref,$source,$method,$start,$stop,
			 $score,$strand,$phase,
			 $gclass,$gname,$tstart,$tstop,$notes);
  }

  $self->finish_load();
}

=head2 setup_load

 Title   : setup_load
 Usage   : $db->setup_load
 Function: called before load_gff_line()
 Returns : void
 Args    : none
 Status  : protected

This method gives subclasses a chance to do any schema-specific
initialization prior to loading a set of GFF records.

=cut

sub setup_load {
  shift->throw("setup_load(): must be implemented by an adaptor");
}

=head2 finish_load

 Title   : finish_load
 Usage   : $db->finish_load
 Function: called after load_gff_line()
 Returns : number of records loaded
 Args    : none
 Status  : protected

This method gives subclasses a chance to do any schema-specific
cleanup after loading a set of GFF records.

=cut

sub finish_load {
  shift->throw("finish_load(): must be implemented by an adaptor");
}

=head2 load_gff_line

 Title   : load_gff_line
 Usage   : $db->load_gff_line(@args)
 Function: called to load one parsed line of GFF
 Returns : true if successfully inserted
 Args    : see below
 Status  : protected

This method is called once per line of the GFF and passed a series of
parsed data items.  The items are:

 $ref          reference sequence
 $source       annotation source
 $method       annotation method
 $start        annotation start
 $stop         annotation stop
 $score        annotation score (may be undef)
 $strand       annotation strand (may be undef)
 $phase        annotation phase (may be undef)
 $group_class  class of annotation's group (may be undef)
 $group_name   ID of annotation's group (may be undef)
 $target_start start of target of a similarity hit
 $target_stop  stop of target of a similarity hit
 $notes        array reference of text items to be attached

=cut

sub load_gff_line {
  shift->throw("load_gff_line(): must be implemented by an adaptor");
}


=head2 do_initialize

 Title   : do_initialize
 Usage   : $db->do_initialize([$erase])
 Function: initialize and possibly erase database
 Returns : true if successful
 Args    : optional erase flag
 Status  : protected

This method implements the initialize() method described above, and
takes the same arguments.

=cut

sub do_initialize {
    shift->throw('do_initialize(): must be implemented by an adaptor');
}

=head2 get_dna

 Title   : get_dna
 Usage   : $db->get_dna($id,$class,$start,$stop)
 Function: get DNA for indicated segment
 Returns : the dna string
 Args    : sequence ID, start, stop and class
 Status  : protected

If start > stop and the sequence is nucleotide, then this method
should return the reverse complement.  The sequence class may be
ignored by those databases that do not recognize different object
types.

=cut

sub get_dna {
  my $self = shift;
  my ($id,$class,$start,$stop) = @_;
  $self->throw("get_dna() must be implemented by an adaptor");
}

=head2 get_features

 Title   : get_features
 Usage   : $db->get_features($isrange,$refseq,$class,$start,$stop,$types,$callback)
 Function: get list of features for a region
 Returns : count of number of features retrieved
 Args    : see below
 Status  : protected

Arguments are as follows:

   $isrange   Flag indicating that a range query is desired, in which 
              case only features that are completely contained within
              start->stop (inclusive) are retrieved.  Otherwise, an
              overlap retrieval is performed to find those 

   $refseq    ID of the landmark that establishes the absolute 
              coordinate system.

   $class     Class of this landmark.  Can be ignored by implementations
              that don't recognize such distinctions.

   $start,$stop  Start and stop of the range, inclusive.

   $types     Array reference containing the list of annotation types
              to fetch from the database.  Each annotation type is an
              array reference consisting of [source,method].

   $callback  A code reference.  As the passed features are retrieved
              they are passed to this callback routine for processing.

This routine is responsible for getting arrays of GFF data out of the
database and passing them to the callback subroutine.  The callback
does the work of constructing a Bio::DB::GFF::Feature object out of
that data.  The callback expects a list of 11 fields:

  $srcseq      source sequence
  $start       feature start
  $stop        feature stop
  $source      feature source
  $method      feature method
  $score       feature score
  $strand      feature strand
  $phase       feature phase
  $groupclass  group class (may be undef)
  $groupname   group ID (may be undef)
  $tstart      target start for similarity hits (may be undef)
  $tstop       target stop for similarity hits (may be undef)

These fields are in the same order as the raw GFF file, with the
exception that the group column has been parsed into class and name
fields.

=cut

sub get_features{
  my $self = shift;
  my ($isrange,$srcseq,$class,$start,$stop,$types,$callback) = @_;
  $self->throw("get_features() must be implemented by an adaptor");
}


=head2 get_abscoords

 Title   : get_abscoords
 Usage   : $db->get_abscoords($name,$class)
 Function: get the absolute coordinates of sequence with name & class
 Returns : ($absref,$absstart,$absstop,$absstrand)
 Args    : name and class of the landmark
 Status  : protected

Given the name and class of a genomic landmark, this function returns
a four-element array consisting of:

  $absref      the ID of the reference sequence that contains this landmark
  $absstart    the position at which the landmark starts
  $absstop     the position at which the landmark stops
  $absstrand   the strand of the landmark, relative to the reference sequence

=cut

sub get_abscoords {
  my $self = shift;
  my ($name,$class) = @_;
  $self->throw("get_abscoords() must be implemented by an adaptor");
}

=head2 get_types

 Title   : get_types
 Usage   : $db->get_types($absref,$start,$stop,$count)
 Function: get list of all feature types on the indicated segment
 Returns : list or hash of Bio::DB::GFF::Typename objects
 Args    : see below
 Status  : protected

Arguments are:

  $absref      the ID of the reference sequence
  $class       the class of the reference sequence
  $start       the position to start counting
  $stop        the position to end counting
  $count       a boolean indicating whether to count the number
	       of occurrences of each feature type

If $count is true, then a hash is returned.  The keys of the hash are
feature type names in the format "method:source" and the values are
the number of times a feature of this type overlaps the indicated
segment.  Otherwise, the call returns a set of Bio::DB::GFF::Typename
objects.  If $start or $stop are undef, then all features on the
indicated segment are enumerated.  If $absref is undef, then the call
returns all feature types in the database.

=cut

sub get_types {
  my $self = shift;
  my ($refseq,$class,$start,$stop,$count) = @_;
  $self->throw("get_types() must be implemented by an adaptor");
}


=head2 make_feature

 Title   : make_feature
 Usage   : $db->make_feature(@args)
 Function: Create a Bio::DB::GFF::Feature object from string data
 Returns : a Bio::DB::GFF::Feature object
 Args    : see below
 Status  : internal

 This takes 14 arguments (really!):

  $parent                A Bio::DB::GFF::RelSegment object
  $group_hash            A hashref containing unique list of GFF groups
  $absref                The reference sequence for this feature
  $start                 Start of feature
  $stop                  Stop of feature
  $source                Feature source field
  $method                Feature method field
  $score                 Feature score field
  $strand                Feature strand
  $phase                 Feature phase
  $group_class           Class of feature group
  $group_name            Name of feature group         
  $tstart                For homologies, start of hit on target
  $tstop                 Stop of hit on target

The $parent argument, if present, is used to establish relative
coordinates in the resulting Bio::DB::Feature object.  This allows one
feature to generate a list of other features that are relative to its
coordinate system (for example, finding the coordinates of the second
exon relative to the coordinates of the first).

The $group_hash allows the group_class/group_name strings to be turned
into rich database objects via the make_obect() method (see above).
Because these objects may be expensive to create, $group_hash is used
to uniquefy them.  The index of this hash is the composite key
{$group_class,$group_name,$tstart,$tstop}.  Values are whatever object
is returned by the make_object() method.

The remainder of the fields are taken from the GFF line, with the
exception that "Target" features, which contain information about the
target of a homology search, are parsed into their components.

=cut

# This call is responsible for turning a line of GFF into a
# feature object.
# The $parent argument is a Bio::DB::GFF::Segment object and is used
# to establish the coordinate system for the new feature.
# The $group_hash argument is an hash ref that holds previously-
# generated group objects.
# Other arguments are taken right out of the GFF table.
sub make_feature {
  my $self = shift;
  my ($parent,$group_hash,
      $srcseq,$start,$stop,
      $source,$method,
      $score,$strand,$phase,
      $group_class,$group_name,
      $tstart,$tstop,$db_id) = @_;

  my $group;  # undefined
  if (defined $group_class && defined $group_name) {
    $tstart ||= '';
    $tstop  ||= '';
    $group = $group_hash->{$group_class,$group_name,$tstart,$tstop} 
      ||= $self->make_object($group_class,$group_name,$tstart,$tstop);
  }

  if (ref $parent) { # note that the src sequence is ignored
    return Bio::DB::GFF::Feature->new_from_parent($parent,$start,$stop,
						  $method,$source,
						  $score,$strand,$phase,
						  $group,$db_id);
  } else {
    return Bio::DB::GFF::Feature->new($self,$srcseq,
				      $start,$stop,
				      $method,$source,
				      $score,$strand,$phase,
				      $group,$db_id);
  }
}

=head2 parse_types

 Title   : parse_types
 Usage   : $db->parse_types(@args)
 Function: parses list of types
 Returns : an array ref containing ['method','source'] pairs
 Args    : a list of types in 'method:source' form
 Status  : internal

This method takes an array of type names in the format "method:source"
and returns an array reference of ['method','source'] pairs.  It will
also accept a single argument consisting of an array reference with
the list of type names.

=cut

# turn feature types in the format "method:source" into a list of [method,source] refs
sub parse_types {
  my $self  = shift;
  return [] if !@_ or !defined($_[0]);

  my @types = ref($_[0]) ? @{$_[0]} : @_;
  my @type_list = map { [split(':',$_,2)] } @types;
  return \@type_list;
}

=head2 make_match_sub

 Title   : make_match_sub
 Usage   : $db->make_match_sub($types)
 Function: creates a subroutine used for filtering features
 Returns : a code reference
 Args    : a list of parsed type names
 Status  : protected

This method is used internally to generate a code subroutine that will
accept or reject a feature based on its method and source.  It takes
an array of parsed type names in the format returned by parse_types(),
and generates an anonymous subroutine.  The subroutine takes a single
Bio::DB::GFF::Feature object and returns true if the feature matches
one of the desired feature types, and false otherwise.

=cut

# a subroutine that matches features indicated by list of types
sub make_match_sub {
  my $self = shift;
  my $types = shift;

  return sub { 1 } unless ref $types && @$types;

  my @expr;
  for my $type (@$types) {
    my ($method,$source) = @$type;
    $method ||= '.*';
    $source ||= '.*';
    push @expr,"$method:$source";
  }
  my $expr = join '|',@expr;
  return $self->{match_subs}{$expr} if $self->{match_subs}{$expr};

  my $sub =<<END;
sub {
  my \$feature = shift;
  return \$feature->type =~ /^$expr\$/i;
}
END
  my $compiled_sub = eval $sub;
  $self->throw($@) if $@;
  return $self->{match_subs}{$expr} = $compiled_sub;
}

=head2 make_object

 Title   : make_object
 Usage   : $db->make_object($name,$class,$start,$stop)
 Function: creates a feature object
 Returns : a feature object
 Args    : see below
 Status  : protected

This method is called to make an object from the GFF "group" field.
By default, all Target groups are turned into Bio::DB::GFF::Homol
objects, and everything else becomes a Bio::DB::GFF::Featname.
However, adaptors are free to override this method to generate more
interesting objects, such as true BioPerl objects, or Acedb objects.

Arguments are:

  $name      database ID for object
  $class     class of object
  $start     for similarities, start of match inside object
  $stop      for similarities, stop of match inside object

=cut

# abstract call to turn a feature into an object, given its class and name
sub make_object {
  my $self = shift;
  my ($name,$class,$start,$stop) = @_;
  return Bio::DB::GFF::Homol->new($self,$name,$class,$start,$stop) if defined $start and length $start;
  return Bio::DB::GFF::Featname->new($class,$name);
}

=head1 Internal Methods

The following methods are internal to Bio::DB::GFF and are not
guaranteed to remain the same.

=head2 _features

 Title   : _features
 Usage   : $db->_features(@args)
 Function: internal method
 Returns : a list of Bio::DB::GFF::Feature objects
 Args    : see below
 Status  : internal

This is an internal method that is called by overlapping_features(),
contained_features() and features() to do the actual work.  It takes
nine positional arguments:

  $range_query	 if true, this is a request for contained features, 
		 otherwise for overlapping features.
  $refseq        reference sequence ID
  $class	 reference sequence class
  $start	 start of range
  $stop		 stop of range
  $types	 list of types
  $parent	 parent sequence, for relative coordinates
  $automerge	 if true, invoke aggregators to merge features
  $iterator	 if true, return an iterator

=cut

sub _features {
  my $self = shift;
  my ($range_query,$refseq,$class,$start,$stop,$types,$parent,$automerge,$iterator) = @_;

  ($start,$stop) = ($stop,$start) if defined($start) && $start > $stop;

  $types = $self->parse_types($types);  # parse out list of types
  my $aggregated_types = $types;         # keep a copy

  my %groups;  # cache groups so that we don't create them unecessarily

  if ($iterator) {
    my $callback = sub { $self->make_feature($parent,\%groups,@_) };
    return $self->get_features_iterator($range_query,$refseq,$class,
					$start,$stop,$aggregated_types,$callback) ;
  }

  # allow the aggregators to operate on the original
  if ($automerge) {
    for my $a ($self->aggregators) {
      $aggregated_types = $a->disaggregate($aggregated_types,$self);
    }
  }

  my $features = [];

  my $callback = sub { push @$features,$self->make_feature($parent,\%groups,@_) };
  $self->get_features($range_query,$refseq,$class,
		      $start,$stop,$aggregated_types,$callback) ;

  if ($automerge) {
    warn "aggregating...\n" if $self->debug;
    my @aggregated;
    foreach my $a (reverse $self->aggregators) {  # last aggregator gets first shot
      $features = $a->aggregate($features,$self);
    }
  }

  warn "filtering...\n" if $self->debug;

  # remove anything from the features list that was not specifically requested.
  my $match = $self->make_match_sub($types);
  return grep { $match->($_) } @$features;
}

=head2 _split_group

 Title   : _split_group
 Usage   : $db->_split_group(@groups)
 Function: parse GFF group field
 Returns : ($gclass,$gname,$tstart,$tstop,$notes)
 Args    : a list of group fields from a GFF line
 Status  : internal

This is an internal method that is called by load_gff_line to parse
out the contents of one or more group fields.  It returns the class of
the group, its name, the start and stop of the target, if any, and an
array reference containing any notes that were stuck into the group
field.

=cut

sub _split_group {
  my $self = shift;
  my @groups = @_;

  my ($gclass,$gname,$tstart,$tstop,@notes);

  for (@groups) {

    my ($tag,$value) = /(\S+)(?:\s+\"([^\"]+)\")?/;
    $value =~ s/\\t/\t/g;
    $value =~ s/\\r/\r/g;

    # if the tag is "Note", then we add this to the
    # notes array
   if ($tag eq 'Note') {  # just a note, not a group!
     push @notes,$value;
   }

    # if the tag eq 'Target' then the class name is embedded in the ID
    # (the GFF format is obviously screwed up here)
    elsif ($tag eq 'Target' && /\"([^:\"]+):([^\"]+)\"/) {
      ($gclass,$gname) = ($1,$2);
      ($tstart,$tstop) = /(\d+) (\d+)/;
    }

    elsif (!$value) {
      push @notes,$tag;  # e.g. "Confirmed_by_EST"
    }

    # otherwise, the tag and value correspond to the
    # group class and name
    else {
      ($gclass,$gname) = ($tag,$value);
    }
  }

  return ($gclass,$gname,$tstart,$tstop,\@notes);
}

1;

__END__

=head1 BUGS

Not really Bio::SeqFeatureI compliant yet.

Schemas need some work.

=head1 SEE ALSO

L<bioperl>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

