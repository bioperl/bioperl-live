=head1 NAME

Bio::DB::GFF -- Storage and retrieval of sequence annotation data

=head1 SYNOPSIS

  use Bio::DB::GFF;

  # Open the sequence database
  my $db      = Bio::DB::GFF->new( -adaptor => 'dbi:mysql',
                                   -dsn     => 'dbi:mysql:elegans42');

  # fetch a 1 megabase segment of sequence starting at landmark "ZK909"
  my $segment = $seqfactory->segment('ZK909', 1 => 1000000);

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
database.  It supports multiple types of database (ACeDB, relational),
and multiple schemas through a system of adaptors.  It allows
annotations to be aggregated together into hierarchical structures
(genes, transcripts, exons, introns) using a system of aggregators.

The following types of operation are supported by this module:

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

-mostly to come-

The sequences used to establish the coordinate system for annotations
can correspond to sequenced clones, clone fragments, contigs or
super-contigs.  Thus, this module can be used throughout the lifecycle
of a sequencing project.

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

=cut

sub segment {
  my $self = shift;
  # (see Ace::Sequence::DBI::Segment for all the arguments)
  return $_[0] =~ /^-/ ? Bio::DB::GFF::RelSegment->new(-factory => $self,@_)
                       : Bio::DB::GFF::RelSegment->new($self,@_);
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
  shift->throw("load_gff(): must be implemented by an adaptor");
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
exception that some parsing has been performed on the group field

=cut

sub get_features{
  my $self = shift;
  my ($isrange,$srcseq,$class,$start,$stop,$types,$callback) = @_;
  $self->throw("get_features() must be implemented by an adaptor");
}


sub get_abscoords {
  my $self = shift;
  my ($name,$class) = @_;
  $self->throw("get_abscoords() must be implemented by an adaptor");
}

sub get_types {
  my $self = shift;
  my ($refseq,$start,$stop,$count) = @_;
  $self->throw("get_types() must be implemented by an adaptor");
}


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
      $tstart,$tstop) = @_;

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
						  $group);
  } else {
    return Bio::DB::GFF::Feature->new($self,$srcseq,
				      $start,$stop,
				      $method,$source,
				      $score,$strand,$phase,
				      $group);
  }
}

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

sub _features {
  my $self = shift;
  my ($range_query,$refseq,$class,$start,$stop,$types,$parent,$automerge,$iterator) = @_;

  ($start,$stop) = ($stop,$start) if $start > $stop;

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

# turn feature types in the format "method:source" into a list of [method,source] refs
sub parse_types {
  my $self  = shift;
  return [] if !@_ or !defined($_[0]);

  my @types = ref($_[0]) ? @{$_[0]} : @_;
  my @type_list = map { [split(':',$_,2)] } @types;
  return \@type_list;
}

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

# abstract call to turn a feature into an object, given its class and name
sub make_object {
  my $self = shift;
  my ($class,$name,$start,$stop) = @_;
  return Bio::DB::GFF::Homol->new($self,$name,$class,$start,$stop) if defined $start and length $start;
  return Bio::DB::GFF::Featname->new($class,$name);
}

# given a sequence class and name, return its coordinates in format (reference,start,stop,strand)
sub abscoords {
  my $self = shift;
  my ($name,$class) = @_;
  $class ||= 'Sequence';
  $self->get_abscoords($name,$class);
}

1;

=head1 BUGS

Not completely Bio::SeqFeatureI compliant yet.

Schemas need some work.

=head1 SEE ALSO

L<bioperl>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.  

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

