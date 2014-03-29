
=head1 NAME

Bio::DB::GFF -- Storage and retrieval of sequence annotation data

=head1 SYNOPSIS

  use Bio::DB::GFF;

  # Open the sequence database
  my $db      = Bio::DB::GFF->new( -adaptor => 'dbi::mysqlopt',
                                   -dsn     => 'dbi:mysql:elegans');

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
  my $upstream = $exons[0]->subseq(-1000,0);

  # get its DNA
  my $dna = $upstream->seq;

  # and get all curated polymorphisms inside it
  @polymorphisms = $upstream->contained_features('polymorphism:curated');

  # get all feature types in the database
  my @types = $db->types;

  # count all feature types in the segment
  my %type_counts = $segment->types(-enumerate=>1);

  # get an iterator on all curated features of type 'exon' or 'intron'
  my $iterator = $db->get_seq_stream(-type     => ['exon:curated','intron:curated']);

  while (my $s = $iterator->next_seq) {
      print $s,"\n";
  }

  # find all transcripts annotated as having function 'kinase'
  my $iterator = $db->get_seq_stream(-type=>'transcript',
			             -attributes=>{Function=>'kinase'});
  while (my $s = $iterator->next_seq) {
      print $s,"\n";
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
  - using a sequence segment to create new segments based on relative 
    offsets

The data model used by Bio::DB::GFF is compatible with the GFF flat
file format (L<http://www.sequenceontology.org/gff3.shtml>). The module 
can load a set of GFF files into the database, and serves objects that
have methods corresponding to GFF fields.

The objects returned by Bio::DB::GFF are compatible with the
SeqFeatureI interface, allowing their use by the Bio::Graphics and
Bio::DAS modules.

=head2 Auxiliary Scripts

The bioperl distribution includes several scripts that make it easier
to work with Bio::DB::GFF databases.  They are located in the scripts
directory under a subdirectory named Bio::DB::GFF:

=over 4

=item *

bp_load_gff.pl

This script will load a Bio::DB::GFF database from a flat GFF file of
sequence annotations.  Only the relational database version of
Bio::DB::GFF is supported.  It can be used to create the database from
scratch, as well as to incrementally load new data.

This script takes a --fasta argument to load raw DNA into the database
as well.  However, GFF databases do not require access to the raw DNA
for most of their functionality.

load_gff.pl also has a --upgrade option, which will perform a
non-destructive upgrade of older schemas to newer ones.

=item *

bp_bulk_load_gff.pl

This script will populate a Bio::DB::GFF database from a flat GFF file
of sequence annotations.  Only the MySQL database version of
Bio::DB::GFF is supported.  It uses the "LOAD DATA INFILE" query in
order to accelerate loading considerably; however, it can only be used
for the initial load, and not for updates.

This script takes a --fasta argument to load raw DNA into the database
as well.  However, GFF databases do not require access to the raw DNA
for most of their functionality.

=item *

bp_fast_load_gff.pl

This script is as fast as bp_bulk_load_gff.pl but uses Unix pipe
tricks to allow for incremental updates.  It only supports the MySQL
database version of Bio::DB::GFF and is guaranteed not to work on
non-Unix platforms.

Arguments are the same as bp_load_gff.pl

=item *

gadfly_to_gff.pl

This script will convert the GFF-like format used by the Berkeley
Drosophila Sequencing project into a format suitable for use with this
module.

=item *

sgd_to_gff.pl

This script will convert the tab-delimited feature files used by the
Saccharomyces Genome Database into a format suitable for use with this
module.

=back

=head2 GFF Fundamentals

The GFF format is a flat tab-delimited file, each line of which
corresponds to an annotation, or feature.  Each line has nine columns
and looks like this:

 Chr1  curated  CDS 365647  365963  .  +  1  Transcript "R119.7"

The 9 columns are as follows:

=over 4

=item 1.

reference sequence

This is the ID of the sequence that is used to establish the
coordinate system of the annotation.  In the example above, the
reference sequence is "Chr1".

=item 2.

source

The source of the annotation.  This field describes how the annotation
was derived.  In the example above, the source is "curated" to
indicate that the feature is the result of human curation.  The names
and versions of software programs are often used for the source field,
as in "tRNAScan-SE/1.2".

=item 3.

method

The annotation method.  This field describes the type of the
annotation, such as "CDS".  Together the method and source describe
the annotation type.

=item 4.

start position

The start of the annotation relative to the reference sequence. 

=item 5.

stop position

The stop of the annotation relative to the reference sequence.  Start
is always less than or equal to stop.

=item 6.

score

For annotations that are associated with a numeric score (for example,
a sequence similarity), this field describes the score.  The score
units are completely unspecified, but for sequence similarities, it is
typically percent identity.  Annotations that don't have a score can
use "."

=item 7.

strand

For those annotations which are strand-specific, this field is the
strand on which the annotation resides.  It is "+" for the forward
strand, "-" for the reverse strand, or "." for annotations that are
not stranded.

=item 8.

phase

For annotations that are linked to proteins, this field describes the
phase of the annotation on the codons.  It is a number from 0 to 2, or
"." for features that have no phase.

=item 9.

group

GFF provides a simple way of generating annotation hierarchies ("is
composed of" relationships) by providing a group field.  The group
field contains the class and ID of an annotation which is the logical
parent of the current one.  In the example given above, the group is
the Transcript named "R119.7".

The group field is also used to store information about the target of
sequence similarity hits, and miscellaneous notes.  See the next
section for a description of how to describe similarity targets.

The format of the group fields is "Class ID" with a single space (not
a tab) separating the class from the ID. It is VERY IMPORTANT to
follow this format, or grouping will not work properly.

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

Double-quotes are sometimes used in GFF files around components of the
group field.  Strictly, this is only necessary if the group name or
class contains whitespace.

=head2 Making GFF files work with this module

Some annotations do not need to be individually named.  For example,
it is probably not useful to assign a unique name to each ALU repeat
in a vertebrate genome.  Others, such as predicted genes, correspond
to named biological objects; you probably want to be able to fetch the
positions of these objects by referring to them by name.

To accommodate named annotations, the GFF format places the object
class and name in the group field.  The name identifies the object,
and the class prevents similarly-named objects, for example clones and
sequences, from collding.

A named object is shown in the following excerpt from a GFF file:

 Chr1  curated transcript  939627 942410 . +  . Transcript Y95B8A.2

This object is a predicted transcript named Y95BA.2.  In this case,
the group field is used to identify the class and name of the object,
even though no other annotation belongs to that group.

It now becomes possible to retrieve the region of the genome covered
by transcript Y95B8A.2 using the segment() method:

  $segment = $db->segment(-class=>'Transcript',-name=>'Y95B8A.2');

It is not necessary for the annotation's method to correspond to the
object class, although this is commonly the case.

As explained above, each annotation in a GFF file refers to a
reference sequence.  It is important that each reference sequence also
be identified by a line in the GFF file.  This allows the Bio::DB::GFF
module to determine the length and class of the reference sequence,
and makes it possible to do relative arithmetic.

For example, if "Chr1" is used as a reference sequence, then it should
have an entry in the GFF file similar to this one:

 Chr1 assembly chromosome 1 14972282 . + . Sequence Chr1

This indicates that the reference sequence named "Chr1" has length
14972282 bp, method "chromosome" and source "assembly".  In addition,
as indicated by the group field, Chr1 has class "Sequence" and name
"Chr1".

The object class "Sequence" is used by default when the class is not
specified in the segment() call.  This allows you to use a shortcut
form of the segment() method:

 $segment = $db->segment('Chr1');          # whole chromosome
 $segment = $db->segment('Chr1',1=>1000);  # first 1000 bp

For your convenience, if, during loading a GFF file, Bio::DB::GFF
encounters a line like the following:

  ##sequence-region Chr1 1 14972282

It will automatically generate the following entry:

 Chr1 reference Component 1 14972282 . + . Sequence Chr1

This is sufficient to use Chr1 as a reference point.
The ##sequence-region line is frequently found in the GFF files
distributed by annotation groups.

=head2 Specifying the group tag

A frequent problem with GFF files is the problem distinguishing
which of the several tag/value pairs in the 9th column is the grouping
pair.  Ordinarily the first tag will be used for grouping, but some
GFF manipulating tools do not preserve the order of attributes.  To
eliminate this ambiguity, this module provides two ways of explicitly
specifying which tag to group on:

=over 4

=item *

Using -preferred_groups

When you create a Bio::DB::GFF object, pass it a -preferred_groups=E<gt>
argument.  This specifies a tag that will be used for grouping.  You
can pass an array reference to specify a list of such tags.

=item *

In the GFF header

The GFF file itself can specify which tags are to be used for
grouping.  Insert a comment like the following:

 ##group-tags Accession Locus

This says to use the Accession tag for grouping.  If it is not
available, use the Locus tag.  If neither tag is available, use the
first pair to appear.

=back

These options only apply when B<loading> a GFF file into the database,
and have no effect on existing databases.

The group-tags comment in the GFF file will *override* the preferred
groups set when you create the Bio::DB::GFF object.

For backward compatibility, the tags Sequence and Transcript are
always treated as grouping tags unless preferred_tags are specified.
The "Target" tag is always used for grouping regardless of the
preferred_groups() setting, and the tags "tstart", "tend" and "Note"
cannot be used for grouping.  These are historical artefacts coming
from various interpretations of GFF2, and cannot be changed.

=head2 Sequence alignments

There are two cases in which an annotation indicates the relationship
between two sequences.  The first case is a similarity hit, where the
annotation indicates an alignment.  The second case is a map assembly,
in which the annotation indicates that a portion of a larger sequence
is built up from one or more smaller ones.

Both cases are indicated by using the B<Target> tag in the group
field.  For example, a typical similarity hit will look like this:

 Chr1 BLASTX similarity 76953 77108 132 + 0 Target Protein:SW:ABL_DROME 493 544

The group field contains the Target tag, followed by an identifier for
the biological object referred to.  The GFF format uses the notation
I<Class>:I<Name> for the biological object, and even though this is
stylistically inconsistent, that's the way it's done.  The object
identifier is followed by two integers indicating the start and stop
of the alignment on the target sequence.

Unlike the main start and stop columns, it is possible for the target
start to be greater than the target end.  The previous example
indicates that the the section of Chr1 from 76,953 to 77,108 aligns to
the protein SW:ABL_DROME starting at position 493 and extending to
position 544.

A similar notation is used for sequence assembly information as shown
in this example:

 Chr1        assembly Link   10922906 11177731 . . . Target Sequence:LINK_H06O01 1 254826
 LINK_H06O01 assembly Cosmid 32386    64122    . . . Target Sequence:F49B2       6 31742

This indicates that the region between bases 10922906 and 11177731 of
Chr1 are composed of LINK_H06O01 from bp 1 to bp 254826.  The region
of LINK_H0601 between 32386 and 64122 is, in turn, composed of the
bases 5 to 31742 of cosmid F49B2.

=head2 Attributes

While not intended to serve as a general-purpose sequence database
(see bioperl-db for that), GFF allows you to tag features with
arbitrary attributes.  Attributes appear in the Group field following
the initial class/name pair.  For example:

 Chr1  cur trans  939 942 . +  . Transcript Y95B8A.2 ; Gene sma-3 ; Alias sma3

This line tags the feature named Transcript Y95B8A.2 as being "Gene"
named sma-3 and having the Alias "sma3".  Features having these
attributes can be looked up using the fetch_feature_by_attribute() method.

Two attributes have special meaning: "Note" is for backward
compatibility and is used for unstructured text remarks.  "Alias" is
considered as a synonym for the feature name and will be consulted
when looking up a feature by its name.

=head2 Adaptors and Aggregators

This module uses a system of adaptors and aggregators in order to make
it adaptable to use with a variety of databases.

=over 4

=item *

Adaptors

The core of the module handles the user API, annotation coordinate
arithmetic, and other common issues.  The details of fetching
information from databases is handled by an adaptor, which is
specified during Bio::DB::GFF construction.  The adaptor encapsulates
database-specific information such as the schema, user authentication
and access methods.

There are currently five adaptors recommended for general use:

  Adaptor Name             Description
  ------------             -----------

  memory                   A simple in-memory database suitable for testing
                            and small data sets.

  berkeleydb               An indexed file database based on the DB_File module,
                            suitable for medium-sized read-only data sets.

  dbi::mysql               An interface to a schema implemented in the Mysql
                            relational database management system.

  dbi::oracle              An interface to a schema implemented in the Oracle
                            relational database management system.

  dbi::pg                  An interface to a schema implemented in the PostgreSQL
                            relational database management system.

Check the Bio/DB/GFF/Adaptor directory and subdirectories for other,
more specialized adaptors, as well as experimental ones.

=item *

Aggregators

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

For example, during disaggregation, the standard
"processed_transcript" aggregator generates a list of component
feature types including "UTR", "CDS", and "polyA_site".  Later, it
aggregates these features into a set of annotations of type
"processed_transcript".

During aggregation, the list of aggregators is called in reverse
order.  This allows aggregators to collaborate to create multi-level
structures: the transcript aggregator assembles transcripts from
introns and exons; the gene aggregator then assembles genes from sets
of transcripts.

Three default aggregators are provided:

      transcript   assembles transcripts from features of type
                   exon, CDS, 5'UTR, 3'UTR, TSS, and PolyA
      clone        assembles clones from Clone_left_end, Clone_right_end
                   and Sequence features.
      alignment    assembles gapped alignments from features of type
                   "similarity".

In addition, this module provides the optional "wormbase_gene"
aggregator, which accommodates the WormBase representation of genes.
This aggregator aggregates features of method "exon", "CDS", "5'UTR",
"3'UTR", "polyA" and "TSS" into a single object.  It also expects to
find a single feature of type "Sequence" that spans the entire gene.

The existing aggregators are easily customized.

Note that aggregation will not occur unless you specifically request
the aggregation type.  For example, this call:

  @features = $segment->features('alignment');

will generate an array of aggregated alignment features.  However,
this call:

  @features = $segment->features();

will return a list of unaggregated similarity segments.

For more informnation, see the manual pages for
Bio::DB::GFF::Aggregator::processed_transcript, Bio::DB::GFF::Aggregator::clone,
etc.

=back

=head2 Loading GFF3 Files

This module will accept GFF3 files, as described at
http://song.sourceforge.net/gff3.shtml. However, the implementation
has some limitations.

=over 4

=item GFF version string is required

The GFF file B<must> contain the version comment:

 ##gff-version 3

Unless this version string is present at the top of the GFF file, the
loader will attempt to parse the file in GFF2 format, with
less-than-desirable results.

=item Only one level of nesting allowed

A major restriction is that Bio::DB::GFF only allows one level of
nesting of features.  For nesting, the Target tag will be used
preferentially followed by the ID tag, followed by the Parent tag.
This means that if genes are represented like this:

  XXXX XXXX gene XXXX XXXX XXXX ID=myGene
  XXXX XXXX mRNA XXXX XXXX XXXX ID=myTranscript;Parent=myGene
  XXXX XXXX exon XXXX XXXX XXXX Parent=myTranscript
  XXXX XXXX exon XXXX XXXX XXXX Parent=myTranscript

Then there will be one group called myGene containing the "gene"
feature and one group called myTranscript containing the mRNA, and two
exons.

You can work around this restriction to some extent by using the Alias
attribute literally:

  XXXX XXXX gene XXXX XXXX XXXX ID=myGene
  XXXX XXXX mRNA XXXX XXXX XXXX ID=myTranscript;Parent=myGene;Alias=myGene
  XXXX XXXX exon XXXX XXXX XXXX Parent=myTranscript;Alias=myGene
  XXXX XXXX exon XXXX XXXX XXXX Parent=myTranscript;Alias=myGene

This limitation will be corrected in the next version of Bio::DB::GFF.

=back

=head1 API

The following is the API for Bio::DB::GFF.

=cut

package Bio::DB::GFF;

use strict;

use IO::File;
use File::Glob ':glob';
use Bio::DB::GFF::Util::Rearrange;
use Bio::DB::GFF::RelSegment;
use Bio::DB::GFF::Feature;
use Bio::DB::GFF::Aggregator;

use base qw(Bio::Root::Root Bio::DasI);

my %valid_range_types = (overlaps     => 1,
			 contains     => 1,
			 contained_in => 1);

=head1 Querying GFF Databases

=head2 new

 Title   : new
 Usage   : my $db = Bio::DB::GFF->new(@args);
 Function: create a new Bio::DB::GFF object
 Returns : new Bio::DB::GFF object
 Args    : lists of adaptors and aggregators
 Status  : Public

These are the arguments:

 -adaptor      Name of the adaptor module to use.  If none
               provided, defaults to "dbi::mysqlopt".

 -aggregator   Array reference to a list of aggregators
               to apply to the database.  If none provided,
	       defaults to ['processed_transcript','alignment'].

  -preferred_groups  When interpreteting the 9th column of a GFF2 file,
                 the indicated group names will have preference over
                 other attributes, even if they do not come first in
                 the list of attributes.  This can be a scalar value
                 or an array reference.

  <other>      Any other named argument pairs are passed to
               the adaptor for processing.

The adaptor argument must correspond to a module contained within the
Bio::DB::GFF::Adaptor namespace.  For example, the
Bio::DB::GFF::Adaptor::dbi::mysql adaptor is loaded by specifying
'dbi::mysql'.  By Perl convention, the adaptors names are lower case
because they are loaded at run time.

The aggregator array may contain a list of aggregator names, a list of
initialized aggregator objects, or a string in the form
"aggregator_name{subpart1,subpart2,subpart3/main_method}" (the
"/main_method" part is optional, but if present a feature with the
main_method must be present in order for aggregation to occur).  For
example, if you wish to change the components aggregated by the
transcript aggregator, you could pass it to the GFF constructor this
way:

  my $transcript = 
     Bio::DB::Aggregator::transcript->new(-sub_parts=>[qw(exon intron utr
                                                          polyA spliced_leader)]);

  my $db = Bio::DB::GFF->new(-aggregator=>[$transcript,'clone','alignment],
                             -adaptor   => 'dbi::mysql',
                             -dsn      => 'dbi:mysql:elegans42');

Alternatively, you could create an entirely new transcript aggregator
this way:

  my $new_agg = 'transcript{exon,intron,utr,polyA,spliced_leader}';
  my $db      = Bio::DB::GFF->new(-aggregator=>[$new_agg,'clone','alignment],
                                  -adaptor   => 'dbi::mysql',
                                  -dsn       => 'dbi:mysql:elegans42');

See L<Bio::DB::GFF::Aggregator> for more details.

The B<-preferred_groups> argument is used to change the default
processing of the 9th column of GFF version 2 files.  By default, the
first tag/value pair is used to establish the group class and name.
If you pass -preferred_groups a scalar, the parser will look for a tag
of the indicated type and use it as the group even if it is not first
in the file.  If you pass this argument a list of group classes as an
array ref, then the list will establish the precedence for searching.

The commonly used 'dbi::mysql' adaptor recognizes the following
adaptor-specific arguments:

  Argument       Description
  --------       -----------

  -dsn           the DBI data source, e.g. 'dbi:mysql:ens0040'
                 If a partial name is given, such as "ens0040", the
                 "dbi:mysql:" prefix will be added automatically.

  -user          username for authentication

  -pass          the password for authentication

  -refclass      landmark Class; defaults to "Sequence"


The commonly used 'dbi::mysqlopt' adaptor also recogizes the following
arguments.

  Argument       Description
  --------       -----------

  -fasta         path to a directory containing FASTA files for the DNA
                 contained in this database (e.g. "/usr/local/share/fasta")

  -acedb         an acedb URL to use when converting features into ACEDB
                    objects (e.g. sace://localhost:2005)

=cut

#'

sub new {
  my $package   = shift;
  my ($adaptor,$aggregators,$args,$refclass,$preferred_groups);

  if (@_ == 1) {  # special case, default to dbi::mysqlopt
    $adaptor = 'dbi::mysqlopt';
    $args = {DSN => shift};
  } else {
    ($adaptor,$aggregators,$refclass,$preferred_groups,$args) = rearrange([
									   [qw(ADAPTOR FACTORY)],
									   [qw(AGGREGATOR AGGREGATORS)],
									   'REFCLASS',
									   'PREFERRED_GROUPS'
									  ],@_);
  }

  $adaptor    ||= 'dbi::mysqlopt';
  my $class = "Bio::DB::GFF::Adaptor::\L${adaptor}\E";
  unless ($class->can('new')) {
    eval "require $class;1;" or $package->throw("Unable to load $adaptor adaptor: $@");
  }

  # this hack saves the memory adaptor, which loads the GFF file in new()
  $args->{PREFERRED_GROUPS} = $preferred_groups if defined $preferred_groups;

  my $self = $class->new($args);

  # handle preferred groups
  $self->preferred_groups($preferred_groups) if defined $preferred_groups;
  $self->default_class($refclass || 'Sequence');

  # handle the aggregators.
  # aggregators are responsible for creating complex multi-part features
  # from the GFF "group" field.  If none are provided, then we provide a
  # list of the two used in WormBase.
  # Each aggregator can be a scalar or a ref.  In the former case
  # it is treated as a class name to call new() on.  In the latter
  # the aggreator is treated as a ready made object.
  $aggregators = $self->default_aggregators unless defined $aggregators;
  my @a = ref($aggregators) eq 'ARRAY' ? @$aggregators : $aggregators;
  for my $a (@a) {
    $self->add_aggregator($a);
  }

  # default settings go here.....
  $self->automerge(1);  # set automerge to true

  $self;
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
  my ($refseq,$start,$stop,$enumerate,$refclass,$types) = rearrange ([
								      [qw(REF REFSEQ)],
								      qw(START),
								      [qw(STOP END)],
								      [qw(ENUMERATE COUNT)],
								      [qw(CLASS SEQCLASS)],
								      [qw(TYPE TYPES)],
								     ],@_);
  $types = $self->parse_types($types) if defined $types;
  $self->get_types($refseq,$refclass,$start,$stop,$enumerate,$types);
}

=head2 classes

 Title   : classes
 Usage   : $db->classes
 Function: return list of landmark classes in database
 Returns : a list of classes
 Args    : none
 Status  : public

This routine returns the list of reference classes known to the
database, or empty if classes are not used by the database.  Classes
are distinct from types, being essentially qualifiers on the reference
namespaces.

=cut

sub classes {
  my $self = shift;
  return ();
}

=head2 segment

 Title   : segment
 Usage   : $db->segment(@args);
 Function: create a segment object
 Returns : segment object(s)
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

 -name         ID of the landmark sequence.

 -class        Database object class for the landmark sequence.
               "Sequence" assumed if not specified.  This is
               irrelevant for databases which do not recognize
               object classes.

 -start        Start of the segment relative to landmark.  Positions
               follow standard 1-based sequence rules.  If not specified,
               defaults to the beginning of the landmark.

 -end          Stop of the segment relative to the landmark.  If not specified,
               defaults to the end of the landmark.

 -stop         Same as -end.

 -offset       For those who prefer 0-based indexing, the offset specifies the
               position of the new segment relative to the start of the landmark.

 -length       For those who prefer 0-based indexing, the length specifies the
               length of the new segment.

 -refseq       Specifies the ID of the reference landmark used to establish the
               coordinate system for the newly-created segment.

 -refclass     Specifies the class of the reference landmark, for those databases
               that distinguish different object classes.  Defaults to "Sequence".

 -absolute
               Return features in absolute coordinates rather than relative to the
               parent segment.

 -nocheck      Don't check the database for the coordinates and length of this
               feature.  Construct a segment using the indicated name as the
               reference, a start coordinate of 1, an undefined end coordinate,
               and a strand of +1.

 -force        Same as -nocheck.

 -seq,-sequence,-sourceseq   Aliases for -name.

 -begin,-end   Aliases for -start and -stop

 -off,-len     Aliases for -offset and -length

 -seqclass     Alias for -class

Here's an example to explain how this works:

  my $db = Bio::DB::GFF->new(-dsn => 'dbi:mysql:human',-adaptor=>'dbi::mysql');

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
establish it as the coordinate reference point:

  $segment->refseq($s2);

Now calling start() will return the start of the segment relative to
the beginning of $s2, accounting for differences in strandedness:

  my $rel_start = $segment->start;

IMPORTANT NOTE: This method can be used to return the segment spanned
by an arbitrary named annotation.  However, if the annotation appears
at multiple locations on the genome, for example an EST that maps to
multiple locations, then, provided that all locations reside on the
same physical segment, the method will return a segment that spans the
minimum and maximum positions.  If the reference sequence occupies
ranges on different physical segments, then it returns them all in an
array context, and raises a "multiple segment exception" exception in
a scalar context.

=cut

#'

sub segment {
  my $self = shift;
  my @segments =  Bio::DB::GFF::RelSegment->new(-factory => $self,
						$self->setup_segment_args(@_));
  foreach (@segments) {
    $_->absolute(1) if $self->absolute;
  }

  $self->_multiple_return_args(@segments);
}

sub _multiple_return_args {
  my $self = shift;
  my @args = @_;
  if (@args == 0) {
    return;
  } elsif (@args == 1) {
    return $args[0];
  } elsif (wantarray) { # more than one reference sequence
    return @args;
  } else {
    $self->error($args[0]->name,
		 " has more than one reference sequence in database.  Please call in a list context to retrieve them all.");
    $self->throw('multiple segment exception');
    return;
  }

}

# backward compatibility -- don't use!
# (deliberately undocumented too)
sub abs_segment {
  my $self = shift;
  return $self->segment($self->setup_segment_args(@_),-absolute=>1);
}

sub setup_segment_args {
  my $self = shift;
  return @_ if defined $_[0] && $_[0] =~ /^-/;
  return (-name=>$_[0],-start=>$_[1],-stop=>$_[2]) if @_ == 3;
  return (-class=>$_[0],-name=>$_[1])              if @_ == 2;
  return (-name=>$_[0])                            if @_ == 1;
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

  -rare      Turn on optimizations suitable for a relatively rare feature type,
             where it makes more sense to filter by feature type first,
             and then by position.

  -attributes A hash reference containing attributes to match.

  -iterator  Whether to return an iterator across the features.

  -binsize   A true value will create a set of artificial features whose
             start and stop positions indicate bins of the given size, and
             whose scores are the number of features in the bin.  The
             class and method of the feature will be set to "bin",
             its source to "method:source", and its group to "bin:method:source".
             This is a handy way of generating histograms of feature density.

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

The -attributes argument is a hashref containing one or more attributes
to match against:

  -attributes => { Gene => 'abc-1',
                   Note => 'confirmed' }

Attribute matching is simple string matching, and multiple attributes
are ANDed together.

=cut

sub features {
  my $self = shift;
  my ($types,$automerge,$sparse,$iterator,$refseq,$start,$end,$other);
  if (defined $_[0] && 
      $_[0] =~ /^-/) {
    ($types,$automerge,$sparse,$iterator,
     $refseq,$start,$end,
     $other) = rearrange([
	[qw(TYPE TYPES)],
	[qw(MERGE AUTOMERGE)],
	[qw(RARE SPARSE)],
	'ITERATOR',
	[qw(REFSEQ SEQ_ID)],
	'START',
	[qw(STOP END)],
			 ],@_);
  } else {
    $types = \@_;
  }

  # for whole database retrievals, we probably don't want to automerge!
  $automerge = $self->automerge unless defined $automerge;
  $other ||= {};
  $self->_features({
		    rangetype => $refseq ? 'overlaps' : 'contains',
		    types     => $types,
		    refseq    => $refseq,
		    start     => $start,
		    stop      => $end,
		   },
		   { sparse    => $sparse,
		     automerge => $automerge,
		     iterator  =>$iterator,
		     %$other,
		   }
		   );
}

=head2 get_seq_stream

 Title   : get_seq_stream
 Usage   : my $seqio = $self->get_seq_sream(@args)
 Function: Performs a query and returns an iterator over it
 Returns : a Bio::SeqIO stream capable of producing sequence
 Args    : As in features()
 Status  : public

This routine takes the same arguments as features(), but returns a
Bio::SeqIO::Stream-compliant object.  Use it like this:

  $stream = $db->get_seq_stream('exon');
  while (my $exon = $stream->next_seq) {
     print $exon,"\n";
  }

NOTE: This is also called get_feature_stream(), since that's what it
really does.

=cut

sub get_seq_stream {
  my $self = shift;
  my @args = !defined($_[0]) || $_[0] =~ /^-/ ? (@_,-iterator=>1)
                                              : (-types=>\@_,-iterator=>1);
  $self->features(@args);
}

*get_feature_stream = \&get_seq_stream;

=head2 get_feature_by_name

 Title   : get_feature_by_name
 Usage   : $db->get_feature_by_name($class => $name)
 Function: fetch features by their name
 Returns : a list of Bio::DB::GFF::Feature objects
 Args    : the class and name of the desired feature
 Status  : public

This method can be used to fetch a named feature from the database.
GFF annotations are named using the group class and name fields, so
for features that belong to a group of size one, this method can be
used to retrieve that group (and is equivalent to the segment()
method).  Any Alias attributes are also searched for matching names.

An alternative syntax allows you to search for features by name within
a circumscribed region:

  @f = $db->get_feature_by_name(-class => $class,-name=>$name,
                                -ref   => $sequence_name,
                                -start => $start,
                                -end   => $end);

This method may return zero, one, or several Bio::DB::GFF::Feature
objects.

Aggregation is performed on features as usual.

NOTE: At various times, this function was called fetch_group(),
fetch_feature(), fetch_feature_by_name() and segments().  These names
are preserved for backward compatibility.

=cut

sub get_feature_by_name {
  my $self = shift;
  my ($gclass,$gname,$automerge,$ref,$start,$end);
  if (@_ == 1) {
    $gclass = $self->default_class;
    $gname  = shift;
  } else  {
    ($gclass,$gname,$automerge,$ref,$start,$end) = rearrange(['CLASS','NAME','AUTOMERGE',
							      ['REF','REFSEQ'],
							      'START',['STOP','END']
							     ],@_);
    $gclass ||= $self->default_class;
  }
  $automerge = $self->automerge unless defined $automerge;

  # we need to refactor this... It's repeated code (see below)...
  my @aggregators;
  if ($automerge) {
    for my $a ($self->aggregators) {
      push @aggregators,$a if $a->disaggregate([],$self);
    }
  }

  my %groups;         # cache the groups we create to avoid consuming too much unecessary memory
  my $features = [];
  my $callback = sub { push @$features,$self->make_feature(undef,\%groups,@_) };
  my $location = [$ref,$start,$end] if defined $ref;
  $self->_feature_by_name($gclass,$gname,$location,$callback);

  warn "aggregating...\n" if $self->debug;
  foreach my $a (@aggregators) {  # last aggregator gets first shot
      $a->aggregate($features,$self) or next;
  }

  @$features;
}

# horrible indecision regarding proper names!
*fetch_group   = *fetch_feature = *fetch_feature_by_name = \&get_feature_by_name;
*segments      = \&segment;

=head2 get_feature_by_target

 Title   : get_feature_by_target
 Usage   : $db->get_feature_by_target($class => $name)
 Function: fetch features by their similarity target
 Returns : a list of Bio::DB::GFF::Feature objects
 Args    : the class and name of the desired feature
 Status  : public

This method can be used to fetch a named feature from the database
based on its similarity hit.

=cut

sub get_feature_by_target {
  shift->get_feature_by_name(@_);
}

=head2 get_feature_by_attribute

 Title   : get_feature_by_attribute
 Usage   : $db->get_feature_by_attribute(attribute1=>value1,attribute2=>value2)
 Function: fetch segments by combinations of attribute values
 Returns : a list of Bio::DB::GFF::Feature objects
 Args    : the class and name of the desired feature
 Status  : public

This method can be used to fetch a set of features from the database.
Attributes are a list of name=E<gt>value pairs.  They will be logically
ANDED together.

=cut

sub get_feature_by_attribute {
  my $self = shift;
  my %attributes = ref($_[0]) ? %{$_[0]} : @_;

  # we need to refactor this... It's repeated code (see above)...
  my @aggregators;
  if ($self->automerge) {
    for my $a ($self->aggregators) {
      unshift @aggregators,$a if $a->disaggregate([],$self);
    }
  }

  my %groups;         # cache the groups we create to avoid consuming too much unecessary memory
  my $features = [];
  my $callback = sub { push @$features,$self->make_feature(undef,\%groups,@_) };
  $self->_feature_by_attribute(\%attributes,$callback);

  warn "aggregating...\n" if $self->debug;
  foreach my $a (@aggregators) {  # last aggregator gets first shot
      $a->aggregate($features,$self) or next;
  }

  @$features;
}

# more indecision...
*fetch_feature_by_attribute = \&get_feature_by_attribute;

=head2 get_feature_by_id

 Title   : get_feature_by_id
 Usage   : $db->get_feature_by_id($id)
 Function: fetch segments by feature ID
 Returns : a Bio::DB::GFF::Feature object
 Args    : the feature ID
 Status  : public

This method can be used to fetch a feature from the database using its
ID.  Not all GFF databases support IDs, so be careful with this.

=cut

sub get_feature_by_id {
  my $self = shift;
  my $id   = ref($_[0]) eq 'ARRAY' ? $_[0] : \@_;
  my %groups;         # cache the groups we create to avoid consuming too much unecessary memory
  my $features = [];
  my $callback = sub { push @$features,$self->make_feature(undef,\%groups,@_) };
  $self->_feature_by_id($id,'feature',$callback);
  return wantarray ? @$features : $features->[0];
}
*fetch_feature_by_id = \&get_feature_by_id;

=head2 get_feature_by_gid

 Title   : get_feature_by_gid
 Usage   : $db->get_feature_by_gid($id)
 Function: fetch segments by feature ID
 Returns : a Bio::DB::GFF::Feature object
 Args    : the feature ID
 Status  : public

This method can be used to fetch a feature from the database using its
group ID.  Not all GFF databases support IDs, so be careful with this.

The group ID is often more interesting than the feature ID, since
groups can be complex objects containing subobjects.

=cut

sub get_feature_by_gid {
  my $self = shift;
  my $id   = ref($_[0]) eq 'ARRAY' ? $_[0] : \@_;
  my %groups;         # cache the groups we create to avoid consuming too much unecessary memory
  my $features = [];
  my $callback = sub { push @$features,$self->make_feature(undef,\%groups,@_) };
  $self->_feature_by_id($id,'group',$callback);
  return wantarray ? @$features : $features->[0];
}
*fetch_feature_by_gid = \&get_feature_by_gid;

=head2 delete_fattribute_to_features

 Title   : delete_fattribute_to_features
 Usage   : $db->delete_fattribute_to_features(@ids_or_features)
 Function: delete one or more fattribute_to_features
 Returns : count of fattribute_to_features deleted
 Args    : list of features or feature ids
 Status  : public

Pass this method a list of numeric feature ids or a set of features.  
It will attempt to remove the fattribute_to_features rows of those features
from the database and return a count of the rows removed.  

NOTE: This method is also called delete_fattribute_to_feature().  Also see
delete_groups() and delete_features().

=cut

*delete_fattribute_to_feature = \&delete_fattribute_to_features;

sub delete_fattribute_to_features {
  my $self = shift;
  my @features_or_ids = @_;
  my @ids = map {UNIVERSAL::isa($_,'Bio::DB::GFF::Feature') ? $_->id : $_} @features_or_ids;
  return unless @ids;
  $self->_delete_fattribute_to_features(@ids);
}

=head2 delete_features

 Title   : delete_features
 Usage   : $db->delete_features(@ids_or_features)
 Function: delete one or more features
 Returns : count of features deleted
 Args    : list of features or feature ids
 Status  : public

Pass this method a list of numeric feature ids or a set of features.
It will attempt to remove the features from the database and return a
count of the features removed.  

NOTE: This method is also called delete_feature().  Also see
delete_groups().

=cut

*delete_feature = \&delete_features;

sub delete_features {
  my $self = shift;
  my @features_or_ids = @_;
  my @ids = map {UNIVERSAL::isa($_,'Bio::DB::GFF::Feature') ? $_->id : $_} @features_or_ids;
  return unless @ids;
  $self->_delete_features(@ids);
}

=head2 delete_groups

 Title   : delete_groups
 Usage   : $db->delete_groups(@ids_or_features)
 Function: delete one or more feature groups
 Returns : count of features deleted
 Args    : list of features or feature group ids
 Status  : public

Pass this method a list of numeric group ids or a set of features.  It
will attempt to recursively remove the features and ALL members of
their group from the database.  It returns a count of the number of
features (not groups) returned.

NOTE: This method is also called delete_group().  Also see
delete_features().

=cut

*delete_group = \&delete_groupss;

sub delete_groups {
  my $self = shift;
  my @features_or_ids = @_;
  my @ids = map {UNIVERSAL::isa($_,'Bio::DB::GFF::Feature') ? $_->group_id : $_} @features_or_ids;
  return unless @ids;
  $self->_delete_groups(@ids);
}

=head2 delete

 Title   : delete
 Usage   : $db->delete(@args)
 Function: delete features
 Returns : count of features deleted -- if available
 Args    : numerous, see below
 Status  : public

This method deletes all features that overlap the specified region or
are of a particular type.  If no arguments are provided and the -force
argument is true, then deletes ALL features.

Arguments:

 -name         ID of the landmark sequence.

 -ref          ID of the landmark sequence (synonym for -name).

 -class        Database object class for the landmark sequence.
               "Sequence" assumed if not specified.  This is
               irrelevant for databases which do not recognize
               object classes.

 -start        Start of the segment relative to landmark.  Positions
               follow standard 1-based sequence rules.  If not specified,
               defaults to the beginning of the landmark.

 -end          Stop of the segment relative to the landmark.  If not specified,
               defaults to the end of the landmark.

 -offset       Zero-based addressing

 -length       Length of region

 -type,-types  Either a single scalar type to be deleted, or an
               reference to an array of types.

 -force        Force operation to be performed even if it would delete
               entire feature table.

 -range_type   Control the range type of the deletion.  One of "overlaps" (default)
               "contains" or "contained_in"

Examples:

  $db->delete(-type=>['intron','repeat:repeatMasker']);  # remove all introns & repeats
  $db->delete(-name=>'chr3',-start=>1,-end=>1000);       # remove annotations on chr3 from 1 to 1000
  $db->delete(-name=>'chr3',-type=>'exon');              # remove all exons on chr3

The short form of this call, as described in segment() is also allowed:

  $db->delete("chr3",1=>1000);
  $db->delete("chr3");

IMPORTANT NOTE: This method only deletes features.  It does *NOT*
delete the names of groups that contain the deleted features.  Group
IDs will be reused if you later load a feature with the same group
name as one that was previously deleted.

NOTE ON FEATURE COUNTS: The DBI-based versions of this call return the
result code from the SQL DELETE operation.  Some dbd drivers return the
count of rows deleted, while others return 0E0.  Caveat emptor.

=cut

sub delete {
  my $self = shift;
  my @args = $self->setup_segment_args(@_);
  my ($name,$class,$start,$end,$offset,$length,$type,$force,$range_type) =
    rearrange([['NAME','REF'],'CLASS','START',[qw(END STOP)],'OFFSET',
	       'LENGTH',[qw(TYPE TYPES)],'FORCE','RANGE_TYPE'],@args);
  $offset = 0 unless defined $offset;
  $start = $offset+1 unless defined $start;
  $end   = $start+$length-1 if !defined $end and $length;
  $class ||= $self->default_class;

  my $types = $self->parse_types($type);  # parse out list of types

  $range_type ||= 'overlaps';
  $self->throw("range type must be one of {".
	       join(',',keys %valid_range_types).
	       "}\n")
    unless $valid_range_types{lc $range_type};


  my @segments;
  if (defined $name && $name ne '') {
    my @args = (-name=>$name,-class=>$class);
    push @args,(-start=>$start) if defined $start;
    push @args,(-end  =>$end)   if defined $end;
    @segments = $self->segment(@args);
    return unless @segments;
  }
  $self->_delete({segments   => \@segments,
		  types      => $types,
		  range_type => $range_type,
		  force      => $force}
		);
}

=head2 absolute

 Title   : absolute
 Usage   : $abs = $db->absolute([$abs]);
 Function: gets/sets absolute mode
 Returns : current setting of absolute mode boolean
 Args    : new setting for absolute mode boolean
 Status  : public

$db-E<gt>absolute(1) will turn on absolute mode for the entire database.
All segments retrieved will use absolute coordinates by default,
rather than relative coordinates.  You can still set them to use
relative coordinates by calling $segment-E<gt>absolute(0).

Note that this is not the same as calling abs_segment(); it continues
to allow you to look up groups that are not used directly as reference
sequences.

=cut

sub absolute {
  my $self = shift;
  my $d = $self->{absolute};
  $self->{absolute} = shift if @_;
  $d;
}

=head2 strict_bounds_checking

 Title   : strict_bounds_checking
 Usage   : $flag = $db->strict_bounds_checking([$flag])
 Function: gets/sets strict bounds checking
 Returns : current setting of bounds checking flag
 Args    : new setting for bounds checking flag
 Status  : public

This flag enables extra checks for segment requests that go beyond the
ends of their reference sequences.  If bounds checking is enabled,
then retrieved segments will be truncated to their physical length,
and their truncated() methods will return true.

If the flag is off (the default), then the module will return segments
that appear to extend beyond their physical boundaries.  Requests for
features beyond the end of the segment will, however, return empty.

=cut

sub strict_bounds_checking {
  my $self = shift;
  my $d = $self->{strict};
  $self->{strict} = shift if @_;
  $d;
}

=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : $seq = $db->get_Seq_by_id('ROA1_HUMAN')
 Function: Gets a Bio::Seq object by its name
 Returns : a Bio::Seq object
 Args    : the id (as a string) of a sequence
 Throws  : "id does not exist" exception

NOTE: Bio::DB::RandomAccessI compliant method

=cut

sub  get_Seq_by_id {
  my $self = shift;
  $self->get_feature_by_name(@_);
}


=head2 get_Seq_by_accession

 Title   : get_Seq_by_accession
 Usage   : $seq = $db->get_Seq_by_accession('AL12234')
 Function: Gets a Bio::Seq object by its accession
 Returns : a Bio::Seq object
 Args    : the id (as a string) of a sequence
 Throws  : "id does not exist" exception

NOTE: Bio::DB::RandomAccessI compliant method

=cut

sub  get_Seq_by_accession {
  my $self = shift;
  $self->get_feature_by_name(@_);
}

=head2 get_Seq_by_acc

 Title   : get_Seq_by_acc
 Usage   : $seq = $db->get_Seq_by_acc('X77802');
 Function: Gets a Bio::Seq object by accession number
 Returns : A Bio::Seq object
 Args    : accession number (as a string)
 Throws  : "acc does not exist" exception

NOTE: Bio::DB::RandomAccessI compliant method

=cut

sub  get_Seq_by_acc {
  my $self = shift;
  $self->get_feature_by_name(@_);
}

=head2 get_Stream_by_name

  Title   : get_Stream_by_name
  Usage   : $seq = $db->get_Stream_by_name(@ids);
  Function: Retrieves a stream of Seq objects given their names
  Returns : a Bio::SeqIO stream object
  Args    : an array of unique ids/accession numbers, or 
            an array reference

NOTE: This is also called get_Stream_by_batch()

=cut

sub get_Stream_by_name {
  my $self = shift;
  my @ids  = @_;
  my $id = ref($ids[0]) ? $ids[0] : \@ids;
  Bio::DB::GFF::ID_Iterator->new($self,$id,'name');
}

=head2 get_Stream_by_id

  Title   : get_Stream_by_id
  Usage   : $seq = $db->get_Stream_by_id(@ids);
  Function: Retrieves a stream of Seq objects given their ids
  Returns : a Bio::SeqIO stream object
  Args    : an array of unique ids/accession numbers, or 
            an array reference

NOTE: This is also called get_Stream_by_batch()

=cut

sub get_Stream_by_id {
  my $self = shift;
  my @ids  = @_;
  my $id = ref($ids[0]) ? $ids[0] : \@ids;
  Bio::DB::GFF::ID_Iterator->new($self,$id,'feature');
}

=head2 get_Stream_by_batch ()

  Title   : get_Stream_by_batch
  Usage   : $seq = $db->get_Stream_by_batch(@ids);
  Function: Retrieves a stream of Seq objects given their ids
  Returns : a Bio::SeqIO stream object
  Args    : an array of unique ids/accession numbers, or 
            an array reference

NOTE: This is the same as get_Stream_by_id().

=cut

*get_Stream_by_batch = \&get_Stream_by_id;


=head2 get_Stream_by_group ()

Bioperl compatibility.

=cut

sub get_Stream_by_group {
  my $self = shift;
  my @ids  = @_;
  my $id = ref($ids[0]) ? $ids[0] : \@ids;
  Bio::DB::GFF::ID_Iterator->new($self,$id,'group');
}

=head2 all_seqfeatures

 Title   : all_seqfeatures
 Usage   : @features = $db->all_seqfeatures(@args)
 Function: fetch all the features in the database
 Returns : an array of features, or an iterator
 Args    : See below
 Status  : public

This is equivalent to calling $db-E<gt>features() without any types, and
will return all the features in the database.  The -merge and
-iterator arguments are recognized, and behave the same as described
for features().

=cut

sub all_seqfeatures {
  my $self = shift;
  my ($automerge,$iterator)= rearrange([
					[qw(MERGE AUTOMERGE)],
					'ITERATOR'
				       ],@_);
  my @args;
  push @args,(-merge=>$automerge)   if defined $automerge;
  push @args,(-iterator=>$iterator) if defined $iterator;
  $self->features(@args);
}

=head1 Creating and Loading GFF Databases

=head2 initialize

 Title   : initialize
 Usage   : $db->initialize(-erase=>$erase,-option1=>value1,-option2=>value2);
 Function: initialize a GFF database
 Returns : true if initialization successful
 Args    : a set of named parameters
 Status  : Public

This method can be used to initialize an empty database.  It takes the following
named arguments:

  -erase     A boolean value.  If true the database will be wiped clean if it
             already contains data.

Other named arguments may be recognized by subclasses.  They become database
meta values that control various settable options.

As a shortcut (and for backward compatibility) a single true argument
is the same as initialize(-erase=E<gt>1).

=cut

sub initialize {
  my $self = shift;

  my ($erase,$meta) = rearrange(['ERASE'],@_);
  $meta ||= {};

  # initialize (possibly erasing)
  return unless $self->do_initialize($erase);
  my @default = $self->default_meta_values;

  # this is an awkward way of uppercasing the 
  # even-numbered values (necessary for case-insensitive SQL databases)
  for (my $i=0; $i<@default; $i++) {
    $default[$i] = uc $default[$i] if !($i % 2);
  }

  my %values = (@default,%$meta);
  foreach (keys %values) {
    $self->meta($_ => $values{$_});
  }
  1;
}


=head2 load_gff

 Title   : load_gff
 Usage   : $db->load_gff($file|$directory|$filehandle [,$verbose]);
 Function: load GFF data into database
 Returns : count of records loaded
 Args    : a directory, a file, a list of files, 
           or a filehandle
 Status  : Public

This method takes a single overloaded argument, which can be any of:

=over 4

=item *

a scalar corresponding to a GFF file on the system

A pathname to a local GFF file.  Any files ending with the .gz, .Z, or
.bz2 suffixes will be transparently decompressed with the appropriate
command-line utility.

=item *

an array reference containing a list of GFF files on the system

For example ['/home/gff/gff1.gz','/home/gff/gff2.gz']

=item *

directory path

The indicated directory will be searched for all files ending in the
suffixes .gff, .gff.gz, .gff.Z or .gff.bz2.

=item *

filehandle

An open filehandle from which to read the GFF data.  Tied filehandles
now work as well.

=item *

a pipe expression

A pipe expression will also work. For example, a GFF file on a remote
web server can be loaded with an expression like this:

  $db->load_gff("lynx -dump -source http://stein.cshl.org/gff_test |");

=back

The optional second argument, if true, will turn on verbose status
reports that indicate the progress.

If successful, the method will return the number of GFF lines
successfully loaded.

NOTE:this method used to be called load(), but has been changed.  The
old method name is also recognized.

=cut

sub load_gff {
  my $self              = shift;
  my $file_or_directory = shift || '.';
  my $verbose           = shift;

  local $self->{__verbose__} = $verbose;
  return $self->do_load_gff($file_or_directory) if ref($file_or_directory) 
                                                   && tied *$file_or_directory;

  my $tied_stdin = tied(*STDIN);
  open my $SAVEIN, "<&STDIN" unless $tied_stdin;
  local @ARGV = $self->setup_argv($file_or_directory,'gff','gff3') or return;  # to play tricks with reader
  my $result = $self->do_load_gff('ARGV');
  open STDIN, '<', $SAVEIN unless $tied_stdin;  # restore STDIN
  return $result;
}

*load = \&load_gff;

=head2 load_gff_file

 Title   : load_gff_file
 Usage   : $db->load_gff_file($file [,$verbose]);
 Function: load GFF data into database
 Returns : count of records loaded
 Args    : a path to a file
 Status  : Public

This is provided as an alternative to load_gff_file. It doesn't munge
STDIN or play tricks with ARGV.

=cut

sub load_gff_file {
  my $self     = shift;
  my $file     = shift;
  my $verbose  = shift;
  my $fh = IO::File->new($file) or return;
  return $self->do_load_gff($fh);
}

=head2 load_fasta

 Title   : load_fasta
 Usage   : $db->load_fasta($file|$directory|$filehandle);
 Function: load FASTA data into database
 Returns : count of records loaded
 Args    : a directory, a file, a list of files, 
           or a filehandle
 Status  : Public

This method takes a single overloaded argument, which can be any of:

=over 4

=item *

scalar corresponding to a FASTA file on the system

A pathname to a local FASTA file.  Any files ending with the .gz, .Z, or
.bz2 suffixes will be transparently decompressed with the appropriate
command-line utility.

=item *

array reference containing a list of FASTA files on the
system

For example ['/home/fasta/genomic.fa.gz','/home/fasta/genomic.fa.gz']

=item *

path to a directory

The indicated directory will be searched for all files ending in the
suffixes .fa, .fa.gz, .fa.Z or .fa.bz2.

=item *

filehandle

An open filehandle from which to read the FASTA data.

=item *

pipe expression

A pipe expression will also work. For example, a FASTA file on a remote
web server can be loaded with an expression like this:

  $db->load_gff("lynx -dump -source http://stein.cshl.org/fasta_test.fa |");

=back

=cut

sub load_fasta {
  my $self              = shift;
  my $file_or_directory = shift || '.';
  my $verbose           = shift;

  local $self->{__verbose__} = $verbose;
  return $self->load_sequence($file_or_directory) if ref($file_or_directory)
                                                     && tied *$file_or_directory;

  my $tied = tied(*STDIN);
  open my $SAVEIN, "<&STDIN" unless $tied;
  local @ARGV = $self->setup_argv($file_or_directory,'fa','dna','fasta') or return;  # to play tricks with reader
  my $result = $self->load_sequence('ARGV');
  open STDIN, '<', $SAVEIN unless $tied;  # restore STDIN
  return $result;
}


=head2 load_fasta_file

 Title   : load_fasta_file
 Usage   : $db->load_fasta_file($file [,$verbose]);
 Function: load FASTA data into database
 Returns : count of records loaded
 Args    : a path to a file
 Status  : Public

This is provided as an alternative to load_fasta. It doesn't munge
STDIN or play tricks with ARGV.

=cut

sub load_fasta_file {
  my $self     = shift;
  my $file     = shift;
  my $verbose  = shift;
  my $fh = IO::File->new($file) or return;
  return $self->do_load_fasta($fh);
}


=head2 load_sequence_string

 Title   : load_sequence_string
 Usage   : $db->load_sequence_string($id,$dna)
 Function: load a single DNA entry
 Returns : true if successfully loaded
 Args    : a raw sequence string (DNA, RNA, protein)
 Status  : Public

=cut

sub load_sequence_string {
  my $self = shift;
  my ($acc,$seq)  = @_;
  my $offset = 0;
  $self->insert_sequence_chunk($acc,\$offset,\$seq) or return;
  $self->insert_sequence($acc,$offset,$seq) or return;
  1;
}

sub setup_argv {
  my $self = shift;
  my $file_or_directory = shift;
  my @suffixes          = @_;
  no strict 'refs';  # so that we can call fileno() on the argument

  my @argv;

  if (-d $file_or_directory) {
    # Because glob() is broken with long file names that contain spaces
    $file_or_directory = Win32::GetShortPathName($file_or_directory)
      if $^O =~ /^MSWin/i && eval 'use Win32; 1';
    @argv = map { glob("$file_or_directory/*.{$_,$_.gz,$_.Z,$_.bz2}")} @suffixes;
  }elsif (my $fd = fileno($file_or_directory)) {
    open STDIN,"<&=$fd" or $self->throw("Can't dup STDIN");
    @argv = '-';
  } elsif (ref $file_or_directory) {
    @argv = @$file_or_directory;
  } else {
    @argv = $file_or_directory;
  }

  foreach (@argv) {
    if (/\.gz$/) {
      $_ = "gunzip -c $_ |";
    } elsif (/\.Z$/) {
      $_ = "uncompress -c $_ |";
    } elsif (/\.bz2$/) {
      $_ = "bunzip2 -c $_ |";
    }
  }
  @argv;
}

=head2 lock_on_load

 Title   : lock_on_load
 Usage   : $lock = $db->lock_on_load([$lock])
 Function: set write locking during load
 Returns : current value of lock-on-load flag
 Args    : new value of lock-on-load-flag
 Status  : Public

This method is honored by some of the adaptors.  If the value is true,
the tables used by the GFF modules will be locked for writing during
loads and inaccessible to other processes.

=cut

sub lock_on_load {
  my $self = shift;
  my $d = $self->{lock};
  $self->{lock} = shift if @_;
  $d;
}

=head2 meta

 Title   : meta
 Usage   : $value = $db->meta($name [,$newval])
 Function: get or set a meta variable
 Returns : a string
 Args    : meta variable name and optionally value
 Status  : abstract

Get or set a named metavalues for the database.  Metavalues can be
used for database-specific settings.

By default, this method does nothing!

=cut

sub meta {
  my $self = shift;
  my ($name,$value) = @_;
  return;
}

=head2 default_meta_values

 Title   : default_meta_values
 Usage   : %values = $db->default_meta_values
 Function: empty the database
 Returns : a list of tag=>value pairs
 Args    : none
 Status  : protected

This method returns a list of tag=E<gt>value pairs that contain default
meta information about the database.  It is invoked by initialize() to
write out the default meta values.  The base class version returns an
empty list.

For things to work properly, meta value names must be UPPERCASE.

=cut

sub default_meta_values {
  my $self = shift;
  return ();
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
  $self->{error} = join '',@_ if @_;
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


=head2 automerge

 Title   : automerge
 Usage   : $db->automerge( [$new automerge] );
 Function: get or set automerge value
 Returns : current value (boolean)
 Args    : an optional argument to set the automerge value
 Status  : Public

By default, this module will use the aggregators to merge groups into
single composite objects.  This default can be changed to false by
calling automerge(0).

=cut

sub automerge {
  my $self = shift;
  my $g = $self->{automerge};
  $self->{automerge} = shift if @_;
  $g;
}

=head2 attributes

 Title   : attributes
 Usage   : @attributes = $db->attributes($id,$name)
 Function: get the "attributes" on a particular feature
 Returns : an array of string
 Args    : feature ID
 Status  : public

Some GFF version 2 files use the groups column to store a series of
attribute/value pairs.  In this interpretation of GFF, the first such
pair is treated as the primary group for the feature; subsequent pairs
are treated as attributes.  Two attributes have special meaning:
"Note" is for backward compatibility and is used for unstructured text
remarks.  "Alias" is considered as a synonym for the feature name.

If no name is provided, then attributes() returns a flattened hash, of
attribute=E<gt>value pairs.  This lets you do:

  %attributes = $db->attributes($id);

If no arguments are provided, attributes() will return the list of
all attribute names:

  @attribute_names = $db->attributes();

Normally, however, attributes() will be called by the feature:

  @notes = $feature->attributes('Note');

In a scalar context, attributes() returns the first value of the
attribute if a tag is present, otherwise a hash reference in which the
keys are attribute names and the values are anonymous arrays
containing the values.

=cut

sub attributes {
  my $self = shift;
  my ($id,$tag) = @_;
  my @result = $self->do_attributes(@_) or return;
  return @result if wantarray;

  # what to do in an array context
  return $result[0] if $tag;
  my %result;
  while (my($key,$value) = splice(@result,0,2)) {
     push @{$result{$key}},$value;
  }
  return \%result;
}

=head2 fast_queries

 Title   : fast_queries
 Usage   : $flag = $db->fast_queries([$flag])
 Function: turn on and off the "fast queries" option
 Returns : a boolean
 Args    : a boolean flag (optional)
 Status  : public

The mysql database driver (and possibly others) support a "fast" query
mode that caches results on the server side.  This makes queries come
back faster, particularly when creating iterators.  The downside is
that while iterating, new queries will die with a "command synch"
error.  This method turns the feature on and off.

For databases that do not support a fast query, this method has no
effect.

=cut

# override this method in order to set the mysql_use_result attribute, which is an obscure
# but extremely powerful optimization for both performance and memory.
sub fast_queries {
  my $self = shift;
  my $d = $self->{fast_queries};
  $self->{fast_queries} = shift if @_;
  $d;
}

=head2 add_aggregator

 Title   : add_aggregator
 Usage   : $db->add_aggregator($aggregator)
 Function: add an aggregator to the list
 Returns : nothing
 Args    : an aggregator
 Status  : public

This method will append an aggregator to the end of the list of
registered aggregators.  Three different argument types are accepted:

  1) a Bio::DB::GFF::Aggregator object -- will be added
  2) a string in the form "aggregator_name{subpart1,subpart2,subpart3/main_method}"
         -- will be turned into a Bio::DB::GFF::Aggregator object (the /main_method
        part is optional).
  3) a valid Perl token -- will be turned into a Bio::DB::GFF::Aggregator
        subclass, where the token corresponds to the subclass name.

=cut

sub add_aggregator {
  my $self       = shift;
  my $aggregator = shift;
  my $list = $self->{aggregators} ||= [];
  if (ref $aggregator) { # an object
    @$list = grep {$_->get_method ne $aggregator->get_method} @$list;
    push @$list,$aggregator;
  }

  elsif ($aggregator =~ /^(\w+)\{([^\/\}]+)\/?(.*)\}$/) {
    my($agg_name,$subparts,$mainpart) = ($1,$2,$3);
    my @subparts = split /,\s*/,$subparts;
    my @args = (-method      => $agg_name,
		-sub_parts   => \@subparts);
    if ($mainpart) {
      push @args,(-main_method => $mainpart,
		  -whole_object => 1);
    }
    warn "making an aggregator with (@args), subparts = @subparts" if $self->debug;
    push @$list,Bio::DB::GFF::Aggregator->new(@args);
  }

  else {
    my $class = "Bio::DB::GFF::Aggregator::\L${aggregator}\E";
    eval "require $class; 1" or  $self->throw("Unable to load $aggregator aggregator: $@");
    push @$list,$class->new();
  }
}

=head2 aggregators

 Title   : aggregators
 Usage   : $db->aggregators([@new_aggregators]);
 Function: retrieve list of aggregators
 Returns : list of aggregators
 Args    : a list of aggregators to set (optional)
 Status  : public

This method will get or set the list of aggregators assigned to
the database.  If 1 or more arguments are passed, the existing
set will be cleared.

=cut

sub aggregators {
  my $self = shift;
  my $d = $self->{aggregators};
  if (@_) {
    $self->clear_aggregators;
    $self->add_aggregator($_) foreach @_;
  }
  return unless $d;
  return @$d;
}

=head2 clear_aggregators

 Title   : clear_aggregators
 Usage   : $db->clear_aggregators
 Function: clears list of aggregators
 Returns : nothing
 Args    : none
 Status  : public

This method will clear the aggregators stored in the database object.
Use aggregators() or add_aggregator() to add some back.

=cut

sub clear_aggregators { shift->{aggregators} = [] }

=head2 preferred_groups

 Title   : preferred_groups
 Usage   : $db->preferred_groups([$group_name_or_arrayref])
 Function: get/set list of groups for altering GFF2 parsing
 Returns : a list of classes
 Args    : new list (scalar or array ref)
 Status  : public

=cut

sub preferred_groups {
  my $self = shift;
  my $d    = $self->{preferred_groups};
  if (@_) {
    my @v = map {ref($_) eq 'ARRAY' ? @$_ : $_} @_;
    $self->{preferred_groups} = \@v;
    delete $self->{preferred_groups_hash};
  }
  return unless $d;
  return @$d;
}

sub _preferred_groups_hash {
  my $self = shift;
  my $gff3 = shift;
  return $self->{preferred_groups_hash} if exists $self->{preferred_groups_hash};
  my $count = 0;

  my @preferred = $self->preferred_groups;

  # defaults
  if (!@preferred) {
    @preferred = $gff3 || $self->{load_data}{gff3_flag} ? qw(Target Parent ID) : qw(Target Sequence Transcript);
  }

  my %preferred = map {lc($_) => @preferred-$count++} @preferred;
  return $self->{preferred_groups_hash} = \%preferred;
}

=head1 Methods for use by Subclasses

The following methods are chiefly of interest to subclasses and are
not intended for use by end programmers.

=head2 abscoords

 Title   : abscoords
 Usage   : $db->abscoords($name,$class,$refseq)
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

If $refseq is present in the argument list, it forces the query to
search for the landmark in a particular reference sequence.

=cut

sub abscoords {
  my $self = shift;
  my ($name,$class,$refseq) = @_;
  $class ||= $self->{default_class};
  $self->get_abscoords($name,$class,$refseq);
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
  return ['processed_transcript','alignment'];
}

=head2 do_load_gff

 Title   : do_load_gff
 Usage   : $db->do_load_gff($handle)
 Function: load a GFF input stream
 Returns : number of features loaded
 Args    : A filehandle.
 Status  : protected

This method is called to load a GFF data stream.  The method will read
GFF features from E<lt>E<gt> and load them into the database.  On exit the
method must return the number of features loaded.

Note that the method is responsible for parsing the GFF lines.  This
is to allow for differences in the interpretation of the "group"
field, which are legion.

You probably want to use load_gff() instead.  It is more flexible
about the arguments it accepts.

=cut

sub do_load_gff {
  my $self      = shift;
  my $io_handle = shift;

  local $self->{load_data} = {
			      lineend => (-t STDERR && !$ENV{EMACS} ? "\r" : "\n"),
			      count   => 0
			     };

  $self->setup_load();
  my $mode = 'gff';

  while (<$io_handle>) {
    chomp;
    if ($mode eq 'gff') {
      if (/^>/) {    # Sequence coming
	$mode = 'fasta';
	$self->_load_sequence_start;
	$self->_load_sequence_line($_);
      } else {
	$self->_load_gff_line($_);
      }
    }
    elsif ($mode eq 'fasta') {
      if (/^##|\t/) {    # Back to GFF mode
	$self->_load_sequence_finish;
	$mode = 'gff';
	$self->_load_gff_line($_);
      } else {
	$self->_load_sequence_line($_);
      }
    }
  }
  $self->finish_load();
  $self->_load_sequence_finish;

  return $self->{load_data}{count};
}

sub _load_gff_line {
  my $self = shift;
  my $line = shift;
  my $lineend = $self->{load_data}{lineend};

  $self->{load_data}{gff3_flag}++           if $line =~ /^\#\#\s*gff-version\s+3/;

  if (defined $self->{load_data}{gff3_flag} and !defined $self->{load_data}{gff3_warning}) {
    $self->print_gff3_warning();
    $self->{load_data}{gff3_warning}=1;
  }

  $self->preferred_groups(split(/\s+/,$1))  if $line =~ /^\#\#\s*group-tags?\s+(.+)/;

  if ($line =~ /^\#\#\s*sequence-region\s+(\S+)\s+(-?\d+)\s+(-?\d+)/i) { # header line
    $self->load_gff_line(
			 {
			  ref    => $1,
			  class  => 'Sequence',
			  source => 'reference',
			  method => 'Component',
			  start  => $2,
			  stop   => $3,
			  score  => undef,
			  strand => undef,
			  phase  => undef,
			  gclass => 'Sequence',
			  gname  => $1,
			  tstart => undef,
			  tstop  => undef,
			  attributes  => [],
			 }
			);
    return $self->{load_data}{count}++;
  }

  return if /^#/;

  my ($ref,$source,$method,$start,$stop,$score,$strand,$phase,$group) = split "\t",$line;
  return unless defined($ref) && defined($method) && defined($start) && defined($stop);
  foreach (\$score,\$strand,\$phase) {
    undef $$_ if $$_ eq '.';
  }

  my ($gclass,$gname,$tstart,$tstop,$attributes) = $self->split_group($group,$self->{load_data}{gff3_flag});

  # no standard way in the GFF file to denote the class of the reference sequence -- drat!
  # so we invoke the factory to do it
  my $class = $self->refclass($ref);

  # call subclass to do the dirty work
  if ($start > $stop) {
    ($start,$stop) = ($stop,$start);
    if ($strand eq '+') {
      $strand = '-';
    } elsif ($strand eq '-') {
      $strand = '+';
    }
  }
  # GFF2/3 transition stuff
  $gclass = [$gclass] unless ref $gclass;
  $gname  = [$gname]  unless ref $gname;
  for (my $i=0; $i<@$gname;$i++) {
    $self->load_gff_line({ref    => $ref,
			  class  => $class,
			  source => $source,
			  method => $method,
			  start  => $start,
			  stop   => $stop,
			  score  => $score,
			  strand => $strand,
			  phase  => $phase,
			  gclass => $gclass->[$i],
			  gname  => $gname->[$i],
			  tstart => $tstart,
			  tstop  => $tstop,
			  attributes  => $attributes}
			);
    $self->{load_data}{count}++;
  }
}

sub _load_sequence_start {
  my $self = shift;
  my $ld   = $self->{load_data};
  undef $ld->{id};
  $ld->{offset} = 0;
  $ld->{seq}    = '';
}
sub _load_sequence_finish {
  my $self = shift;
  my $ld   = $self->{load_data};
  $self->insert_sequence($ld->{id},$ld->{offset},$ld->{seq}) if defined $ld->{id};
}

sub _load_sequence_line {
  my $self = shift;
  my $line = shift;
  my $ld   = $self->{load_data};
  my $lineend = $ld->{lineend};

  if (/^>(\S+)/) {
    $self->insert_sequence($ld->{id},$ld->{offset},$ld->{seq}) if defined $ld->{id};
    $ld->{id}     = $1;
    $ld->{offset} = 0;
    $ld->{seq}    = '';
    $ld->{count}++;
    print STDERR $ld->{count}," sequences loaded$lineend" if $self->{__verbose__} && $ld->{count} % 1000 == 0;
  } else {
    $ld->{seq} .= $_;
    $self->insert_sequence_chunk($ld->{id},\$ld->{offset},\$ld->{seq});
  }

}

=head2 load_sequence

 Title   : load_sequence
 Usage   : $db->load_sequence($handle)
 Function: load a FASTA data stream
 Returns : number of sequences
 Args    : a filehandle to the FASTA file
 Status  : protected

You probably want to use load_fasta() instead.

=cut

# note - there is some repeated code here
sub load_sequence {
  my $self = shift;
  my $io_handle = shift;

  local $self->{load_data} = {
			      lineend => (-t STDERR && !$ENV{EMACS} ? "\r" : "\n"),
			      count   => 0
			     };

  $self->_load_sequence_start;
  while (<$io_handle>) {
    chomp;
    $self->_load_sequence_line($_);
  }
  $self->_load_sequence_finish;
  return $self->{load_data}{count};
}

sub insert_sequence_chunk {
  my $self = shift;
  my ($id,$offsetp,$seqp) = @_;
  if (my $cs = $self->dna_chunk_size) {
    while (length($$seqp) >= $cs) {
      my $chunk = substr($$seqp,0,$cs);
      $self->insert_sequence($id,$$offsetp,$chunk);
      $$offsetp += length($chunk);
      substr($$seqp,0,$cs) = '';
    }
  }
  return 1;  # the calling routine may expect success or failure
}

# used to store big pieces of DNA in itty bitty pieces
sub dna_chunk_size {
  return 0;
}

sub insert_sequence {
  my $self = shift;
  my($id,$offset,$seq) = @_;
  $self->throw('insert_sequence(): must be defined in subclass');
}

# This is the default class for reference points.  Defaults to Sequence.
sub default_class {
   my $self = shift;
   return 'Sequence' unless ref $self;
   my $d = $self->{default_class};
   $self->{default_class} = shift if @_;
   $d;
}

# gets name of the reference sequence, and returns its class
# currently just calls default_class
sub refclass {
  my $self = shift;
  my $name = shift;
  return $self->default_class;
}

=head2 setup_load

 Title   : setup_load
 Usage   : $db->setup_load
 Function: called before load_gff_line()
 Returns : void
 Args    : none
 Status  : abstract

This abstract method gives subclasses a chance to do any
schema-specific initialization prior to loading a set of GFF records.
It must be implemented by a subclass.

=cut

sub setup_load {
  # default, do nothing
}

=head2 finish_load

 Title   : finish_load
 Usage   : $db->finish_load
 Function: called after load_gff_line()
 Returns : number of records loaded
 Args    : none
 Status  :abstract

This method gives subclasses a chance to do any schema-specific
cleanup after loading a set of GFF records.

=cut

sub finish_load {
  # default, do nothing
}

=head2 load_gff_line

 Title   : load_gff_line
 Usage   : $db->load_gff_line(@args)
 Function: called to load one parsed line of GFF
 Returns : true if successfully inserted
 Args    : see below
 Status  : abstract

This abstract method is called once per line of the GFF and passed a
hashref containing parsed GFF fields.  The fields are:

 {ref    => $ref,
  class  => $class,
  source => $source,
  method => $method,
  start  => $start,
  stop   => $stop,
  score  => $score,
  strand => $strand,
  phase  => $phase,
  gclass => $gclass,
  gname  => $gname,
  tstart => $tstart,
  tstop  => $tstop,
  attributes  => $attributes}

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

=head2 dna

 Title   : dna
 Usage   : $db->dna($id,$start,$stop,$class)
 Function: return the raw DNA string for a segment
 Returns : a raw DNA string
 Args    : id of the sequence, its class, start and stop positions
 Status  : public

This method is invoked by Bio::DB::GFF::Segment to fetch the raw DNA
sequence.

Arguments: -name          sequence name
           -start         start position
           -stop          stop position
           -class         sequence class

If start and stop are both undef, then the entire DNA is retrieved.
So to fetch the whole dna, call like this:

  $db->dna($name_of_sequence);

or like this:

   $db->dna(-name=>$name_of_sequence,-class=>$class_of_sequence);

NOTE: you will probably prefer to create a Segment and then invoke its
dna() method.

=cut

# call to return the DNA string for the indicated region
# real work is done by get_dna()
sub dna {
  my $self = shift;
  my ($id,$start,$stop,$class)  = rearrange([
					     [qw(NAME ID REF REFSEQ)],
					     qw(START),
					     [qw(STOP END)],
    					    'CLASS',
					   ],@_);
# return unless defined $start && defined $stop;
  $self->get_dna($id,$start,$stop,$class);
}

sub fetch_sequence { shift->dna(@_) } 

sub features_in_range {
  my $self = shift;
  my ($range_type,$refseq,$class,$start,$stop,$types,$parent,$sparse,$automerge,$iterator,$other) =
    rearrange([
	       [qw(RANGE_TYPE)],
	       [qw(REF REFSEQ)],
	       qw(CLASS),
	       qw(START),
	       [qw(STOP END)],
	       [qw(TYPE TYPES)],
	       qw(PARENT),
	       [qw(RARE SPARSE)],
	       [qw(MERGE AUTOMERGE)],
	       'ITERATOR'
	      ],@_);
  $other ||= {};
  # $automerge = $types && $self->automerge unless defined $automerge;
  $automerge = $self->automerge unless defined $automerge;
  $self->throw("range type must be one of {".
	       join(',',keys %valid_range_types).
	       "}\n")
    unless $valid_range_types{lc $range_type};
  $self->_features({
		    rangetype => lc $range_type,
		    refseq    => $refseq,
		    refclass  => $class,
		    start     => $start,
		    stop      => $stop,
		    types     => $types },
		   {
		    sparse    => $sparse,
		    automerge => $automerge,
		    iterator  => $iterator,
		    %$other,
		   },
		   $parent);
}

=head2 get_dna

 Title   : get_dna
 Usage   : $db->get_dna($id,$start,$stop,$class)
 Function: get DNA for indicated segment
 Returns : the dna string
 Args    : sequence ID, start, stop and class
 Status  : protected

If start E<gt> stop and the sequence is nucleotide, then this method
should return the reverse complement.  The sequence class may be
ignored by those databases that do not recognize different object
types.

=cut

sub get_dna {
  my $self = shift;
  my ($id,$start,$stop,$class,) = @_;
  $self->throw("get_dna() must be implemented by an adaptor");
}

=head2 get_features

 Title   : get_features
 Usage   : $db->get_features($search,$options,$callback)
 Function: get list of features for a region
 Returns : count of number of features retrieved
 Args    : see below
 Status  : protected

The first argument is a hash reference containing search criteria for
retrieving features.  It contains the following keys:

   rangetype One of "overlaps", "contains" or "contained_in".  Indicates
              the type of range query requested.

   refseq    ID of the landmark that establishes the absolute 
              coordinate system.

   refclass  Class of this landmark.  Can be ignored by implementations
              that don't recognize such distinctions.

   start     Start of the range, inclusive.

   stop      Stop of the range, inclusive.

   types     Array reference containing the list of annotation types
              to fetch from the database.  Each annotation type is an
              array reference consisting of [source,method].

The second argument is a hash reference containing certain options
that affect the way information is retrieved:

   sort_by_group
             A flag.  If true, means that the returned features should be
             sorted by the group that they're in.

   sparse    A flag.  If true, means that the expected density of the 
             features is such that it will be more efficient to search
             by type rather than by range.  If it is taking a long
             time to fetch features, give this a try.

   binsize   A true value will create a set of artificial features whose
             start and stop positions indicate bins of the given size, and
             whose scores are the number of features in the bin.  The
             class of the feature will be set to "bin", and its name to
             "method:source".  This is a handy way of generating histograms
             of feature density.

The third argument, the $callback, is a code reference to which
retrieved features are passed.  It is described in more detail below.

This routine is responsible for getting arrays of GFF data out of the
database and passing them to the callback subroutine.  The callback
does the work of constructing a Bio::DB::GFF::Feature object out of
that data.  The callback expects a list of 13 fields:

  $refseq      The reference sequence
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
  $feature_id  A unique feature ID (may be undef)

These fields are in the same order as the raw GFF file, with the
exception that the group column has been parsed into group class and
group name fields.

The feature ID, if provided, is a unique identifier of the feature
line.  The module does not depend on this ID in any way, but it is
available via Bio::DB::GFF-E<gt>id() if wanted.  In the dbi::mysql and
dbi::mysqlopt adaptor, the ID is a unique row ID.  In the acedb
adaptor it is not used.

=cut

=head2 feature_summary(), coverage_array()

The DBI adaptors provide methods for rapidly fetching coverage
statistics across a region of interest. Please see
L<Bio::DB::GFF::Adaptor::dbi> for more information about these
methods.

=cut

sub get_features{
  my $self = shift;
  my ($search,$options,$callback) = @_;
  $self->throw("get_features() must be implemented by an adaptor");
}


=head2 _feature_by_name

 Title   : _feature_by_name
 Usage   : $db->_feature_by_name($class,$name,$location,$callback)
 Function: get a list of features by name and class
 Returns : count of number of features retrieved
 Args    : name of feature, class of feature, and a callback
 Status  : abstract

This method is used internally.  The callback arguments are the same
as those used by make_feature().  This method must be overidden by
subclasses.

=cut

sub _feature_by_name {
  my $self = shift;
  my ($class,$name,$location,$callback) = @_;
  $self->throw("_feature_by_name() must be implemented by an adaptor");
}

sub _feature_by_attribute {
  my $self = shift;
  my ($attributes,$callback) = @_;
  $self->throw("_feature_by_name() must be implemented by an adaptor");
}

=head2 _feature_by_id

 Title   : _feature_by_id
 Usage   : $db->_feature_by_id($ids,$type,$callback)
 Function: get a feature based
 Returns : count of number of features retrieved
 Args    : arrayref to feature IDs to fetch
 Status  : abstract

This method is used internally to fetch features either by their ID or
their group ID.  $ids is a arrayref containing a list of IDs, $type is
one of "feature" or "group", and $callback is a callback.  The
callback arguments are the same as those used by make_feature().  This
method must be overidden by subclasses.

=cut

sub _feature_by_id {
  my $self = shift;
  my ($ids,$type,$callback) = @_;
  $self->throw("_feature_by_id() must be implemented by an adaptor");
}

=head2 overlapping_features

 Title   : overlapping_features
 Usage   : $db->overlapping_features(@args)
 Function: get features that overlap the indicated range
 Returns : a list of Bio::DB::GFF::Feature objects
 Args    : see below
 Status  : public

This method is invoked by Bio::DB::GFF::Segment-E<gt>features() to find
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
  -rare      Turn on an optimization suitable for a relatively rare feature type,
             where it will be faster to filter by feature type first
             and then by position, rather than vice versa.
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
  $self->features_in_range(-range_type=>'overlaps',@_);
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
  $self->features_in_range(-range_type=>'contains',@_);
}

=head2 contained_in

 Title   : contained_in
 Usage   : @features = $s->contained_in(@args)
 Function: get features that contain this segment
 Returns : a list of Bio::DB::GFF::Feature objects
 Args    : see features()
 Status  : Public

This is identical in behavior to features() except that it returns
only those features that completely contain the segment.

=cut

sub contained_in {
  my $self = shift;
  $self->features_in_range(-range_type=>'contained_in',@_);
}

=head2 get_abscoords

 Title   : get_abscoords
 Usage   : $db->get_abscoords($name,$class,$refseq)
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

If $refseq is provided, the function searches only within the
specified reference sequence.

=cut

sub get_abscoords {
  my $self = shift;
  my ($name,$class,$refseq) = @_;
  $self->throw("get_abscoords() must be implemented by an adaptor");
}

=head2 get_types

 Title   : get_types
 Usage   : $db->get_types($absref,$class,$start,$stop,$count)
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
  my ($refseq,$class,$start,$stop,$count,$types) = @_;
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
  $refname               The name of the reference sequence for this feature
  $refclass              The class of the reference sequence for this feature
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
  my ($parent,$group_hash,          # these arguments provided by generic mechanisms
      $srcseq,                      # the rest is provided by adaptor
      $start,$stop,
      $source,$method,
      $score,$strand,$phase,
      $group_class,$group_name,
      $tstart,$tstop,
      $db_id,$group_id) = @_;

  return unless $srcseq;            # return undef if called with no arguments.  This behavior is used for
                                    # on-the-fly aggregation.

  my $group;  # undefined
  if (defined $group_class && defined $group_name) {
    $tstart ||= '';
    $tstop  ||= '';
    if ($group_hash) {
      $group = $group_hash->{$group_class,$group_name,$tstart,$tstop}
	||= $self->make_object($group_class,$group_name,$tstart,$tstop);
    } else {
      $group = $self->make_object($group_class,$group_name,$tstart,$tstop);
    }
  }

# fix for some broken GFF files
# unfortunately - has undesired side effects
#  if (defined $tstart && defined $tstop && !defined $strand) {
#    $strand = $tstart <= $tstop ? '+' : '-';
#  }

  if (ref $parent) { # note that the src sequence is ignored
    return Bio::DB::GFF::Feature->new_from_parent($parent,$start,$stop,
						  $method,$source,
						  $score,$strand,$phase,
						  $group,$db_id,$group_id,
						  $tstart,$tstop);
  } else {
    return Bio::DB::GFF::Feature->new($self,$srcseq,
				      $start,$stop,
				      $method,$source,
				      $score,$strand,$phase,
				      $group,$db_id,$group_id,
				      $tstart,$tstop);
  }
}

sub make_aggregated_feature {
  my $self                 = shift;
  my ($accumulated_features,$parent,$aggregators) = splice(@_,0,3);
  my $feature = $self->make_feature($parent,undef,@_);
  return [$feature] if $feature && !$feature->group;

  # if we have accumulated features and either: 
  # (1) make_feature() returned undef, indicated very end or
  # (2) the current group is different from the previous one

  local $^W = 0;  # irritating uninitialized value warning in next statement
  if (@$accumulated_features &&
      (!defined($feature) || ($accumulated_features->[-1]->group ne $feature->group))) {
    foreach my $a (@$aggregators) {  # last aggregator gets first shot
      $a->aggregate($accumulated_features,$self) or next;
    }
    my @result = @$accumulated_features;
    @$accumulated_features = $feature ? ($feature) : ();
    return unless @result;
    return \@result ;
  }
  push @$accumulated_features,$feature;
  return;
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
    $method = $method ? "\\Q$method\\E"  : ".*";
    $source = $source ? ":\\Q$source\\E" : "(?::.+)?";
    push @expr,"${method}${source}";
  }
  my $expr = join '|',@expr;
  return $self->{match_subs}{$expr} if $self->{match_subs}{$expr};

  my $sub =<<END;
sub {
  my \$feature = shift or return;
  return \$feature->type =~ /^($expr)\$/i;
}
END
  warn "match sub: $sub\n" if $self->debug;
  undef $@;
  my $compiled_sub = eval $sub;
  $self->throw($@) if $@;
  return $self->{match_subs}{$expr} = $compiled_sub;
}

=head2 make_object

 Title   : make_object
 Usage   : $db->make_object($class,$name,$start,$stop)
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
  my ($class,$name,$start,$stop) = @_;
  return Bio::DB::GFF::Homol->new($self,$class,$name,$start,$stop)
    if defined $start and length $start;
  return Bio::DB::GFF::Featname->new($class,$name);
}


=head2 do_attributes

 Title   : do_attributes
 Usage   : $db->do_attributes($id [,$tag]);
 Function: internal method to retrieve attributes given an id and tag
 Returns : a list of Bio::DB::GFF::Feature objects
 Args    : a feature id and a attribute tag (optional)
 Status  : protected

This method is overridden by subclasses in order to return a list of
attributes.  If called with a tag, returns the value of attributes of
that tag type.  If called without a tag, returns a flattened array of
(tag=E<gt>value) pairs.  A particular tag can be present multiple times.

=cut

sub do_attributes {
  my $self = shift;
  my ($id,$tag) = @_;
  return ();
}

=head2 clone

The clone() method should be used when you want to pass the
Bio::DB::GFF object to a child process across a fork(). The child must
call clone() before making any queries.

The default behavior is to do nothing, but adaptors that use the DBI
interface may need to implement this in order to avoid database handle
errors. See the dbi adaptor for an example.

=cut

sub clone { }


=head1 Internal Methods

The following methods are internal to Bio::DB::GFF and are not
guaranteed to remain the same.

=head2 _features

 Title   : _features
 Usage   : $db->_features($search,$options,$parent)
 Function: internal method
 Returns : a list of Bio::DB::GFF::Feature objects
 Args    : see below
 Status  : internal

This is an internal method that is called by overlapping_features(),
contained_features() and features() to create features based on a
parent segment's coordinate system.  It takes three arguments, a
search options hashref, an options hashref, and a parent segment.

The search hashref contains the following keys:

  rangetype     One of "overlaps", "contains" or "contained_in".  Indicates
                the type of range query requested.
  refseq        reference sequence ID
  refclass      reference sequence class
  start	        start of range
  stop		stop of range
  types	        arrayref containing list of types in "method:source" form

The options hashref contains zero or more of the following keys:

  sparse	turn on optimizations for a rare feature
  automerge	if true, invoke aggregators to merge features
  iterator	if true, return an iterator

The $parent argument is a scalar object containing a
Bio::DB::GFF::RelSegment object or descendent.

=cut

#'

sub _features {
  my $self = shift;
  my ($search,$options,$parent) = @_;
  (@{$search}{qw(start stop)}) = (@{$search}{qw(stop start)})
    if defined($search->{start}) && $search->{start} > $search->{stop};
  $search->{refseq} = $search->{seq_id} if exists $search->{seq_id};

  my $types = $self->parse_types($search->{types});  # parse out list of types
  my @aggregated_types = @$types;         # keep a copy

  # allow the aggregators to operate on the original
  my @aggregators;
  if ($options->{automerge}) {
    for my $a ($self->aggregators) {
      $a = $a->clone if $options->{iterator};
      unshift @aggregators,$a
	if $a->disaggregate(\@aggregated_types,$self);
    }
  }

  if ($options->{iterator}) {
    my @accumulated_features;
    my $callback = $options->{automerge} ? sub { $self->make_aggregated_feature(\@accumulated_features,$parent,\@aggregators,@_) }
                                         : sub { [$self->make_feature($parent,undef,@_)] };
    return $self->get_features_iterator({ %$search, 
					  types => \@aggregated_types  },
					{ %$options,
					 sort_by_group => $options->{automerge}  },
					$callback
				       );
  }

  my %groups;         # cache the groups we create to avoid consuming too much unecessary memory
  my $features = [];

  my $callback = sub { push @$features,$self->make_feature($parent,\%groups,@_) };
  $self->get_features({ %$search,
			types  => \@aggregated_types },
		        $options,
		        $callback);

  if ($options->{automerge}) {
    warn "aggregating...\n" if $self->debug;
    foreach my $a (@aggregators) {  # last aggregator gets first shot
      warn "Aggregator $a:\n" if $self->debug;
      $a->aggregate($features,$self);
    }
  }

  @$features;
}

=head2 get_features_iterator

 Title   : get_features_iterator
 Usage   : $db->get_features_iterator($search,$options,$callback)
 Function: get an iterator on a features query
 Returns : a Bio::SeqIO object
 Args    : as per get_features()
 Status  : Public

This method takes the same arguments as get_features(), but returns an
iterator that can be used to fetch features sequentially, as per
Bio::SeqIO.

Internally, this method is simply a front end to range_query().
The latter method constructs and executes the query, returning a
statement handle. This routine passes the statement handle to the
constructor for the iterator, along with the callback.

=cut

sub get_features_iterator {
  my $self = shift;
  my ($search,$options,$callback) = @_;
  $self->throw('feature iteration is not implemented in this adaptor');
}

=head2 split_group

 Title   : split_group
 Usage   : $db->split_group($group_field,$gff3_flag)
 Function: parse GFF group field
 Returns : ($gclass,$gname,$tstart,$tstop,$attributes)
 Args    : the gff group column and a flag indicating gff3 compatibility
 Status  : internal

This is a method that is called by load_gff_line to parse out the
contents of one or more group fields.  It returns the class of the
group, its name, the start and stop of the target, if any, and an
array reference containing any attributes that were stuck into the
group field, in [attribute_name,attribute_value] format.

=cut

sub split_group {
  my $self = shift;
  my ($group,$gff3) = @_;
  if ($gff3) {
    my @groups = split /[;&]/,$group;  # so easy!
    return $self->_split_gff3_group(@groups);
  } else {
    # handle group parsing
    # protect embedded semicolons in the group; there must be faster/more elegant way
    # to do this.
    $group =~ s/\\;/$;/g;
    while ($group =~ s/( \"[^\"]*);([^\"]*\")/$1$;$2/) { 1 }
    my @groups = split(/\s*;\s*/,$group);
    foreach (@groups) { s/$;/;/g }
    return $self->_split_gff2_group(@groups);
  }
}

=head2 _split_gff2_group

This is an internal method called by split_group().

=cut

# this has gotten quite nasty due to transition from GFF2 to GFF2.5
# (artemis) to GFF3.

sub _split_gff2_group {
  my $self = shift;
  my @groups = @_;
  my $target_found;

  my ($gclass,$gname,$tstart,$tstop,@attributes,@notes);

  for (@groups) {

    my ($tag,$value) = /^(\S+)(?:\s+(.+))?/;
    $value = '' unless defined $value;
    if ($value =~ /^\"(.+)\"$/) {  #remove quotes
      $value = $1;
    }
    $value =~ s/\\t/\t/g;
    $value =~ s/\\r/\r/g;
    $value =~ s/\s+$//;

    # Any additional groups become part of the attributes hash
    # For historical reasons, the tag "Note" is treated as an
    # attribute, even if it is the only group.
    $tag ||= '';
    if ($tag eq 'tstart' && $target_found) {
      $tstart = $value;
    }

    elsif ($tag eq 'tend' && $target_found) {
      $tstop = $value;
    }

    elsif (ucfirst $tag eq 'Note') {
      push @notes, [$tag => $value];
    }

    elsif ($tag eq 'Target' && /([^:\"\s]+):([^\"\s]+)/) { # major disagreement in implementors of GFF2 here
      $target_found++;
      ($gclass,$gname) = ($1,$2);
      ($tstart,$tstop) = / (\d+) (\d+)/;
    }

    elsif (!defined($value)) {
      push @notes, [Note => $tag];  # e.g. "Confirmed_by_EST"
    }

    else {
      push @attributes, [$tag => $value];
    }
  }

  # group assignment
  if (@attributes && !($gclass && $gname) ) {

    my $preferred = ref($self) ? $self->_preferred_groups_hash : {};

    for my $pair (@attributes) {
      my ($c,$n) = @$pair;
      ($gclass,$gname) = ($c,$n) 
	if !$gclass # pick up first one
	  ||
	    ($preferred->{lc $gclass}||0) < ($preferred->{lc $c}||0); # pick up higher priority one
    }

    @attributes = grep {$gclass ne $_->[0]} @attributes;
  }

  push @attributes, @notes;

  return ($gclass,$gname,$tstart,$tstop,\@attributes);
}


=head2 gff3_name_munging

 Title   : gff3_name_munging
 Usage   : $db->gff3_name_munging($boolean)
 Function: get/set gff3_name_munging flag
 Returns : $current value of flag
 Args    : new value of flag (optional)
 Status  : utility

If this is set to true (default false), then features identified in
gff3 files with an ID in the format foo:bar will be parsed so that
"foo" is the class and "bar" is the name.  This is mostly for backward
compatibility with GFF2.

=cut

sub gff3_name_munging {
  my $self = shift;
  my $d = $self->{gff3_name_munging};
  $self->{gff3_name_munging} = shift if @_;
  $d;
}

=head2 _split_gff3_group

This is called internally from split_group().

=cut

sub _split_gff3_group {
  my $self   = shift;
  my @groups = @_;
  my $dc     = $self->default_class;
  my (%id,@attributes);

  for my $group (@groups) {
    my ($tag,$value) = split /=/,$group;
    $tag             = unescape($tag);
    my @values       = map {unescape($_)} split /,/,$value;

    # GFF2 traditionally did not distinguish between a feature's name
    # and the group it belonged to.  This code is a transition between
    # gff2 and the new parent/ID dichotomy in gff3.
    if ($tag eq 'Parent') {
      my (@names,@classes);
      for (@values) {
	my ($name,$class) = $self->_gff3_name_munging($_,$dc);
	push @names,$name;
	push @classes,$class;
      }
      $id{$tag} = @names > 1 ? [\@names,\@classes] : [$names[0],$classes[0]];
    }
    elsif ($tag eq 'ID' || $tag eq 'Name') {
      $id{$tag} = [$self->_gff3_name_munging(shift(@values),$dc)];
    }
    elsif ($tag eq 'Target') {
      my ($gname,$tstart,$tstop) = split /\s+/,shift @values;
      $id{$tag} = [$self->_gff3_name_munging($gname,$dc),$tstart,$tstop];
    }
    elsif ($tag =~ /synonym/i) {
      $tag = 'Alias';
    }
    push @attributes,[$tag=>$_] foreach @values;
  }

  my $priorities = $self->_preferred_groups_hash(1);
  my ($gclass,$gname,$tstart,$tstop);
  for my $preferred (sort {$priorities->{lc $b}<=>$priorities->{lc $a}}
		     keys %id) {
    unless (defined $gname) {
      ($gname,$gclass,$tstart,$tstop) = @{$id{$preferred}};
    }
  }

  # set null gclass to empty string to preserve compatibility with
  # programs that expect a defined gclass if no gname
  $gclass ||= '' if defined $gname;

  return ($gclass,$gname,$tstart,$tstop,\@attributes);
}

# accomodation for wormbase style of class:name naming
sub _gff3_name_munging {
  my $self = shift;
  my ($name,$default_class) = @_;
  return ($name,$default_class) unless $self->gff3_name_munging;

  if ($name =~ /^(\w+):(.+)/) {
    return ($2,$1);
  } else {
    return ($name,$default_class);
  }
}

=head2 _delete_features(), _delete_groups(),_delete(),_delete_fattribute_to_features()

 Title   : _delete_features(), _delete_groups(),_delete(),_delete_fattribute_to_features()
 Usage   : $count = $db->_delete_features(@feature_ids)
           $count = $db->_delete_groups(@group_ids)
           $count = $db->_delete(\%delete_spec)
           $count = $db->_delete_fattribute_to_features(@feature_ids)
 Function: low-level feature/group deleter
 Returns : count of groups removed
 Args    : list of feature or group ids removed
 Status  : for implementation by subclasses

These methods need to be implemented in adaptors.  For _delete_features,
_delete_groups and _delete_fattribute_to_features, the arguments are a list of
feature or group IDs to remove.  For _delete(), the argument is a hashref with
the three keys 'segments', 'types' and 'force'.  The first contains an arrayref
of Bio::DB::GFF::RelSegment objects to delete (all FEATURES within the segment
are deleted).  The second contains an arrayref of [method,source] feature types
to delete.  The two are ANDed together.  If 'force' has a true value, this
forces the operation to continue even if it would delete all features.

=cut

sub _delete_features {
  my $self = shift;
  my @feature_ids = @_;
  $self->throw('_delete_features is not implemented in this adaptor');
}

sub _delete_groups {
  my $self = shift;
  my @group_ids = @_;
  $self->throw('_delete_groups is not implemented in this adaptor');
}

sub _delete {
  my $self = shift;
  my $delete_options = shift;
  $self->throw('_delete is not implemented in this adaptor');
}

sub _delete_fattribute_to_features {
  my $self = shift;
  my @feature_ids = @_;
  $self->throw('_delete_fattribute_to_features is not implemented in this adaptor');
}


sub unescape {
  my $v = shift;
  $v =~ tr/+/ /;
  $v =~ s/%([0-9a-fA-F]{2})/chr hex($1)/ge;
  return $v;
}

sub print_gff3_warning {
  my $self = shift;
  print STDERR <<END

You are loading a Bio::DB::GFF database with GFF3 formatted data.
While this will likely work fine, the Bio::DB::GFF schema does not
always faithfully capture the complexity represented in GFF3 files.
Unless you have a specific reason for using Bio::DB::GFF, we suggest
that you use a Bio::DB::SeqFeature::Store database and its corresponding
loader, bp_seqfeature_load.pl.

END
;

  return;
}


package Bio::DB::GFF::ID_Iterator;
use strict;

use base qw(Bio::Root::Root);

sub new {
  my $class            = shift;
  my ($db,$ids,$type)  = @_;
  return bless {ids=>$ids,db=>$db,type=>$type},$class;
}

sub next_seq {
  my $self = shift;
  my $next = shift @{$self->{ids}};
  return unless $next;
  my $name = ref($next) eq 'ARRAY' ? Bio::DB::GFF::Featname->new(@$next) : $next;
  my $segment = $self->{type} eq 'name'      ? $self->{db}->segment($name)
                : $self->{type} eq 'feature' ? $self->{db}->fetch_feature_by_id($name)
                : $self->{type} eq 'group'   ? $self->{db}->fetch_feature_by_gid($name)
                : $self->throw("Bio::DB::GFF::ID_Iterator called to fetch an unknown type of identifier");
  $self->throw("id does not exist") unless $segment;
  return $segment;
}

package Bio::DB::GFF::FeatureIterator;

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


1;

__END__

=head1 BUGS

Features can only belong to a single group at a time.  This must be
addressed soon.

Start coordinate can be greater than stop coordinate for relative
addressing.  This breaks strict BioPerl compatibility and must be
fixed.

=head1 SEE ALSO

L<Bio::DB::GFF::RelSegment>,
L<Bio::DB::GFF::Aggregator>,
L<Bio::DB::GFF::Feature>,
L<Bio::DB::GFF::Adaptor::dbi::mysqlopt>,
L<Bio::DB::GFF::Adaptor::dbi::oracle>,
L<Bio::DB::GFF::Adaptor::memory>
L<Bio::DB::GFF::Adaptor::berkeleydb>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut
