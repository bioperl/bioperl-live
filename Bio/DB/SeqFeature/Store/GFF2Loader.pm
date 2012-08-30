package Bio::DB::SeqFeature::Store::GFF2Loader;

# $Id: GFF2Loader.pm 11755 2007-11-08 02:19:29Z cjfields $

=head1 NAME

Bio::DB::SeqFeature::Store::GFF2Loader -- GFF2 file loader for Bio::DB::SeqFeature::Store

=head1 SYNOPSIS

  use Bio::DB::SeqFeature::Store;
  use Bio::DB::SeqFeature::Store::GFF2Loader;

  # Open the sequence database
  my $db      = Bio::DB::SeqFeature::Store->new( -adaptor => 'DBI::mysql',
                                                 -dsn     => 'dbi:mysql:test',
                                                 -write   => 1 );

  my $loader = Bio::DB::SeqFeature::Store::GFF2Loader->new(-store    => $db,
							   -verbose  => 1,
							   -fast     => 1);

  $loader->load('./my_genome.gff');


=head1 DESCRIPTION

The Bio::DB::SeqFeature::Store::GFF2Loader object parsers GFF2-format
sequence annotation files and loads Bio::DB::SeqFeature::Store
databases. For certain combinations of SeqFeature classes and
SeqFeature::Store databases it features a "fast load" mode which will
greatly accelerate the loading of GFF2 databases by a factor of 5-10.

The GFF2 file format has been extended very slightly to accommodate
Bio::DB::SeqFeature::Store. First, the loader recognizes is a new
directive:

  # #index-subfeatures [0|1]

Note that you can place a space between the two #'s in order to
prevent GFF2 validators from complaining.

If this is true, then subfeatures are indexed (the default) so that
they can be retrieved with a query. See L<Bio::DB::SeqFeature::Store>
for an explanation of this. If false, then subfeatures can only be
accessed through their parent feature. The default is to index all
subfeatures.

Second, the loader recognizes a new attribute tag called index, which
if present, controls indexing of the current feature. Example:

 ctg123	. TF_binding_site 1000 1012 . + . ID=tfbs00001;index=1

You can use this to turn indexing on and off, overriding the default
for a particular feature.

=cut


# load utility - incrementally load the store based on GFF2 file
#
# two modes:
#   slow mode -- features can occur in any order in the GFF2 file
#   fast mode -- all features with same ID must be contiguous in GFF2 file

use strict;
use Carp 'croak';
use Bio::DB::GFF::Util::Rearrange;
use Text::ParseWords 'quotewords';
use base 'Bio::DB::SeqFeature::Store::GFF3Loader';

my %Special_attributes =(
			 Gap    => 1, Target => 1,
			 Parent => 1, Name   => 1,
			 Alias  => 1, ID     => 1,
			 index  => 1, Index  => 1,
			);

=head2 new

 Title   : new
 Usage   : $loader = Bio::DB::SeqFeature::Store::GFF2Loader->new(@options)
 Function: create a new parser
 Returns : a Bio::DB::SeqFeature::Store::GFF2Loader gff2 parser and loader
 Args    : several - see below
 Status  : public

This method creates a new GFF2 loader and establishes its connection
with a Bio::DB::SeqFeature::Store database. Arguments are -name=E<gt>$value
pairs as described in this table:

 Name               Value
 ----               -----

 -store             A writable Bio::DB::SeqFeature::Store database handle.

 -seqfeature_class  The name of the type of Bio::SeqFeatureI object to create
                      and store in the database (Bio::DB::SeqFeature by default)

 -sf_class          A shorter alias for -seqfeature_class

 -verbose           Send progress information to standard error.

 -fast              If true, activate fast loading (see below)

 -chunk_size        Set the storage chunk size for nucleotide/protein sequences
                       (default 2000 bytes)

 -tmp               Indicate a temporary directory to use when loading non-normalized
                       features.

When you call new(), a connection to a Bio::DB::SeqFeature::Store
database should already have been established and the database
initialized (if appropriate).

Some combinations of Bio::SeqFeatures and Bio::DB::SeqFeature::Store
databases support a fast loading mode. Currently the only reliable
implementation of fast loading is the combination of DBI::mysql with
Bio::DB::SeqFeature. The other important restriction on fast loading
is the requirement that a feature that contains subfeatures must occur
in the GFF2 file before any of its subfeatures. Otherwise the
subfeatures that occurred before the parent feature will not be
attached to the parent correctly. This restriction does not apply to
normal (slow) loading.

If you use an unnormalized feature class, such as
Bio::SeqFeature::Generic, then the loader needs to create a temporary
database in which to cache features until all their parts and subparts
have been seen. This temporary databases uses the "berkeleydb" adaptor. The
-tmp option specifies the directory in which that database will be
created. If not present, it defaults to the system default tmp
directory specified by File::Spec-E<gt>tmpdir().

The -chunk_size option allows you to tune the representation of
DNA/Protein sequence in the Store database. By default, sequences are
split into 2000 base/residue chunks and then reassembled as
needed. This avoids the problem of pulling a whole chromosome into
memory in order to fetch a short subsequence from somewhere in the
middle. Depending on your usage patterns, you may wish to tune this
parameter using a chunk size that is larger or smaller than the
default.

=cut

# sub new { } inherited

=head2 load

 Title   : load
 Usage   : $count = $loader->load(@ARGV)
 Function: load the indicated files or filehandles
 Returns : number of feature lines loaded
 Args    : list of files or filehandles
 Status  : public

Once the loader is created, invoke its load() method with a list of
GFF2 or FASTA file paths or previously-opened filehandles in order to
load them into the database. Compressed files ending with .gz, .Z and
.bz2 are automatically recognized and uncompressed on the fly. Paths
beginning with http: or ftp: are treated as URLs and opened using the
LWP GET program (which must be on your path).

FASTA files are recognized by their initial "E<gt>" character. Do not feed
the loader a file that is neither GFF2 nor FASTA; I don't know what
will happen, but it will probably not be what you expect.

=cut

# sub load { } inherited

=head2 accessors

The following read-only accessors return values passed or created during new():

 store()          the long-term Bio::DB::SeqFeature::Store object

 tmp_store()      the temporary Bio::DB::SeqFeature::Store object used
                    during loading

 sfclass()        the Bio::SeqFeatureI class

 fast()           whether fast loading is active

 seq_chunk_size() the sequence chunk size

 verbose()        verbose progress messages

=cut

# sub store           inherited
# sub tmp_store       inherited
# sub sfclass         inherited
# sub fast            inherited
# sub seq_chunk_size  inherited
# sub verbose         inherited

=head2 Internal Methods

The following methods are used internally and may be overidden by
subclasses.

=over 4

=item default_seqfeature_class

  $class = $loader->default_seqfeature_class

Return the default SeqFeatureI class (Bio::DB::SeqFeature).

=cut

# sub default_seqfeature_class { } inherited

=item subfeatures_normalized

  $flag = $loader->subfeatures_normalized([$new_flag])

Get or set a flag that indicates that the subfeatures are
normalized. This is deduced from the SeqFeature class information.

=cut

# sub subfeatures_normalized { } inherited

=item subfeatures_in_table

  $flag = $loader->subfeatures_in_table([$new_flag])

Get or set a flag that indicates that feature/subfeature relationships
are stored in a table. This is deduced from the SeqFeature class and
Store information.

=cut

# sub subfeatures_in_table { } inherited

=item load_fh

  $count = $loader->load_fh($filehandle)

Load the GFF2 data at the other end of the filehandle and return true
if successful. Internally, load_fh() invokes:

  start_load();
  do_load($filehandle);
  finish_load();

=cut

# sub load_fh { } inherited

=item start_load, finish_load

These methods are called at the start and end of a filehandle load.

=cut

# sub create_load_data { } #inherited

# sub finish_load { } #inherite

=item do_load

  $count = $loader->do_load($fh)

This is called by load_fh() to load the GFF2 file's filehandle and
return the number of lines loaded.

=cut

# sub do_load { } inherited

=item load_line

    $loader->load_line($data);

Load a line of a GFF2 file. You must bracket this with calls to
start_load() and finish_load()!

    $loader->start_load();
    $loader->load_line($_) while <FH>;
    $loader->finish_load();

=cut

# sub load_line { } # inherited

=item handle_meta

  $loader->handle_meta($meta_directive)

This method is called to handle meta-directives such as
##sequence-region. The method will receive the directive with the
initial ## stripped off.

=cut

# sub handle_meta {} # inherited

=item handle_feature

  $loader->handle_feature($gff2_line)

This method is called to process a single GFF2 line. It manipulates
information stored a data structure called $self-E<gt>{load_data}.

=cut

# sub handle_feature { } # inherited

=item store_current_feature

  $loader->store_current_feature()

This method is called to store the currently active feature in the
database. It uses a data structure stored in $self-E<gt>{load_data}.

=cut

# sub store_current_feature { } inherited

=item build_object_tree

 $loader->build_object_tree()

This method gathers together features and subfeatures and builds the graph that connects them.

=cut

# sub build_object_tree { } # inherited

=item build_object_tree_in_tables

 $loader->build_object_tree_in_tables()

This method gathers together features and subfeatures and builds the
graph that connects them, assuming that parent/child relationships
will be stored in a database table.

=cut

# sub build_object_tree_in_tables { } # inherited

=item build_object_tree_in_features

 $loader->build_object_tree_in_features()

This method gathers together features and subfeatures and builds the
graph that connects them, assuming that parent/child relationships are
stored in the seqfeature objects themselves.

=cut

# sub build_object_tree_in_features { } # inherited

=item attach_children

 $loader->attach_children($store,$load_data,$load_id,$feature)

This recursively adds children to features and their subfeatures. It
is called when subfeatures are directly contained within other
features, rather than stored in a relational table.

=cut

# sub attach_children { } # inherited

=item fetch

 my $feature = $loader->fetch($load_id)

Given a load ID (from the ID= attribute) this method returns the
feature from the temporary database or the permanent one, depending on
where it is stored.

=cut

# sub fetch { } # inherited

=item add_segment

 $loader->add_segment($parent,$child)

This method is used to add a split location to the parent.

=cut

# sub add_segment { } # inherited

=item parse_attributes

 ($reserved,$unreserved) = $loader->parse_attributes($attribute_line)

This method parses the information contained in the $attribute_line
into two hashrefs, one containing the values of reserved attribute
tags (e.g. ID) and the other containing the values of unreserved ones.

=cut

sub parse_attributes { # overridden
  my $self = shift;
  my $att  = shift;
  my @groups = quotewords('\s*;\s*',0,$att);
  my (%reserved,%unreserved);

  my $found_name;

  for (@groups) {
      my ($tag,$value);
      if (/^(\S+)\s+(.+)/) { # Tag value pair
	  ($tag,$value) = ($1,$2);
      } else {
	  $tag    = 'Note';
	  $value  = $_;
      }
      if ($tag eq 'Target') {
	  my ($target,$start,$end) = split /\s+/,$value;
	  push @{$reserved{ID}},$target;
	  $found_name++;
	  if ($start <= $end) { $value .= ' +' }
	  else                { $value .= ' -' }
      }

      if (!$found_name++) {
	  push @{$reserved{Alias}},$value;
	  $value = "$tag:$value";
	  push @{$reserved{ID}},$value;
	  $tag   = 'Name';
      }

      if ($Special_attributes{$tag}) {  # reserved attribute
	  push @{$reserved{$tag}},$value;
      } else {
	  push @{$unreserved{$tag}},$value;
      }
  }
  return (\%reserved,\%unreserved);
}

=item start_or_finish_sequence

  $loader->start_or_finish_sequence('Chr9')

This method is called at the beginning and end of a fasta section.

=cut

# sub start_or_finish_sequence { } inherited

=item load_sequence

  $loader->load_sequence('gatttcccaaa')

This method is called to load some amount of sequence after
start_or_finish_sequence() is first called.

=cut

# sub load_sequence { } inherited

=item open_fh

 my $io_file = $loader->open_fh($filehandle_or_path)

This method opens up the indicated file or pipe, using some
intelligence to recognized compressed files and URLs and doing the
right thing.

=cut

# sub open_fh { } inherited

# sub msg { } inherited

=item time

 my $time = $loader->time

This method returns the current time in seconds, using Time::HiRes if available.

=cut

# sub time { } inherited

=item unescape

 my $unescaped = GFF2Loader::unescape($escaped)

This is an internal utility.  It is the same as CGI::Util::unescape,
but doesn't change pluses into spaces and ignores unicode escapes.

=cut

# sub unescape { } inherited

1;
__END__

=back

=head1 BUGS

This is an early version, so there are certainly some bugs. Please
use the BioPerl bug tracking system to report bugs.

=head1 SEE ALSO

L<bioperl>,
L<Bio::DB::SeqFeature::Store>,
L<Bio::DB::SeqFeature::Segment>,
L<Bio::DB::SeqFeature::NormalizedFeature>,
L<Bio::DB::SeqFeature::GFF3Loader>,
L<Bio::DB::SeqFeature::Store::DBI::mysql>,
L<Bio::DB::SeqFeature::Store::berkeleydb>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2006 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut


