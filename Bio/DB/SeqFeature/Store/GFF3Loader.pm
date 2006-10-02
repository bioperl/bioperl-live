package Bio::DB::SeqFeature::Store::GFF3Loader;

# $Id$

=head1 NAME

Bio::DB::SeqFeature::Store::GFF3Loader -- GFF3 file loader for Bio::DB::SeqFeature::Store

=head1 SYNOPSIS

  use Bio::DB::SeqFeature::Store;

  # Open the sequence database
  my $db      = Bio::DB::SeqFeature::Store->new( -adaptor => 'DBI::mysql',
                                                 -dsn     => 'dbi:mysql:test',
                                                 -write   => 1 );

  my $loader = Bio::DB::SeqFeature::Store::GFF3Loader->new(-store    => $db,
							   -verbose  => 1,
							   -fast     => 1);

  $loader->load('./my_genome.gff3');


=head1 DESCRIPTION

The Bio::DB::SeqFeature::Store::GFF3Loader object parsers GFF3-format
sequence annotation files and loads Bio::DB::SeqFeature::Store
databases. For certain combinations of SeqFeature classes and
SeqFeature::Store databases it features a "fast load" mode which will
greatly accelerate the loading of GFF3 databases by a factor of 5-10.

The GFF3 file format has been extended very slightly to accomodate
Bio::DB::SeqFeature::Store. First, the loader recognizes is a new
directive:

  # #index-subfeatures [0|1]

Note that you can place a space between the two #'s in order to
prevent GFF3 validators from complaining.

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


# load utility - incrementally load the store based on GFF3 file
#
# two modes:
#   slow mode -- features can occur in any order in the GFF3 file
#   fast mode -- all features with same ID must be contiguous in GFF3 file

use strict;
use Carp 'croak';
use IO::File;
use Bio::DB::GFF::Util::Rearrange;
use Bio::DB::SeqFeature::Store;
use File::Spec;
use base 'Bio::Root::Root';

use constant DEFAULT_SEQ_CHUNK_SIZE => 2000;

my %Special_attributes =(
			 Gap    => 1, Target => 1,
			 Parent => 1, Name   => 1,
			 Alias  => 1, ID     => 1,
			 index  => 1, Index  => 1,
			);
my %Strandedness = ( '+'  => 1,
		     '-'  => -1,
		     '.'  => 0,
		     ''   => 0,
		     0    => 0,
		     1    => 1,
		     -1   => -1,
		     +1   => 1,
		     undef => 0,
		   );

=head2 new

 Title   : new
 Usage   : $loader = Bio::DB::SeqFeature::Store::GFF3Loader->new(@options)
 Function: create a new parser
 Returns : a Bio::DB::SeqFeature::Store::GFF3Loader gff3 parser and loader
 Args    : several - see below
 Status  : public

This method creates a new GFF3 loader and establishes its connection
with a Bio::DB::SeqFeature::Store database. Arguments are -name=E<gt>$value
pairs as described in this table:

 Name               Value
 ----               -----

 -store             A writeable Bio::DB::SeqFeature::Store database handle.

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
in the GFF3 file before any of its subfeatures. Otherwise the
subfeatures that occurred before the parent feature will not be
attached to the parent correctly. This restriction does not apply to
normal (slow) loading.

If you use an unnormalized feature class, such as
Bio::SeqFeature::Generic, then the loader needs to create a temporary
database in which to cache features until all their parts and subparts
have been seen. This temporary databases uses the "bdb" adaptor. The
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

sub new {
  my $self = shift;
  my ($store,$seqfeature_class,$tmpdir,$verbose,$fast,$seq_chunk_size) = rearrange(['STORE',
										    ['SF_CLASS','SEQFEATURE_CLASS'],
										    ['TMP','TMPDIR'],
										    'VERBOSE',
										    'FAST',
										    'CHUNK_SIZE',
										   ],@_);

  $seqfeature_class ||= $self->default_seqfeature_class;
  eval "require $seqfeature_class" unless $seqfeature_class->can('new');
  $self->throw($@) if $@;

  my $normalized = $seqfeature_class->can('subfeatures_are_normalized')
    && $seqfeature_class->subfeatures_are_normalized;

  my $in_table = $seqfeature_class->can('subfeatures_are_stored_in_a_table')
    && $seqfeature_class->subfeatures_are_stored_in_a_table;

  if ($fast) {
    my $canfast = $normalized && $in_table;
    warn <<END unless $canfast;
Only features that support the Bio::DB::SeqFeature::NormalizedTableFeature interface
can be loaded using the -fast method. Reverting to slower feature-by-feature method.
END
    $fast &&= $canfast;
  }

  # try to bring in highres time() function
  eval "require Time::HiRes";

  $tmpdir ||= File::Spec->tmpdir();

  my $tmp_store = Bio::DB::SeqFeature::Store->new(-adaptor  => 'berkeleydb',
						  -temporary=> 1,
						  -dsn      => $tmpdir,
						  -cache    => 1,
						  -write    => 1)
    unless $normalized;

  return bless {
		store            => $store,
		tmp_store        => $tmp_store,
		seqfeature_class => $seqfeature_class,
		fast             => $fast,
		seq_chunk_size   => $seq_chunk_size || DEFAULT_SEQ_CHUNK_SIZE,
		verbose          => $verbose,
		load_data        => {},
		subfeatures_normalized => $normalized,
		subfeatures_in_table   => $in_table,
	       },ref($self) || $self;
}

=head2 load

 Title   : load
 Usage   : $count = $loader->load(@ARGV)
 Function: load the indicated files or filehandles
 Returns : number of feature lines loaded
 Args    : list of files or filehandles
 Status  : public

Once the loader is created, invoke its load() method with a list of
GFF3 or FASTA file paths or previously-opened filehandles in order to
load them into the database. Compressed files ending with .gz, .Z and
.bz2 are automatically recognized and uncompressed on the fly. Paths
beginning with http: or ftp: are treated as URLs and opened using the
LWP GET program (which must be on your path).

FASTA files are recognized by their initial "E<gt>" character. Do not feed
the loader a file that is neither GFF3 nor FASTA; I don't know what
will happen, but it will probably not be what you expect.

=cut

sub load {
  my $self       = shift;
  my $start      = $self->time();
  my $count = 0;

  for my $file_or_fh (@_) {
    $self->msg("loading $file_or_fh...\n");
    my $fh = $self->open_fh($file_or_fh) or $self->throw("Couldn't open $file_or_fh: $!");
    $count += $self->load_fh($fh);
    $self->msg(sprintf "load time: %5.2fs\n",$self->time()-$start);
  }
  $count;
}

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

sub store          { shift->{store}            }
sub tmp_store      { shift->{tmp_store}        }
sub sfclass        { shift->{seqfeature_class} }
sub fast           { shift->{fast}             }
sub seq_chunk_size { shift->{seq_chunk_size}             }
sub verbose        { shift->{verbose}          }

=head2 Internal Methods

The following methods are used internally and may be overidden by
subclasses.

=over 4

=item default_seqfeature_class

  $class = $loader->default_seqfeature_class

Return the default SeqFeatureI class (Bio::DB::SeqFeature).

=cut

sub default_seqfeature_class {
  my $self = shift;
  return 'Bio::DB::SeqFeature';
}

=item subfeatures_normalized

  $flag = $loader->subfeatures_normalized([$new_flag])

Get or set a flag that indicates that the subfeatures are
normalized. This is deduced from the SeqFeature class information.

=cut

sub subfeatures_normalized {
  my $self = shift;
  my $d    = $self->{subfeatures_normalized};
  $self->{subfeatures_normalized} = shift if @_;
  $d;
}

=item subfeatures_in_table

  $flag = $loader->subfeatures_in_table([$new_flag])

Get or set a flag that indicates that feature/subfeature relationships
are stored in a table. This is deduced from the SeqFeature class and
Store information.

=cut

sub subfeatures_in_table {
  my $self = shift;
  my $d    = $self->{subfeatures_in_table};
  $self->{subfeatures_in_table} = shift if @_;
  $d;
}

=item load_fh

  $count = $loader->load_fh($filehandle)

Load the GFF3 data at the other end of the filehandle and return true
if successful. Internally, load_fh() invokes:

  start_load();
  do_load($filehandle);
  finish_load();

=cut

sub load_fh {
  my $self = shift;
  my $fh   = shift;
  $self->start_load();
  my $count = $self->do_load($fh);
  $self->finish_load();
  $count;
}


=item start_load, finish_load

These methods are called at the start and end of a filehandle load.

=cut

sub start_load {
  my $self = shift;
  $self->{load_data}{Parent2Child}     = {};
  $self->{load_data}{Local2GlobalID}   = {};
  $self->{load_data}{TemporaryID}      = "GFFLoad0000000";
  $self->{load_data}{IndexSubfeatures} = 1;
  $self->{load_data}{CurrentFeature}   = undef;
  $self->{load_data}{CurrentID}        = undef;
  $self->store->start_bulk_update() if $self->fast;
}

sub finish_load {
  my $self  = shift;

  $self->msg("Building object tree...");
  my $start = $self->time();
  $self->build_object_tree;
  $self->msg(sprintf "%5.2fs\n",$self->time()-$start);

  if ($self->fast) {
    $self->msg("Loading bulk data into database...");
    $start = $self->time();
    $self->store->finish_bulk_update;
    $self->msg(sprintf "%5.2fs\n",$self->time()-$start);
  }
  eval {$self->store->commit};
  delete $self->{load_data};
}

=item do_load

  $count = $loader->do_load($fh)

This is called by load_fh() to load the GFF3 file's filehandle and
return the number of lines loaded.

=cut

sub do_load {
  my $self = shift;
  my $fh   = shift;

  my $start = $self->time();
  my $count = 0;
  my $mode  = 'gff';  # or 'fasta'

  while (<$fh>) {
    chomp;

    next unless /^\S/;     # blank line
    $mode = 'gff' if /\t/;  # if it has a tab in it, switch to gff mode

    if (/^\#\s?\#\s*(.+)/) {  ## meta instruction
      $mode = 'gff';
      $self->handle_meta($1);

    } elsif (/^\#/) {
      $mode = 'gff';  # just to be safe
      next;  # comment
    }

    elsif (/^>\s*(\S+)/) { # FASTA lines are coming
      $mode = 'fasta';
      $self->start_or_finish_sequence($1);
    }

    elsif ($mode eq 'fasta') {
      $self->load_sequence($_);
    }

    elsif ($mode eq 'gff') {
      $self->handle_feature($_);
      if (++$count % 1000 == 0) {
	my $now = $self->time();
	my $nl = -t STDOUT && !$ENV{EMACS} ? "\r" : "\n";
	$self->msg(sprintf("%d features loaded in %5.2fs...$nl",$count,$now - $start));
	$start = $now;
      }
    }

    else {
      $self->throw("I don't know what to do with this line:\n$_");
    }
  }
  $self->store_current_feature();      # during fast loading, we will have a feature left at the very end
  $self->start_or_finish_sequence();   # finish any half-loaded sequences
  $self->msg(' 'x80,"\n"); #clear screen
  $count;
}

=item handle_meta

  $loader->handle_meta($meta_directive)

This method is called to handle meta-directives such as
##sequence-region. The method will receive the directive with the
initial ## stripped off.

=cut

sub handle_meta {
  my $self = shift;
  my $instruction = shift;

  if ($instruction =~ /sequence-region\s+(.+)\s+(-?\d+)\s+(-?\d+)/i) {
    my $feature = $self->sfclass->new(-name        => $1,
				      -seq_id      => $1,
				      -start       => $2,
				      -end         => $3,
				      -primary_tag => 'region');
    $self->store->store($feature);
    return;
  }

  if ($instruction =~/index-subfeatures\s+(\S+)/i) {
    $self->{load_data}{IndexSubfeatures} = $1;
    $self->store->index_subfeatures($1);
    return;
  }
}

=item handle_feature

  $loader->handle_feature($gff3_line)

This method is called to process a single GFF3 line. It manipulates
information stored a data structure called $self-E<gt>{load_data}.

=cut

sub handle_feature {
  my $self     = shift;
  my $gff_line = shift;
  my $ld       = $self->{load_data};

  my @columns = map {$_ eq '.' ? undef : $_ } split /\t/,$gff_line;
  return unless @columns >= 8;
  my ($refname,$source,$method,$start,$end, $score,$strand,$phase,$attributes)      = @columns;
  $strand = $Strandedness{$strand||0};

  my ($reserved,$unreserved) = $self->parse_attributes($attributes);

  my $name        = ($reserved->{Name}   && $reserved->{Name}[0]);

  my $has_loadid  = defined $reserved->{ID}[0];

  my $feature_id  = $reserved->{ID}[0] || $ld->{TemporaryID}++;
  my @parent_ids  = @{$reserved->{Parent}} if $reserved->{Parent};

  my $index_it = $ld->{IndexSubfeatures};
  if (exists $reserved->{Index} || exists $reserved->{index}) {
    $index_it = $reserved->{Index}[0] || $reserved->{index}[0];
  }

  # Everything in the unreserved hash becomes an attribute, so we copy
  # some attributes over
  $unreserved->{Note}   = $reserved->{Note}   if exists $reserved->{Note};
  $unreserved->{Alias}  = $reserved->{Alias}  if exists $reserved->{Alias};
  $unreserved->{Target} = $reserved->{Target} if exists $reserved->{Target};
  $unreserved->{Gap}    = $reserved->{Gap}    if exists $reserved->{Gap};
  $unreserved->{load_id}= $reserved->{ID}     if exists $reserved->{ID};

  # TEMPORARY HACKS TO SIMPLIFY DEBUGGING
  push @{$unreserved->{Alias}},$feature_id  if $has_loadid;
  $unreserved->{parent_id} = \@parent_ids   if @parent_ids;

  # POSSIBLY A PERMANENT HACK -- TARGETS BECOME ALIASES
  # THIS IS TO ALLOW FOR TARGET-BASED LOOKUPS
  if (exists $reserved->{Target}) {
    my %aliases = map {$_=>1} @{$unreserved->{Alias}};
    for my $t (@{$reserved->{Target}}) {
      (my $tc = $t) =~ s/\s+.*$//;  # get rid of coordinates
      $name ||= $tc;
      push @{$unreserved->{Alias}},$tc unless $name eq $tc || $aliases{$tc};
    }
  }

  my @args = (-display_name => $name,
	      -seq_id       => $refname,
	      -start        => $start,
	      -end          => $end,
	      -strand       => $strand || 0,
	      -score        => $score,
	      -phase        => $phase,
	      -primary_tag  => $method || 'feature',
	      -source       => $source,
	      -tag          => $unreserved,
	      -attributes   => $unreserved,
	     );

  # Here's where we handle feature lines that have the same ID (multiple locations, not
  # parent/child relationships)

  my $old_feat;

  # Current feature is the same as the previous feature, which hasn't yet been loaded
  if (defined $ld->{CurrentID} && $ld->{CurrentID} eq $feature_id) {
    $old_feat = $ld->{CurrentFeature};
  }

  # Current feature is the same as a feature that was loaded earlier
  elsif (my $id = $self->{load_data}{Local2GlobalID}{$feature_id}) {
    $old_feat = $self->fetch($feature_id)
      or $self->warn(<<END);
ID=$feature_id has been used more than once, but it cannot be found in the database.
This can happen if you have specified fast loading, but features sharing the same ID
are not contiguous in the GFF file. This will be loaded as a separate feature.
Line $.: "$_"
END
  }

  # contiguous feature, so add a segment
  if (defined $old_feat) {
    $self->add_segment($old_feat,$self->sfclass->new(@args));
    return;
  }

  # we get here if this is a new feature
  # first of all, store the current feature if it is there
  $self->store_current_feature() if defined $ld->{CurrentID};

  # now create the new feature
  # (index top-level features only if policy asks us to)
  my $feature = $self->sfclass->new(@args);
  $feature->object_store($self->store) if $feature->can('object_store');  # for lazy table features
  $ld->{CurrentFeature} = $feature;
  $ld->{CurrentID}      = $feature_id;

  my $top_level = !@parent_ids;
  my $has_id    = defined $reserved->{ID}[0];
  $index_it   ||= $top_level;

  $ld->{IndexIt}{$feature_id}++    if $index_it;
  $ld->{TopLevel}{$feature_id}++   if !$self->{fast} && $top_level;  # need to track top level features

  # remember parentage
  for my $parent (@parent_ids) {
    push @{$ld->{Parent2Child}{$parent}},$feature_id;
  }

}

=item store_current_feature

  $loader->store_current_feature()

This method is called to store the currently active feature in the
database. It uses a data structure stored in $self-E<gt>{load_data}.

=cut


sub store_current_feature {
  my $self    = shift;

  my $ld   = $self->{load_data};
  defined $ld->{CurrentFeature} or return;
  my $f    = $ld->{CurrentFeature};

  my $normalized = $self->subfeatures_normalized;
  my $indexed    = $ld->{IndexIt}{$ld->{CurrentID}};

  # logic is as follows:
  # 1. If the feature is an indexed feature, then we store it into the main database
  #    so that it can be searched. It doesn't matter whether it is a top-level feature
  #    or a subfeature.
  # 2. If the feature class is normalized, but not indexed, then we store it into the
  #    main database using the "no_index" method. This will make it accessible to
  #    queries on the top level parent, but it won't come up by itself in range or
  #    attribute searches.
  # 3. Otherwise, this is an unindexed subfeature; we store it in the temporary database
  #    until the object build step, at which point it gets integrated into its object tree
  #    and copied into the main database.

  if ($indexed) {
    $self->store->store($f);
  }

  elsif ($normalized) {
    $self->store->store_noindex($f)
  }

  else {
    $self->tmp_store->store_noindex($f)
  }
	
  my $id        = $f->primary_id;    # assigned by store()
  $ld->{Local2GlobalID}{$ld->{CurrentID}} = $id;

  undef $ld->{IndexIt}{$ld->{CurrentID}} if $normalized;  # no need to remember this
  undef $ld->{CurrentID};
  undef $ld->{CurrentFeature};
}

=item build_object_tree

 $loader->build_object_tree()

This method gathers together features and subfeatures and builds the graph that connects them.

=cut

###
# put objects together
#
sub build_object_tree {
  my $self = shift;
  $self->subfeatures_in_table ? $self->build_object_tree_in_tables : $self->build_object_tree_in_features;
}

=item build_object_tree_in_tables

 $loader->build_object_tree_in_tables()

This method gathers together features and subfeatures and builds the
graph that connects them, assuming that parent/child relationships
will be stored in a database table.

=cut

sub build_object_tree_in_tables {
  my $self = shift;
  my $store = $self->store;
  my $ld    = $self->{load_data};

  while (my ($load_id,$children) = each %{$ld->{Parent2Child}}) {
    my $parent_id = $ld->{Local2GlobalID}{$load_id} or $self->throw("$load_id doesn't have a primary id");
    my @children  = map {$ld->{Local2GlobalID}{$_}} @$children;

    # this updates the table that keeps track of parent/child relationships,
    # but does not update the parent object -- so (start,end) had better be right!!!
    $store->add_SeqFeature($parent_id,@children);
  }

}

=item build_object_tree_in_features

 $loader->build_object_tree_in_features()

This method gathers together features and subfeatures and builds the
graph that connects them, assuming that parent/child relationships are
stored in the seqfeature objects themselves.

=cut

sub build_object_tree_in_features {
  my $self  = shift;
  my $store      = $self->store;
  my $tmp        = $self->tmp_store;
  my $ld         = $self->{load_data};
  my $normalized = $self->subfeatures_normalized;

  while (my ($load_id) = each %{$ld->{TopLevel}}) {
    my $feature  = $self->fetch($load_id)
      or $self->throw("$load_id (id=$ld->{Local2GlobalID}{$load_id}) should have a database entry, but doesn't");
    $self->attach_children($store,$ld,$load_id,$feature);
    $feature->primary_id(undef) unless $ld->{IndexIt}{$load_id};  # Indexed objects are updated, not created anew
    $store->store($feature);
  }

}

=item attach_children

 $loader->attach_children($store,$load_data,$load_id,$feature)

This recursively adds children to features and their subfeatures. It
is called when subfeatures are directly contained within other
features, rather than stored in a relational table.

=cut

sub attach_children {
  my $self = shift;
  my ($store,$ld,$load_id,$feature)  = @_;

  my $children   = $ld->{Parent2Child}{$load_id} or return;
  for my $child_id (@$children) {
    my $child = $self->fetch($child_id)
      or $self->throw("$child_id should have a database entry, but doesn't");
    $self->attach_children($store,$ld,$child_id,$child);   # recursive call
    $feature->add_SeqFeature($child);
  }
}

=item fetch

 my $feature = $loader->fetch($load_id)

Given a load ID (from the ID= attribute) this method returns the
feature from the temporary database or the permanent one, depending on
where it is stored.

=cut

sub fetch {
  my $self    = shift;
  my $load_id = shift;
  my $ld      = $self->{load_data};
  my $id      = $ld->{Local2GlobalID}{$load_id};

  return
    $self->subfeatures_normalized || $ld->{IndexIt}{$load_id}
      ? $self->store->fetch($id)
      : $self->tmp_store->fetch($id);
}

=item add_segment

 $loader->add_segment($parent,$child)

This method is used to add a split location to the parent.

=cut

sub add_segment {
  my $self = shift;
  my ($parent,$child) = @_;

  if ($parent->can('add_segment')) { # probably a lazy table feature
    my $segment_count =  $parent->can('denormalized_segment_count') ? $parent->denormalized_segment_count
                       : $parent->can('denormalized_segments ')     ? $parent->denormalized_segments
		       : $parent->can('segments')                   ? $parent->segments
		       : 0;
    unless ($segment_count) {  # convert into a segmented object
      my $segment;
      if ($parent->can('clone')) {
	$segment = $parent->clone;
      } else {
	my %clone   = %$parent;
	$segment = bless \%clone,ref $parent;
      }
      delete $segment->{segments};
      eval {$segment->object_store(undef) };
      $segment->primary_id(undef);

      # this updates the object and expands its start and end positions without writing
      # the segments into the database as individual objects
      $parent->add_segment($segment);
    }
    $parent->add_segment($child);
    1; # for debugging
  }

  # a conventional Bio::SeqFeature::Generic object - create a split location
  else {
    my $current_location = $parent->location;
    if ($current_location->can('add_sub_Location')) {
      $current_location->add_sub_Location($child->location);
    } else {
      eval "require Bio::Location::Split" unless Bio::Location::Split->can('add_sub_Location');
      my $new_location = Bio::Location::Split->new();
      $new_location->add_sub_Location($current_location);
      $new_location->add_sub_Location($child->location);
      $parent->location($new_location);
    }
  }
}

=item parse_attributes

 ($reserved,$unreserved) = $loader->parse_attributes($attribute_line)

This method parses the information contained in the $attribute_line
into two hashrefs, one containing the values of reserved attribute
tags (e.g. ID) and the other containing the values of unreserved ones.

=cut

sub parse_attributes {
  my $self = shift;
  my $att  = shift;
  my @pairs =  map {my ($name,$value) = split /=/; [unescape($name) => unescape($value)] } split /;/,$att;
  my (%reserved,%unreserved);
  foreach (@pairs) {
    my $tag    = $_->[0];
    my @values = split /,/,$_->[1];

    if ($Special_attributes{$tag}) {  # reserved attribute
      push @{$reserved{$tag}},@values;
    } else {
      push @{$unreserved{$tag}},@values
    }
  }
  return (\%reserved,\%unreserved);
}

=item start_or_finish_sequence

  $loader->start_or_finish_sequence('Chr9')

This method is called at the beginning and end of a fasta section.

=cut


# this gets called at the beginning and end of a fasta section
sub start_or_finish_sequence {
  my $self  = shift;
  my $seqid = shift;
  if (my $sl    = $self->{fasta_load}) {
    if (defined $sl->{seqid}) {
      $self->store->insert_sequence($sl->{seqid},$sl->{sequence},$sl->{offset});
      delete $self->{fasta_load};
    }
  }
  if (defined $seqid) {
    $self->{fasta_load} = {seqid  => $seqid,
			   offset => 0,
			   sequence => ''};
  }
}

=item load_sequence

  $loader->load_sequence('gatttcccaaa')

This method is called to load some amount of sequence after
start_or_finish_sequence() is first called.

=cut

sub load_sequence {
  my $self = shift;
  my $seq  = shift;
  my $sl   = $self->{fasta_load} or return;
  my $cs   = $self->seq_chunk_size;
  $sl->{sequence} .= $seq;
  while (length $sl->{sequence} >= $cs) {
    my $chunk = substr($sl->{sequence},0,$cs);
    $self->store->insert_sequence($sl->{seqid},$chunk,$sl->{offset});
    $sl->{offset} += length $chunk;
    substr($sl->{sequence},0,$cs) = '';
  }
}

=item open_fh

 my $io_file = $loader->open_fh($filehandle_or_path)

This method opens up the indicated file or pipe, using some intelligence to recognized compressed files and URLs and doing 
the right thing.

=cut


sub open_fh {
  my $self  = shift;
  my $thing = shift;

  no strict 'refs';

  return $thing                                  if defined fileno($thing);
  return IO::File->new("gunzip -c $thing |")     if $thing =~ /\.gz$/;
  return IO::File->new("uncompress -c $thing |") if $thing =~ /\.Z$/;
  return IO::File->new("bunzip2 -c $thing |")    if $thing =~ /\.bz2$/;
  return IO::File->new("GET $thing |")           if $thing =~ /^(http|ftp):/;
  return IO::File->new($thing);
}

sub msg {
  my $self = shift;
  my @msg  = @_;
  return unless $self->verbose;
  print STDERR @msg;
}

=item time

 my $time = $loader->time

This method returns the current time in seconds, using Time::HiRes if available.

=cut

sub time {
  return Time::HiRes::time() if Time::HiRes->can('time');
  return time();
}

=item escape

 my $unescaped = GFF3Loader::unescape($escaped)

This is an internal utility.  It is the same as CGI::Util::unescape,
but don't change pluses into spaces and ignores unicode escapes.

=cut

sub unescape {
  my $todecode = shift;
  $todecode =~ s/%([0-9a-fA-F]{2})/chr hex($1)/ge;
  return $todecode;
}

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
L<Bio::DB::SeqFeature::Store::bdb>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2006 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut


