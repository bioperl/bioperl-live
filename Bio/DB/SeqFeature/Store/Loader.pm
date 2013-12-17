package Bio::DB::SeqFeature::Store::Loader;


=head1 NAME

Bio::DB::SeqFeature::Store::Loader -- Loader

=head1 SYNOPSIS

 # non-instantiable base class

=head1 DESCRIPTION

This is the base class for Bio::DB::SeqFeature::Loader::GFF3Loader,
Bio::DB::SeqFeature::Loader::GFFLoader, and
Bio::DB::SeqFeature::FeatureFileLoader. Please see the manual pages
for these modules.

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
use File::Temp 'tempdir';
use base 'Bio::Root::Root';

use constant DEFAULT_SEQ_CHUNK_SIZE => 2000;

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

 -map_coords        A code ref that will transform a list of ($ref,[$start1,$end1]...)
                       coordinates into a list of ($newref,[$newstart1,$newend1]...)

 -index_subfeatures Indicate true if subfeatures should be indexed. Default is true.

 -summary_stats     Rebuild summary stats at the end of loading (not incremental,
                     so takes a long time)

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

sub new {
  my $self = shift;
  my ($store,$seqfeature_class,$tmpdir,$verbose,$fast,
      $seq_chunk_size,$coordinate_mapper,$index_subfeatures,$summary_stats) = 
      rearrange(['STORE',
		 ['SF_CLASS','SEQFEATURE_CLASS'],
		 ['TMP','TMPDIR'],
		 'VERBOSE',
		 'FAST',
		 'CHUNK_SIZE',
		 'MAP_COORDS',
		 'INDEX_SUBFEATURES',
		 'SUMMARY_STATS'
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

  $tmpdir      ||= File::Spec->tmpdir();

  my ($tmp_store,$temp_load);
  unless ($normalized) {

      # remember the temporary directory in order to delete it on exit
      $temp_load = tempdir(
	  'BioDBSeqFeature_XXXXXXX',
	  DIR=>$tmpdir,
	  CLEANUP=>1
	  );

      $tmp_store = Bio::DB::SeqFeature::Store->new(-adaptor  => 'berkeleydb',
						      -temporary=> 1,
						      -dsn      => $temp_load,
						      -cache    => 1,
						      -write    => 1)
	  unless $normalized;
  }

  $index_subfeatures = 1 unless defined $index_subfeatures;

  return bless {
		store            => $store,
		tmp_store        => $tmp_store,
		seqfeature_class => $seqfeature_class,
		fast             => $fast,
		seq_chunk_size   => $seq_chunk_size || DEFAULT_SEQ_CHUNK_SIZE,
		verbose          => $verbose,
		load_data        => {},
		tmpdir           => $tmpdir,
		temp_load        => $temp_load,
		subfeatures_normalized => $normalized,
		subfeatures_in_table   => $in_table,
		coordinate_mapper      => $coordinate_mapper,
		index_subfeatures      => $index_subfeatures,
		summary_stats          => $summary_stats,
	       },ref($self) || $self;
}

sub coordinate_mapper {
    my $self = shift;
    my $d    = $self->{coordinate_mapper};
    $self->{coordinate_mapper} = shift if @_;
    $d;
}

sub index_subfeatures {
    my $self = shift;
    my $d    = $self->{index_subfeatures};
    $self->{index_subfeatures} = shift if @_;
    $d;
}


sub summary_stats {
    my $self = shift;
    my $d    = $self->{summary_stats};
    $self->{summary_stats} = shift if @_;
    $d;
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
  
  if ($self->summary_stats) {
      $self->msg("Building summary statistics for coverage graphs...");
      my $start = $self->time();
      $self->build_summary;
      $self->msg(sprintf "coverage graph build time: %5.2fs\n",$self->time()-$start);
  }
  $self->msg(sprintf "total load time: %5.2fs\n",$self->time()-$start);
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
    $self->create_load_data;
    $self->store->start_bulk_update() if $self->fast;
}

sub create_load_data {
    my $self = shift;
    $self->{load_data}{CurrentFeature}   = undef;
    $self->{load_data}{CurrentID}        = undef;
    $self->{load_data}{IndexIt}          = {};
    $self->{load_data}{Local2GlobalID}   = {};
    $self->{load_data}{count}            = 0;
    $self->{load_data}{mode}             = undef;
    $self->{load_data}{start_time}       = 0;
}

sub delete_load_data {
    my $self = shift;
    delete $self->{load_data};
}

sub finish_load {
  my $self  = shift;

  $self->store_current_feature();      # during fast loading, we will have a feature left at the very end
  $self->start_or_finish_sequence();   # finish any half-loaded sequences

  if ($self->fast) {
    $self->{load_data}{start_time} = $self->time();
    $self->store->finish_bulk_update;
  }
  $self->msg(sprintf "%5.2fs\n",$self->time()-$self->{load_data}{start_time});
  eval {$self->store->commit};

  # don't delete load data so that caller can ask for the loaded IDs
  # $self->delete_load_data;
}

=item build_summary

  $loader->build_summary

Call this to rebuild the summary coverage statistics. This is done automatically
if new() was passed a true value for -summary_stats at create time.

=cut

sub build_summary {
    my $self = shift;
    $self->store->build_summary_statistics;
}

=item do_load

  $count = $loader->do_load($fh)

This is called by load_fh() to load the GFF3 file's filehandle and
return the number of lines loaded.

=cut

sub do_load {
  my $self = shift;
  my $fh   = shift;

  $self->{load_data}{start_time}       = $self->time();
  $self->{load_data}->{millenium_time} = $self->{load_data}{start_time};
  $self->load_line($_) while <$fh>;
  $self->msg(sprintf "%d features loaded in %5.2fs%s\r",
	     $self->{load_data}->{count},
	     $self->time()-$self->{load_data}{start_time},
	     ' 'x80
      );
  $self->{load_data}{count};
}

=item load_line

    $loader->load_line($data);

Load a line of a GFF3 file. You must bracket this with calls to
start_load() and finish_load()!

    $loader->start_load();
    $loader->load_line($_) while <FH>;
    $loader->finish_load();

=cut

sub load_line {
    my $self = shift;
    my $line = shift;
    # don't do anything
}


=item handle_feature

  $loader->handle_feature($data_line)

This method is called to process a single data line. It manipulates
information stored a data structure called $self-E<gt>{load_data}.

=cut

sub handle_feature {
  my $self     = shift;
  my $line = shift;
  # do nothing 
}

=item handle_meta

  $loader->handle_meta($data_line)

This method is called to process a single data line. It manipulates
information stored a data structure called $self-E<gt>{load_data}.

=cut

sub handle_meta {
  my $self     = shift;
  my $line = shift;
  # do nothing 
}

sub _indexit {
    my $self      = shift;
    my $id        = shift;
    $id         ||= '';     # avoid uninit warnings
    my $indexhash = $self->{load_data}{IndexIt};
    $indexhash->{$id} = shift if @_;
    return $indexhash->{$id};
}

sub _local2global {
    my $self      = shift;
    my $id        = shift;
    $id         ||= '';  # avoid uninit warnings
    my $indexhash = $self->{load_data}{Local2GlobalID};
    $indexhash->{$id} = shift if @_;
    return $indexhash->{$id};
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
  my $indexed    = $self->_indexit($ld->{CurrentID});

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

  $self->_local2global($ld->{CurrentID} => $id);
  $self->_indexit($ld->{CurrentID} => 0)if $normalized;  # no need to remember this
  undef $ld->{CurrentID};
  undef $ld->{CurrentFeature};
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
  # do nothing
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

This method opens up the indicated file or pipe, using some
intelligence to recognized compressed files and URLs and doing the
right thing.

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
  return $thing                                  if ref $thing && $thing->isa('IO::String');
  return IO::File->new($thing);
}

sub msg {
  my $self = shift;
  my @msg  = @_;
  return unless $self->verbose;
  print STDERR @msg;
}

=item loaded_ids

 my $ids    = $loader->loaded_ids;
 my $id_cnt = @$ids;

After performing a load, this returns an array ref containing all the
feature primary ids that were created during the load.

=cut

sub loaded_ids {
    my $self = shift;
    my @ids  = values %{$self->{load_data}{Local2GlobalID}}
                     if $self->{load_data};
    return \@ids;
}

=item local_ids

 my $ids    = $self->local_ids;
 my $id_cnt = @$ids;

After performing a load, this returns an array ref containing all the
load file IDs that were contained within the file just loaded.

=cut

sub local_ids {
    my $self = shift;
    my @ids  = keys %{$self->{load_data}{Local2GlobalID}}
                   if $self->{load_data};
    return \@ids;
}

=item time

 my $time = $loader->time

This method returns the current time in seconds, using Time::HiRes if available.

=cut

sub time {
  return Time::HiRes::time() if Time::HiRes->can('time');
  return time();
}

=item unescape

 my $unescaped = GFF3Loader::unescape($escaped)

This is an internal utility.  It is the same as CGI::Util::unescape,
but doesn't change pluses into spaces and ignores unicode escapes.

=cut

sub unescape {
    my $self = shift;
    my $todecode = shift;
    $todecode =~ s/%([0-9a-fA-F]{2})/chr hex($1)/ge;
    return $todecode;
}

sub DESTROY {
    my $self = shift;
    # Close filehandles, so temporal files can be properly deleted
    my $store = $self->store;
    if (   $store->isa('Bio::DB::SeqFeature::Store::memory')
	or $store->isa('Bio::DB::SeqFeature::Store::berkeleydb3')
	) {
      $store->private_fasta_file->close;

      if ($store->{fasta_db}) {
	while (my ($file, $fh) = each %{ $store->{fasta_db}->{fhcache} }) {
	  $fh->close;
	}
	$store->{fasta_db}->_close_index($store->{fasta_db}->{offsets});
      }
    }
    elsif ($store->isa('Bio::DB::SeqFeature::Store::DBI::SQLite')) {
      if (%DBI::installed_drh) {
	DBI->disconnect_all;
	%DBI::installed_drh = ();
      }
      undef $store->{dbh};
    }

    if (my $ld = $self->{temp_load}) {
	unlink $ld;
    }
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
L<Bio::DB::SeqFeature::Store::GFF3Loader>,
L<Bio::DB::SeqFeature::Store::DBI::mysql>,
L<Bio::DB::SeqFeature::Store::berkeleydb>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2006 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut
