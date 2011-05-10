package Bio::DB::SeqFeature::Store::GFF3Loader;


=head1 NAME

Bio::DB::SeqFeature::Store::GFF3Loader -- GFF3 file loader for Bio::DB::SeqFeature::Store

=head1 SYNOPSIS

  use Bio::DB::SeqFeature::Store;
  use Bio::DB::SeqFeature::Store::GFF3Loader;

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

The GFF3 file format has been extended very slightly to accommodate
Bio::DB::SeqFeature::Store. First, the loader recognizes is a new
directive:

  # #index-subfeatures [0|1]

Note that you can place a space between the two #'s in order to
prevent GFF3 validators from complaining.

If this is true, then subfeatures are indexed (the default) so that
they can be retrieved with a query. See L<Bio::DB::SeqFeature::Store>
for an explanation of this. If false, then subfeatures can only be
accessed through their parent feature.

Second, the loader recognizes a new attribute tag called index, which
if present, controls indexing of the current feature. Example:

 ctg123	. TF_binding_site 1000 1012 . + . ID=tfbs00001;index=1

You can use this to turn indexing on and off, overriding the default
for a particular feature.

Note that the loader keeps a record -- in memory -- of each feature
that it has processed. If you find the loader running out of memory on
particularly large GFF3 files, please split the input file into
smaller pieces and do the load in steps.

=cut


# load utility - incrementally load the store based on GFF3 file
#
# two modes:
#   slow mode -- features can occur in any order in the GFF3 file
#   fast mode -- all features with same ID must be contiguous in GFF3 file

use strict;
use Carp 'croak';
use Bio::DB::GFF::Util::Rearrange;
use Bio::DB::SeqFeature::Store::LoadHelper;

use base 'Bio::DB::SeqFeature::Store::Loader';


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

 -ignore_seqregion  Ignore ##sequence-region directives. The default is to create a
                       feature corresponding to the directive.

 -noalias_target    Don't create an Alias attribute for a target_id named in a 
                    Target attribute. The default is to create an Alias
                    attribute containing the target_id found in a Target 
                    attribute.

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
have been seen. This temporary databases uses the "berkeleydb"
adaptor. The -tmp option specifies the directory in which that
database will be created. If not present, it defaults to the system
default tmp directory specified by File::Spec-E<gt>tmpdir().

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
    my $class = shift;
    my $self  = $class->SUPER::new(@_);
    my ($ignore_seqregion) = rearrange(['IGNORE_SEQREGION'],@_);
    $self->ignore_seqregion($ignore_seqregion);
    my ($noalias_target) = rearrange(['NOALIAS_TARGET'],@_);
    $self->noalias_target($noalias_target);
    $self;
}

=head2 ignore_seqregion

  $ignore_it = $loader->ignore_seqregion([$new_flag])

Get or set the ignore_seqregion flag, which if true, will cause 
GFF3 ##sequence-region directives to be ignored. The default behavior
is to create a feature corresponding to the region.

=cut

sub ignore_seqregion {
    my $self = shift;
    my $d    = $self->{ignore_seqregion};
    $self->{ignore_seqregion} = shift if @_;
    $d;
}

=head2 noalias_target

  $noalias_target = $loader->noalias_target([$new_flag])

Get or set the noalias_target flag, which if true, will disable the creation of
an Alias attribute for a target_id named in a Target attribute. The default is 
to create an Alias attribute containing the target_id found in a Target 
attribute.

=cut

sub noalias_target {
    my $self = shift;
    my $d    = $self->{noalias_target};
    $self->{noalias_target} = shift if @_;
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

Load the GFF3 data at the other end of the filehandle and return true
if successful. Internally, load_fh() invokes:

  start_load();
  do_load($filehandle);
  finish_load();

=cut

# sub load_fh { } inherited

=item start_load, finish_load

These methods are called at the start and end of a filehandle load.

=cut

sub create_load_data { #overridden
  my $self = shift;
  $self->SUPER::create_load_data;
  $self->{load_data}{TemporaryID}      = "GFFLoad0000000";
  $self->{load_data}{IndexSubfeatures} = $self->index_subfeatures();
  $self->{load_data}{mode}             = 'gff';

  $self->{load_data}{Helper}           = 
      Bio::DB::SeqFeature::Store::LoadHelper->new($self->{tmpdir});
}

sub finish_load { #overridden
  my $self  = shift;

  $self->store_current_feature();      # during fast loading, we will have a feature left at the very end
  $self->start_or_finish_sequence();   # finish any half-loaded sequences

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

  # don't delete load data so that caller can ask for the loaded IDs
  # $self->delete_load_data;
}

=item do_load

  $count = $loader->do_load($fh)

This is called by load_fh() to load the GFF3 file's filehandle and
return the number of lines loaded.

=cut

# sub do_load { } inherited

=item load_line

    $loader->load_line($data);

Load a line of a GFF3 file. You must bracket this with calls to
start_load() and finish_load()!

    $loader->start_load();
    $loader->load_line($_) while <FH>;
    $loader->finish_load();

=cut

sub load_line { #overridden
    my $self = shift;
    my $line = shift;

    chomp($line);
    my $load_data = $self->{load_data};
    $load_data->{line}++;

    return unless $line =~ /^\S/;     # blank line

    # if it has a tab in it or looks like a chrom.sizes file, switch to gff mode
    $load_data->{mode} = 'gff' if $line =~ /\t/
	or $line =~ /^\w+\s+\d+\s*$/;

    if ($line =~ /^\#\s?\#\s*(.+)/) {  ## meta instruction
      $load_data->{mode} = 'gff';
      $self->handle_meta($1);

    } elsif ($line =~ /^\#/) {
      $load_data->{mode} = 'gff';  # just to be safe
      return; # comment
    }

    elsif ($line =~ /^>\s*(\S+)/) { # FASTA lines are coming
	$load_data->{mode} = 'fasta';
	$self->start_or_finish_sequence($1);
    }

    elsif ($load_data->{mode} eq 'fasta') {
      $self->load_sequence($line);
    }

    elsif ($load_data->{mode} eq 'gff') {
      $self->handle_feature($line);
      if (++$load_data->{count} % 1000 == 0) {
	my $now = $self->time();
	my $nl = -t STDOUT && !$ENV{EMACS} ? "\r" : "\n";
	local $^W = 0; # kill uninit variable warning
	$self->msg(sprintf("%d features loaded in %5.2fs (%5.2fs/1000 features)...%s$nl",
			   $load_data->{count},$now - $load_data->{start_time},
			   $now - $load_data->{millenium_time},
			   ' ' x 80
		   ));
	$load_data->{millenium_time} = $now;
      }
    }

    else {
      $self->throw("I don't know what to do with this line:\n$line");
    }
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

  if ( $instruction =~ /^#$/ ) {
    $self->store_current_feature() ;   # during fast loading, we will have a feature left at the very end
    $self->start_or_finish_sequence();    # finish any half-loaded sequences
    if ( $self->store->can('handle_resolution_meta') ) {
      $self->store->handle_resolution_meta($instruction);
    }
    return;
  }

  if ($instruction =~ /sequence-region\s+(.+)\s+(-?\d+)\s+(-?\d+)/i 
      && !$self->ignore_seqregion()) {
      my($ref,$start,$end,$strand)    = $self->_remap($1,$2,$3,+1);
      my $feature = $self->sfclass->new(-name        => $ref,
					-seq_id      => $ref,
					-start       => $start,
					-end         => $end,
					-strand      => $strand,
					-primary_tag => 'region');
    $self->store->store($feature);
    return;
  }

  if ($instruction =~/index-subfeatures\s+(\S+)/i) {
    $self->{load_data}{IndexSubfeatures} = $1;
    $self->store->index_subfeatures($1);
    return;
  }

  if ( $self->store->can('handle_unrecognized_meta') ) {
    $self->store->handle_unrecognized_meta($instruction);
    return;
  }
}

=item handle_feature

  $loader->handle_feature($gff3_line)

This method is called to process a single GFF3 line. It manipulates
information stored a data structure called $self-E<gt>{load_data}.

=cut

sub handle_feature { #overridden
  my $self     = shift;
  my $gff_line = shift;
  my $ld       = $self->{load_data};

  my $allow_whitespace = $self->allow_whitespace;

  # special case for a chrom.sizes-style line
  my @columns;
  if ($gff_line =~ /^(\w+)\s+(\d+)\s*$/) {
      @columns = ($1,undef,'chromosome',1,$2,undef,undef,undef,"Name=$1");
  } else {
      $gff_line    =~ s/\s+/\t/g if $allow_whitespace;
      @columns = map {$_ eq '.' ? undef : $_ } split /\t/,$gff_line;
  }

  $self->invalid_gff($gff_line) if @columns < 4;
  $self->invalid_gff($gff_line) if @columns > 9 && $allow_whitespace;

  {
      local $^W = 0;
      if (@columns > 9) { #oops, split too much due to whitespace
	  $columns[8] = join(' ',@columns[8..$#columns]);
      }
  }

  my ($refname,$source,$method,$start,$end,$score,$strand,$phase,$attributes) = @columns;
  
  $self->invalid_gff($gff_line) unless defined $refname;
  $self->invalid_gff($gff_line) unless !defined $start || $start =~ /^[\d.-]+$/;
  $self->invalid_gff($gff_line) unless !defined $end   || $end   =~ /^[\d.-]+$/;
  $self->invalid_gff($gff_line) unless defined $method;

  $strand = $Strandedness{$strand||0};
  my ($reserved,$unreserved) = $attributes ? $self->parse_attributes($attributes) : ();

  my $name        = ($reserved->{Name}   && $reserved->{Name}[0]);

  my $has_loadid  = defined $reserved->{ID}[0];

  my $feature_id  = defined $reserved->{ID}[0] ? $reserved->{ID}[0] : $ld->{TemporaryID}++;
  my @parent_ids  = @{$reserved->{Parent}}     if defined $reserved->{Parent};

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

  # mec@stowers-institute.org, wondering why not all attributes are
  # carried forward, adds ID tag in particular service of
  # round-tripping ID, which, though present in database as load_id
  # attribute, was getting lost as itself
  # $unreserved->{ID}= $reserved->{ID}     if exists $reserved->{ID}; 

  # TEMPORARY HACKS TO SIMPLIFY DEBUGGING
  $feature_id = '' unless defined $feature_id;
  $name       = '' unless defined $name;  # prevent uninit variable warnings
  # push @{$unreserved->{Alias}},$feature_id  if $has_loadid && $feature_id ne $name;
  $unreserved->{parent_id} = \@parent_ids   if @parent_ids;

  # POSSIBLY A PERMANENT HACK -- TARGETS BECOME ALIASES
  # THIS IS TO ALLOW FOR TARGET-BASED LOOKUPS
  if (exists $reserved->{Target} && !$self->{noalias_target}) {
    my %aliases = map {$_=>1} @{$unreserved->{Alias}};
    for my $t (@{$reserved->{Target}}) {
      (my $tc = $t) =~ s/\s+.*$//;  # get rid of coordinates
      $name ||= $tc;
      push @{$unreserved->{Alias}},$tc unless $name eq $tc || $aliases{$tc};
    }
  }

  ($refname,$start,$end,$strand) = $self->_remap($refname,$start,$end,$strand) or return;

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
  elsif (defined(my $id = $self->{load_data}{Helper}->local2global($feature_id))) {
    $old_feat = $self->fetch($feature_id)
      or $self->warn(<<END);
ID=$feature_id has been used more than once, but it cannot be found in the database.
This can happen if you have specified fast loading, but features sharing the same ID
are not contiguous in the GFF file. This will be loaded as a separate feature.
Line $.: "$_"
END
  }

  # contiguous feature, so add a segment
  warn $old_feat if defined $old_feat and !ref $old_feat;
  if (defined $old_feat) {
      # set this to 1 to disable split-location behavior
      if (0 && @parent_ids) {                  # If multiple features are held together by the same ID
	  $feature_id = $ld->{TemporaryID}++;  # AND they have a Parent attribute, this causes an undesirable
      }                                        # additional layer of aggregation. Changing the ID fixes this.
      elsif       (
	  $old_feat->seq_id ne $refname || 
	  $old_feat->start  != $start || 
	  $old_feat->end    != $end # make sure endpoints are distinct
	  )
      {
	  $self->add_segment($old_feat,$self->sfclass->new(@args));
	  return;
      }
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

  my $helper = $ld->{Helper};
  $helper->indexit($feature_id=>1)  if $index_it;
  $helper->toplevel($feature_id=>1) if !$self->{fast} 
                                       && $top_level;  # need to track top level features


  # remember parentage
  for my $parent (@parent_ids) {
      $helper->add_children($parent=>$feature_id);
  }

}

sub invalid_gff {
    my $self = shift;
    my $line = shift;
    $self->throw("invalid GFF line at line $self->{load_data}{line}.\n".$line);
}

=item allow_whitespace

   $allow_it = $loader->allow_whitespace([$newvalue]);

Get or set the allow_whitespace flag. If true, then GFF3 files are
allowed to be delimited with whitespace in addition to tabs.

=cut

sub allow_whitespace {
    my $self = shift;
    my $d    = $self->{allow_whitespace};
    $self->{allow_whitespace} = shift if @_;
    $d;
}

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
  my $store  = $self->store;
  my $helper = $self->{load_data}{Helper};

  while (my ($load_id,$children) = $helper->each_family()) {

      my $parent_id = $helper->local2global($load_id);
      die $self->throw("$load_id doesn't have a primary id") 
	  unless defined $parent_id;

      my @children  = map {$helper->local2global($_)} @$children;
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

  my $helper     = $ld->{Helper};

  while (my $load_id = $helper->each_toplevel) {
    my $feature  = $self->fetch($load_id)
      or $self->throw("$load_id (id="
		      .$helper->local2global($load_id)
		      ." should have a database entry, but doesn't");
    $self->attach_children($store,$ld,$load_id,$feature);
    # Indexed objects are updated, not created anew
    $feature->primary_id(undef) unless $helper->indexit($load_id);
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

  my $children   = $ld->{Helper}->children() or return;
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
  my $helper  = $self->{load_data}{Helper};
  my $id      = $helper->local2global($load_id);

  return
      ($self->subfeatures_normalized || $helper->indexit($load_id)
       ? $self->store->fetch($id)
       : $self->tmp_store->fetch($id)
      );
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

  unless ($att =~ /=/) {  # ouch! must be a GFF line
      require Bio::DB::SeqFeature::Store::GFF2Loader
	  unless Bio::DB::SeqFeature::Store::GFF2Loader->can('parse_attributes');
      return $self->Bio::DB::SeqFeature::Store::GFF2Loader::parse_attributes($att);
  }

  my @pairs =  map { my ($name,$value) = split '=';
                    [$self->unescape($name) => $value];  
                   } split ';',$att;
  my (%reserved,%unreserved);
  foreach (@pairs) {
    my $tag    = $_->[0];

    unless (defined $_->[1]) {
      warn "$tag does not have a value at GFF3 file line $.\n";
      next;
    }

    my @values = split ',',$_->[1];
    map {$_ = $self->unescape($_);} @values;
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

 my $unescaped = GFF3Loader::unescape($escaped)

This is an internal utility.  It is the same as CGI::Util::unescape,
but doesn't change pluses into spaces and ignores unicode escapes.

=cut

# sub unescape { } inherited

sub _remap {
    my $self = shift;
    my ($ref,$start,$end,$strand) = @_;
    my $mapper = $self->coordinate_mapper;
    return ($ref,$start,$end,$strand) unless $mapper;

    my ($newref,$coords) = $mapper->($ref,[$start,$end]);
    return unless defined $coords->[0];
    if ($coords->[0] > $coords->[1]) {
	@{$coords} = reverse(@{$coords}); 
	$strand *= -1;
    }
    return ($newref,@{$coords},$strand);
}

sub _indexit { # override
    my $self      = shift;
    return $self->{load_data}{Helper}->indexit(@_);
}

sub _local2global { # override
    my $self      = shift;
    return $self->{load_data}{Helper}->local2global(@_);
}

=item local_ids

 my $ids    = $self->local_ids;
 my $id_cnt = @$ids;

After performing a load, this returns an array ref containing all the
load file IDs that were contained within the file just loaded.

=cut

sub local_ids { # override
    my $self = shift;
    return $self->{load_data}{Helper}->local_ids(@_);
}

=item loaded_ids

 my $ids    = $loader->loaded_ids;
 my $id_cnt = @$ids;

After performing a load, this returns an array ref containing all the
feature primary ids that were created during the load.

=cut

sub loaded_ids { # override
    my $self = shift;
    return $self->{load_data}{Helper}->loaded_ids(@_);
}


1;
__END__

=back

=head1 BUGS

This is an early version, so there are certainly some bugs. Please
use the BioPerl bug tracking system to report bugs.

=head1 SEE ALSO

L<Bio::DB::SeqFeature::Store>,
L<Bio::DB::SeqFeature::Segment>,
L<Bio::DB::SeqFeature::NormalizedFeature>,
L<Bio::DB::SeqFeature::GFF2Loader>,
L<Bio::DB::SeqFeature::Store::DBI::mysql>,
L<Bio::DB::SeqFeature::Store::berkeleydb>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2006 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut


