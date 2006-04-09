package Bio::DB::SeqFeature::Store::GFF3Loader;

# $Id$

# load utility - incrementally load the store based on GFF3 file
#
# two modes:
#   slow mode -- features can occur in any order in the GFF3 file
#   fast mode -- all features with same ID must be contiguous in GFF3 file

use strict;
use Carp 'croak';
use IO::File;
use Bio::DB::GFF::Util::Rearrange;
use Bio::DB::SeqFeature::Store::Cacher;
use CGI::Util 'unescape';
use base 'Bio::Root::Root';

use constant DEFAULT_SEQ_CHUNK_SIZE => 2000;

my %Special_attributes =(
			 Gap    => 1, Target => 1,
			 Parent => 1, Name   => 1,
			 Alias  => 1, ID     => 1,
			 Index  => 1,
			);
my %Strandedness = ( '+' => 1,
		     '-' => -1,
		     '.' => 0,
		     ''  => 0,
		   );

sub new {
  my $self = shift;
  my ($store,$seqfeature_class,$tmpdir,$verbose,$fast,$seq_chunk_size) = rearrange(['STORE',
										    ['SF_CLASS','SEQFEATURE_CLASS'],
										    'TMP',
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

  my $tmp_store = Bio::DB::SeqFeature::Store::Cacher->new(-adaptor=>'bdb',-tmp=>1,-dir=>$tmpdir) unless $normalized;

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

sub default_seqfeature_class {
  my $self = shift;
  return 'Bio::DB::SeqFeature::LazyTableFeature';
}

sub store          { shift->{store}            }
sub tmp_store      { shift->{tmp_store}        }
sub sfclass        { shift->{seqfeature_class} }
sub fast           { shift->{fast}             }
sub seq_chunk_size { shift->{seq_chunk_size}             }
sub verbose        { shift->{verbose}          }

sub subfeatures_normalized {
  my $self = shift;
  my $d    = $self->{subfeatures_normalized};
  $self->{subfeatures_normalized} = shift if @_;
  $d;
}

sub subfeatures_in_table {
  my $self = shift;
  my $d    = $self->{subfeatures_in_table};
  $self->{subfeatures_in_table} = shift if @_;
  $d;
}

sub load {
  my $self       = shift;
  my $start      = $self->time();

  for my $file_or_fh (@_) {
    $self->msg("loading $file_or_fh...\n");
    my $fh = $self->open_fh($file_or_fh) or $self->throw("Couldn't open $file_or_fh: $!");
    $self->load_fh($fh);
    $self->msg(sprintf "load time: %5.2fs\n",$self->time()-$start);
  }
}

sub load_fh {
  my $self = shift;
  my $fh   = shift;
  $self->start_load();
  $self->do_load($fh);
  $self->finish_load();
}

sub start_load {
  my $self = shift;
  $self->{load_data}{Parent2Child}     = {};
  $self->{load_data}{Local2GlobalID}   = {};
  $self->{load_data}{TemporaryID}      = 1;
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
  delete $self->{load_data};
}

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

    if (/^\#\#\s*(.+)/) {  ## meta instruction
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
}

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

  if (/^\#\#\s*index-subfeatures\s+(\S+)/) {
    $self->{load_data}{IndexSubfeatures} = $1;
    $self->store->index_subfeatures($1);
    return;
  }
}

sub handle_feature {
  my $self     = shift;
  my $gff_line = shift;
  my $ld       = $self->{load_data};

  my @columns = map {$_ eq '.' ? undef : $_ } split /\t/,$gff_line;
  return unless @columns >= 8;
  my ($refname,$source,$method,$start,$end, $score,$strand,$phase,$attributes)      = @columns;
  $strand = $Strandedness{$strand};

  my ($reserved,$unreserved) = $self->parse_attributes($attributes);

  my $name        = ($reserved->{Name}   && $reserved->{Name}[0]);

  my $feature_id  = $reserved->{ID}[0] || $ld->{TemporaryID}++;
  my @parent_ids  = @{$reserved->{Parent}} if $reserved->{Parent};

  my $index_it = $ld->{IndexSubfeatures};
  if (exists $reserved->{Index}) {
    $index_it = $reserved->{Index}[0];
  }

  # Everything in the unreserved hash becomes an attribute, so we copy
  # some attributes over
  $unreserved->{Note}  = $reserved->{Note}   if exists $reserved->{Note};
  $unreserved->{Alias} = $reserved->{Alias}  if exists $reserved->{Alias};
  $unreserved->{Target}= $reserved->{Target} if exists $reserved->{Target};
  $unreserved->{Gap}   = $reserved->{Gap}    if exists $reserved->{Gap};

  # TEMPORARY HACKS TO SIMPLIFY DEBUGGING
  $unreserved->{load_id}   = $feature_id    if defined $feature_id;
  push @{$unreserved->{Alias}},$feature_id  if defined $feature_id;
  $unreserved->{parent_id} = \@parent_ids   if @parent_ids;

  my @args = (-display_name => $name || undef,
	      -seq_id       => $refname,
	      -start        => $start,
	      -end          => $end,
	      -strand       => $strand || 0,
	      -score        => $score  || undef,
	      -phase        => $phase  || undef,
	      -primary_tag  => $method || 'feature',
	      -source       => $source || undef,
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
END
  }

  # contiguous feature, so add a segment
  if ($old_feat) {
    $self->add_segment($old_feat,$self->sfclass->new(@args));
    return;
  }

  # we get here if this is a new feature
  # first of all, store the current feature if it is there
  $self->store_current_feature() if defined $ld->{CurrentID};

  # now create the new feature
  # (index top-level features only if policy asks us to)
  my $feature = $self->sfclass->new(@args);
  eval {$feature->object_store($self->store)};  # for lazy table features
  $ld->{CurrentFeature} = $feature;
  $ld->{CurrentID}      = $feature_id;

  $ld->{IndexIt}{$feature_id}++    if $index_it || !@parent_ids || !(defined $reserved->{ID}[0]);
  $ld->{TopLevel}{$feature_id}++   if !$self->{fast} && !@parent_ids;  # need to track top level features

  # remember parentage
  for my $parent (@parent_ids) {
    push @{$ld->{Parent2Child}{$parent}},$feature_id;
  }

}

sub store_current_feature {
  my $self    = shift;

  my $ld   = $self->{load_data};
  my $f    = $ld->{CurrentFeature} or return;

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
    $self->tmp_store->store($f)
  }
	
  my $id        = $f->primary_id;    # assigned by store()
  $ld->{Local2GlobalID}{$ld->{CurrentID}} = $id;

  undef $ld->{CurrentID};
  undef $ld->{CurrentFeature};
  undef $ld->{IndexIt}{$ld->{CurrentID}} if $normalized;  # no need to remember this
}

###
# put objects together
#
sub build_object_tree {
  my $self = shift;
  $self->subfeatures_in_table ? $self->build_object_tree_in_tables : $self->build_object_tree_in_features;
}

sub build_object_tree_in_tables {
  my $self = shift;
  my $store = $self->store;
  my $ld    = $self->{load_data};

  while (my ($load_id,$children) = each %{$ld->{Parent2Child}}) {
    my $parent_id = $ld->{Local2GlobalID}{$load_id} or die "$load_id doesn't have a primary id";
    my @children  = map {$ld->{Local2GlobalID}{$_}} @$children;

    # this updates the table that keeps track of parent/child relationships,
    # but does not update the parent object -- so (start,end) had better be right!!!
    $store->add_SeqFeature($parent_id,@children);
  }

}

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

sub add_segment {
  my $self = shift;
  my ($parent,$child) = @_;

  if ($parent->can('add_segment')) { # probably a lazy table feature
    my @segments = $parent->can('denormalized_segments')
      ? $parent->denormalized_segments 
      : $parent->segments;
    unless (@segments) {  # convert into a segmented object
      my %clone   = %$parent;
      my $segment = bless \%clone,ref $parent;
      delete $segment->{segments};
      eval {$segment->object_store(undef) };
      $segment->primary_id(undef);

      # this updates the object and expands its start and end positions without writing
      # the segments into the database as individual objects
      $parent->add_segment($segment);
    }
    $parent->add_segment($child);
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

# this gets called at the beginning and end of a fasta section
sub start_or_finish_sequence {
  my $self  = shift;
  my $seqid = shift;
  if (my $sl    = $self->{fasta_load}) {
    if (defined $sl->{seqid}) {
      $self->store->insert_sequence($sl->{seqid},$sl->{offset},$sl->{sequence});
      delete $self->{fasta_load};
    }
  }
  if (defined $seqid) {
    $self->{fasta_load} = {seqid  => $seqid,
			   offset => 0,
			   sequence => ''};
  }
}

sub load_sequence {
  my $self = shift;
  my $seq  = shift;
  my $sl   = $self->{fasta_load} or return;
  my $cs   = $self->seq_chunk_size;
  $sl->{sequence} .= $seq;
  while (length $sl->{sequence} >= $cs) {
    my $chunk = substr($sl->{sequence},0,$cs);
    $self->store->insert_sequence($sl->{seqid},$sl->{offset},$chunk);
    $sl->{offset} += length $chunk;
    substr($sl->{sequence},0,$cs) = '';
  }
}

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

sub time {
  return Time::HiRes::time() if Time::HiRes->can('time');
  return time();
}

1;
