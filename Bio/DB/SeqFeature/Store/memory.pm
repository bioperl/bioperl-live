package Bio::DB::SeqFeature::Store::memory;

=head1 NAME

Bio::DB::SeqFeature::Store::memory -- In-memory implementation of Bio::DB::SeqFeature::Store

=head1 SYNOPSIS

  use Bio::DB::SeqFeature::Store;

  # Open the sequence database
  my $db      = Bio::DB::SeqFeature::Store->new( -adaptor => 'memory',
                                                 -dsn     => '/var/databases/test');
  # search... by id
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
  my @features = $segment->features(-type=>['mRNA','match']);

  # getting & storing sequence information
  # Warning: this returns a string, and not a PrimarySeq object
  $db->insert_sequence('Chr1','GATCCCCCGGGATTCCAAAA...');
  my $sequence = $db->fetch_sequence('Chr1',5000=>6000);

  # what feature types are defined in the database?
  my @types    = $db->types;

  # create a new feature in the database
  my $feature = $db->new_feature(-primary_tag => 'mRNA',
                                 -seq_id      => 'chr3',
                                 -start      => 10000,
                                 -end        => 11000);

=head1 DESCRIPTION

Bio::DB::SeqFeature::Store::memory is the in-memory adaptor for
Bio::DB::SeqFeature::Store. You will not create it directly, but
instead use Bio::DB::SeqFeature::Store-E<gt>new() to do so.

See L<Bio::DB::SeqFeature::Store> for complete usage instructions.

=head2 Using the memory adaptor

Before using the memory adaptor, populate a readable-directory on the
file system with annotation and/or sequence files. The annotation
files must be in GFF3 format, and sholud end in the extension .gff or
.gff3. They may be compressed with "compress", "gzip" or "bzip2" (in
which case the appropriate compression extension must be present as
well.)

You may include sequence data inline in the GFF3 files, or put the
sequence data in one or more separate FASTA-format files. These files
must end with .fa or .fasta and may be compressed. Because of the way
the adaptor works, you will get much better performance if you keep
the sequence data in separate FASTA files.

Initialize the database using the -dsn option. This should point to
the directory creating the annotation and sequence files, or to a
single GFF3 file. Examples:


  # load all GFF3 and FASTA files located in /var/databases/test directory
  $db  = Bio::DB::SeqFeature::Store->new( -adaptor => 'memory',
                                          -dsn     => '/var/databases/test');


  # load the data in a single compressed GFF3 file located at
  # /usr/annotations/worm.gf33.gz
  $db  = Bio::DB::SeqFeature::Store->new( -adaptor => 'memory',
                                          -dsn     => '/usr/annotations/worm.gff3.gz');

For compatibility with the Bio::DB::GFF memory adaptor, -gff is
recognized as an alias for -dsn.

See L<Bio::DB::SeqFeature::Store> for all the access methods supported
by this adaptor. The various methods for storing and updating features
and sequences into the database are supported, including GFF3 loading
support, but since this is an in-memory adaptor all changes you make
will be lost when the script exits.

=cut

use strict;
use base 'Bio::DB::SeqFeature::Store';
use Bio::DB::SeqFeature::Store::GFF3Loader;
use Bio::DB::GFF::Util::Rearrange 'rearrange';
use File::Temp 'tempdir';
use IO::File;
use Bio::DB::Fasta;
use File::Glob ':glob';

use constant BINSIZE => 10_000;

###
# object initialization
#
sub init {
  my ($self, $args) = @_;
  $self->SUPER::init($args);
  $self->{_data}     = {};
  $self->{_children} = {};
  $self->{_index}    = {};
  $self;
}

sub post_init {
  my $self = shift;
  my ($file_or_dir) = rearrange([['DIR','DSN','FILE','GFF']],@_);
  return unless $file_or_dir;

  my $loader = Bio::DB::SeqFeature::Store::GFF3Loader->new(-store    => $self,
							   -sf_class => $self->seqfeature_class) 
    or $self->throw("Couldn't create GFF3Loader");
  my @argv;
  if (-d $file_or_dir) {
    @argv = (
	     bsd_glob("$file_or_dir/*.size*"),
	     bsd_glob("$file_or_dir/*.gff"),            bsd_glob("$file_or_dir/*.gff3"),
	     bsd_glob("$file_or_dir/*.gff.{gz,Z,bz2}"), bsd_glob("$file_or_dir/*.gff3.{gz,Z,bz2}")
	     );
  } else {
    @argv = $file_or_dir;
  }
  local $self->{file_or_dir} = $file_or_dir;
  $loader->load(@argv);
  warn $@ if $@;
}

sub commit { # reindex fasta files
  my $self = shift;

  my $db;
  if (my $fh = $self->{fasta_fh}) {
      $fh->close;
      $db = Bio::DB::Fasta->new($self->{fasta_file});
  } elsif (exists $self->{file_or_dir} && -d $self->{file_or_dir}) {
      $db = Bio::DB::Fasta->new($self->{file_or_dir});
  }
  $self->{fasta_db} = $db if $db;
}

sub can_store_parentage { 1 }

# return a hash ref in which each key is primary id
sub data {
  shift->{_data};
}

sub _init_database { shift->init }

sub _store {
  my $self    = shift;
  my $indexed = shift;
  my @objs    = @_;
  my $data = $self->data;
  my $count = 0;
  for my $obj (@objs) {
    # Add unique ID to feature if needed
    my $primary_id = $self->_autoid($obj);
    # Store feature (overwriting any existing feature with the same primary ID
    # as required by Bio::DB::SF::Store)
    $data->{$primary_id} = $obj;
    if ($indexed) {
      $self->{_index}{ids}{$primary_id} = undef;
      $self->_update_indexes($obj);
    }
    $count++;
  }
  return $count;
}


sub _autoid {
  # If a feature has no ID, assign it a unique ID
  my ($self, $obj) = @_;
  my $data = $self->data;
  my $primary_id = $obj->primary_id;
  if (not defined $primary_id) {
    # Create a unique ID
    $primary_id = 1 + scalar keys %{$data};
    while (exists $data->{$primary_id}) {
      $primary_id++;
    }
    $obj->primary_id($primary_id);
  }
  return $primary_id;
}


sub _deleteid {
  my ($self, $id) = @_;
  if (exists $self->{_index}{ids}{$id}) {
    # $indexed was true
    $self->_update_indexes( $self->fetch($id), 1 );
    delete $self->{_index}{ids}{$id};
  }
  delete $self->data->{$id};
  return 1;
}

sub _fetch {
  my ($self, $id) = @_;
  return $self->data->{$id};
}

sub _add_SeqFeature {
  my ($self, $parent, @children) = @_;
  my $count = 0;
  my $parent_id = ref $parent ? $parent->primary_id : $parent;
  defined $parent_id or $self->throw("Parent $parent should have a primary ID");
  for my $child (@children) {
    my $child_id = ref $child ? $child->primary_id : $child;
    defined $child_id or $self->throw("Child $child should have a primary ID");
    $self->{_children}{$parent_id}{$child_id}++;
    $count++;
  }
  return $count;
}

sub _fetch_SeqFeatures {
  my ($self, $parent, @types) = @_;
  my $parent_id = $parent->primary_id;
  defined $parent_id or $self->throw("Parent $parent should have a primary ID");
  my @children_ids  = keys %{$self->{_children}{$parent_id}};
  my @children      = map {$self->fetch($_)} @children_ids;
  
  if (@types) {
    my $data;
    for my $c (@children) {
      push @{$$data{$c->primary_tag}{$c->source_tag||''}}, $c;
    }
    @children = ();
    for my $type (@types) {
      $type .= ':' if (not $type =~ m/:/);
      my ($primary_tag, undef, $source_tag) = ($type =~ m/^(.*?)(:(.*?))$/);
      $source_tag ||= '';
      if ($source_tag eq '') {
        for my $source (keys %{$$data{$primary_tag}}) {
          if (exists $$data{$primary_tag}{$source_tag}) {
            push @children, @{$$data{$primary_tag}{$source_tag}};
          }
        }
      } else {
        if (exists $$data{$primary_tag}{$source_tag}) {
          push @children, @{$$data{$primary_tag}{$source_tag}};
        }
      }
    }
  }

  return @children;
}

sub _update_indexes {
  my ($self, $obj, $del) = @_;
  defined (my $id = $obj->primary_id) or return;
  $del ||= 0;
  $self->_update_name_index($obj,$id, $del);
  $self->_update_type_index($obj,$id, $del);
  $self->_update_location_index($obj, $id, $del);
  $self->_update_attribute_index($obj,$id, $del);
}

sub _update_name_index {
  my ($self, $obj, $id, $del) = @_;
  my ($names, $aliases) = $self->feature_names($obj);
  foreach (@$names) {
    if (not $del) {
      $self->{_index}{name}{lc $_}{$id}   = 1;
    } else {
      delete $self->{_index}{name}{lc $_}{$id};
      if (scalar keys %{ $self->{_index}{name}{lc $_} } == 0) {
        delete $self->{_index}{name}{lc $_};
      }
    };
  }
  foreach (@$aliases) {
    if (not $del) {
      $self->{_index}{name}{lc $_}{$id} ||= 2;
    } else {
      delete $self->{_index}{name}{lc $_}{$id};
      if (scalar keys %{ $self->{_index}{name}{lc $_} } == 0) {
        delete $self->{_index}{name}{lc $_};
      }
    }
  }
}

sub _update_type_index {
  my ($self, $obj, $id, $del) = @_;
  my $primary_tag = lc($obj->primary_tag) || return;
  my $source_tag  = lc($obj->source_tag || '');
  if (not $del) {
    $self->{_index}{type}{$primary_tag}{$source_tag}{$id} = undef;
  } else {
    delete $self->{_index}{type}{$primary_tag}{$source_tag}{$id};
    if ( scalar keys %{$self->{_index}{type}{$primary_tag}{$source_tag}} == 0 ) {
      delete $self->{_index}{type}{$primary_tag}{$source_tag};
      if (scalar keys %{$self->{_index}{type}{$primary_tag}} == 0 ) {
        delete $self->{_index}{type}{$primary_tag};
      }
    }
  }
}

sub _update_location_index {
  my ($self, $obj, $id, $del) = @_;
  my $seq_id      = $obj->seq_id || '';
  my $start       = $obj->start  || 0;
  my $end         = $obj->end    || 0;
  my $strand      = $obj->strand;
  my $bin_min     = int $start/BINSIZE;
  my $bin_max     = int $end/BINSIZE;
  for (my $bin = $bin_min; $bin <= $bin_max; $bin++ ) {
    if (not $del) {
      $self->{_index}{location}{lc $seq_id}{$bin}{$id} = undef;
    } else {
      delete $self->{_index}{location}{lc $seq_id}{$bin}{$id};
      if (scalar keys %{$self->{_index}{location}{lc $seq_id}{$bin}{$id}} == 0) {
        delete $self->{_index}{location}{lc $seq_id}{$bin}{$id};
      }
      if (scalar keys %{$self->{_index}{location}{lc $seq_id}{$bin}} == 0) {
        delete $self->{_index}{location}{lc $seq_id}{$bin};
      }
      if (scalar keys %{$self->{_index}{location}{lc $seq_id}} == 0) {
        delete $self->{_index}{location}{lc $seq_id};
      }
    }
  }
}

sub _update_attribute_index {
  my ($self, $obj, $id, $del) = @_;
  for my $tag ($obj->get_all_tags) {
    for my $value ($obj->get_tag_values($tag)) {
      if (not $del) {
        $self->{_index}{attribute}{lc $tag}{lc $value}{$id} = undef;
      } else {
        delete $self->{_index}{attribute}{lc $tag}{lc $value}{$id};
        if ( scalar keys %{$self->{_index}{attribute}{lc $tag}{lc $value}} == 0) {
          delete $self->{_index}{attribute}{lc $tag}{lc $value};
        }
        if ( scalar keys %{$self->{_index}{attribute}{lc $tag}} == 0) {
          delete $self->{_index}{attribute}{lc $tag};
        }
        if ( scalar keys %{$self->{_index}{attribute}} == 0) {
          delete $self->{_index}{attribute};
        }
      }
    }
  }
}

sub _features {
  my $self = shift;
  my ($seq_id,$start,$end,$strand,
      $name,$class,$allow_aliases,
      $types,
      $attributes,
      $range_type,
      $iterator
     ) = rearrange([['SEQID','SEQ_ID','REF'],'START',['STOP','END'],'STRAND',
		    'NAME','CLASS','ALIASES',
		    ['TYPES','TYPE','PRIMARY_TAG'],
		    ['ATTRIBUTES','ATTRIBUTE'],
		    'RANGE_TYPE',
		    'ITERATOR',
		   ],@_);

  my (@from,@where,@args,@group);
  $range_type ||= 'overlaps';

  my @result;
  unless (defined $name or defined $seq_id or defined $types or defined $attributes) {
    @result = keys %{$self->{_index}{ids}};
  }

  my %found  = ();
  my $result = 1;

  if (defined($name)) {
    # hacky backward compatibility workaround
    undef $class if $class && $class eq 'Sequence';
    $name     = "$class:$name" if defined $class && length $class > 0;
    $result &&= $self->filter_by_name($name,$allow_aliases,\%found);
  }

  if (defined $seq_id) {
    $result &&= $self->filter_by_location($seq_id,$start,$end,$strand,$range_type,\%found);
  }

  if (defined $types) {
    $result &&= $self->filter_by_type($types,\%found);
  }

  if (defined $attributes) {
    $result &&= $self->filter_by_attribute($attributes,\%found);
  }

  push @result,keys %found if $result;
  return $iterator ? Bio::DB::SeqFeature::Store::memory::Iterator->new($self,\@result)
                   : map {$self->fetch($_)} @result;
}


sub filter_by_type {
  my ($self, $types_req, $filter) = @_;
  my @types_req = ref $types_req eq 'ARRAY' ?  @$types_req : $types_req;

  my $types = $self->{_index}{type};
  my @types_found = $self->find_types(\@types_req);

  my @results;
  for my $type_found (@types_found) {
    my ($primary_tag, undef, $source_tag) = ($type_found =~ m/^(.*?)(:(.*?))$/);
    next unless exists $types->{$primary_tag}{$source_tag};
    push @results, keys %{$types->{$primary_tag}{$source_tag}};
  }

  $self->update_filter($filter,\@results);
}

sub find_types {
  my ($self, $types_req) = @_;
  my @types_found;

  my $types = $self->{_index}{type};

  for my $type_req (@$types_req) {

    # Type is the primary tag and an optional source tag
    my ($primary_tag, $source_tag);
    if (ref $type_req && $type_req->isa('Bio::DB::GFF::Typename')) {
      $primary_tag = $type_req->method;
      $source_tag  = $type_req->source;
    } else {
      ($primary_tag, undef, $source_tag) = ($type_req =~ m/^(.*?)(:(.*))?$/); 
    }
    ($primary_tag, $source_tag) = (lc $primary_tag, lc($source_tag || ''));

    next if not exists $$types{$primary_tag};

    if ($source_tag eq '') {
      # Match all sources for this primary_tag
      push @types_found, map {"$primary_tag:$_"} (keys %{ $$types{$primary_tag} });
    } else {
      # Match only the requested source
      push @types_found, "$primary_tag:$source_tag";
    }

  }
  return @types_found;
}

sub attributes {
  my $self = shift;
  return keys %{$self->{_index}{attribute}};
}

sub filter_by_attribute {
  my ($self, $attributes, $filter) = @_;

  my $index = $self->{_index}{attribute};
  my $result;

  for my $att_name (keys %$attributes) {
    my @result;
    my @matching_values;
    my @search_terms = ref($attributes->{$att_name}) && ref($attributes->{$att_name}) eq 'ARRAY'
                           ? @{$attributes->{$att_name}} : $attributes->{$att_name};
    my @regexp_terms;
    my @terms;

    for my $v (@search_terms) {
      if (my $regexp = $self->glob_match($v)) {
	@regexp_terms      = keys %{$index->{lc $att_name}} unless @regexp_terms;
	push @terms,grep {/^$v$/i} @regexp_terms;
      } else {
	push @terms,lc $v;
      }
    }

    for my $v (@terms) {
      push @result,keys %{$index->{lc $att_name}{$v}};
    }

    $result ||= $self->update_filter($filter,\@result);
  }

  $result;
}

sub filter_by_location {
  my ($self, $seq_id, $start, $end, $strand, $range_type, $filter) = @_;
  $strand ||= 0;

  my $index = $self->{_index}{location}{lc $seq_id};
  my @bins;

  if (!defined $start or !defined $end or $range_type eq 'contained_in') {
    @bins = sort {$a<=>$b} keys %{$index};
    $start = $bins[0]  * BINSIZE  unless defined $start;
    $end   = (($bins[-1] + 1) * BINSIZE) - 1 unless defined $end;
  }
  my %seenit;
  my $bin_min       = int $start/BINSIZE;
  my $bin_max       = int $end/BINSIZE;
  my @bins_in_range = $range_type eq 'contained_in' ? ($bins[0]..$bin_min,$bin_max..$bins[-1])
                                                    : ($bin_min..$bin_max);

  my @results;
  for my $bin (@bins_in_range) {
    next unless exists $index->{$bin};
    my @found = keys %{$index->{$bin}};
    for my $f (@found) {
      next if $seenit{$f}++;
      my $feature = $self->_fetch($f) or next;
      next if $strand && $feature->strand != $strand;

      if ($range_type eq 'overlaps') {
	next unless $feature->end >= $start && $feature->start <= $end;
      }
      elsif ($range_type eq 'contains') {
	next unless $feature->start >= $start && $feature->end <= $end;
      }
      elsif ($range_type eq 'contained_in') {
	next unless $feature->start <= $start && $feature->end >= $end;
      }

      push @results,$f;
    }
  }
  $self->update_filter($filter,\@results);
}


sub filter_by_name {
  my ($self, $name, $allow_aliases, $filter) = @_;

  my $index = $self->{_index}{name};

  my @names_to_fetch;
  if (my $regexp = $self->glob_match($name)) {
    @names_to_fetch = grep {/^$regexp$/i} keys %{$index};
  } else {
    @names_to_fetch = lc $name;
  }

  my @results;
  for my $n (@names_to_fetch) {
    if ($allow_aliases) {
      push @results,keys %{$index->{$n}};
    } else {
      push @results,grep {$index->{$n}{$_} == 1} keys %{$index->{$n}};
    }
  }
  $self->update_filter($filter,\@results);
}

sub glob_match {
  my ($self, $term) = @_;
  return unless $term =~ /(?:^|[^\\])[*?]/;
  $term =~ s/(^|[^\\])([+\[\]^{}\$|\(\).])/$1\\$2/g;
  $term =~ s/(^|[^\\])\*/$1.*/g;
  $term =~ s/(^|[^\\])\?/$1./g;
  return $term;
}


sub update_filter {
  my ($self, $filter, $results) = @_;
  return unless @$results;

  if (%$filter) {
    my @filtered = grep {$filter->{$_}} @$results;
    %$filter     = map {$_=>1} @filtered;
  } else {
    %$filter     = map {$_=>1} @$results;
  }

}

sub _search_attributes {
  my ($self, $search_string, $attribute_array, $limit) = @_;

  $search_string =~ tr/*?//d;

  my @words = map {quotemeta($_)} $search_string =~ /(\w+)/g;
  my $search = join '|',@words;

  my (%results,%notes);

  my $index  = $self->{_index}{attribute};
  for my $tag (@$attribute_array) {
    my $attributes = $index->{lc $tag};
    for my $value (keys %{$attributes}) {
      next unless $value =~ /$search/i;
      my @ids = keys %{$attributes->{$value}};
      for my $w (@words) {
	my @hits = $value =~ /($w)/ig or next;
	$results{$_} += @hits foreach @ids;
      }
      $notes{$_} .= "$value " foreach @ids;
    }
  }

  my @results;
  for my $id (keys %results) {
    my $hits = $results{$id};
    my $note = $notes{$id};
    $note =~ s/\s+$//;
    my $relevance = 10 * $hits;
    my $feature   = $self->fetch($id);
    my $name      = $feature->display_name or next;
    my $type      = $feature->type;
    push @results,[$name,$note,$relevance,$type,$id];
  }

  return @results;
}

=head2 types

 Title   : types
 Usage   : @type_list = $db->types
 Function: Get all the types in the database
 Returns : array of Bio::DB::GFF::Typename objects (arrayref in scalar context)
 Args    : none
 Status  : public

=cut

sub types {
  my $self = shift;
  eval "require Bio::DB::GFF::Typename" unless Bio::DB::GFF::Typename->can('new');
  my @types;
  for my $primary_tag ( keys %{$$self{_index}{type}} ) {
    for my $source_tag ( keys %{$$self{_index}{type}{$primary_tag}} ) {
	push @types, Bio::DB::GFF::Typename->new($primary_tag,$source_tag);
    }
  }
  return @types;
}

# this is ugly
sub _insert_sequence {
  my ($self, $seqid, $seq, $offset) = @_;
  my $dna_fh = $self->private_fasta_file or return;
  if ($offset == 0) { # start of the sequence
    print $dna_fh ">$seqid\n";
  }
  print $dna_fh $seq,"\n";
}

sub _fetch_sequence {
  my ($self, $seqid, $start, $end) = @_;
  my $db = $self->{fasta_db} or return;
  $db->seq($seqid,$start,$end);
}

sub private_fasta_file {
  my $self = shift;
  return $self->{fasta_fh} if exists $self->{fasta_fh};
  my $dir = tempdir (CLEANUP => 1);
  $self->{fasta_file}   = "$dir/sequence.$$.fasta";
  return $self->{fasta_fh} = IO::File->new($self->{fasta_file},">");
}

# summary support
sub coverage_array {
    my $self = shift;

    my ($seq_name,$start,$end,$types,$bins) = 
	rearrange([['SEQID','SEQ_ID','REF'],'START',['STOP','END'],
		   ['TYPES','TYPE','PRIMARY_TAG'],'BINS'],@_);

    my @features = $self->_features(-seq_id=> $seq_name,
				    -start => $start,
				    -end   => $end,
				    -types => $types);

    my $binsize = ($end-$start+1)/$bins;
    my $report_tag;
    my @coverage_array = (0) x $bins;
    
    for my $f (@features) {
	$report_tag ||= $f->primary_tag;
	my $fs        = $f->start;
	my $fe        = $f->end;
	my $start_bin = int(($fs-$start)/$binsize);
	my $end_bin   = int(($fe-$start)/$binsize);
	$start_bin    = 0       if $start_bin < 0;
	$end_bin      = $bins-1 if $end_bin  >= $bins;
	$coverage_array[$_]++ for ($start_bin..$end_bin);
    }
    return wantarray ? (\@coverage_array,$report_tag) : \@coverage_array;
}

sub _seq_ids {
    my $self = shift;

    if (my $fa = $self->{fasta_db}) {
	if (my @s = eval {$fa->ids}) {
	    return @s;
	}
    } 
    
    my $l    = $self->{_index}{location} or return;
    return keys %$l;
}

package Bio::DB::SeqFeature::Store::memory::Iterator;

sub new {
  my ($class, $store, $ids) = @_;
  return bless {store => $store,
		ids   => $ids},ref($class) || $class;
}

sub next_seq {
  my $self  = shift;
  my $store = $self->{store} or return;
  my $id    = shift @{$self->{ids}};
  defined $id or return;
  return $store->fetch($id);
}

1;

__END__

=head1 BUGS

This is an early version, so there are certainly some bugs. Please
use the BioPerl bug tracking system to report bugs.

=head1 SEE ALSO

L<bioperl>,
L<Bio::DB::SeqFeature>,
L<Bio::DB::SeqFeature::Store>,
L<Bio::DB::SeqFeature::GFF3Loader>,
L<Bio::DB::SeqFeature::Segment>,
L<Bio::DB::SeqFeature::Store::berkeleydb>,
L<Bio::DB::SeqFeature::Store::DBI::mysql>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2006 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

