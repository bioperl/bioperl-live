package Bio::DB::GFF::Adaptor::memory;

=head1 NAME

Bio::DB::GFF::Adaptor::dbi::mysql -- Database adaptor for a specific mysql schema

=head1 SYNOPSIS

  use Bio::DB::GFF;
  my $db = Bio::DB::GFF->new(-adaptor=> 'memory',
                             -gff    => 'my_features.gff',
                             -fasta  => 'my_dna.fa'
                            );

See L<Bio::DB::GFF> for other methods.

=head1 DESCRIPTION

This adaptor implements an in-memory version of Bio::DB::GFF.  It can be used to
store and retrieve SHORT GFF files. It inherits from Bio::DB::GFF.

=head1 CONSTRUCTOR

Use Bio::DB::GFF-E<gt>new() to construct new instances of this class.
Three named arguments are recommended:

   Argument         Description

   -adaptor         Set to "memory" to create an instance of this class.
   -gff             Read the indicated file or directory of .gff file.
   -fasta           Read the indicated file or directory of fasta files.
   -dir             Indicates a directory containing .gff and .fa files

If you use the -dsn option and the indicated directory is writable by
the current process, then this library will create a FASTA file index
that greatly diminishes the memory usage of this module.

=head1 METHODS

See L<Bio::DB::GFF> for inherited methods.

=head1 BUGS

none ;-)

=head1 SEE ALSO

L<Bio::DB::GFF>, L<bioperl>

=head1 AUTHOR

Shuly Avraham E<lt>avraham@cshl.orgE<gt>.

Copyright (c) 2002 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

use strict;
# $Id$
# AUTHOR: Shulamit Avraham
# This module needs to be cleaned up and documented

# Bio::DB::GFF::Adaptor::memory --  in-memory db adaptor
# implements the low level handling of data which stored in memory.
# This adaptor implements a specific in memory schema that is compatible with Bio::DB::GFF.
# Inherits from Bio::DB::GFF.


#use lib './blib/lib';
#use lib '/u/swiss/shuly/bioperl-live';
# use lib '/a/swiss/export/home/shuly/bioperl-live';
use Bio::DB::GFF;
use Bio::DB::GFF::Util::Rearrange; # for rearrange()
use Bio::DB::GFF::Adaptor::memory_iterator;
use File::Basename 'dirname';

use vars qw(@ISA $VERSION);

use constant MAX_SEGMENT => 100_000_000;  # the largest a segment can get

@ISA =  qw(Bio::DB::GFF);
$VERSION = 1.00;

sub new {
  my $class = shift ;
  my ($file,$fasta,$dbdir) = rearrange([
					[qw(GFF FILE)],
					'FASTA',
					[qw(DSN DB DIR DIRECTORY)],
				],@_);

  # fill in object
  my $self = bless{ data => [] },$class;
  $file  ||= $dbdir;
  $fasta ||= $dbdir;
  $self->load_gff($file)             if $file;
  $self->load_or_store_fasta($fasta) if $fasta;
  return $self;
}

sub load_or_store_fasta {
  my $self  = shift;
  my $fasta = shift;
  if ((-f $fasta && -w dirname($fasta))
      or
      (-d $fasta && -w $fasta)) {
    require Bio::DB::Fasta;
    my $dna_db = eval {Bio::DB::Fasta->new($fasta)} 
      or warn "No sequence available. Use -gff instead of -dir if you wish to load features without sequence.\n";
    $dna_db && $self->dna_db($dna_db);
  } else {
    $self->load_fasta($fasta);
  }
}

sub dna_db {
  my $self = shift;
  my $d    = $self->{dna_db};
  $self->{dna_db} = shift if @_;
  $d;
}

sub insert_sequence {
  my $self = shift;
  my($id,$offset,$seq) = @_;
  $self->{dna}{$id} .= $seq;
}

# low-level fetch of a DNA substring given its
# name, class and the desired range.
sub get_dna {
  my $self = shift;
  my ($id,$start,$stop,$class) = @_;
  if (my $dna_db = $self->dna_db) {
    return $dna_db->seq($id,$start=>$stop);
  }

  return $self->{dna}{$id} if !defined $start || !defined $stop;
  $start = 1 if !defined $start;

  my $reversed = 0;
  if ($start > $stop) {
    $reversed++;
    ($start,$stop) = ($stop,$start);
  }
  my $dna = substr($self->{dna}{$id},$start-1,$stop-$start+1);
  if ($reversed) {
    $dna =~ tr/gatcGATC/ctagCTAG/;
    $dna = reverse $dna;
  }

  $dna;
}

# this method loads the feature as a hash into memory -
# keeps an array of features-hashes as an in-memory db
sub load_gff_line {
  my $self = shift;
  my $feature_hash  = shift;
  $feature_hash->{strand} = '' if $feature_hash->{strand} && $feature_hash->{strand} eq '.';
  $feature_hash->{phase} = ''  if $feature_hash->{phase}  && $feature_hash->{phase} eq '.';
  push @{$self->{data}},$feature_hash;
}

# given sequence name, return (reference,start,stop,strand)
sub get_abscoords {
  my $self = shift;
  my ($name,$class,$refseq) = @_;
  my %refs;
  my $regexp;

  if ($name =~ /[*?]/) {  # uh oh regexp time
    $name =~ quotemeta($name);
    $name =~ s/\\\*/.*/g;
    $name =~ s/\\\?/.?/g;
    $regexp++;
  }

  # Find all features that have the requested name and class.
  # Sort them by reference point.
  for my $feature (@{$self->{data}}) {

    my $no_match_class_name;
    my $empty_class_name;
    if (defined $feature->{gname} and defined $feature->{gclass}){
      my $matches = $feature->{gclass} eq $class
	&& ($regexp ? $feature->{gname} =~ /$name/i : $feature->{gname} eq $name);
      $no_match_class_name = !$matches;  # to accomodate Shuly's interesting logic
    }

    else{
      $empty_class_name = 1;
    }

    if ($no_match_class_name || $empty_class_name){

      my $feature_attributes = $feature->{attributes};
      my $attributes = {Alias => $name};
      if (!_matching_attributes($feature_attributes,$attributes)){
	next;
      }
    }
    push @{$refs{$feature->{ref}}},$feature;
  }

  # find out how many reference points we recovered

  if (! %refs) {
    $self->error("$name not found in database");
    return;
  } elsif (keys %refs > 1) {
    $self->error("$name has more than one reference sequence in database");
    return;
  }

  # compute min and max
  my ($ref) = keys %refs;
  my @found = @{$refs{$ref}};
  my ($strand,$start,$stop);

  my @found_segments;
  foreach my $ref (keys %refs) {
    next if defined($refseq) and $ref ne $refseq;
    my @found = @{$refs{$ref}};
    my ($strand,$start,$stop);
    foreach (@found) {
      $strand ||= $_->{strand};
      $strand = '+' if $strand && $strand eq '.'; 
      $start  = $_->{start} if !defined($start) || $start > $_->{start};
      $stop   = $_->{stop}  if !defined($stop)  || $stop  < $_->{stop};
    }
    push @found_segments,[$ref,$class,$start,$stop,$strand];

  }
  return \@found_segments;
}

sub search_notes {
  my $self = shift;
  my ($search_string,$limit) = @_;
  my @results;
  my @words = map {quotemeta($_)} $search_string =~ /(\w+)/g;

  for my $feature (@{$self->{data}}) {
    next unless defined $feature->{gclass} && defined $feature->{gname}; # ignore NULL objects
    next unless $feature->{attributes};
    my @attributes = @{$feature->{attributes}};
    my @values     = map {$_->[1]} @attributes;
    my $value      = "@values";
    my $matches    = 0;
    my $note;
    for my $w (@words) {
      my @hits = $value =~ /($w)/g;
      $note ||= $value if @hits;
      $matches += @hits;
    }
    next unless $matches;

    my $relevance = 10 * $matches;
    my $featname = Bio::DB::GFF::Featname->new($feature->{gclass}=>$feature->{gname});
    push @results,[$featname,$note,$relevance];
    last if @results >= $limit;
  }
  @results;
}

sub _delete_features {
  my $self        = shift;
  my @feature_ids = sort {$b<=>$a} @_;
  my $removed = 0;
  foreach (@feature_ids) {
    next unless $_ >= 0 && $_ < @{$self->{data}};
    $removed += defined splice(@{$self->{data}},$_,1);
  }
  $removed;
}

sub _delete {
  my $self = shift;
    my $delete_spec = shift;
  my $ranges      = $delete_spec->{segments} || [];
  my $types       = $delete_spec->{types}    || [];
  my $force       = $delete_spec->{force};
  my $range_type  = $delete_spec->{range_type};

  my $deleted = 0;
  if (@$ranges) {
    my @args = @$types ? (-type=>$types) : ();
    push @args,(-range_type => $range_type);
    my %ids_to_remove = map {$_->id => 1} map {$_->features(@args)} @$ranges;
    $deleted = $self->delete_features(keys %ids_to_remove);
  } elsif (@$types) {
    my %ids_to_remove = map {$_->id => 1} $self->features(-type=>$types);
    $deleted = $self->delete_features(keys %ids_to_remove);
  } else {
    $self->throw("This operation would delete all feature data and -force not specified")
      unless $force;
    $deleted = @{$self->{data}};
    @{$self->{data}} = ();
  }
  $deleted;
}

# attributes -

# Some GFF version 2 files use the groups column to store a series of
# attribute/value pairs.  In this interpretation of GFF, the first such
# pair is treated as the primary group for the feature; subsequent pairs
# are treated as attributes.  Two attributes have special meaning:
# "Note" is for backward compatibility and is used for unstructured text
# remarks.  "Alias" is considered as a synonym for the feature name.
# If no name is provided, then attributes() returns a flattened hash, of
# attribute=>value pairs.

sub do_attributes{
  my $self = shift;
  my ($feature_id,$tag) = @_;
  my $attr ;

  my $feature = ${$self->{data}}[$feature_id];
  my @result;
  for my $attr (@{$feature->{attributes}}) {
    my ($attr_name,$attr_value) = @$attr ;
    if (defined($tag) && $attr_name eq $tag){push @result,$attr_value;}
    elsif (!defined($tag)) {push @result,($attr_name,$attr_value);}
  }
  return @result;
}


#sub get_feature_by_attribute{
sub _feature_by_attribute{
  my $self = shift;
  my ($attributes,$callback) = @_;
  $callback || $self->throw('must provide a callback argument');
  my $count = 0;
  my $feature_id = -1;
  my $feature_group_id = undef;

  for my $feature (@{$self->{data}}) {

    $feature_id++;
    for my $attr (@{$feature->{attributes}}) {
      my ($attr_name,$attr_value) = @$attr ;
      #there could be more than one set of attributes......
      foreach (keys %$attributes) {
	if ($_ eq $attr_name && $attributes->{$_} eq $attr_value){
           $callback->($feature->{ref},
	        $feature->{start},
	        $feature->{stop},
	        $feature->{source},
	        $feature->{method},
	        $feature->{score},
	        $feature->{strand},
	        $feature->{phase},
	        $feature->{gclass},
	        $feature->{gname},
		$feature->{tstart},
		$feature->{tstop},
	        $feature_id,
		$feature_group_id);
	   $count++;
        }
      }
    }
  }

}



# This is the low-level method that is called to retrieve GFF lines from
# the database.  It is responsible for retrieving features that satisfy
# range and feature type criteria, and passing the GFF fields to a
# callback subroutine.

sub get_features{
  my $self = shift;
  my $count = 0;
  my ($search,$options,$callback) = @_;
  my $data = \@{$self->{data}};

  my $found_features;

  $found_features = _get_features_by_search_options($data,$search,$options);

  # only true if the sort by group option was specified
  @{$found_features} = sort {"$a->{gclass}:$a->{gname}" cmp "$b->{gclass}:$b->{gname}"} 
    @{$found_features} if $options->{sort_by_group} ;

  for my $feature (@{$found_features}) {  # only true if the sort by group option was specified
    $count++;
    $callback->(
		@{$feature}{qw(ref start stop source method score strand phase gclass gname tstart tstop feature_id feature_group_id)}
	       );
  }

  return $count;
}


# Low level implementation of fetching a named feature.
# GFF annotations are named using the group class and name fields.
# May return zero, one, or several Bio::DB::GFF::Feature objects.

=head2 _feature_by_name

 Title   : _feature_by_name
 Usage   : $db->get_features_by_name($name,$class,$callback)
 Function: get a list of features by name and class
 Returns : count of number of features retrieved
 Args    : name of feature, class of feature, and a callback
 Status  : protected

This method is used internally.  The callback arguments are those used
by make_feature().

=cut

sub _feature_by_name {
  my $self = shift;
  my ($class,$name,$location,$callback) = @_;
  $callback || $self->throw('must provide a callback argument');
  my $count = 0;
  my $id    = -1;
  my $regexp;

  if ($name =~ /[*?]/) {  # uh oh regexp time
    $name =~ quotemeta($name);
    $name =~ s/\\\*/.*/g;
    $name =~ s/\\\?/.?/g;
    $regexp++;
  }


  for my $feature (@{$self->{data}}) {
    $id++;
    next unless ($regexp && $feature->{gname} =~ /$name/i) || $feature->{gname}  eq $name;
    next unless $feature->{gclass} eq $class;
    if ($location) {
      next if $location->[0] ne $feature->{ref};
      next if $location->[1] && $location->[1] > $feature->{stop};
      next if $location->[2] && $location->[2] < $feature->{start};
    }
    $count++;
    $callback->(@{$feature}{qw(
			       ref
			       start
			       stop
			       source
			       method
			       score
			       strand
			       phase
			       gclass
			       gname
			       tstart
			       tstop
			      )},$id,0
	       );
  }
  return $count;
}

# Low level implementation of fetching a feature by it's id. 
# The id of the feature as implemented in the in-memory db, is the location of the 
# feature in the features hash array.
sub _feature_by_id{
  my $self = shift;
  my ($ids,$type,$callback) = @_;
  $callback || $self->throw('must provide a callback argument');

  my $feature_group_id = undef;

  my $count = 0;
  if ($type eq 'feature'){
    for my $feature_id (@$ids){
       my $feature = ${$self->{data}}[$feature_id];
       
       $callback->($feature->{ref},
	        $feature->{start},
	        $feature->{stop},
	        $feature->{source},
	        $feature->{method},
	        $feature->{score},
	        $feature->{strand},
	        $feature->{phase},
	        $feature->{gclass},
	        $feature->{gname},
		$feature->{tstart},
		$feature->{tstop},
	        $feature_id,
		$feature_group_id);
	   $count++;			
    
    }
  }
}


# This method is similar to get_features(), except that it returns an
# iterator across the query.  
# See Bio::DB::GFF::Adaptor::memory_iterator.

sub get_features_iterator {
  my $self = shift;
  my ($search,$options,$callback) = @_;
  $callback || $self->throw('must provide a callback argument');

  my $data = \@{$self->{data}};
  my $results = _get_features_by_search_options($data,$search,$options);
  my $results_array = _convert_feature_hash_to_array($results);

  return Bio::DB::GFF::Adaptor::memory_iterator->new($results_array,$callback);
}




# This method is responsible for fetching the list of feature type names.
# The query may be limited to a particular range, in
# which case the range is indicated by a landmark sequence name and
# class and its subrange, if any.  These arguments may be undef if it is
# desired to retrieve all feature types.

# If the count flag is false, the method returns a simple list of
# Bio::DB::GFF::Typename objects.  If $count is true, the method returns
# a list of $name=>$count pairs, where $count indicates the number of
# times this feature occurs in the range.

sub get_types {
  my $self = shift;
  my ($srcseq,$class,$start,$stop,$want_count,$typelist) = @_;
  my(%result,%obj);

  for my $feature (@{$self->{data}}) {
    my $feature_start = $feature->{start};
    my $feature_stop  = $feature->{stop};
    my $feature_ref   = $feature->{ref};
    my $feature_class = $feature->{class};
    my $feature_method = $feature->{method};
    my $feature_source = $feature->{source};
   
    if (defined $srcseq){
      next unless $feature_ref eq $srcseq ;
    }
    
    if (defined $class){ 
      next unless $feature_class eq $class ;
    }
    
     # the requested range should OVERLAP the retrieved features
     if (defined $start or defined $stop) {
      $start = 1           unless defined $start;
      $stop  = MAX_SEGMENT unless defined $stop;
      next unless $feature_stop >= $start && $feature_start <= $stop;
    }
    
    if (defined $typelist && @$typelist){
      next unless _matching_typelist($feature_method,$feature_source,$typelist);
    }

    my $type = Bio::DB::GFF::Typename->new($feature_method,$feature_source);
    $result{$type}++;
    $obj{$type} = $type;

  }   #end features loop
  
  return $want_count ? %result : values %obj;
 
}




# Internal method that performs a search on the features array, 
# sequentialy retrieves the features, and performs a check on each feature
# according to the search options.
sub _get_features_by_search_options{
  my $count = 0;
  my ($data,$search,$options) = @_;
  my ($rangetype,$refseq,$class,$start,$stop,$types,$sparse,$order_by_group,$attributes) = 
    (@{$search}{qw(rangetype refseq refclass start stop types)},
    @{$options}{qw(sparse sort_by_group ATTRIBUTES)}) ;
					       
  my @found_features;

  my $feature_id = -1 ;
  my $feature_group_id = undef;

  for my $feature (@{$data}) {

    $feature_id++;
    
    my $feature_start = $feature->{start};
    my $feature_stop  = $feature->{stop};
    my $feature_ref   = $feature->{ref};
    
    if (defined $refseq){
      next unless $feature_ref eq $refseq;
    }

     if (defined $start or defined $stop) {
      $start = 0               unless defined($start);
      $stop  = MAX_SEGMENT     unless defined($stop);
    
      if ($rangetype eq 'overlaps') {
	next unless $feature_stop >= $start && $feature_start <= $stop;
      } elsif ($rangetype eq 'contains') {
	next unless $feature_start >= $start && $feature_stop <= $stop;
      } elsif ($rangetype eq 'contained_in') {
	next unless $feature_start <= $start && $feature_stop >= $stop;
      } else {
	next unless $feature_start == $start && $feature_stop == $stop;
      }

    }
    
    my $feature_source = $feature->{source};
    my $feature_method = $feature->{method};

    if (defined $types && @$types){
      next unless _matching_typelist($feature_method,$feature_source,$types);
    } 

    my $feature_attributes = $feature->{attributes};
    if (defined $attributes){
      next unless _matching_attributes($feature_attributes,$attributes);
    } 
    
    # if we get here, then we have a feature that meets the criteria.
    # Then we just push onto an array
    # of found features and continue. 
   
    my $found_feature = $feature ;
    $found_feature->{feature_id} = $feature_id;
    $found_feature->{group_id} = $feature_group_id;
    push @found_features,$found_feature;
   
  }

  return \@found_features; 
}





# this subroutine is needed for convertion of the feature from hash to array in order to 
# pass it to the callback subroutine
sub _convert_feature_hash_to_array{
  my @features_hash_array = @_;

  use constant FREF    => 0;
  use constant FSTART  => 1;
  use constant FSTOP   => 2;
  use constant FSOURCE => 3;
  use constant FMETHOD => 4;
  use constant FSCORE  => 5;
  use constant FSTRAND => 6;
  use constant FPHASE  => 7;
  use constant GCLASS  => 8;
  use constant GNAME   => 9;
  use constant TSTART  => 10;
  use constant TSTOP   => 11;
  use constant FID     => 12;
  use constant GID     => 13;

  my @features_array_array;
  my $feature_count = 0;
   
  for my $feature_hash (@{$features_hash_array[0]}){
    my @feature_array;

    $feature_array[FREF]    = $feature_hash->{ref};
    $feature_array[FSTART]  = $feature_hash->{start};
    $feature_array[FSTOP]   = $feature_hash->{stop};  
    $feature_array[FSOURCE] = $feature_hash->{source};
    $feature_array[FMETHOD] = $feature_hash->{method};
    $feature_array[FSCORE]  = $feature_hash->{score};
    $feature_array[FSTRAND] = $feature_hash->{strand};  
    $feature_array[FPHASE ] = $feature_hash->{phase};
    $feature_array[GCLASS]  = $feature_hash->{gclass};  
    $feature_array[GNAME]   = $feature_hash->{gname};
    $feature_array[TSTART]  = $feature_hash->{tstart};
    $feature_array[TSTOP]   = $feature_hash->{tstop};
    $feature_array[FID]     = $feature_hash->{feature_id};  
    $feature_array[GID]     = $feature_hash->{group_id};

    $features_array_array[$feature_count] = \@feature_array;
    $feature_count++;
  }
  return \@features_array_array;
}

sub _matching_typelist{ 
  my ($feature_method,$feature_source,$typelist) = @_; 
  foreach (@$typelist) {
	 my ($search_method,$search_source) = @$_;
	 next if $search_method ne $feature_method;
	 next if defined($search_source) && $search_source ne $feature_source;
	 return 1;
  }
  return 0;
}

sub _matching_attributes{
  my ($feature_attributes,$attributes) = @_ ;
  foreach (keys %$attributes) {
    return 0 if !_match_all_attr_in_feature($_,$attributes->{$_},$feature_attributes)
   
  }
  return 1;
}

sub _match_all_attr_in_feature{
  my ($attr_name,$attr_value,$feature_attributes) = @_;
  for my $attr (@$feature_attributes) {
      my ($feature_attr_name,$feature_attr_value) = @$attr ;
      next if ($attr_name ne $feature_attr_name || $attr_value ne $feature_attr_value);
      return 1;
  }
  return 0;
}


sub do_initialize { 1; }
sub setup_load { }
sub finish_load { 1; }
sub get_feature_by_group_id{ 1; }

1;

