package Bio::SeqFeature::Lite;

=head1 NAME

Bio::SeqFeature::Lite - Lightweight Bio::SeqFeatureI class

=head1 SYNOPSIS

 # create a simple feature with no internal structure
 $f = Bio::SeqFeature::Lite->new(-start => 1000,
                                  -stop  => 2000,
                                  -type  => 'transcript',
                                  -name  => 'alpha-1 antitrypsin',
                                  -desc  => 'an enzyme inhibitor',
                                 );

 # create a feature composed of multiple segments, all of type "similarity"
 $f = Bio::SeqFeature::Lite->new(-segments => [[1000,1100],[1500,1550],[1800,2000]],
                                  -name     => 'ABC-3',
                                  -type     => 'gapped_alignment',
                                  -subtype  => 'similarity');

 # build up a gene exon by exon
 $e1 = Bio::SeqFeature::Lite->new(-start=>1,-stop=>100,-type=>'exon');
 $e2 = Bio::SeqFeature::Lite->new(-start=>150,-stop=>200,-type=>'exon');
 $e3 = Bio::SeqFeature::Lite->new(-start=>300,-stop=>500,-type=>'exon');
 $f  = Bio::SeqFeature::Lite->new(-segments=>[$e1,$e2,$e3],-type=>'gene');

=head1 DESCRIPTION

This is a simple Bio::SeqFeatureI-compliant object that is compatible
with Bio::Graphics::Panel.  With it you can create lightweight feature
objects for drawing.

All methods are as described in L<Bio::SeqFeatureI> with the following additions:

=head2 The new() Constructor

 $feature = Bio::SeqFeature::Lite->new(@args);

This method creates a new feature object.  You can create a simple
feature that contains no subfeatures, or a hierarchically nested object.

Arguments are as follows:

  -seq_id      the reference sequence
  -start       the start position of the feature
  -end         the stop position of the feature
  -stop        an alias for end
  -name        the feature name (returned by seqname())
  -type        the feature type (returned by primary_tag())
  -primary_tag the same as -type
  -source      the source tag
  -score       the feature score (for GFF compatibility)
  -desc        a description of the feature
  -segments    a list of subfeatures (see below)
  -subtype     the type to use when creating subfeatures
  -strand      the strand of the feature (one of -1, 0 or +1)
  -phase       the phase of the feature (0..2)
  -seq         a dna or protein sequence string to attach to feature
  -id          an alias for -name
  -seqname     an alias for -name
  -display_id  an alias for -name
  -display_name an alias for -name  (do you get the idea the API has changed?)
  -primary_id  unique database ID
  -url         a URL to link to when rendered with Bio::Graphics
  -attributes  a hashref of tag value attributes, in which the key is the tag
               and the value is an array reference of values
  -factory     a reference to a feature factory, used for compatibility with
               more obscure parts of Bio::DB::GFF

The subfeatures passed in -segments may be an array of
Bio::SeqFeature::Lite objects, or an array of [$start,$stop]
pairs. Each pair should be a two-element array reference.  In the
latter case, the feature type passed in -subtype will be used when
creating the subfeatures.

If no feature type is passed, then it defaults to "feature".

=head2 Non-SeqFeatureI methods

A number of new methods are provided for compatibility with
Ace::Sequence, which has a slightly different API from SeqFeatureI:

=over 4

=item url()

Get/set the URL that the graphical rendering of this feature will link to.

=item add_segment(@segments)

Add one or more segments (a subfeature).  Segments can either be
Feature objects, or [start,stop] arrays, as in the -segments argument
to new().  The feature endpoints are automatically adjusted.

=item segments()

An alias for sub_SeqFeature().

=item get_SeqFeatures()

Alias for sub_SeqFeature()

=item get_all_SeqFeatures()

Alias for sub_SeqFeature()

=item merged_segments()

Another alias for sub_SeqFeature().

=item stop()

An alias for end().

=item name()

An alias for seqname().

=item exons()

An alias for sub_SeqFeature() (you don't want to know why!)

=back

=cut

use strict;

use base qw(Bio::Root::Root Bio::SeqFeatureI Bio::LocationI Bio::SeqI);

*stop        = \&end;
*info        = \&name;
*seqname     = \&name;
*exons       = *sub_SeqFeature = *merged_segments = \&segments;
*get_all_SeqFeatures = *get_SeqFeatures = \&segments;
*method         = \&primary_tag;
*source         = \&source_tag;
*get_tag_values = \&each_tag_value;
*add_SeqFeature = \&add_segment;
*get_all_tags   = \&all_tags;
*abs_ref        = \&ref;

# implement Bio::SeqI and FeatureHolderI interface

sub primary_seq { return $_[0] }
sub annotation { 
    my ($obj,$value) = @_;
    if( defined $value ) {
	$obj->throw("object of class ".ref($value)." does not implement ".
		    "Bio::AnnotationCollectionI. Too bad.")
	    unless $value->isa("Bio::AnnotationCollectionI");
	$obj->{'_annotation'} = $value;
    } elsif( ! defined $obj->{'_annotation'}) {
	$obj->{'_annotation'} = Bio::Annotation::Collection->new();
    }
    return $obj->{'_annotation'};
}
sub species {
    my ($self, $species) = @_;
    if ($species) {
        $self->{'species'} = $species;
    } else {
        return $self->{'species'};
    }
}
sub is_remote { return }
sub feature_count { return scalar @{shift->{segments} || []} }

sub target { return; }
sub hit    { shift->target }

sub type {
  my $self = shift;
  my $method = $self->primary_tag;
  my $source = $self->source_tag;
  return $source ne '' ? "$method:$source" : $method;
}

# usage:
# Bio::SeqFeature::Lite->new(
#                         -start => 1,
#                         -end   => 100,
#                         -name  => 'fred feature',
#                         -strand => +1);
#
# Alternatively, use -segments => [ [start,stop],[start,stop]...]
# to create a multisegmented feature.
sub new {
  my $class= shift;
  $class = ref($class) if ref $class;
  my %arg = @_;

  my $self = bless {},$class;

  $arg{-strand} ||= 0;
  if ($arg{-strand} =~ /^[\+\-\.]$/){
	($arg{-strand} eq "+") && ($self->{strand} = '1');
	($arg{-strand} eq "-") && ($self->{strand} = '-1');
	($arg{-strand} eq ".") && ($self->{strand} = '0');
  } else {
	  $self->{strand}  = $arg{-strand} ? ($arg{-strand} >= 0 ? +1 : -1) : 0;
  }
  $self->{name}    = $arg{-name}   || $arg{-seqname} || $arg{-display_id} 
    || $arg{-display_name} || $arg{-id};
  $self->{type}    = $arg{-type}   || $arg{-primary_tag} || 'feature';
  $self->{subtype} = $arg{-subtype} if exists $arg{-subtype};
  $self->{source}  = $arg{-source} || $arg{-source_tag} || '';
  $self->{score}   = $arg{-score}   if exists $arg{-score};
  $self->{start}   = $arg{-start};
  $self->{stop}    = exists $arg{-end} ? $arg{-end} : $arg{-stop};
  $self->{ref}     = $arg{-seq_id} || $arg{-ref};
  $self->{attributes}     = $arg{-attributes} || $arg{-tag};
  for my $option (qw(class url seq phase desc primary_id)) {
    $self->{$option} = $arg{"-$option"} if exists $arg{"-$option"};
  }

  # is_circular is needed for Bio::PrimarySeqI compliance
  $self->{is_circular} = $arg{-is_circular} || 0;

  # fix start, stop
  if (defined $self->{stop} && defined $self->{start}
      && $self->{stop} < $self->{start}) {
    @{$self}{'start','stop'} = @{$self}{'stop','start'};
    $self->{strand} *= -1;
  }

  my @segments;
  if (my $s = $arg{-segments}) {
    # NB: when $self ISA Bio::DB::SeqFeature the following invokes
    # Bio::DB::SeqFeature::add_segment and not
    # Bio::DB::SeqFeature::add_segment (as might be expected?)
    $self->add_segment(@$s);
  }

  $self;
}

sub add_segment {
  my $self        = shift;
  my $type = $self->{subtype} || $self->{type};
  $self->{segments} ||= [];
  my $ref   = $self->seq_id;
  my $name  = $self->name;
  my $class = $self->class;
  my $source_tag = $self->source_tag;

  my $min_start = $self->start ||  999_999_999_999;
  my $max_stop  = $self->end   || -999_999_999_999;

  my @segments = @{$self->{segments}};

  for my $seg (@_) {
    if (ref($seg) eq 'ARRAY') {
      my ($start,$stop) = @{$seg};
      next unless defined $start && defined $stop;  # fixes an obscure bug somewhere above us
      my $strand = $self->{strand};

      if ($start > $stop) {
	($start,$stop) = ($stop,$start);
	$strand = -1;
      }

      push @segments,$self->new(-start  => $start,
				-stop   => $stop,
				-strand => $strand,
				-ref    => $ref,
				-type   => $type,
			        -name   => $name,
			        -class  => $class,
				-phase  => $self->{phase},
				-score  => $self->{score},
				-source_tag  => $source_tag,
				-attributes  => $self->{attributes},
			       );
      $min_start = $start if $start < $min_start;
      $max_stop  = $stop  if $stop  > $max_stop;

    } elsif (ref $seg) {
      push @segments,$seg;

      $min_start = $seg->start if ($seg->start && $seg->start < $min_start);
      $max_stop  = $seg->end   if ($seg->end && $seg->end > $max_stop);
    }
  }
  if (@segments) {
    local $^W = 0;  # some warning of an uninitialized variable...
    $self->{segments} = \@segments;
    $self->{ref}    ||= $self->{segments}[0]->seq_id;
    $self->{start}    = $min_start;
    $self->{stop}     = $max_stop;
  }
}

sub segments {
  my $self = shift;
  my $s = $self->{segments} or return wantarray ? () : 0;
  @$s;
}
sub score    {
  my $self = shift;
  my $d = $self->{score};
  $self->{score} = shift if @_;
  $d;
}
sub primary_tag     { 
    my $self = shift;
    my $d    = $self->{type};
    $self->{type} = shift if @_;
    $d;
}
sub name            {
  my $self = shift;
  my $d    = $self->{name};
  $self->{name} = shift if @_;
  $d;
}
sub seq_id          { shift->ref(@_)         }
sub ref {
  my $self = shift;
  my $d = $self->{ref};
  $self->{ref} = shift if @_;
  $d;
}
sub start    {
  my $self = shift;
  my $d = $self->{start};
  $self->{start} = shift if @_;
  if (my $rs = $self->{refseq}) {
    my $strand = $rs->strand || 1;
    return $strand >= 0 ? ($d - $rs->start + 1) : ($rs->end - $d + 1);
  } else {
    return $d;
  }
}
sub end    {
  my $self = shift;
  my $d = $self->{stop};
  $self->{stop} = shift if @_;
  if (my $rs = $self->{refseq}) {
    my $strand = $rs->strand || 1;
    return $strand >= 0 ? ($d - $rs->start + 1) : ($rs->end - $d + 1);
  }
  $d;
}
sub strand {
  my $self = shift;
  my $d = $self->{strand};
  $self->{strand} = shift if @_;
  if (my $rs = $self->{refseq}) {
    my $rstrand = $rs->strand;
    return  0 unless $d;
    return  1 if $rstrand == $d;
    return -1 if $rstrand != $d;
  }
  $d;
}

# this does nothing, but it is here for compatibility reasons
sub absolute {
  my $self = shift;
  my $d = $self->{absolute};
  $self->{absolute} = shift if @_;
  $d;
}

sub abs_start {
  my $self = shift;
  local $self->{refseq} = undef;
  $self->start(@_);
}
sub abs_end {
  my $self = shift;
  local $self->{refseq} = undef;
  $self->end(@_);
}
sub abs_strand {
  my $self = shift;
  local $self->{refseq} = undef;
  $self->strand(@_);
}

sub length {
  my $self = shift;
  return $self->end - $self->start + 1;
}

#is_circular is needed for Bio::PrimarySeqI
sub is_circular {
  my $self = shift;
  my $d = $self->{is_circular};
  $self->{is_circular} = shift if @_;
  $d;
}


sub seq {
  my $self = shift;
  my $seq =  exists $self->{seq} ? $self->{seq} : '';
  return $seq;
}

sub dna {
  my $seq = shift->seq;
  $seq    = $seq->seq if CORE::ref($seq);
  return $seq;
}

=head2 display_name

 Title   : display_name
 Usage   : $id = $obj->display_name or $obj->display_name($newid);
 Function: Gets or sets the display id, also known as the common name of
           the Seq object.

           The semantics of this is that it is the most likely string
           to be used as an identifier of the sequence, and likely to
           have "human" readability.  The id is equivalent to the LOCUS
           field of the GenBank/EMBL databanks and the ID field of the
           Swissprot/sptrembl database. In fasta format, the >(\S+) is
           presumed to be the id, though some people overload the id
           to embed other information. Bioperl does not use any
           embedded information in the ID field, and people are
           encouraged to use other mechanisms (accession field for
           example, or extending the sequence object) to solve this.

           Notice that $seq->id() maps to this function, mainly for
           legacy/convenience issues.
 Returns : A string
 Args    : None or a new id


=cut

sub display_name { shift->name(@_) }

*display_id = \&display_name;

=head2 accession_number

 Title   : accession_number
 Usage   : $unique_biological_key = $obj->accession_number;
 Function: Returns the unique biological id for a sequence, commonly
           called the accession_number. For sequences from established
           databases, the implementors should try to use the correct
           accession number. Notice that primary_id() provides the
           unique id for the implemetation, allowing multiple objects
           to have the same accession number in a particular implementation.

           For sequences with no accession number, this method should return
           "unknown".
 Returns : A string
 Args    : None


=cut

sub accession_number {
    return 'unknown';
}

=head2 alphabet

 Title   : alphabet
 Usage   : if( $obj->alphabet eq 'dna' ) { /Do Something/ }
 Function: Returns the type of sequence being one of
           'dna', 'rna' or 'protein'. This is case sensitive.

           This is not called <type> because this would cause
           upgrade problems from the 0.5 and earlier Seq objects.

 Returns : a string either 'dna','rna','protein'. NB - the object must
           make a call of the type - if there is no type specified it
           has to guess.
 Args    : none
 Status  : Virtual


=cut

sub alphabet{
    return 'dna'; # no way this will be anything other than dna!
}



=head2 desc

 Title   : desc
 Usage   : $seqobj->desc($string) or $seqobj->desc()
 Function: Sets or gets the description of the sequence
 Example :
 Returns : The description
 Args    : The description or none


=cut

sub desc {
  my $self = shift;
  my ($d)    = $self->notes;
  $self->{desc} = shift if @_;
  $d;
}

sub attributes {
  my $self = shift;
  if (@_) {
    return $self->get_tag_values(@_);
  } else {
    return $self->{attributes} ? %{$self->{attributes}} : ();
  }
}

sub primary_id {
  my $self = shift;
  my $d = $self->{primary_id};
  $self->{primary_id} = shift if @_;
  return $d;
#  return $d if defined $d;
#  return (overload::StrVal($self) =~ /0x([a-f0-9]+)/)[0];
}

sub notes {
    my $self  = shift;
    my $notes = $self->{desc};
    return $notes if defined $notes;
    return $self->attributes('Note');
}

sub aliases {
    my $self  = shift;
    return $self->attributes('Alias');
}

sub low {
  my $self = shift;
  return $self->start < $self->end ? $self->start : $self->end;
}

sub high {
  my $self = shift;
  return $self->start > $self->end ? $self->start : $self->end;
}

=head2 location

 Title   : location
 Usage   : my $location = $seqfeature->location()
 Function: returns a location object suitable for identifying location
	   of feature on sequence or parent feature
 Returns : Bio::LocationI object
 Args    : none

=cut

sub location {
   my $self = shift;
   require Bio::Location::Split unless Bio::Location::Split->can('new');
   my $location;
   if (my @segments = $self->segments) {
       $location = Bio::Location::Split->new();
       foreach (@segments) {
	 $location->add_sub_Location($_);
       }
   } else {
       $location = $self;
   }
   $location;
}

sub each_Location {
  my $self = shift;
  require Bio::Location::Simple unless Bio::Location::Simple->can('new');
  if (my @segments = $self->segments) {
    return map {
	Bio::Location::Simple->new(-start  => $_->start,
				   -end    => $_->end,
				   -strand => $_->strand);
      } @segments;
  } else {
    return Bio::Location::Simple->new(-start  => $self->start,
				      -end    => $self->end,
				      -strand => $self->strand);
  }
}

=head2 location_string

 Title   : location_string
 Usage   : my $string = $seqfeature->location_string()
 Function: Returns a location string in a format recognized by gbrowse
 Returns : a string
 Args    : none

This is a convenience function used by the generic genome browser. It
returns the location of the feature and its subfeatures in the compact
form "start1..end1,start2..end2,...".  Use
$seqfeature-E<gt>location()-E<gt>toFTString() to obtain a standard
GenBank/EMBL location representation.

=cut

sub location_string {
  my $self = shift;
  my @segments = $self->segments or return $self->to_FTstring;
  join ',',map {$_->to_FTstring} @segments;
}

sub coordinate_policy {
   require Bio::Location::WidestCoordPolicy unless Bio::Location::WidestCoordPolicy->can('new');
   return Bio::Location::WidestCoordPolicy->new();
}

sub min_start { shift->low }
sub max_start { shift->low }
sub min_end   { shift->high }
sub max_end   { shift->high}
sub start_pos_type { 'EXACT' }
sub end_pos_type   { 'EXACT' }
sub to_FTstring {
  my $self = shift;
  my $low  = $self->min_start;
  my $high = $self->max_end;
  my $strand = $self->strand;
  my $str = defined $strand && $strand<0 ? "complement($low..$high)" : "$low..$high";
  if (my $id = $self->seq_id()) {
    $str = $id . ":" . $str;
  }
  $str;
}
sub phase {
    my $self = shift;
    my $d    = $self->{phase};
    $self->{phase} = shift if @_;
    $d;
}

sub class {
  my $self = shift;
  my $d = $self->{class};
  $self->{class} = shift if @_;
  return defined($d) ? $d : 'Sequence';  # acedb is still haunting me - LS
}

# set GFF dumping version
sub version {
  my $self = shift;
  my $d    = $self->{gff3_version} || 2;
  $self->{gff3_version} = shift if @_;
  $d;
}

sub gff_string {
  my $self    = shift;
  
  if ($self->version == 3) {
    return $self->gff3_string(@_);
  }
  
  my $recurse = shift;
  my $name  = $self->name;
  my $class = $self->class;
  my $group = "$class $name" if $name;
  my $strand = ('-','.','+')[$self->strand+1];
  my $string;
  $string .= join("\t",
		  $self->ref||'.',$self->source||'.',$self->method||'.',
		  $self->start||'.',$self->stop||'.',
		  defined($self->score) ? $self->score : '.',
		  $strand||'.',
		  defined($self->phase) ? $self->phase : '.',
		  $group||''
		 );
  $string .= "\n";
  if ($recurse) {
    foreach ($self->sub_SeqFeature) {
      $string .= $_->gff_string($recurse);
    }
  }
  $string;
}

# Suggested strategy for dealing with the multiple parentage issue.
# First recurse through object tree and record parent tree.
# Then recurse again, skipping objects we've seen before.
sub gff3_string {
    my ($self,$recurse,$parent_tree,$seenit,$force_id) = @_;
    $parent_tree ||= {};
    $seenit      ||= {};
    my @rsf      =   ();
    my @parent_ids;

    if ($recurse) {
	$self->_traverse($parent_tree) unless %$parent_tree;  # this will record parents of all children
	my $primary_id = defined $force_id ? $force_id : $self->_real_or_dummy_id;

	return if $seenit->{$primary_id}++;

	@rsf = $self->get_SeqFeatures;
	if (@rsf) {
	    # Detect case in which we have a split location feature. In this case we 
	    # skip to the grandchildren and trick them into thinking that our parent is theirs.
	    my %types = map {$_->primary_tag=>1} @rsf;
	    my @types = keys %types;
	    if (@types == 1 && $types[0] eq $self->primary_tag) {
		return join ("\n",map {$_->gff3_string(1,$parent_tree,{},$primary_id)} @rsf);
	    }
	}

	@parent_ids = keys %{$parent_tree->{$primary_id}};
    }

    my $group      = $self->format_attributes(\@parent_ids,$force_id);
    my $name       = $self->name;

    my $class = $self->class;
    my $strand = ('-','.','+')[$self->strand+1];
    my $p = join("\t",
		 $self->seq_id||'.',
		 $self->source||'.',
		 $self->method||'.',
		 $self->start||'.',
		 $self->stop||'.',
		 defined($self->score) ? $self->score : '.',
		 $strand||'.',
		 defined($self->phase) ? $self->phase : '.',
		 $group||'');
    return join("\n",
		$p,
		map {$_->gff3_string(1,$parent_tree,$seenit)} @rsf);
}

sub _real_or_dummy_id {
    my $self = shift;
    my $id   = $self->primary_id;
    return $id if defined $id;
    return return (overload::StrVal($self) =~ /0x([a-f0-9]+)/)[0];
}

sub _traverse {
    my $self   = shift;
    my $tree   = shift;  # tree => {$child}{$parent} = 1
    my $parent = shift;
    my $id     = $self->_real_or_dummy_id;
    defined $id or return;
    $tree->{$id}{$parent->_real_or_dummy_id}++ if $parent;
    $_->_traverse($tree,$self) foreach $self->get_SeqFeatures;
}

sub db { return }

sub source_tag {
  my $self = shift;
  my $d = $self->{source};
  $self->{source} = shift if @_;
  $d;
}

# This probably should be deleted.  Not sure why it's here, but might
# have been added for Ace::Sequence::Feature-compliance.
sub introns {
  my $self = shift;
  return;
}

sub has_tag { 
    my $self = shift;
    my $tag  = shift;
    return exists $self->{attributes}{$tag};
}

sub escape {
  my $self     = shift;
  my $toencode = shift;
  $toencode    =~ s/([^a-zA-Z0-9_.:?^*\(\)\[\]@!+-])/uc sprintf("%%%02x",ord($1))/eg;
  $toencode;
}

sub all_tags {
  my $self = shift;
  return keys %{$self->{attributes}};
}

sub add_tag_value {
  my $self = shift;
  my ($tag_name,@tag_values) = @_;
  push @{$self->{attributes}{$tag_name}},@tag_values;
}

sub remove_tag {
  my $self = shift;
  my $tag_name = shift;
  delete $self->{attributes}{$tag_name};
}

sub each_tag_value {
  my $self = shift;
  my $tag  = shift;
  my $value = $self->{attributes}{$tag} or return;
  my $ref = CORE::ref $value;
  return $ref && $ref eq 'ARRAY' ? @{$self->{attributes}{$tag}}
                                 : $self->{attributes}{$tag};
}

sub get_Annotations {
  my $self = shift;
  my $tag  = shift;
  my @values = $self->get_tag_values($tag);
  return $values[0] if @values == 1;
  return @values;
}

sub format_attributes {
  my $self        = shift;
  my $parent      = shift;
  my $fallback_id = shift;

  my @tags = $self->get_all_tags;
  my @result;
  for my $t (@tags) {
    my @values = $self->get_tag_values($t);
    push @result,join '=',$self->escape($t),join(',', map {$self->escape($_)} @values) if @values;
  }
  #my $id        = $self->escape($self->_real_or_dummy_id) || $fallback_id;
  my $id        = $fallback_id || $self->escape($self->_real_or_dummy_id);

  my $parent_id;
  if (@$parent) {
      $parent_id  = join (',',map {$self->escape($_)} @$parent);
  }

  my $name = $self->display_name;
  unshift @result,"ID=".$id                                    if defined $id;
  unshift @result,"Parent=".$parent_id                         if defined $parent_id;
  unshift @result,"Name=".$self->escape($name)                 if defined $name;
  return join ';',@result;
}

=head2 clone

 Title   : clone
 Usage   : my $feature = $seqfeature->clone
 Function: Create a deep copy of the feature
 Returns : A copy of the feature
 Args    : none

=cut

sub clone {
  my $self  = shift;
  my %clone  = %$self;
  # overwrite attributes
  my $clone = bless \%clone,CORE::ref($self);
  $clone{attributes} = {};
  for my $k (keys %{$self->{attributes}}) {
    @{$clone{attributes}{$k}} = @{$self->{attributes}{$k}};
  }
  return $clone;
}

=head2 refseq

 Title   : refseq
 Usage   : $ref = $s->refseq([$newseq] [,$newseqclass])
 Function: get/set reference sequence
 Returns : current reference sequence
 Args    : new reference sequence and class (optional)
 Status  : Public

This method will get or set the reference sequence.  Called with no
arguments, it returns the current reference sequence.  Called with any
Bio::SeqFeatureI object that provides the seq_id(), start(), end() and
strand() methods.

The method will generate an exception if you attempt to set the
reference sequence to a sequence that has a different seq_id from the
current feature.

=cut

sub refseq {
  my $self = shift;
  my $d    = $self->{refseq};
  if (@_) {
    my $newref = shift;
    $self->throw("attempt to set refseq using a feature that does not share the same seq_id")
      unless $newref->seq_id eq $self->seq_id;
    $self->{refseq} = $newref;
  }
  return $d;
}

sub DESTROY { }

1;

__END__

=head1 SEE ALSO

L<Bio::Graphics::Feature> 

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.eduE<gt>.

Copyright (c) 2006 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
