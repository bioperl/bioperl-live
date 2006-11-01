package Bio::Graphics::FeatureBase;

=head1 NAME

Bio::Graphics::FeatureBase - Base class for Bio::Graphics::Feature

=head1 SYNOPSIS

 See Bio::Graphics::Feature for full synopsis.

=head1 DESCRIPTION

This is the base class for Bio::Graphics::Feature. It has all the
methods of Bio::Graphics::Feature except for those that are required
to interface with Bio::Graphics::FeatureFile, namely factory(),
configurator(), url(), and make_link().  Please see
L<Bio::Graphics::Feature> for full documentation.

=cut

use strict;

use base qw(Bio::Root::Root Bio::SeqFeatureI Bio::LocationI Bio::SeqI Bio::RangeI);

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
*abs_start      = \&start;
*abs_end        = \&end;
*abs_strand     = \&strand;

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
	$obj->{'_annotation'} = new Bio::Annotation::Collection;
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
# Bio::Graphics::Feature->new(
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
	$arg{-strand} = "+" && $self->{strand} ='1';
	$arg{-strand} = "-" && $self->{strand} = '-1';
	$arg{-strand} = "." && $self->{strand} = '0';
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
  $self->{stop}    = $arg{-end} || $arg{-stop};
  $self->{ref}     = $arg{-seq_id} || $arg{-ref};
  for my $option (qw(class url seq phase desc attributes primary_id)) {
    $self->{$option} = $arg{"-$option"} if exists $arg{"-$option"};
  }

  # fix start, stop
  if (defined $self->{stop} && defined $self->{start}
      && $self->{stop} < $self->{start}) {
    @{$self}{'start','stop'} = @{$self}{'stop','start'};
    $self->{strand} *= -1;
  }

  my @segments;
  if (my $s = $arg{-segments}) {
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
			        -class  => $class);
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
    # this was killing performance!
    #  $self->{segments} = [ sort {$a->start <=> $b->start } @segments ];
    # this seems much faster and seems to work still
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
sub primary_tag     { shift->{type}        }
sub name            {
  my $self = shift;
  my $d    = $self->{name};
  $self->{name} = shift if @_;
  $d;
}
sub seq_id          { shift->ref()         }
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
  $d;
}
sub end    {
  my $self = shift;
  my $d = $self->{stop};
  $self->{stop} = shift if @_;
  $d;
}
sub strand {
  my $self = shift;
  my $d = $self->{strand};
  $self->{strand} = shift if @_;
  $d;
}
sub length {
  my $self = shift;
  return $self->end - $self->start + 1;
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

sub display_name { shift->name }

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
  my $d    = $self->{desc};
  $self->{desc} = shift if @_;
  $d;
}

sub attributes {
  my $self = shift;
  if (@_) {
    return $self->each_tag_value(@_);
  } else {
    return $self->{attributes} ? %{$self->{attributes}} : ();
  }
}

sub primary_id {
  my $self = shift;
  my $d = $self->{primary_id};
  $self->{primary_id} = shift if @_;
  $d;
}

sub notes {
  return shift->desc;
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
  return "$low..$high";
}
sub phase { shift->{phase} }
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
  my $recurse = shift;

  if ($self->version == 3) {
    return $self->gff3_string(@_);
  }

  my $name  = $self->name;
  my $class = $self->class;
  my $group = "$class $name" if $name;
  my $strand = ('-','.','+')[$self->strand+1];
  my $string;
  $string .= join("\t",$self->ref||'.',$self->source||'.',$self->method||'.',
                       $self->start||'.',$self->stop||'.',
                       $self->score||'.',$strand||'.',$self->phase||'.',
                       $group||'');
  $string .= "\n";
  if ($recurse) {
    foreach ($self->sub_SeqFeature) {
      $string .= $_->gff_string($recurse);
    }
  }
  $string;
}

sub gff3_string {
  my $self              = shift;
  my ($recurse,$parent) = @_;

  my $name  = $self->name;
  my $class = $self->class;
  my $group = $self->format_attributes($parent);
  my $strand = ('-','.','+')[$self->strand+1];
  my $p      = join("\t",$self->ref||'.',$self->source||'.',$self->method||'.',
		    $self->start||'.',$self->stop||'.',
		    $self->score||'.',$strand||'.',$self->phase||'.',
		    $group||'');

  # the "homogeneous" flag will be true if the parent and children are all of the same type,
  # meaning that they can be collapsed into a set of children with all the same ID
  my ($parent_type,$homogeneous);
  $homogeneous = 1;
  my @children;
  if ($recurse) {
    foreach ($self->sub_SeqFeature) {
      push @children,$_->gff3_string(1,$self);
      $parent_type   ||= $self->type;
      $homogeneous   &&= $_->type eq $parent_type && !defined $_->primary_id;
    }
  }

  # if we get here we're dealing with a homogeneous set of Parent[child,child...]
  # where parent and child all have the same type. In this case, we omit the Parent
  # and give the children the same ID. This removes an extraneous level of parentage.

  if (@children && $homogeneous) {
    foreach (@children) { 
      s/Parent=/ID=/g; 
    } # replace Parent tag with ID
    return join "\n",@children;
  }

  return join("\n",$p,@children);
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

sub has_tag { exists shift->{attributes}{shift()} }

sub escape {
  my $self    = shift;
  my $toencode = shift;
  $toencode    =~ s/([^a-zA-Z0-9_. :?^*\(\)\[\]@!+-])/uc sprintf("%%%02x",ord($1))/eg;
#  $toencode    =~ tr/ /+/;  # not needed in GFF3
  $toencode;
}

sub all_tags {
  my $self = shift;
  return keys %{$self->{attributes}};
}
sub each_tag_value {
  my $self = shift;
  my $tag  = shift;
  my $value = $self->{attributes}{$tag} or return;
  return CORE::ref $value ? @{$self->{attributes}{$tag}}
                          : $self->{attributes}{$tag};
}

sub format_attributes {
  my $self   = shift;
  my $parent = shift;
  my @tags = $self->all_tags;
  my @result;
  for my $t (@tags) {
    my @values = $self->each_tag_value($t);
    push @result,join '=',$self->escape($t),$self->escape($_) foreach @values;
  }
  my $id   = $self->primary_id;
  my $name = $self->display_name;
  push @result,"ID=".$self->escape($id)                     if defined $id;
  push @result,"Parent=".$self->escape($parent->primary_id) if defined $parent;
  push @result,"Name=".$self->escape($name)                 if defined $name;
  return join ';',@result;
}

sub DESTROY { }

1;

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


__END__

=head1 SEE ALSO

L<Bio::Graphics::Feature>, L<Bio::Graphics::FeatureFile>,
L<Bio::Graphics::Panel>,L<Bio::Graphics::Glyph>, L<GD>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2006 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
