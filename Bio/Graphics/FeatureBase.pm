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
=back

=cut

use strict;
use Bio::Root::Root;
use Bio::SeqFeatureI;
use Bio::SeqI;
use Bio::LocationI;

use vars '@ISA';
@ISA  = qw(Bio::Root::Root Bio::SeqFeatureI Bio::LocationI Bio::SeqI);

*stop        = \&end;
*info        = \&name;
*seqname     = \&name;
*type        = \&primary_tag;
*exons       = *sub_SeqFeature = *merged_segments = \&segments;
*get_all_SeqFeatures = *get_SeqFeatures = \&segments;
*method      = \&type;
*source      = \&source_tag;
*get_tag_values = \&each_tag_value;
*add_SeqFeature = \&add_segment;
*get_all_tags   = \&all_tags;
*abs_ref        = \&ref;
*abs_start      = \&start;
*abs_end        = \&end;

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
sub hit    { return; }

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

      $min_start = $seg->start if $seg->start < $min_start;
      $max_stop  = $seg->end   if $seg->end   > $max_stop;
    }
  }
  if (@segments) {
    local $^W = 0;  # some warning of an uninitialized variable...
    $self->{segments} = [ sort {$a->start <=> $b->start } @segments ];
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
  my $dna =  exists $self->{seq} ? $self->{seq} : '';
  # $dna .= 'n' x ($self->length - CORE::length($dna));
  return $dna;
}
*dna = \&seq;

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
    return $self->{attributes};
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
  return defined($d) ? $d : ucfirst $self->method;
}

sub gff_string {
  my $self    = shift;
  my $recurse = shift;

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
  my $string;
  $string .= join("\t",$self->ref||'.',$self->source||'.',$self->method||'.',
                       $self->start||'.',$self->stop||'.',
                       $self->score||'.',$strand||'.',$self->phase||'.',
                       $group||'');
  $string .= "\n";
  if ($recurse) {
    foreach ($self->sub_SeqFeature) {
      $string .= $_->gff3_string(1,$self->name);
    }
  }
  $string;
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

sub has_tag { }

sub escape {
  my $toencode = shift;
  $toencode    =~ s/([^a-zA-Z0-9_. :?^*\(\)\[\]@!-])/uc sprintf("%%%02x",ord($1))/eg;
  $toencode    =~ tr/ /+/;
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
    push @result,join '=',escape($t),escape($_) foreach @values;
  }
  push @result,"ID=".escape($self->name) if defined $self->name;
  push @result,"Parent=".escape($parent) if defined $parent;
  return join ';',@result;
}

sub DESTROY { }

1;

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
