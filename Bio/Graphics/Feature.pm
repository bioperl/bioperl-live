package Bio::Graphics::Feature;

=head1 NAME

Bio::Graphics::Feature - A simple feature object for use with Bio::Graphics::Panel

=head1 SYNOPSIS

 use Bio::Graphics::Feature;

 # create a simple feature with no internal structure
 $f = Bio::Graphics::Feature->new(-start => 1000,
                                  -stop  => 2000,
                                  -type  => 'transcript',
                                  -name  => 'alpha-1 antitrypsin',
				  -desc  => 'an enzyme inhibitor',
                                 );

 # create a feature composed of multiple segments, all of type "similarity"
 $f = Bio::Graphics::Feature->new(-segments => [[1000,1100],[1500,1550],[1800,2000]],
                                  -name     => 'ABC-3',
                                  -type     => 'gapped_alignment',
                                  -subtype  => 'similarity');

 # build up a gene exon by exon
 $e1 = Bio::Graphics::Feature->new(-start=>1,-stop=>100,-type=>'exon');
 $e2 = Bio::Graphics::Feature->new(-start=>150,-stop=>200,-type=>'exon');
 $e3 = Bio::Graphics::Feature->new(-start=>300,-stop=>500,-type=>'exon');
 $f  = Bio::Graphics::Feature->new(-segments=>[$e1,$e2,$e3],-type=>'gene');

=head1 DESCRIPTION

This is a simple Bio::SeqFeatureI-compliant object that is compatible
with Bio::Graphics::Panel.  With it you can create lightweight feature
objects for drawing.

All methods are as described in L<Bio::SeqFeatureI> with the following additions:

=head2 The new() Constructor

 $feature = Bio::Graphics::Feature->new(@args);

This method creates a new feature object.  You can create a simple
feature that contains no subfeatures, or a hierarchically nested object.

Arguments are as follows:

  -start       the start position of the feature
  -end         the stop position of the feature
  -stop        an alias for end
  -name        the feature name (returned by seqname())
  -type        the feature type (returned by primary_tag())
  -source      the source tag
  -desc        a description of the feature
  -segments    a list of subfeatures (see below)
  -subtype     the type to use when creating subfeatures
  -strand      the strand of the feature (one of -1, 0 or +1)
  -id          an alias for -name
  -seqname     an alias for -name
  -primary_id  an alias for -name
  -display_id  an alias for -name
  -display_name an alias for -name  (do you get the idea the API has changed?)
  -url         a URL to link to when rendered with Bio::Graphics
  -configurator an object (like a Bio::Graphics::FeatureFile) that knows how 
               to configure the graphical representation of the object based
               on its type.
  -attributes  a hashref of tag value attributes, in which the key is the tag
               and the value is an array reference of values
  -factory     a reference to a feature factory, used for compatibility with
               more obscure parts of Bio::DB::GFF

The subfeatures passed in -segments may be an array of
Bio::Graphics::Feature objects, or an array of [$start,$stop]
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

=item configurator()

Get/set the configurator that knows how to adjust the graphical
representation of this feature based on its type.  Currently the only
configurator that will work is Bio::Graphics::FeatureFile.

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
*method      = \&type;
*source      = \&source_tag;

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

*get_all_SeqFeatures = *get_SeqFeatures = \&segments;
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
    || $arg{-display_name} || $arg{-id} || $arg{-primary_id};
  $self->{type}    = $arg{-type}   || 'feature';
  $self->{subtype} = $arg{-subtype} if exists $arg{-subtype};
  $self->{source}  = $arg{-source} || $arg{-source_tag} || '';
  $self->{score}   = $arg{-score}   if exists $arg{-score};
  $self->{start}   = $arg{-start};
  $self->{stop}    = $arg{-end} || $arg{-stop};
  $self->{ref}     = $arg{-ref};
  for my $option (qw(class url seq phase desc attributes factory configurator)) {
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
    } else {
      push @segments,$seg;
    }
  }
  if (@segments) {
    local $^W = 0;  # some warning of an uninitialized variable...
    $self->{segments} = [ sort {$a->start <=> $b->start } @segments ];
    $self->{start}    = $self->{segments}[0]->start;
    ($self->{stop})   = sort { $b <=> $a } map { $_->end } @segments;
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


=head2 factory

 Title   : factory
 Usage   : $factory = $obj->factory([$new_factory])
 Function: Returns the feature factory from which this feature was generated.
           Mostly for compatibility with weird dependencies in gbrowse.
 Returns : A feature factory
 Args    : None

=cut

sub factory {
  my $self = shift;
  my $d = $self->{factory};
  $self->{factory} = shift if @_;
  $d;
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

# get/set the configurator (Bio::Graphics::FeatureFile) for this feature
sub configurator {
  my $self = shift;
  my $d = $self->{configurator};
  $self->{configurator} = shift if @_;
  $d;
}

# get/set the url for this feature
sub url {
  my $self = shift;
  my $d = $self->{url};
  $self->{url} = shift if @_;
  $d;
}

# make a link
sub make_link {
  my $self = shift;

  if (my $url = $self->url) {
    return $url;
  }

  elsif (my $configurator = $self->configurator) {
    return $configurator->make_link($self);
  }

  else {
    return;
  }
}

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

L<Bio::Graphics::Panel>,L<Bio::Graphics::Glyph>,
L<GD>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
