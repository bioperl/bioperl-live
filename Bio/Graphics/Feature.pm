package Bio::Graphics::Feature;
use strict;
use Bio::SeqFeatureI;

use vars '$VERSION','@ISA';
$VERSION = '1.31';
@ISA  = 'Bio::SeqFeatureI';

*stop        = \&end;
*info        = \&name;
*seqname     = \&name;
*type        = \&primary_tag;
*exons       = *sub_SeqFeature = *merged_segments = \&segments;
*class       = *method = \&type;
*source      = \&source_tag;

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
  $self->{strand}  = $arg{-strand} >= 0 ? +1 : -1;
  $self->{name}    = $arg{-name};
  $self->{type}    = $arg{-type}   || 'feature';
  $self->{subtype} = $arg{-subtype} if exists $arg{-subtype};
  $self->{source}  = $arg{-source} || $arg{-source_tag} || '';
  $self->{score}   = $arg{-score}  || 0;
  $self->{start}   = $arg{-start};
  $self->{stop}    = $arg{-end} || $arg{-stop};
  $self->{ref}     = $arg{-ref};
  $self->{url}     = $arg{-url} if $arg{-url};

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

  my @segments = @{$self->{segments}};

  for my $seg (@_) {
    if (ref($seg) eq 'ARRAY') {
      my ($start,$stop) = @{$seg};
      next unless defined $start && defined $stop;  # fixes an obscure bug somewhere above us
      my $strand = $self->{strand};

      if ($start > $stop) {
	($start,$stop) = ($stop,$start);
	$strand *= -1;
      }
      push @segments,$self->new(-start=>$start,
				-stop=>$stop,
				-strand=>$strand,
				-type  => $type);
    } else {
      push @segments,$seg;
    }
  }
  if (@segments) {
    $self->{segments} = [ sort {$a->start <=> $b->start } @segments ];
    $self->{start}    = $self->{segments}[0]->start;
    ($self->{stop})   = sort { $b <=> $a } map { $_->end } @segments;
  }
}

sub location {
  my $self = shift;

  require Bio::Location::Split;
  my @segments = $self->segments;
  if (@segments) {
    my $split = Bio::Location::Split->new;
    foreach (@segments) {
      $split->add_sub_Location(Bio::Location::Simple->new(-start  => $_->start,
							  -end    => $_->end,
							  -strand => $_->strand
							 ));
    }
    return $split;
  }
  return Bio::Location::Simple->new(-start  => $self->start,
				    -end    => $self->end,
				    -strand => $self->strand);
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
sub name            { shift->{name}        }
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
  return scalar('n' x $self->length);
}
*dna = \&seq;

sub low {
  my $self = shift;
  return $self->start < $self->end ? $self->start : $self->end;
}

sub high {
  my $self = shift;
  return $self->start > $self->end ? $self->start : $self->end;
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

sub DESTROY { }

1;

__END__

=head1 NAME

Bio::Graphics::Feature - A simple feature object for use with Bio::Graphics::Panel

=head1 SYNOPSIS

 use Bio::Graphics::Feature;

 # create a simple feature with no internal structure
 $f = Bio::Graphics::Feature->new(-start => 1000,
                                  -stop  => 2000,
                                  -type  => 'transcript',
                                  -name  => 'alpha-1 antitrypsin'
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
  -stop        the stop position of the feature
  -end         an alias for stop
  -name        the feature name (returned by seqname())
  -type        the feature type (returned by primary_tag())
  -source      the source tag
  -segments    a list of subfeatures (see below)
  -subtype     the type to use when creating subfeatures

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

=item add_segment(@segments)

Add one or more segments (a subfeature).  Segments can either be
Feature objects, or [start,stop] arrays, as in the -segments argument
to new().  The feature endpoints are automatically adjusted.

=item segments()

An alias for sub_SeqFeatures().

=item merged_segments()

Another alias for sub_SeqFeatures().

=item stop()

An alias for end().

=item name()

An alias for seqname().

=item exons()

An alias for sub_SeqFeatures() (you don't want to know why!)

=back

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
