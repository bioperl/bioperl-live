package Bio::Graphics::Glyph::Factory;

use strict;
use Carp qw(:DEFAULT cluck);
use GD;

my %LOADED_GLYPHS = ();
my %GENERIC_OPTIONS = (
		       bgcolor    => 'turquoise',
		       fgcolor    => 'black',
		       fontcolor  => 'black',
		       font2color => 'turquoise',
		       height     => 8,
		       font       => gdSmallFont,
		       bump       => +1,       # bump by default (perhaps a mistake?)
		       );

sub new {
  my $class = shift;
  my $panel = shift;
  my %args = @_;
  my $stylesheet = $args{-stylesheet};   # optional, for Bio::Das compatibility
  my $map        = $args{-map};          # map type name to glyph name
  my $options    = $args{-options};      # map type name to glyph options
  return bless {
		stylesheet => $stylesheet,
		glyph_map  => $map,
		options    => $options,
		panel      => $panel,
		},$class;
}
sub clone {
  my $self = shift;
  my %new = %$self;
  my $new = bless \%new,ref($self);
  $new;
}
sub stylesheet { shift->{stylesheet}  }
sub glyph_map  { shift->{glyph_map}   }
sub option_map { shift->{options}     }
sub global_opts{ shift->{global_opts} }
sub panel      { shift->{panel}       }
sub scale      { shift->{panel}->scale }
sub font       {
  my $self = shift;
  my $glyph = shift;
  $self->option($glyph,'font') || $self->{font};
}

sub map_pt {
  my $self = shift;
  my @result = $self->panel->map_pt(@_);
  return wantarray ? @result : $result[0];
}

sub map_no_trunc {
  my $self = shift;
  my @result = $self->panel->map_no_trunc(@_);
  return wantarray ? @result : $result[0];
}

sub translate_color {
  my $self = shift;
  my $color_name = shift;
  $self->panel->translate_color($color_name);
}

# create a glyph
sub make_glyph {
  my $self  = shift;
  my $level = shift;
  my @result;
  my $panel = $self->panel;
  my ($leftmost,$rightmost) = ($panel->left,$panel->right);

  for my $f (@_) {

    my $type = $self->feature_to_glyph($f);
    my $glyphclass = 'Bio::Graphics::Glyph';
    $type ||= 'generic';
    $glyphclass .= "\:\:\L$type";

    unless ($LOADED_GLYPHS{$glyphclass}++) {
      carp("the requested glyph class, ``$type'' is not available: $@")
	unless (eval "require $glyphclass");
    }
    my $glyph = $glyphclass->new(-feature  => $f,
				 -factory  => $self,
				 -level    => $level);

    # this is removing glyphs that are not onscreen at all.
    # But never remove tracks!
    push @result,$glyph if $type eq 'track'
	|| ($glyph->{left} + $glyph->{width} > $leftmost && $glyph->{left} < $rightmost);

  }
  return wantarray ? @result : $result[0];
}

sub feature_to_glyph {
  my $self    = shift;
  my $feature = shift;

  return scalar $self->{stylesheet}->glyph($feature) if $self->{stylesheet};
  my $map = $self->glyph_map    or return 'generic';
  return $map->($feature)              || 'generic' if ref($map) eq 'CODE';
  return $map->{$feature->primary_tag} || 'generic';
}

sub set_option {
  my $self = shift;
  my ($option_name,$option_value) = @_;
  $self->{overriding_options}{lc $option_name} = $option_value;
}

# options:
#    the overriding_options hash has precedence
#    ...followed by the option_map
#    ...followed by the stylesheet
#    ...followed by generic options
sub option {
  my $self = shift;
  my ($glyph,$option_name,$partno,$total_parts) = @_;
  return unless defined $option_name;
  $option_name = lc $option_name;   # canonicalize

  return $self->{overriding_options}{$option_name} 
    if exists $self->{overriding_options} && exists $self->{overriding_options}{$option_name};

  if (my $map    = $self->option_map) {
    if (defined(my $value  = $map->{$option_name})) {
      my $feature = $glyph->feature;
      return $value unless ref $value eq 'CODE';
      return unless $feature->isa('Bio::SeqFeatureI');
      my $val = $value->($feature,$option_name,$partno,$total_parts,$glyph);
      return defined $val && $val eq '*default*' ? $GENERIC_OPTIONS{$option_name} : $val;
    }
  }

  if (my $ss = $self->stylesheet) {
    my($glyph,%options) = $ss->glyph($glyph->feature);
    my $value = $options{$option_name};
    return $value if defined $value;
  }

  return $GENERIC_OPTIONS{$option_name};
}

# return names of all the options in the option hashes
sub options {
  my $self = shift;
  my %options;
  if (my $map    = $self->option_map) {
    $options{lc($_)}++ foreach keys %$map;
  }
  $options{lc($_)}++ foreach keys %GENERIC_OPTIONS;
  return keys %options;
}

1;
