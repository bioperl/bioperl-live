=head1 NAME

Bio::Graphics::Glyph::Factory - Factory for Bio::Graphics::Glyph objects

=head1 SYNOPSIS

See L<Bio::Graphics::Panel>.

=head1 DESCRIPTION

This class is used internally by Bio::Graphics to generate new Glyph
objects by combining a list of features with the user's desired
configuration.  It is intended to be used internally by Bio::Graphics.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org             - General discussion
  http://bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via
email or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Lincoln Stein

Email - lstein@cshl.org

=head1 SEE ALSO

L<Bio::Graphics::Panel>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with an "_"
(underscore).

=cut

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

=head2 new

  Title   : new
  Usage   : $f = Bio::Graphics::Glyph::Factory->new(
                     -stylesheet => $stylesheet,
		     -glyph_map  => $glyph_map,
		     -options    => $options);
  Function : create a new Bio::Graphics::Glyph::Factory object
  Returns  : the new object
  Args     : $stylesheet is a Bio::Das::Stylesheet object that can
                 convert Bio::Das feature objects into glyph names and
                 associated options.
             $glyph_map is a hash that maps primary tags to glyph names.
             $options is a hash that maps option names to their values.
  Status   : Internal to Bio::Graphics

=cut

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

=head2 clone

  Title    : clone
  Usage    : $f2 = $f->clone
  Function : Deep copy of a factory object
  Returns  : a deep copy of the factory object
  Args     : None
  Status   : Internal to Bio::Graphics

=cut

sub clone {
  my $self = shift;
  my %new = %$self;
  my $new = bless \%new,ref($self);
  $new;
}

=head2 stylesheet

  Title    : stylesheet
  Usage    : $stylesheet = $f->stylesheet
  Function : accessor for stylesheet
  Returns  : a Bio::Das::Stylesheet object
  Args     : None
  Status   : Internal to Bio::Graphics

=cut

sub stylesheet { shift->{stylesheet}  }

=head2 glyph_map

  Title    : glyph_map
  Usage    : $map = $f->glyph_map
  Function : accessor for the glyph map
  Returns  : a hash mapping primary tags to glyphs
  Args     : None
  Status   : Internal to Bio::Graphics

=cut

sub glyph_map  { shift->{glyph_map}   }

=head2 option_map

  Title    : option_map
  Usage    : $map = $f->option_map
  Function : accessor for the option map
  Returns  : a hash mapping option names to values
  Args     : None
  Status   : Internal to Bio::Graphics

=cut

sub option_map { shift->{options}     }

=head2 global_opts

  Title    : global_opts
  Usage    : $map = $f->global_opts
  Function : accessor for global options
  Returns  : a hash mapping option names to values
  Args     : None
  Status   : Internal to Bio::Graphics

This returns a set of defaults for option values.

=cut

sub global_opts{ shift->{global_opts} }

=head2 panel

  Title    : panel
  Usage    : $panel = $f->panel
  Function : accessor for Bio::Graphics::Panel
  Returns  : a Bio::Graphics::Panel
  Args     : None
  Status   : Internal to Bio::Graphics

This returns the panel with which the factory is associated.

=cut

sub panel      { shift->{panel}       }

=head2 scale

  Title    : scale
  Usage    : $scale = $f->scale
  Function : accessor for the scale
  Returns  : a floating point number
  Args     : None
  Status   : Internal to Bio::Graphics

This returns the scale, in pixels/bp for glyphs constructed by this
factory.

=cut

sub scale      { shift->{panel}->scale }

=head2 font

  Title    : font
  Usage    : $font = $f->font
  Function : accessor for the font
  Returns  : a font name
  Args     : None
  Status   : Internal to Bio::Graphics

This returns a GD font name.

=cut

sub font       {
  my $self = shift;
  my $glyph = shift;
  $self->option($glyph,'font') || $self->{font};
}

=head2 map_pt

  Title    : map_pt
  Usage    : @pixel_positions = $f->map_pt(@bp_positions)
  Function : map bp positions to pixel positions
  Returns  : a list of pixel positions
  Args     : a list of bp positions
  Status   : Internal to Bio::Graphics

The real work is done by the panel, but factory subclasses can
override if desired.

=cut

sub map_pt {
  my $self = shift;
  my @result = $self->panel->map_pt(@_);
  return wantarray ? @result : $result[0];
}

=head2 map_no_trunc

  Title    : map_no_trunc
  Usage    : @pixel_positions = $f->map_no_trunc(@bp_positions)
  Function : map bp positions to pixel positions
  Returns  : a list of pixel positions
  Args     : a list of bp positions
  Status   : Internal to Bio::Graphics

Same as map_pt(), but it will NOT clip pixel positions to be within
the drawing frame.

=cut

sub map_no_trunc {
  my $self = shift;
  my @result = $self->panel->map_no_trunc(@_);
  return wantarray ? @result : $result[0];
}

=head2 translate_color

  Title    : translate_color
  Usage    : $index = $f->translate_color($color_name)
  Function : translate symbolic color names into GD indexes
  Returns  : an integer
  Args     : a color name in format "green" or "#00FF00"
  Status   : Internal to Bio::Graphics

The real work is done by the panel, but factory subclasses can
override if desired.

=cut

sub translate_color {
  my $self = shift;
  my $color_name = shift;
  $self->panel->translate_color($color_name);
}

=head2 glyph

  Title    : glyph
  Usage    : @glyphs = $f->glyph($level,$feature1,$feature2...)
  Function : transform features into glyphs.
  Returns  : a list of Bio::Graphics::Glyph objects
  Args     : a feature "level", followed by a list of FeatureI objects.
  Status   : Internal to Bio::Graphics

The level is used to track the level of nesting of features that have
subfeatures.

=cut

# create a glyph
sub make_glyph {
  my $self  = shift;
  my $level = shift;
  my @result;
  my $panel = $self->panel;
  my ($leftmost,$rightmost) = ($panel->left,$panel->right);
  my $flip   = $panel->flip;

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
				 -flip     => $flip,
				 -level    => $level);

    # this is removing glyphs that are not onscreen at all.
    # But never remove tracks!
    push @result,$glyph if $type eq 'track'
	|| ($glyph->{left} + $glyph->{width} > $leftmost && $glyph->{left} < $rightmost);

  }
  return wantarray ? @result : $result[0];
}

=head2 feature_to_glyph

  Title    : feature_to_glyph
  Usage    : $glyph_name = $f->feature_to_glyph($feature)
  Function : choose the glyph name given a feature
  Returns  : a glyph name
  Args     : a Bio::Seq::FeatureI object
  Status   : Internal to Bio::Graphics

=cut

sub feature_to_glyph {
  my $self    = shift;
  my $feature = shift;

  return scalar $self->{stylesheet}->glyph($feature) if $self->{stylesheet};
  my $map = $self->glyph_map    or return 'generic';
  if (ref($map) eq 'CODE') {
    my $val = eval {$map->($feature)};
    warn $@ if $@;
    return $val || 'generic';
  }
  return $map->{$feature->primary_tag} || 'generic';
}


=head2 set_option

  Title    : set_option
  Usage    : $f->set_option($option_name=>$option_value)
  Function : set or change an option
  Returns  : nothing
  Args     : a name/value pair
  Status   : Internal to Bio::Graphics

=cut

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
      my $val = eval { $value->($feature,$option_name,$partno,$total_parts,$glyph)};
      warn $@ if $@;
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


=head2 options

  Title    : options
  Usage    : @option_names = $f->options
  Function : return all configured option names
  Returns  : a list of option names
  Args     : none
  Status   : Internal to Bio::Graphics

=cut

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
