### TODO ###

## From FeatureFile: ##
# render our features onto a panel using configuration data
# return the number of tracks inserted

=over 4

=item ($rendered,$panel) = $features-E<gt>render([$panel])

Render features in the data set onto the indicated
Bio::Graphics::Panel.  If no panel is specified, creates one.

In a scalar context returns the number of tracks rendered.  In a list
context, returns a two-element list containing the number of features
rendered and the panel.  Use this form if you want the panel created
for you.

=back

=cut

sub render {
  my $self = shift;
  my $panel = shift;
  my ($position_to_insert,$options,$max_bump,$max_label) = @_;

  $panel ||= $self->new_panel;

  # count up number of tracks inserted
  my $tracks = 0;
  my $color;
  my %types = map {$_=>1} $self->configured_types;

  my @configured_types   = grep {exists $self->features->{$_}} $self->configured_types;
  my @unconfigured_types = sort grep {!exists $types{$_}}      $self->types;

  my @base_config = $self->style('general');

  my @override = ();
  if ($options && ref $options eq 'HASH') {
    @override = %$options;
  } else {
    $options ||= 0;
    if ($options == 1) {  # compact
      push @override,(-bump => 0,-label=>0);
    } elsif ($options == 2) { #expanded
      push @override,(-bump=>1);
    } elsif ($options == 3) { #expand and label
      push @override,(-bump=>1,-label=>1);
    } elsif ($options == 4) { #hyperexpand
      push @override,(-bump => 2);
    } elsif ($options == 5) { #hyperexpand and label
      push @override,(-bump => 2,-label=>1);
    }
  }

  for my $type (@configured_types,@unconfigured_types) {
    my $features = $self->features($type);
    my @auto_bump;
    push @auto_bump,(-bump  => @$features < $max_bump)  if defined $max_bump;
    push @auto_bump,(-label => @$features < $max_label) if defined $max_label;

    my @config = ( -glyph   => 'segments',         # really generic
		   -bgcolor => $COLORS[$color++ % @COLORS],
		   -label   => 1,
		   -key     => $type,
		   @auto_bump,
		   @base_config,         # global
		   $self->style($type),  # feature-specific
		   @override,
		 );
    if (defined($position_to_insert)) {
      $panel->insert_track($position_to_insert++,$features,@config);
    } else {
      $panel->add_track($features,@config);
    }
    $tracks++;
  }
  return wantarray ? ($tracks,$panel) : $tracks;
}

