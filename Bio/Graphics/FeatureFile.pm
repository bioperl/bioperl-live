package Bio::Graphics::FeatureFile;

# $Id$
# This package parses and renders a simple tab-delimited format for features.
# It is simpler than GFF, but still has a lot of expressive power.
# See __END__ for the file format

=head1 NAME

Bio::Graphics::FeatureFile -- A set of Bio::Graphics features, stored in a file

=head1 SYNOPSIS

 use Bio::Graphics::FeatureFile;
 my $data  = Bio::Graphics::FeatureFile->new(-file => 'features.txt');


 # create a new panel and render contents of the file onto it
 my $panel = $data->new_panel;
 my $tracks_rendered = $data->render($panel);

 # or do it all in one step
 my ($tracks_rendered,$panel) = $data->render;

 # for more control, render tracks individually
 my @feature_types = $data->types;
 for my $type (@feature_types) {
    my $features = $data->features($type);
    my %options  = $data->style($type);
    $panel->add_track($features,%options);  # assuming we have a Bio::Graphics::Panel
 }

 # get individual settings
 my $est_fg_color = $data->setting(EST => 'fgcolor');

 # or create the FeatureFile by hand

 # add a type
 $data->add_type(EST => {fgcolor=>'blue',height=>12});

 # add a feature
 my $feature = Bio::Graphics::Feature->new(....); # or some other SeqI
 $data->add_feature($feature=>'EST');

=head1 DESCRIPTION

The Bio::Graphics::FeatureFile module reads and parses files that
describe sequence features and their renderings.  It accepts both GFF
format and a more human-friendly file format described below.  Once a
FeatureFile object has been initialized, you can interrogate it for
its consistuent features and their settings, or render the entire file
onto a Bio::Graphics::Panel.

This moduel is a precursor of Jason Stajich's
Bio::Annotation::Collection class, and fulfills a similar function of
storing a collection of sequence features.  However, it also stores
rendering information about the features, and does not currently
follow the CollectionI interface.

=head2 The File Format

There are two types of entry in the file format: feature entries, and
formatting entries.  They can occur in any order.  See the Appendix
for a full example.

Feature entries can take several forms.  At their simplest, they look
like this:

 Gene	B0511.1	516-11208

This means that a feature of type "Gene" and name "B0511.1" occupies
the range between bases 516 and 11208.  A range can be specified
equally well using a hyphen, or two dots as in 516..11208.  Negative
coordinates are allowed, such as -187..1000.

A discontinuous range ("split location") uses commas to separate the
ranges.  For example:

 Gene B0511.1  516-619,3185-3294,10946-11208

Alternatively, the locations can be split by repeating the features
type and name on multiple adjacent lines:

 Gene	B0511.1	516-619
 Gene	B0511.1	3185-3294
 Gene	B0511.1	10946-11208

A comment can be added to features by adding a fourth column.  These
comments will be rendered as under-the-glyph descriptions by those
glyphs that honor descriptions:

 Gene  B0511.1  516-619,3185-3294,10946-11208 "Putative primase"

Columns are separated using whitespace, not (necessarily) tabs.
Embedded whitespace can be escaped using quote marks or backslashes in
the same way as in the shell:

 'Putative Gene' my\ favorite\ gene 516-11208

Features can be grouped so that they are rendered by the "group" glyph
(so far this has only been used to relate 5' and 3' ESTs).  To start a
group, create a two-column feature entry showing the group type and a
name for the group.  Follow this with a list of feature entries with a
blank type.  For example:

 EST	yk53c10
 	yk53c10.3	15000-15500,15700-15800
 	yk53c10.5	18892-19154

This example is declaring that the ESTs named yk53c10.3 and yk53c10.5
belong to the same group named yk53c10.  

=cut

use strict;
use Bio::Graphics::Feature;
use Bio::Graphics::Panel;
use Carp;
use IO::File;
use Text::Shellwords;
use vars '$VERSION';
$VERSION = '1.04';

# default colors for unconfigured features
my @COLORS = qw(cyan blue red yellow green wheat turquoise orange);
use constant WIDTH => 600;

=head2 METHODS

=over 4

=item $features = Bio::Graphics::FeatureFile-E<gt>new(@args)

Create a new Bio::Graphics::FeatureFile using @args to initialize the
object.  Arguments are -name=E<gt>value pairs:

  Argument         Value
  --------         -----

   -file           Read data from a file path or filehandle.  Use
                   "-" to read from standard input.

   -text           Read data from a text scalar.

   -map_coords     Coderef containing a subroutine to use for remapping
                   all coordinates.

   -smart_features Flag indicating that the features created by this
                   module should be made aware of the FeatureFile
		   object by calling their configurator() method.

   -safe           Indicates that the contents of this file is trusted.
                   Any option value that begins with the string "sub {"
                   will be evaluated as a code reference.

The -file and -text arguments are mutually exclusive, and -file will
supersede the other if both are present.

-map_coords points to a coderef with the following signature:

  ($newref,[$start1,$end1],[$start2,$end2]....)
            = coderef($ref,[$start1,$end1],[$start2,$end2]...)

See the Bio::Graphics::Browser (part of the generic genome browser
package) for an illustration of how to use this to do wonderful stuff.

The -smart_features flag is used by the generic genome browser to
provide features with a way to access the link-generation code.  See
gbrowse for how this works.

=back

=cut

# args array:
# -file => parse from a file (- allowed for ARGV)
# -text => parse from a text scalar
# -map_coords => code ref to do coordinate mapping
#                called with ($ref,[$start1,$stop1],[$start2,$stop2]...)
#                returns     ($newref,$new_coord1,$new_coord2...)

sub new {
  my $class = shift;
  my %args  = @_;
  my $self = bless {
		    config   => {},
		    features => {},
		    groups   => {},
		    seenit   => {},
		    types    => [],
		    max      => undef,
		    min      => undef,
		    stat     => [],
		    refs     => {},
                    safe     => undef,
		   },$class;
  $self->{coordinate_mapper} = $args{-map_coords} 
    if exists $args{-map_coords} && ref($args{-map_coords}) eq 'CODE';
  $self->{smart_features}    = $args{-smart_features} if exists $args{-smart_features};
  $self->{safe}              = $args{-safe}           if exists $args{-safe};

  # call with
  #   -file
  #   -text
  my $fh;
  if (my $file = $args{-file}) {
    no strict 'refs';
    if (defined fileno($file)) {
      $fh = $file;
    } elsif ($file eq '-') {
      $self->parse_argv();
    } else {
      $fh = IO::File->new($file) or croak("Can't open $file: $!\n");
    }
    $self->parse_file($fh);
  } elsif (my $text = $args{-text}) {
    $self->parse_text($text);
  }
  close($fh) or warn "Error closing file: $!" if $fh;
  $self;
}

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

#"

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

sub _stat {
  my $self = shift;
  my $fh   = shift;
  $self->{stat} = [stat($fh)];
}

=over 4

=item $error = $features-E<gt>error([$error])

Get/set the current error message.

=back

=cut

sub error {
  my $self = shift;
  my $d = $self->{error};
  $self->{error} = shift if @_;
  $d;
}

=over 4

=item $smart_features = $features-E<gt>smart_features([$flag]

Get/set the "smart_features" flag.  If this is set, then any features
added to the featurefile object will have their configurator() method
called using the featurefile object as the argument.

=back

=cut

sub smart_features {
  my $self = shift;
  my $d = $self->{smart_features};
  $self->{smart_features} = shift if @_;
  $d;
}

sub parse_argv {
  my $self = shift;

  $self->init_parse;
  while (<>) {
    chomp;
    $self->parse_line($_);
  }
  $self->finish_parse;
}

sub parse_file {
  my $self = shift;
  my $fh   = shift or return;
  $self->_stat($fh);

  $self->init_parse;
  while (<$fh>) {
    chomp;
    $self->parse_line($_);
  }
  $self->finish_parse;
}

sub parse_text {
  my $self = shift;
  my $text = shift;

  $self->init_parse;
  foreach (split /\r?\n|\r\n?/,$text) {
    $self->parse_line($_);
  }
  $self->finish_parse;
}

sub parse_line {
  my $self = shift;
  local $_ = shift;

  s/\r//g;  # get rid of carriage returns left over by MS-DOS/Windows systems

  return if /^\s*[\#]/;

  if (/^\s+(.+)/ && $self->{current_tag}) { # continuation line
    my $value = $1;
    my $cc = $self->{current_config} ||= 'general';       # in case no configuration named
    $self->{config}{$cc}{$self->{current_tag}} .= ' ' . $value;
    # respect newlines in code subs
    $self->{config}{$cc}{$self->{current_tag}} .= "\n"
      if $self->{config}{$cc}{$self->{current_tag}}=~ /^sub\s*\{/;
    return;
  }

  if (/^\s*\[([^\]]+)\]/) {  # beginning of a configuration section
    my $label = $1;
    my $cc = $label =~ /^(general|default)$/i ? 'general' : $label;  # normalize
    push @{$self->{types}},$cc unless $cc eq 'general';
    $self->{current_config} = $cc;
    return;
  }

  if (/^([\w: -]+?)\s*=\s*(.*)/) {   # key value pair within a configuration section
    my $tag = lc $1;
    my $cc = $self->{current_config} ||= 'general';       # in case no configuration named
    my $value = defined $2 ? $2 : '';
    $self->{config}{$cc}{$tag} = $value;
    $self->{current_tag} = $tag;
    return;
  }


  if (/^$/) { # empty line
    undef $self->{current_tag};
    return;
  }

  # parse data lines
  my @tokens = eval { shellwords($_||'') };
  unshift @tokens,'' if /^\s+/;

  # close any open group
  undef $self->{grouptype} if length $tokens[0] > 0;

  if (@tokens < 3) {      # short line; assume a group identifier
    $self->{grouptype}     = shift @tokens;
    $self->{groupname}     = shift @tokens;
    return;
  }

  my($ref,$type,$name,$strand,$bounds,$description,$url);

  if (@tokens >= 8) { # conventional GFF file
    my ($r,$source,$method,$start,$stop,$score,$s,$phase,@rest) = @tokens;
    my $group = join ' ',@rest;
    $type   = join(':',$method,$source);
    $bounds = join '..',$start,$stop;
    $strand = $s;
    if ($group) {
      my ($notes,@notes);
      (undef,$self->{groupname},undef,undef,$notes) = split_group($group);
      foreach (@$notes) {
	if (m!^(http|ftp)://!) { $url = $_ } else { push @notes,$_ }
      }
      $description = join '; ',@notes if @notes;
    }
    $name ||= $self->{groupname};
    $ref = $r;
  }

  elsif ($tokens[2] =~ /^([+-.]|[+-]?[01])$/) { # old simplified version
    ($type,$name,$strand,$bounds,$description,$url) = @tokens;
  } else {                              # new simplified version
    ($type,$name,$bounds,$description,$url) = @tokens;
  }

  $type ||= $self->{grouptype} || '';
  $type =~ s/\s+$//;  # get rid of excess whitespace

  # the reference is specified by the GFF reference line first,
  # the last reference line we saw second,
  # or the reference line in the "general" section.
  {
    local $^W = 0;
    $ref  ||= $self->{config}{$self->{current_config}}{'reference'}
      || $self->{config}{general}{reference};
  }
  $self->{refs}{$ref}++ if defined $ref;

  my @parts = map { [/(-?\d+)(?:-|\.\.)(-?\d+)/]} split /(?:,| )\s*/,$bounds;

  foreach (@parts) { # max and min calculation, sigh...
    $self->{min} = $_->[0] if !defined $self->{min} || $_->[0] < $self->{min};
    $self->{max} = $_->[1] if !defined $self->{max} || $_->[1] > $self->{max};
  }

  if ($self->{coordinate_mapper} && $ref) {
    ($ref,@parts) = $self->{coordinate_mapper}->($ref,@parts);
    return unless $ref;
  }

  $type = '' unless defined $type;
  $name = '' unless defined $name;

  # attribute handling
  my %attributes;
  my $score;
  if (defined $description && $description =~ /\w+=\w+/) { # attribute line
    my @attributes = split /;\s*/,$description;
    foreach (@attributes) {
      my ($name,$value) = split /=/,$_,2;
      Bio::Root::Root->throw(qq("$_" is not a valid attribute=value pair)) unless defined $value;
      _unescape($name);
      my @values = split /,/,$value;
      _unescape(@values);
      if ($name =~ /^(note|description)/) {
	$description = "@values";
      } elsif ($name eq 'url') {
	$url = $value;
      } elsif ($name eq 'score') {
	$score = $value;
      } else {
	push @{$attributes{$name}},@values;
      }
    }
  }

  # either create a new feature or add a segment to it
  if (my $feature = $self->{seenit}{$type,$name}) {
    $feature->add_segment(@parts);
  } else {
    $feature = $self->{seenit}{$type,$name} = 
      Bio::Graphics::Feature->new(-name       => $name,
				  -type       => $type,
				  $strand ? (-strand   => make_strand($strand)) : (),
				  defined $score ? (-score=>$score) : (),
				  -segments   => \@parts,
				  -desc       => $description,
				  -ref        => $ref,
				  -attributes => \%attributes,
				  defined($url) ? (-url      => $url) : (),
				 );
    $feature->configurator($self) if $self->smart_features;
    if ($self->{grouptype}) {
      push @{$self->{groups}{$self->{grouptype}}{$self->{groupname}}},$feature;
    } else {
      push @{$self->{features}{$type}},$feature;  # for speed; should use add_feature() instead
    }
  }
}

sub _unescape {
  foreach (@_) {
    tr/+/ /;       # pluses become spaces
    s/%([0-9a-fA-F]{2})/chr hex($1)/g;
  }
  @_;
}

=over 4

=item $features-E<gt>add_feature($feature [=E<gt>$type])

Add a new Bio::FeatureI object to the set.  If $type is specified, the
object will be added with the indicated type.  Otherwise, the
feature's primary_tag() method will be invoked to get the type.

=back

=cut

# add a feature of given type to our list
# we use the primary_tag() method
sub add_feature {
  my $self = shift;
  my ($feature,$type) = @_;
  $type = $feature->primary_tag unless defined $type;
  push @{$self->{features}{$type}},$feature;
}


=over 4

=item $features-E<gt>add_type($type=E<gt>$hashref)

Add a new feature type to the set.  The type is a string, such as
"EST".  The hashref is a set of key=E<gt>value pairs indicating options to
set on the type.  Example:

  $features->add_type(EST => { glyph => 'generic', fgcolor => 'blue'})

When a feature of type "EST" is rendered, it will use the generic
glyph and have a foreground color of blue.

=back

=cut

# Add a type to the list.  Hash values are used for key/value pairs
# in the configuration.  Call as add_type($type,$configuration) where
# $configuration is a hashref.
sub add_type {
  my $self = shift;
  my ($type,$type_configuration) = @_;
  my $cc = $type =~ /^(general|default)$/i ? 'general' : $type;  # normalize
  push @{$self->{types}},$cc unless $cc eq 'general' or $self->{config}{$cc};
  if (defined $type_configuration) {
    for my $tag (keys %$type_configuration) {
      $self->{config}{$cc}{lc $tag} = $type_configuration->{$tag};
    }
  }
}



=over 4

=item $features-E<gt>set($type,$tag,$value)

Change an individual option for a particular type.  For example, this
will change the foreground color of EST features to my favorite color:

  $features->set('EST',fgcolor=>'chartreuse')

=back

=cut

# change configuration of a type.  Call as set($type,$tag,$value)
# $type will be added if not already there.
sub set {
  my $self = shift;
  croak("Usage: \$featurefile->set(\$type,\$tag,\$value\n")
    unless @_ == 3;
  my ($type,$tag,$value) = @_;
  unless ($self->{config}{$type}) {
    return $self->add_type($type,{$tag=>$value});
  } else {
    $self->{config}{$type}{lc $tag} = $value;
  }
}

# break circular references
sub destroy {
  my $self = shift;
  delete $self->{features};
}

sub DESTROY { shift->destroy(@_) }

=over 4

=item $value = $features-E<gt>setting($stanza =E<gt> $option)

In the two-element form, the setting() method returns the value of an
option in the configuration stanza indicated by $stanza.  For example:

  $value = $features->setting(general => 'height')

will return the value of the "height" option in the [general] stanza.

Call with one element to retrieve all the option names in a stanza:

  @options = $features->setting('general');

Call with no elements to retrieve all stanza names:

  @stanzas = $features->setting;

=back

=cut

sub setting {
  my $self = shift;
  if ($self->safe) {
     $self->code_setting(@_);
  } else {
     $self->_setting(@_);
  }
}

# return configuration information
# arguments are ($type) => returns tags for type
#               ($type=>$tag) => returns values of tag on type
sub _setting {
  my $self = shift;
  my $config = $self->{config} or return;
  return keys %{$config} unless @_;
  return keys %{$config->{$_[0]}} if @_ == 1;
  return $config->{$_[0]}{$_[1]}  if @_ > 1;
}


=over 4

=item $value = $features-E<gt>code_setting($stanza=E<gt>$option);

This works like setting() except that it is also able to evaluate code
references.  These are options whose values begin with the characters
"sub {".  In this case the value will be passed to an eval() and the
resulting codereference returned.  Use this with care!

=back

=cut

sub code_setting {
  my $self = shift;
  my $section = shift;
  my $option  = shift;

  my $setting = $self->_setting($section=>$option);
  return unless defined $setting;
  return $setting if ref($setting) eq 'CODE';
  return $setting unless $setting =~ /^sub\s*\{/;
  my $coderef = eval $setting;
  warn $@ if $@;
  $self->set($section,$option,$coderef);
  return $coderef;
}

=over 4

=item $flag = $features-E<gt>safe([$flag]);

This gets or sets and "safe" flag.  If the safe flag is set, then
calls to setting() will invoke code_setting(), allowing values that
begin with the string "sub {" to be interpreted as anonymous
subroutines.  This is a potential security risk when used with
untrusted files of features, so use it with care.

=back

=cut

sub safe {
   my $self = shift;
   my $d = $self->{safe};
   $self->{safe} = shift if @_;
   $self->evaluate_coderefs if $self->{safe} && !$d;
   $d;
}


=over 4

=item @args = $features-E<gt>style($type)

Given a feature type, returns a list of track configuration arguments
suitable for suitable for passing to the
Bio::Graphics::Panel-E<gt>add_track() method.

=back

=cut

# turn configuration into a set of -name=>value pairs suitable for add_track()
sub style {
  my $self = shift;
  my $type = shift;

  my $config  = $self->{config}  or return;
  my $hashref = $config->{$type} or return;

  return map {("-$_" => $hashref->{$_})} keys %$hashref;
}


=over 4

=item $glyph = $features-E<gt>glyph($type);

Return the name of the glyph corresponding to the given type (same as
$features-E<gt>setting($type=E<gt>'glyph')).

=back

=cut

# retrieve just the glyph part of the configuration
sub glyph {
  my $self = shift;
  my $type = shift;
  my $config  = $self->{config}  or return;
  my $hashref = $config->{$type} or return;
  return $hashref->{glyph};
}


=over 4

=item @types = $features-E<gt>configured_types()

Return a list of all the feature types currently known to the feature
file set.  Roughly equivalent to:

  @types = grep {$_ ne 'general'} $features->setting;

=back

=cut

# return list of configured types, in proper order
sub configured_types {
  my $self = shift;
  my $types = $self->{types} or return;
  return @{$types};
}

=over 4

=item  @types = $features-E<gt>types()

This is similar to the previous method, but will return *all* feature
types, including those that are not configured with a stanza.

=back

=cut

sub types {
  my $self = shift;
  my $features = $self->{features} or return;
  return keys %{$features};
}


=over 4

=item @features = $features-E<gt>features($type)

Return a list of all the feature types of type "$type".  If the
featurefile object was created by parsing a file or text scalar, then
the features will be of type Bio::Graphics::Feature (which follow the
Bio::FeatureI interface).  Otherwise the list will contain objects of
whatever type you added with calls to add_feature().

=back

=cut

# return features
sub features {
  my $self = shift;
  return $self->{features}{shift()} if @_;
  return $self->{features};
}



sub make_strand {
  local $^W = 0;
  return +1 if $_[0] =~ /^\+/ || $_[0] > 0;
  return -1 if $_[0] =~ /^\-/ || $_[0] < 0;
  return 0;
}

=over 4

=item @refs = $features-E<gt>refs

Return the list of reference sequences referred to by this data file.

=back

=cut

sub refs {
  my $self = shift;
  my $refs = $self->{refs} or return;
  keys %$refs;
}

=over 4

=item  $min = $features-E<gt>min

Return the minimum coordinate of the leftmost feature in the data set.

=back

=cut

sub min { shift->{min} }

=over 4

=item $max = $features-E<gt>max

Return the maximum coordinate of the rightmost feature in the data set.

=back

=cut

sub max { shift->{max} }

sub init_parse {
  my $s = shift;

  $s->{seenit} = {}; 
  $s->{max}      = $s->{min} = undef;
  $s->{types}    = [];
  $s->{groups}   = {};
  $s->{features} = {};
  $s->{config}   = {}
}

sub finish_parse {
  my $s = shift;
  $s->consolidate_groups;
  $s->evaluate_coderefs if $s->safe;
  $s->{seenit} = {};
  $s->{groups} = {};
}

sub consolidate_groups {
  my $self = shift;
  my $groups = $self->{groups} or return;

  for my $type (keys %$groups) {
    my @groups = values %{$groups->{$type}};
    push @{$self->{features}{$type}},@groups;
  }
}

sub evaluate_coderefs {
  my $self = shift;
  for my $s ($self->_setting) {
    for my $o ($self->_setting($s)) {
      $self->code_setting($s,$o);
    }
  }
}

sub split_group {
  my $group = shift;

  $group =~ s/\\;/$;/g;  # protect embedded semicolons in the group
  $group =~ s/( \"[^\"]*);([^\"]*\")/$1$;$2/g;
  my @groups = split(/\s*;\s*/,$group);
  foreach (@groups) { s/$;/;/g }

  my ($gclass,$gname,$tstart,$tstop,@notes);

  foreach (@groups) {

    my ($tag,$value) = /^(\S+)\s*(.*)/;
    $value =~ s/\\t/\t/g;
    $value =~ s/\\r/\r/g;
    $value =~ s/^"//;
    $value =~ s/"$//;

    # if the tag is "Note", then we add this to the
    # notes array
   if ($tag eq 'Note') {  # just a note, not a group!
     push @notes,$value;
   }

    # if the tag eq 'Target' then the class name is embedded in the ID
    # (the GFF format is obviously screwed up here)
    elsif ($tag eq 'Target' && $value =~ /([^:\"]+):([^\"]+)/) {
      ($gclass,$gname) = ($1,$2);
      ($tstart,$tstop) = /(\d+) (\d+)/;
    }

    elsif (!$value) {
      push @notes,$tag;  # e.g. "Confirmed_by_EST"
    }

    # otherwise, the tag and value correspond to the
    # group class and name
    else {
      ($gclass,$gname) = ($tag,$value);
    }
  }

  return ($gclass,$gname,$tstart,$tstop,\@notes);
}

# create a panel if needed
sub new_panel {
  my $self = shift;

  # general configuration of the image here
  my $width         = $self->setting(general => 'pixels')
                      || $self->setting(general => 'width')
			|| WIDTH;

  my ($start,$stop);
  my $range_expr = '(-?\d+)(?:-|\.\.)(-?\d+)';

  if (my $bases = $self->setting(general => 'bases')) {
    ($start,$stop) =  $bases =~ /([\d-]+)(?:-|\.\.)([\d-]+)/;
  }

  if (!defined $start || !defined $stop) {
    $start = $self->min unless defined $start;
    $stop  = $self->max unless defined $stop;
  }

  my $new_segment = Bio::Graphics::Feature->new(-start=>$start,-stop=>$stop);
  my $panel = Bio::Graphics::Panel->new(-segment   => $new_segment,
					-width     => $width,
					-key_style => 'between');
  $panel;
}

=over 4

=item $mtime = $features-E<gt>mtime

=item $atime = $features-E<gt>atime

=item $ctime = $features-E<gt>ctime

=item $size = $features-E<gt>size

Returns stat() information about the data file, for featurefile
objects created using the -file option.  Size is in bytes.  mtime,
atime, and ctime are in seconds since the epoch.

=back

=cut

sub mtime {
  my $self = shift;
  my $d = $self->{m_time} || $self->{stat}->[9];
  $self->{m_time} = shift if @_;
  $d;
}
sub atime { shift->{stat}->[8];  }
sub ctime { shift->{stat}->[10]; }
sub size  { shift->{stat}->[7];  }

=over 4

=item $label = $features-E<gt>feature2label($feature)

Given a feature, determines the configuration stanza that bests
describes it.  Uses the feature's type() method if it has it (DasI
interface) or its primary_tag() method otherwise.

=back

=cut

sub feature2label {
  my $self = shift;
  my $feature = shift;
  my $type  = eval {$feature->type} || $feature->primary_tag or return;
  (my $basetype = $type) =~ s/:.+$//;
  my $label = $self->type2label($type) || $self->type2label($basetype) || $type;
  $label;
}

=over 4

=item $link = $features-E<gt>make_link($feature)

Given a feature, tries to generate a URL to link out from it.  This
uses the 'link' option, if one is present.  This method is a
convenience for the generic genome browser.

=back

=cut

sub make_link {
  my $self     = shift;
  my $feature  = shift;
  my $label    = $self->feature2label($feature) or return;
  my $link     = $self->setting($label,'link');
  $link        = $self->setting(general=>'link') unless defined $link;
  return unless $link;
  return $self->link_pattern($link,$feature);
}

sub link_pattern {
  my $self = shift;
  my ($pattern,$feature,$panel) = @_;
  $pattern =~ s/\$(\w+)/
    $1 eq 'ref'           ? $feature->location->seq_id
      : $1 eq 'name'      ? $feature->display_name
      : $1 eq 'class'     ? $feature->class
      : $1 eq 'type'      ? $feature->method
      : $1 eq 'method'    ? $feature->method
      : $1 eq 'source'    ? $feature->source
      : $1 eq 'start'     ? $feature->start
      : $1 eq 'end'       ? $feature->end
      : $1 eq 'segstart'  ? $panel->start
      : $1 eq 'segend'    ? $panel->end
      : $1
       /exg;
  return $pattern;
}

# given a feature type, return its label
sub type2label {
  my $self = shift;
  my $type = shift;
  $self->{_type2label} ||= $self->invert_types;
  $self->{_type2label}{$type};
}

sub invert_types {
  my $self = shift;
  my $config  = $self->{config} or return;
  my %inverted;
  for my $label (keys %{$config}) {
    my $feature = $config->{$label}{feature} or next;
    foreach (shellwords($feature||'')) {
      $inverted{$_} = $label;
    }
  }
  \%inverted;
}

=over 4

=item $citation = $features-E<gt>citation($feature)

Given a feature, tries to generate a citation for it, using the
"citation" option if one is present.  This method is a convenience for
the generic genome browser.

=back

=cut

# This routine returns the "citation" field.  It is here in order to simplify the logic
# a bit in the generic browser
sub citation {
  my $self = shift;
  my $feature = shift || 'general';
  return $self->setting($feature=>'citation');
}

=over 4

=item $name = $features-E<gt>name([$feature])

Get/set the name of this feature set.  This is a convenience method
useful for keeping track of multiple feature sets.

=back

=cut

# give this feature file a nickname
sub name {
  my $self = shift;
  my $d = $self->{name};
  $self->{name} = shift if @_;
  $d;
}

1;

__END__

=head1 Appendix -- Sample Feature File

 # file begins
 [general]
 pixels = 1024
 bases = 1-20000
 reference = Contig41
 height = 12

 [Cosmid]
 glyph = segments
 fgcolor = blue
 key = C. elegans conserved regions

 [EST]
 glyph = segments
 bgcolor= yellow
 connector = dashed
 height = 5;

 [FGENESH]
 glyph = transcript2
 bgcolor = green
 description = 1

 Cosmid	B0511	516-619
 Cosmid	B0511	3185-3294
 Cosmid	B0511	10946-11208
 Cosmid	B0511	13126-13511
 Cosmid	B0511	11394-11539
 EST	yk260e10.5	15569-15724
 EST	yk672a12.5	537-618,3187-3294
 EST	yk595e6.5	552-618
 EST	yk595e6.5	3187-3294
 EST	yk846e07.3	11015-11208
 EST	yk53c10
 	yk53c10.3	15000-15500,15700-15800
 	yk53c10.5	18892-19154
 EST	yk53c10.5	16032-16105
 SwissProt	PECANEX	13153-13656	Swedish fish
 FGENESH	Predicted gene 1	1-205,518-616,661-735,3187-3365,3436-3846	Pfam domain
 FGENESH	Predicted gene 2	5513-6497,7968-8136,8278-8383,8651-8839,9462-9515,10032-10705,10949-11340,11387-11524,11765-12067,12876-13577,13882-14121,14169-14535,15006-15209,15259-15462,15513-15753,15853-16219	Mysterious
 FGENESH	Predicted gene 3	16626-17396,17451-17597
 FGENESH	Predicted gene 4	18459-18722,18882-19176,19221-19513,19572-19835	Transmembrane protein
 # file ends

=head1 SEE ALSO

L<Bio::Graphics::Panel>,
L<Bio::Graphics::Glyph>,
L<Bio::Graphics::Feature>,
L<Bio::Graphics::FeatureFile>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut



