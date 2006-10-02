package Bio::Graphics::Glyph::image;

# $Id$

use strict;
use GD;
use base 'Bio::Graphics::Glyph::generic';
our @ISA;

#
#       |--------------------| true position  ('height' high)
#       .                    .
#      .                      .     diagonal  (vertical spacing high)
#     .                        .
#    +--------------------------+
#    |                          |
#    |                          |
#    |                          |   image
#    |                          |
#    |                          |
#    |                          |
#    +--------------------------+

use constant VERTICAL_SPACING => 20;

sub new {
  my $self  = shift->SUPER::new(@_);
  $self->{image} = $self->get_image();
  return $self;
}

sub get_image {
  my $self    = shift;
  my ($format,$image)   = eval { $self->image_data };
  unless ($image) {
    warn $@ if $@;
    return;
  }
  my $gd      =   $format eq 'image/png'  ? GD::Image->newFromPngData($image,1)
                : $format eq 'image/jpeg' ? GD::Image->newFromJpegData($image,1)
		: $format eq 'image/gif'  ? GD::Image->newFromGifData($image)
		: $format eq 'image/gd'   ? GD::Image->newFromGdData($image)
		: $format eq 'image/gd2'  ? GD::Image->newFromGd2Data($image)
		: $self->throw("This module cannot handle images of type $format");
  return $gd;
}

sub _guess_format {
  my $self = shift;
  my $path = shift;
  return 'image/png'   if $path =~ /\.png$/i;
  return 'image/jpeg'  if $path =~ /\.jpe?g$/i;
  return 'image/gif'   if $path =~ /\.gif(87)?$/i;
  return 'image/gd'    if $path =~ /\.gd$/i;
  return 'image/gd2'   if $path =~ /\.gd2$/i;
  my ($extension) = $path =~ /\.(\w+)$/;  #cop-out
  return $extension;
}

sub image_path {
  my $self = shift;
  my $feature  = $self->feature  or $self->throw("no feature!");
  my $dirname  = $self->image_dir;
  my $basename = $self->option('image');

  # can't get it from callback, so try looking for an 'image' attribute
  if (!$basename && $feature->can('has_tag') && $feature->has_tag('image')) {
    ($basename)  = $feature->get_tag_values('image');
  }

  return unless $basename;
  return $basename             if $basename =~ m!^\w+:/!;  # looks like a URL
  return $basename             if $basename =~ m!^/!;      # looks like an abs path
  return "$dirname/$basename";
}

sub image_data {
  my $self = shift;
  my $path = $self->image_path;

  if ($path =~ m!^\w+:/!) { # looks like a URL
    require LWP::UserAgent;
    my $ua = LWP::UserAgent->new(env_proxy => 1);
    my $response = $ua->get($path);
    if ($response->is_success) {
      return ($response->content_type,$response->content);
    } else {
      $self->throw($response->status_line);
    }


  } else {
    my $content_type = $self->_guess_format($path);
    open F,$path or $self->throw("Can't open $path: $!");
    binmode F;
    my $data;
    $data .= $_ while read(F,$_,1024);
    close F;
    return ($content_type,$data);
  }
}

sub pad_left {
  my $self = shift;
  my $pad          = $self->SUPER::pad_left;
  my $image        = $self->{image} or return $pad;
  my $width_needed = ($image->width - $self->width)/2;
  return $pad > $width_needed ? $pad : $width_needed;
}

sub pad_right {
  my $self = shift;
  my $pad          = $self->SUPER::pad_right;
  my $image        = $self->{image} or return $pad;
  my $width_needed = ($image->width - $self->width)/2;
  return $pad > $width_needed ? $pad : $width_needed;
}

sub pad_bottom {
  my $self   = shift;
  my $pb     = $self->SUPER::pad_bottom;
  my $image  = $self->{image} or return $pb;
  $pb       += $self->vertical_spacing;
  $pb       += $image->height;
  return $pb;
}

sub vertical_spacing {
  my $self  = shift;
  my $vs    = $self->option('vertical_spacing');
  return $vs if defined $vs;
  return VERTICAL_SPACING;
}

sub draw_description {
  my $self = shift;
  my ($gd,$left,$top,$partno,$total_parts) = @_;
  $top += $self->{image}->height+$self->vertical_spacing if $self->{image};
  $self->SUPER::draw_description($gd,$left,$top,$partno,$total_parts);
}

sub image_dir {
  my $self = shift;
  return $self->option('image_prefix');
}

sub draw_component {
  my $self  = shift;
  my $gd    = shift;
  my($x1,$y1,$x2,$y2) = $self->bounds(@_);

  my $delegate = $self->option('glyph_delegate') || 'generic';
  if ($delegate eq 'generic') {
    $self->SUPER::draw_component($gd,@_);
  } else {
    eval "require Bio::Graphics::Glyph::$delegate";
    local @ISA = ("Bio::Graphics::Glyph::$delegate");
    my $method = "Bio::Graphics::Glyph::${delegate}::draw_component";
    $self->$method($gd,@_);
  }

  my $image  = $self->{image} or return;

  my $fgcolor = $self->fgcolor;
  my $bgcolor = $self->bgcolor;
  my $height  = $self->option('height');
  my $half    = 4;
  my $vs      = $self->vertical_spacing;

  my $delta = (($x2-$x1) - $image->width)/2;
  my($x,$y) = ($x1+$delta,$y1+$vs+$self->height);
  if ($gd->can('copy')) {
    $gd->copy($image,$x,$y,0,0,$image->width,$image->height) ;
  } else {
    my $gray = $self->panel->translate_color('gray');
    $gd->filledRectangle($x,$y,$x+$image->width,$y+$image->height,$gray);
  }

  if ($vs > 0) {
    $gd->line($x1,$y2+2,$x1,$y2+$half,$fgcolor);
    $gd->line($x2,$y2+2,$x2,$y2+$half,$fgcolor);
    $gd->line($x1,$y2+$half,$x,$y-$half,$fgcolor);
    $gd->line($x2,$y2+$half,$x+$image->width-1,$y-$half,$fgcolor);
    $gd->line($x,$y-$half,$x,$y-2,$fgcolor);
    $gd->line($x+$image->width-1,$y-$half,$x+$image->width-1,$y-2,$fgcolor);
  }
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::image - A glyph that draws photographs & other images

=head1 SYNOPSIS

 use Bio::Graphics;
 use Bio::Seq;
 use Bio::SeqFeature::Generic;

 my $bsg = 'Bio::SeqFeature::Generic';

 my $seq    = Bio::Seq->new(-length=>1000);

 my $whole  = $bsg->new(-display_name => 'Clone82',
 		        -start        => 1,
		        -end          => $seq->length);

 my $image1 = $bsg->new(-start        => 100,
		        -end          => 300,
		        -display_name => 'Excretory System',
		        -tag=>{
			      image=>"http://www.flybase.org/anatomy/image-browser_files/excretory-system.gif"
			      }
		       );

 my $image2 = $bsg->new(-start        => 500,
		        -end          => 800,
		        -display_name => 'Expression Pattern',
		        -tag=>{
			      image=>"http://www.flybase.org/anatomy/image-browser_files/embryonic-expression-pattern.gif"
			      }
		       );

 my $panel = Bio::Graphics::Panel->new(-length    => $seq->length,
				       -width     => 800,
				       -truecolor => 1,
				       -key_style => 'between',
				       -pad_left  => 10,
				       -pad_right => 10,
				      );

 $panel->add_track($whole,
		   -glyph    => 'arrow',
		   -double   => 1,
		   -tick     => 2,
		   -label    => 1,
		   );

 $panel->add_track([$image1,$image2],
		   -glyph    => 'image',
		   -label    => 1,
		   -key       => 'Example images');

 binmode STDOUT;
 print $panel->png;

=head1 DESCRIPTION

This glyph inserts an image into the track at the indicated feature
coordinates. The image can be in PNG, JPEG, GIF or GD format, and can
be either 8-bit or 24-bit ("truecolor"). The image can be located on
the local filesystem or located at a remote URL (provided that you
have the LWP module installed).

When working with photographic images, you may wish to have
Bio::Graphics::Panel create 24-bit (truecolor) images in order to
avoid running out of colors. The symptom of this is that images appear
posterized. To turn on truecolor images, pass the -truecolor option to
Bio::Graphics::Panel as shown in the synopsis.

=head2 OPTIONS

The following options are standard among all Glyphs.  See
L<Bio::Graphics::Glyph> for a full explanation.

  Option      Description                      Default
  ------      -----------                      -------

  -fgcolor      Foreground color	       black

  -outlinecolor	Synonym for -fgcolor

  -bgcolor      Background color               turquoise

  -fillcolor    Synonym for -bgcolor

  -linewidth    Line width                     1

  -height       Height of glyph		       10

  -font         Glyph font		       gdSmallFont

  -connector    Connector type                 0 (false)

  -connector_color
                Connector color                black

  -label        Whether to draw a label	       0 (false)

  -description  Whether to draw a description  0 (false)

  -hilite       Highlight color                undef (no color)

The following additional options are available to the "image" glyph:

  Option            Description                       Default
  ------            -----------                       -------

  -image            Specify the image path or URL     none
                    to use for this feature.

  -image_prefix     String to prepend to              none
                    each image path. You may prepend
                    a directory or a partial URL.

  -vertical_spacing Vertical distance from the box    20
                    that shows the physical span of
                    of the feature to the top of
                    the picture (in pixels).

  -glyph_delegate   Glyph to use for the part of      'generic'
                    the glyph that shows the physical
                    span of the feature.

Set B<-vertical_spacing> to 0 to completely suppress the diagonal
lines that connect the physical span of the feature to the image.

=head2 Specifying the Image

The path to the image can be specified in two ways. First, you can
place it in the feature itself using a tag named "image". Second, you
can specify it as a track option using a callback:

  $panel->add_track(\@features,
                    -glyph=>'image',
                    -image => sub { my $feature = shift;
                                    my $image_path = do_something();
                                    return $image }
                    );

You can of course give -image a constant string, in which case each
feature will show the same image.

The image can be a file on the local operating system or a
URL. However, URL fetching will only work if the LWP module is
installed on your system. Otherwise the glyph will fail with an error
message.

If the image is a relative path (it does not begin with a slash or a
URL protocol), then the contents of -image_prefix will be prepended to
it. This allows you to specify images that are relative to a
particular directory or a partial URL. Example:

  $panel->add_track(\@features,
                    -glyph => 'image',
                    -image_prefix => 'http://www.flybase.org/anatomy/image-browser_files',
                   );

This specifies that each feature's "image" tag is to be appended to
the partial FlyBase URL, thereby saving space.

=head2 Glyph Delegation

The image glyph consists of two parts: an upper part that shows the
extent of the feature in base pair coordinates, and a lower part that
shows the image. No scaling of the image is done; its height and width
are fixed.

By default the upper part uses the "generic" glyph, which is a simple
rectangle filled with the bgcolor and outlined with the fgcolor. To
use a different glyph in the upper part, specify the -glyph_delegate
option, giving the name of the glyph you wish to use. For instance, to
use the "span" glyph:

  $panel->add_track(\@features,
                    -glyph          => 'image',
                    -glyph_delegate => 'span'
                   );

This feature does not work with all glyphs, and in particular requires
a recent CVS checkout of Bio::Perl to work properly with the "arrow",
"span" and "primers" glyphs (support for the feature did not make it
into version 1.5).

=head1 BUGS AND LIMITATIONS

This glyph does not work with GD::SVG. If you try to render it onto a
GD::SVG panel, the image will be shown as a gray box. This will be
fixed in a future version of GD::SVG.

=head1 SEE ALSO

L<Bio::Graphics::Panel>,
L<Bio::Graphics::Glyph>,
L<Bio::Graphics::Glyph::arrow>,
L<Bio::Graphics::Glyph::cds>,
L<Bio::Graphics::Glyph::crossbox>,
L<Bio::Graphics::Glyph::diamond>,
L<Bio::Graphics::Glyph::dna>,
L<Bio::Graphics::Glyph::dot>,
L<Bio::Graphics::Glyph::ellipse>,
L<Bio::Graphics::Glyph::extending_arrow>,
L<Bio::Graphics::Glyph::generic>,
L<Bio::Graphics::Glyph::graded_segments>,
L<Bio::Graphics::Glyph::heterogeneous_segments>,
L<Bio::Graphics::Glyph::line>,
L<Bio::Graphics::Glyph::pinsertion>,
L<Bio::Graphics::Glyph::primers>,
L<Bio::Graphics::Glyph::rndrect>,
L<Bio::Graphics::Glyph::segments>,
L<Bio::Graphics::Glyph::ruler_arrow>,
L<Bio::Graphics::Glyph::toomany>,
L<Bio::Graphics::Glyph::transcript>,
L<Bio::Graphics::Glyph::transcript2>,
L<Bio::Graphics::Glyph::translation>,
L<Bio::Graphics::Glyph::triangle>,
L<Bio::DB::GFF>,
L<Bio::SeqI>,
L<Bio::SeqFeatureI>,
L<Bio::Das>,
L<GD>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>, Todd Harris E<lt>harris@cshl.orgE<gt>

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
