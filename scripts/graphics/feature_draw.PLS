#!/usr/bin/perl -w

use strict;
use lib './blib/lib','../blib/lib';
use Bio::Graphics::Panel;
use Bio::Graphics::Feature;
use Bio::Graphics::FeatureFile;

use Getopt::Long;

use constant WIDTH => 600;
my ($WIDTH,$RANGE,$BOXES);

GetOptions ('width:i'  => \$WIDTH,
            'range:s'  => \$RANGE,
	    'boxes'    => \$BOXES,
	   ) || die <<USAGE;
Usage: $0 [feature file 1] [feature file 2] ...

 Options:
    --width  <pixels>     Set width of image (${\WIDTH} pixels default)
    --range  <start-stop> Set range of region (base pairs, start-stop)
    --boxes               Draw grey boxes around the features (for debugging)

Render a Bio::Graphics feature file and produce a PNG image.
See the manual page for Bio::Graphics::FeatureFile for a
description of the file format.
USAGE

my @COLORS = qw(cyan blue red yellow green wheat turquoise orange);  # default colors
my $color = 0;      # position in color cycle

my $data = Bio::Graphics::FeatureFile->new(-file => '-');

# general configuration of the image here
my $width         = $WIDTH || $data->setting(general => 'pixels')
                           || $data->setting(general => 'width')
                           || WIDTH;

my ($start,$stop);
my $range_expr = '(-?\d+)(?:-|\.\.)(-?\d+)';

if (defined $RANGE) {
   ($start,$stop) = $RANGE =~ /$range_expr/o or die "$RANGE: invalid range specification";
} elsif (my $bases = $data->setting(general => 'bases')) {
   ($start,$stop) =  $bases =~ /([\d-]+)(?:-|\.\.)([\d-]+)/;
}

if (!defined $start || !defined $stop) {
       $start = $data->min unless defined $start;
       $stop  = $data->max unless defined $stop;
}

# Use the order of the stylesheet to determine features.  Whatever is left
# over is presented in alphabetic order
my %types = map {$_=>1} $data->configured_types;

my @configured_types   = grep {exists $data->features->{$_}} $data->configured_types;
my @unconfigured_types = sort grep {!exists $types{$_}}      $data->types;

# create the segment,the panel and the arrow with tickmarks
my $segment = Bio::Graphics::Feature->new(-start=>$start,-stop=>$stop);
my $panel = Bio::Graphics::Panel->new(-segment   => $segment,
				      -width     => $width,
				      -key_style => 'between');
$panel->add_track($segment,-glyph=>'arrow',-tick=>2);

my @base_config = $data->style('general');

for my $type (@configured_types,@unconfigured_types) {
  my @config = ( -glyph   => 'segments',         # really generic
		 -bgcolor => $COLORS[$color++ % @COLORS],
		 -label   => 1,
		 -key     => $type,
		 @base_config,             # global
		 $data->style($type),  # feature-specificp
	       );
  my $features = $data->features($type);
  $panel->add_track($features,@config);
}

my $gd = $panel->gd;

if ($BOXES) {  # debugging code
  my $boxes = $panel->boxes;
  debugging_rectangles($gd,$boxes);
}

print $gd->can('png') ? $gd->png : $gd->gif;

sub debugging_rectangles {
  my ($image,$boxes) = @_;
  my $grey = $image->colorClosest(100,100,100);
  foreach (@$boxes) {
    my @rect = @{$_}[1,2,3,4];
    $image->rectangle(@{$_}[1,2,3,4],$grey);
  }
}

=head1 NAME

feature_draw.pl -- Render a Bio::Graphics Feature File

=head1 SYNOPSIS

 feature_draw.pl [options] file.txt [file2.txt...] > rendering.png
 feature_draw.pl [options] file.txt [file2.txt...] | display -

=head1 DESCRIPTION

The feature_draw.pl script is a thin front end around the
Bio::Graphics module.  It accepts a list of files containing sequence
(protein, nucleotide) feature coordinates from the file(s) listed on
the command line or on standard input, renders them, and produces a
PNG file on standard output.

=head2 Options

This script uses GNU-style long options.  This allows you to specify
the image width option, for example, with any of the following alternative
forms:

 --width=800
 --width 800
 -width 800
 -w 800

=over 4

=item --width

This sets the width of the image, in pixels.  The default is 800
pixels.

=item --range

This sets the range of the region displayed, in base pairs from start
to stop. Any of the following formats are accepted:

  --range 1..1000
  --range 1,1000
  --range 1-1000

Negative ranges are allowed.

=back

=head1 Feature Files Format

This script accepts and processes sequence annotations in a simple
tab-delimited format or in GFF format.

The feature file format has a configuration section and a data
section. The configuration section sets up the size and overall
properties of the image, and the data section gives the feature
data itself.

=head2 Configuration Section

If not provided, this scripts generates a reasonable default
configuration section for you, so you do not need to provide a
configuration section to get a reasonable image. However, to tune the
appearance of the image, you will probably want to tweak the
configuration. Here is an excerpt from the configuration section:


 # example file
 [general]
 bases = -1000..21000
 height = 12

 [EST]
 glyph = segments
 bgcolor= yellow
 connector = dashed
 height = 5

 [FGENES]
 glyph = transcript2
 bgcolor = green
 description = 1


The configuration section is divided into a set of sections, each one
labeled with a [section title]. The [general] section specifies global
options for the entire image. Other sections apply to particular
feature types. In the example above, the configuration in the [EST]
section applies to features labeled as ESTs, while the configuration
in the [FGENES] section applies to features labeled as predictions
from the FGENES gene prediction program.

Inside each section is a series of name=value pairs, where the name is
the name of an option to set. You can put whitespace around the = sign
to make it more readable, or even use a colon (:) if you prefer. The
following option names are recognized:

 Option     Value                                       Example
 ------     -----                                       -------

 bases      Min & max of the sequence range (bp)           1200..60000
 width      width of the image (pixels)                    600
 height     Height of each graphical element (pixels)      10
 glyph      Style of each graphical element (see below)    transcript
 fgcolor    Foreground color of each element               yellow
 bgcolor    Background color of each element               blue
 linewidth  Width of lines                                 3
 label      Print the feature's name                       1
 description Whether to print the feature's description    0
 bump       Elements are not allowed to collide            1
 ticks      Print tick marks on arrows                     1
 connector  Type of group connector (dashed, hat or solid) dashed

The "bases" and "width" options are only relevant in the [general]
section. They are overridden by the like-named command-line options.

The remainder of the options can be located in any section, but if
present in the [general] section will set defaults for the others.

Colors are English-language color names or Web-style #RRGGBB colors
(see a book on HTML for an explanation). True/false values are 1 for
true, and 0 for false. Numeric ranges can be expressed in start..end
fashion with two dots, or as start-end with a hyphen.

The "glyph" option controls how the features are rendered. The
following glyphs are implemented:

  Name                Description
  ----                -----------

  box                 A filled rectangle, nondirectional.
  ellipse             An oval.
  arrow               An arrow; can be unidirectional or
		      bidirectional.  It is also capable of displaying
                      a scale with major and minor tickmarks, and can 
                      be oriented horizontally or vertically. 
  segments            A set of filled rectangles connected by solid
		      lines. Used for interrupted features, such as 
		      gapped alignments and exon groups.
  transcript          Similar to segments, but the connecting line is
		      a "hat" shape, and the direction of
		      transcription is indicated by a small arrow. 
  transcript2         Similar to transcript, but the direction of
		      transcription is indicated by a terminal segment
		      in the shape of an arrow. 
  primers             Two inward pointing arrows connected by a line. Used for STSs. 

The bump option is the most important option for controlling the look
of the image. If set to false (the number 0), then the features are
allowed to overlap. If set to true (the number 1), then the features
will move vertically to avoid colliding. If not specified, bump is
turned on if the number of any given type of sequence feature is
greater than 50.

=head2 Data Section

The data section can follow or proceed the configuration section. The
two sections can also be intermixed. The data section is a tab or
whitespace-delimited file which you can export from a spreadsheet
application or word processor file (be sure to save as text only!)

Here is an example data section:


Cosmid     B0511        .       516-619
Cosmid     B0511        .       3185-3294
Cosmid     B0511        .       10946-11208
Cosmid     B0511        .       13126-13511
Cosmid     B0511        .       66-208
Cosmid     B0511        .       6354-6499
Cosmid     B0511        .       13955-14115
EST        yk595e6.5    +       3187-3294
EST        yk846e07.3   -       11015-11208
EST        yk53c10
           yk53c10.5    +       18892-19154
           yk53c10.3    -       15000-15500,15700-15800
EST        yk53c10.5    +       16032-16105
SwissProt  PECANEX      +       13153-13656     Swedish fish
FGENESH    "Gene 1"     -       1-205,518-616,661-735,3187-3365,3436-3846       Transmembrane domain
FGENESH    "Gene 2"     -       16626-17396,17451-17597 Kinase and sushi domains


Each line of the file contains five columns. The columns are: 

 Column #   Description
 --------   -----------

 1          feature type
 2          feature name
 3          strand
 4          coordinates
 5          description

=over 4

=item Feature type

The feature type should correspond to one of the [feature type]
headings in the configuration section. If it doesn't, the [general]
options will be applied to the feature when rendering it. The feature
name is a name for the feature. Use a "." or "-" if this is not
relevant. If the name contains whitespace, put single or double quotes
("") around the name.

=item Strand

The strand indicates which strand the feature is on. It is one of "+"
for the forward strand, "-" for the reverse strand, or "." for
features that are not stranded.

=item Coordinates

The coordinates column is a set of one or more ranges that the feature
occupies. Ranges are written using ".." as in start..stop, or with
hyphens, as in start-stop. For features that are composed of multiple
ranges &em; for example transcripts that have multiple exons &em; you
can either put the ranges on the same line separated by commas or
spaces, or put the ranges on individual lines and just use the same
feature name and type to group them. In the example above, the Cosmid
B0511 features use the individual line style, while the FGENESH
features use the all-ranges-on-one-line style.

=item Description

The last column contains some descriptive text. If the description
option is set to true, this text will be printed underneath the
feature in the rendering.

=back

Finally, it is possible to group related features together. An example
is the ESTs yk53c10.5 and yk53c10.3, which are related by being reads
from the two ends of the clone yk53c10. To indicate this relationship,
generate a section that looks like this:

 EST        yk53c10
            yk53c10.5    +       18892-19154
            yk53c10.3    -       15000-15500,15700-15800


The group is indicated by a line that contains just two columns
containing the feature type and a unique name for the group. Follow
this line with all the features that form the group, but leave the
first column (the feature type) blank. The group will be rendered by
drawing a dashed line between all the members of the group. You can
change this by specifying a different connector option in the
configuration section for this feature type.

=head1 BUGS

Please report them to the author.

=head1 SEE ALSO

L<Bio::Graphics>

=head1 AUTHOR

Lincoln Stein, lstein@cshl.org

=cut


