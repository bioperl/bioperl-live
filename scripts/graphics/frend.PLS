#!/usr/bin/perl -w

=head1 NAME

frend.pl -- Render a Bio::Graphics Feature File on the web

=head1 SYNOPSIS

 http://your.host.com/cgi-bin/frend.pl

=head1 DESCRIPTION

The frend.pl script is a thin front end around the Bio::Graphics
module.  It accepts a list of files containing sequence (protein,
nucleotide) feature coordinates from the file(s) listed on the command
line or on standard input, renders them, and produces a PNG file on
standard output.

=head1 INSTALLATION

Copy this script into your web site's cgi-bin directory.  Name it
whatever you want.

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
 connector = solid
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

L<Bio::Graphics>, L<feature_draw.pl>

=head1 AUTHOR

Lincoln Stein, lstein@cshl.org

=cut

use strict;
use Bio::Graphics::Panel;
use Bio::Graphics::Feature;
use Bio::Graphics::FeatureFile;
use CGI qw(:standard);
use CGI::Carp;
use File::Temp ':mktemp';
use File::Spec;
use File::Basename 'basename';
use File::Path 'mkpath';
use vars '@COLORS';

use constant WIDTH          => 600;  # default width
use constant BUMP_THRESHOLD => 50;  # if more than this # of features, will stop bumping
@COLORS = qw(cyan blue red yellow green wheat turquoise orange);  # default colors

if (param('cat')) {
  catfile(param('cat'));
  exit 0;
}

print header,start_html('Sequence Feature Renderer');
print h1('Sequence Feature Renderer');

print p('This is a front end to the Bio::Graphics package, a part of the',
	a({-href=>'http://www.bioperl.org'},'BioPerl library.'),
	  'Cut and paste your sequence annotation data into the text field below, or upload it using the',
	'upload button.',
	'The format of the annotation data is explained',a({-href=>'#format'},'below.'));

my $self = url(-relative=>1);
print h3('Instant examples'),
  p('For the impatient, you can paste in an',
    b(a({-href=>"$self?Paste+Example+1"},'example file.')));

read_file() if param('file');

my $example = param('Example 1') 
  ? test_data(0) 
  : param('Example 2')
  ? test_data(1)
  : '';
param(text => $example) if length $example;

render() if param('text') || param('file') =~ /\w/;

print start_multipart_form(),
  table({-border=>0,-width=>300,-cellspacing=>0,-cellpadding=>0},
	TR({-class=>'resultsbody'},
	   td({-colspan=>1},
	      'Cut and Paste the annotation file...'
	     ),
	   td({-colspan=>2},
	      'Image width: ',
	      popup_menu(-name=>'width',-values=>[480,640,800,1024,1280,1600],-default=>800)
	     ),
	   TR({-class=>'resultsbody'},
	      td({-colspan=>3},
		 pre(
		     textarea(-name=>'text',-value=>$example,
			      -cols=>80,-rows=>10,-wrap=>'off',-override=>length $example || param('Clear'))
		    )
		)
	     )
	  ),
	TR({-class=>'resultsbody'},
	   td({-colspan=>1},'Upload it... ',filefield(-name=>'file',-size=>30)),
	   td({-align=>'left',-colspan=>2},
	      'Or paste one of the example files...',
	      submit('Example 1'),
	      submit('Example 2'),
	      submit('Clear'),
	     )
	  ),
	TR({-class=>'resultstitle'},
	   td({-align=>'left',-colspan=>3},
	      "Press",b('Render'),'when ready...',
	      b(submit('Render...'))
	     ),
	   )),
  end_form;

print_format();

print hr(),a({-href=>'http://www.bioperl.org'},'www.bioperl.org'),end_html();

exit 0;

sub read_file {
  my $text;
  my $fh = param('file') or return;
  $text .= $_ while <$fh>;
  param(text => $text);
}

sub render {
  my $text = shift;
  my $color = 0;      # position in color cycle

  $text ||= param('text');
  my $data = $text ? Bio::Graphics::FeatureFile->new(-text => $text) 
                   : Bio::Graphics::FeatureFile->new(-file => param('file'));

  unless ($data->min < $data->max) {
    AceError("This doesn't look like a valid annotation file.  No annotations found.");
    exit 0;
  }

  # adjust the width if requested
  $data->setting(general => 'width',param('width')) if param('width');

  # render the panel
  my $panel = $data->new_panel;
  $data->render($panel);

  # we create the file and write it out
  my $gd = $panel->gd;
  my $suffix = $gd->can('gif') ? '.gif' : '.png';
  my $dir  = tmpdir();
  mkpath($dir) unless -e $dir;
  my($fh,$filename) = mkstemps(tmpfile('XXXXXXXX'),$suffix);

  print $fh ($gd->can('gif') ? $gd->gif : $gd->png);
  close $fh;

  # now we send the link to the user
  my $self = url(-relative=>1);
  my $base = basename($filename);
  my $url  = "$self/features$suffix?cat=$base";
  my ($w,$h) = $gd->getBounds;

  print hr(),h2('Rendering');
  print a({-name=>'rendering'},
	  img({-src=>$url,-alt=>'Right-click and "Save As..." to save this image',
	       -border=>0,-width=>$w,-height=>$h})
	  );
}

sub tmpdir {
  return File::Spec->catfile(File::Spec->tmpdir,'frend');
}

sub tmpfile {
  return File::Spec->catfile(tmpdir(),shift);
}

sub catfile {
  my $file = shift;
  my $path = tmpfile($file);
  print header($path =~ /\.gif$/ ? 'image/gif' : 'image/png');
  open F,$path or die "Couldn't open $file for reading: $!";
  print while <F>;
  close F;
  unlink $path;
}

sub print_format {
  print hr();
  print a({-name=>'format'},h2('Annotation file format'));
  print <<END;
<p>
The annotation file format has a configuration section and a data section.  The configuration section
sets up the size and overall properties of the image, and the data section gives the annotation data 
itself.
<p>
<h3>Configuration Section</h3>
<p>
If not provided, this page generates a reasonable default configuration section for you, so you 
do not need to provide a configuration section to get a reasonable image.  However, to tune the
appearance of the image, you will probably want to tweak the configuration.  Here is an excerpt 
from the configuration section:
<blockquote>
<pre>
# example file
[general]
bases = -1000..21000
height = 12

[EST]
glyph = segments
bgcolor= yellow
connector = solid
height = 5

[FGENES]
glyph = transcript2
bgcolor = green
description = 1
</pre>
</blockquote>

<p>
The configuration section is divided into a set of sections, each one labeled with a [section title].
The [general] section specifies global options for the entire image.  Other sections apply to particular
feature types.  In the example above, the configuration in the [EST] section applies to features labeled
as ESTs, while the configuration in the [FGENES] section applies to features labeled as predictions from
the FGENES gene prediction program.
<p>
Inside each section is a series of <i>name</i>=<i>value</i> pairs, where the name is the name of
an option to set.  You can put whitespace around the = sign to make it more readable, or even use
a colon (:) if you prefer.  The following option names are recognized:
<p>
<table border="1">
<tr>
  <th>Option</th><th>Value</th><th>Example</th>
</tr>
<tr>
  <th>bases</th><td>Min &amp; max of the sequence range (bp)</td><td>1200..60000</td>
</tr>
<tr>
  <th>width</th><td>width of the image (pixels)</td>                 <td>600</td>
</tr>
<tr>
  <th>height</th><td>Height of each graphical element (pixels)</td><td>10</td>
</tr>
<tr>
  <th>glyph</th><td>Style of each graphical element (see below)</td><td>transcript</td>
</tr>
<tr>
  <th>fgcolor</th>      <td>Foreground color of each element</td>            <td>yellow</td>
</tr>
<tr>
  <th>bgcolor</th>      <td>Background color of each element</td>            <td>blue</td>
</tr>
<tr>
  <th>linewidth</th>      <td>Width of lines</td>            <td>3</td>
</tr>
<tr>
  <th>label</th>        <td>Print the feature's name</td>         <td>1</td>
</tr>
<tr>
  <th>description</th>  <td>Whether to print the feature's description </td> <td>0</td>
</tr>
<tr>
  <th>bump</th>         <td>Elements are not allowed to collide</td> <td>1</td>
</tr>
<tr>
  <th>ticks</th>        <td>Print tick marks on arrows</td>       <td>1</td>
</tr>
<tr>
  <th>connector</th>    <td>Type of group connector (dashed, hat or solid)</td>       <td>dashed</td>
</tr>
</table>
<p>

The "bases" and "width" options are only relevant in the [general]
section.  The rest can be located in any section, but if present in
the [general] section will set defaults for the others.

<p>

Colors are English-language color names or Web-style #RRGGBB colors
(see a book on HTML for an explanation).  True/false values are 1 for
true, and 0 for false.  Numeric ranges can be expressed in
<i>start</i>..<i>end</i> fashion with two dots, or as
<i>start</i>-<i>end</i> with a hyphen.

<p>
The "glyph" option controls how the features are rendered.  The
following glyphs are implemented:

<p>

<table border="1">

<tr><th>Name</th><th>Description</th></tr>
<tr>
  <th>
  box
  </th>
  <td>A filled rectangle, nondirectional.</td>
</tr>
<tr>
  <th>ellipse</th><td>An oval.</td>
</tr>
<tr>
<th>arrow</th>
<td>	      An arrow; can be unidirectional or bidirectional.
	      It is also capable of displaying a scale with
	      major and minor tickmarks, and can be oriented
	      horizontally or vertically.
</td>
</tr>
<tr>
  <th>segments</th>
  <td>    A set of filled rectangles connected by solid lines.
  Used for interrupted features, such as gapped
  alignments.
</td>
</tr>
<tr>
  <th>transcript</th>
<td>
  Similar to segments, but the connecting line is
  a "hat" shape, and the direction of transcription
  is indicated by a small arrow.
  </td>
</tr>
<tr>
<th>
  transcript2</th>
<td>  Similar to transcript, but the direction of
  transcription is indicated by a terminal segment
  in the shape of an arrow.
</td>
</tr>
<tr>
<th>
  primers
</th>
<td>     Two inward pointing arrows connected by a line.
	      Used for STSs.
</td>
</tr>
</table>
<p>

The <b>bump</b> option is the most important option for controlling the look
of the image.  If set to false (the number 0), then the features are allowed
to overlap.  If set to true (the number 1), then the features will move
vertically to avoid colliding.  If not specified, bump is turned on
if the number of any given type of sequence feature is greater than
${\BUMP_THRESHOLD}.

<h3>Data Section</h3>
<p>

The data section can follow or proceed the configuration section.  The two sections
can also be intermixed.  The data section is a tab or whitespace-delimited file which you can
export from a spreadsheet application or word processor file (be sure to save as text only!)

<p>

Here is an example data section:

<p>

<blockquote>
<pre>
Cosmid	   B0511	.	516-619
Cosmid	   B0511	.	3185-3294
Cosmid	   B0511	.	10946-11208
Cosmid	   B0511	.	13126-13511
Cosmid	   B0511	.	66-208
Cosmid	   B0511	.	6354-6499
Cosmid	   B0511	.	13955-14115
EST	   yk595e6.5	+	3187-3294
EST	   yk846e07.3	-	11015-11208
EST	   yk53c10
	   yk53c10.5	+	18892-19154
	   yk53c10.3	-	15000-15500,15700-15800
EST	   yk53c10.5	+	16032-16105
SwissProt  PECANEX	+	13153-13656	Swedish fish
FGENESH	   "Gene 1"	-	1-205,518-616,661-735,3187-3365,3436-3846	Transmembrane domain
FGENESH	   "Gene 2"	-	16626-17396,17451-17597	Kinase and sushi domains
</pre>
</blockquote>

<p>

Each line of the file contains five columns.  The columns are:

<p>

<table border="1">
<tr><th>Column #</th><th>Column Description</th></tr>
<tr><td align="right">1</td><td>feature type</td></tr>
<tr><td align="right">2</td><td>feature name</td></tr>
<tr><td align="right">3</td><td>strand</td></tr>
<tr><td align="right">4</td><td>coordinates</td></tr>
<tr><td align="right">5</td><td>description</td></tr>
</table>
<p>

The <b>feature type</b> should correspond to one of the [feature type] headings
in the configuration section.  If it doesn't, the [general] options will
be applied to the feature when rendering it.  The <b>feature name</b> is a
name for the feature.  Use a "." or "-" if this is not relevant.  If
the name contains whitespace, put single or double quotes ("") around
the name.

<p>

The <b>strand</b>
indicates which strand the feature is on.  It is one of "+" for the 
forward strand, "-" for the reverse strand, or "." for features that are not
stranded.  

<p>

The <b>coordinates</b> column is a set of one or more ranges that the
feature occupies.  Ranges are written using ".." as in <i>start</i>..<i>stop</i>,
or with hyphens, as in <i>start</i>-<i>stop</i>. For features that are composed
of multiple ranges &em; for example transcripts that have multiple exons &em;
you can either put the ranges on the same line separated by commas or spaces,
or put the ranges on individual lines and just use the same feature name and
type to group them.  In the example above, the Cosmid B0511 features use
the individual line style, while the FGENESH features use the all-ranges-on-one-line
style.

<p>

The last column contains some descriptive text.  If the <b>description</b> option
is set to true, this text will be printed underneath the feature in the rendering.

<p>

Finally, it is possible to group related features together.  An example is
the ESTs yk53c10.5 and yk53c10.3, which are related by being reads from 
the two ends of the clone yk53c10.  To indicate this relationship, generate
a section that looks like this:

<p>

<blockquote>
<pre>
EST	   yk53c10
	   yk53c10.5	+	18892-19154
	   yk53c10.3	-	15000-15500,15700-15800
</pre>
</blockquote>

<p>

The group is indicated by a line that contains just two columns 
containing the feature type and a unique name for the group.
Follow this line with all 
the features that form the group, but leave the first column 
(the feature type) blank.  The group will be rendered by
drawing a dashed line between all the members of the group.  
You can change this by specifying a different <b>connector</b>
option in the configuration section for this feature type.

END
;

}

sub test_data {
  my $config = shift;
  my $header = <<'END';
[general]
bases = -1000..21000
height = 12
reference = B0511

[Cosmid]
glyph = segments
fgcolor = blue
key = C. elegans conserved regions

[EST]
glyph = segments
bgcolor= yellow
connector = solid
height = 5

[FGENESH]
glyph = transcript2
bgcolor = green
description = 1

[SwissProt]
glyph = arrow
base  = 1
linewidth = 2
fgcolor = red
description = 1

[P-element]
glyph = triangle
orient = S
bgcolor = red
fgcolor = white
label = 1
point = 1

END
;

my $data =<<'END';
Cosmid	B0511	516-619
Cosmid	B0511	3185-3294
Cosmid	B0511	10946-11208
Cosmid	B0511	13126-13511
Cosmid	B0511	11394-11539
Cosmid	B0511	14383-14490
Cosmid	B0511	15569-15755
Cosmid	B0511	18879-19178
Cosmid	B0511	15850-16110
Cosmid	B0511	66-208
Cosmid	B0511	6354-6499
Cosmid	B0511	13955-14115
Cosmid	B0511	7985-8042
Cosmid	B0511	11916-12046
P-element	""	500-500
P-element	MrQ	700-700
P-element	MrR	10000-10000
EST	yk260e10.5	15569-15724
EST	yk672a12.5	537-618,3187-3294
EST	yk595e6.5	552-618
EST	yk595e6.5	3187-3294
EST	yk846e07.3	11015-11208
EST	yk53c10
	yk53c10.3	12876-13577,13882-14121,14169-14535
	yk53c10.5	18892-19154,15853-16219
SwissProt	"PECANEX Protein"	5513-16656	"From SwissProt"
FGENESH	"Predicted gene 1"	-1200--500,518-616,661-735,3187-3365,3436-3846	Pfam domain
FGENESH	"Predicted gene 2"	5513-6497,7968-8136,8278-8383,8651-8839,9462-9515,10032-10705,10949-11340,11387-11524,11765-12067,12876-13577,13882-14121,14169-14535,15006-15209,15259-15462,15513-15753,15853-16219	Mysterious
FGENESH	"Predicted gene 3"	16626-17396,17451-17597
FGENESH	"Predicted gene 4"	18459-18722,18882-19176,19221-19513,19572-30000	"Transmembrane protein"
END

  return $config ? $header . $data : $data;
}

__END__

