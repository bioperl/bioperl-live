#!/usr/bin/perl -w
# Generate a simple display of all glyphs for comparison testing
# T. Harris (harris@cshl.org)

# Usage:
# ./all_glyphs GD > all.png
# ./all_glyphs GD 0 1000 > all.png       # output in png with a wide view
# ./all_glyphs GD::SVG 100 150 > all.svg # output in SVG, zoomed

use lib '.','../..','./blib/lib','../../blib/lib','../..';
use strict;
use Bio::Seq;
use Bio::Graphics::Panel;
use Bio::Graphics::Feature;

chomp (my $CLASS = shift);
$CLASS or die "\nUsage: lots_of_glyphs IMAGE_CLASS
\t- where IMAGE_CLASS is one of GD or GD::SVG
\t- GD generate png output; GD::SVG generates SVG.\n";

chomp (my $start = shift);
chomp (my $end   = shift);

$start ||= -100;
$end   ||= 1000;

my $ftr = 'Bio::Graphics::Feature';
my $segment = $ftr->new(-start=>$start,-end=>$end,-name=>'ZK154',-type=>'clone');
my $panel = Bio::Graphics::Panel->new(
				      -grid => [50,100,150,200,250,300,310,320,330],
				      -gridcolor => 'lightcyan',
				      -grid => 1,
				      -segment => $segment,
				    #  -offset => 300,
				    #  -length  => 1000,
				      -spacing => 15,
				      -width   => 600,
				      -pad_top  => 20,
				      -pad_bottom  => 20,
				      -pad_left => 20,
				      -pad_right=> 20,
				    #  -bgcolor => 'teal',
				    #  -key_style => 'between',
				      -key_style => 'bottom',
				      -image_class => $CLASS,
				     );

my $zk154_1 = $ftr->new(-start=>-50,-end=>800,-name=>'ZK154.1',-type=>'gene',-source=>'predicted');
my $zk154_2 = $ftr->new(-start=>380,-end=>500,-name=>'ZK154.2',-type=>'gene',-source=>'predicted');
my $zk154_3 = $ftr->new(-start=>900,-end=>1200,-name=>'ZK154.3',-type=>'gene',-source=>'confirmed');
my $xyz4 = $ftr->new(-segments=>[[40,80],[100,120],[200,280],[300,320]],
		     -name   =>'xyz4',
		     -source =>'mysterious',
		     -subtype=>'predicted',
		     -type   =>'alignment');

# alignment
add_scores($xyz4);
$panel->add_track([$zk154_1,[$zk154_2,$xyz4]],
		  -glyph => 'alignment',
		  -label => 'alignment',
		  -key   => 'alignment',
		  -height => 10,
		  -font => 'gdSmallFont',
		  -bump => 1,
		  -bgcolor   => sub { shift->primary_tag eq 'predicted' ? 'green' : 'blue'},
		  -connector => sub { my $primary_tag = shift->primary_tag;
				      $primary_tag eq 'transcript' ? 'hat'
				    : $primary_tag eq 'alignment'  ? 'solid'
				    : undef},
		  -connector_color => 'black',
		 );

# anchored_arrow
my $short_segment = $ftr->new(-start=>200,-end=>1000);
$panel->add_track($short_segment,
		  -glyph  => 'anchored_arrow',
		  -label  => 'anchored_arrow',
		  -key    => 'anchored_arrow',
		  -double => 1,
		  -bump   => 0,
		  -height => 10,
		  -linewidth =>1,
		  -arrowstyle=>'regular',
		  -tick      =>1,
		 );

# arrow
$panel->add_track($segment,
		  -glyph  => 'arrow',
		  -label  => 'arrow-minor ticks',
		  -key    => 'arrow',
		  -double => 1,
		  -fgcolor=> 'red',
		  -bump   => 0,
		  -height => 10,
		  -arrowstyle=>'regular',
		  -tick      =>1,
		  -linewidth =>1,
		 );

$panel->add_track($segment,
		  -glyph  => 'arrow',
		  -label  => 'arrow-major ticks',
		  -double => 1,
		  -bump   => 0,
		  -height => 10,
		  -linewidth =>1,
		  -arrowstyle=>'filled',
		  -tick=>2,
		 );

# box
my $box = $ftr->new(-start=>100,-end=>600,-name=>'JC8',-type=>'clone');
$panel->add_track($box,
		  -glyph  => 'box',
		  -label  => 'box',
		  -key => 'box',
		  -bump   => 0,
		  -height => 10,
		  -font   => 'gdLargeFont',
		  -linewidth =>1,
		  -bgcolor => 'turquoise',
		  -fgcolor => 'black',
		 );


# cds
my $cds = $ftr->new(-segments=>[[1,50],[100,150],[222,280],[380,400],[520,599],[801,900]],
			     -name=>'cds',
			     -type=>'gene',
			     -strand=>'+1',
			     -subtype=>'predicted',
			    );

my $cds2 = $ftr->new(-segments =>[[23,90],[157,201],[256,375],[439,502],[600,725]],
		     -name     =>'cds',
		     -strand   => '-1',
		     -subtype  =>'predicted',
		     -type     =>'gene');

$panel->add_track([$cds],
		  -glyph    => 'cds',
		  -label    => 'cds',
		  -key      => 'cds',
		  -bump     => 0,
		  -height   =>30,
		  -linewidth=>1,
		  -frame0f => 'blue',
		  -frame1f => 'green',
		  -frame2f => 'yellow',
		  -frame0r => 'red',
		  -frame1r => 'black',
		  -frame2r => 'purple',
		  #  -sixframe => 1,
		  -require_subparts=>1,
		 );


# crossbox
my $crossbox = $ftr->new(-start=>200,-end=>600);
$panel->add_track($crossbox,
		  -glyph  => 'crossbox',
		  -label  => 'crossbox',
		  -key    => 'crossbox',
		  -bump   => 0,
		  -height => 20,
		  # -font   => 'gdMediumBold',
		  -linewidth =>1,
		  -bgcolor => 'red',
		  -fgcolor => 'black',
		 );

# diamond
my $diamonds = $ftr->new(-segments=>[[10,11],[100,101],[201,202],[214,215],[237,238],
				     [300,301],[350,351],[400,550],[601,602],[775,776]],
			 -name=>'SNPs');
$panel->add_track([$diamonds],
		  -glyph   => 'diamond',
		  -label   => 'diamond',
		  -key     => 'diamond',
		  -height  => 10,
		  -bgcolor => 'aqua',
		 );



# dna
my $string =
'tcgtcaaatgtctattgggtcgaaaagaaggtgaacgagtgctcggtgatgcgttcaaaactcaacacaaatcttcacatttcgctccactagtcgactttatcgattttgattatcatgctcaaatgaagatttccaaagaggcaattgtgcagttgaaaaagaaaatgagcccacatatgacaaagcatggatttttctattcaatgggaaaagaaatagtgaaacgacaaactggagtaattcgaacaaattgtctagattgtctggataggacgaatgccgtacaaacagccatcggacttcaaatgtcacatgatcaagttgcatttctgaatttaaacgcgggaaaagtgaatgtagagcaacgagttgaagagattcttcgtgatttgtggcagaaaaatggagatcagtgtagtacgatctacgcgggaactggagctcttgacggaaagagcaagttgaaagacgcgtcgagatcgcttgcaagaactattcagaataatttgatggatggtgcaaagcaggaatcatttgatttatttttgactggagccgcatatgatccgaggcttttcgatagagcatgtaatatattgccacctagtttgatacaagaatacgctgacgccgtatcgcagcttgtcgagcgaagtcccgaaatcgccgaacctcaatccattaaaatattcgttggaacttggaatgtgaatggaggaaagaatattcataatgtggcattccgtaatgaatcgagtctctcccactggatatttgccaattcaatgacacgtctcgtatctgtagaagatgagcaactagctgatattgtagcaattggagttgaagaacttgttgatttgaatgcaagtaatatggttaaagcaagtaccacaaatcaacgaatgtggtgtgaaagtattcgaaaaactctttctgaaaaagctccatttgtgctcattggctccgagcagctcgtcggtgtttgtctattcctcttcgcaagaccacgtgtctcaccatacctgaaagactttgcagtggcttctgtaaagactggaatgggtggagcaactggaaataagggatccgttgccttccgaatcgtcgtattctccacttctatttgttttatttgttctcactttgcagccgggcaaaacgagattcgagacagaaatgaggattttgcgacgacgttgaaaaagattcgattcccgttgggcagagaaattgactcgcatgacgtcatattttggttgggagatttcaactatcgaattaatttgtcgggggatgaagttaagaatgctgttagaaatggagactatgcgaaattagtcgaaaatgatcaattgacacagcagaaagctcttggacagacatttgttggcttcaacgaaggacagctcacgttcgcaccaacatacaaatacgacacattcagtgatgactatgatacgagtgaaaagtgtcgtgcacccgcatggactgatcgaattctttggaaagatcagagaaagaagggaaaaacgcaacttctcagctatgatagatcagaattaaaaacttctgatcatcgacctgttggagctgttttcaaagtggaaacttttaaagttggcggcagaaaatgtgtggagctcatcgaggatgttgtagaatctatgggtccaccggacggaacaatcattgtcagtattgccggaaaacctcgattcccgccgcaaatgtttccgccgattcatgagaagttgaaggaactcggtgctcaagttcagctgagcaaattcgacgatggcgatctatggattgtactgaatagtggagaaatggcattagccgcattaagtatggatgggctgaagatcggaggaacagatcagattaatgtgaagttgaagtcaccggattgggcttatgctttgaagccacatctttcagattttgatttggaatcgtttgaagtgacggcagaggaagaggcattacttggtggtactgatggtgccgtttttgaatttgcagacgaagacgaggacgcaatcagtgtgtctagtctgacgcttactggttcggctcccgatcgacctcgtccaccatcagcaagaagtgaagcgatcagtgtagccaaacttgaatggccaacagaacaaccaaacgtcctctccacatcaatgccaacacgagcttcatcagcttctcttgccaatagttcttggtatgagcatgtaccaccacttgctccacctcaatcaaacaataataaaagccctccacaagcttgtctattcaatccattcactcaatctgcaccatccccggctccaccaccatccacgattcctcttccaccgactcgtggagcatcagttggaccaggtcctccagcggttcccgtcaggaaggcacccccaccgccacctcggcctgtcattccacctagaccaaaaaatatgtag';

my $fragment = Bio::Seq->new(-seq=>$string);
my $dna = $ftr->new(-seq=>$fragment,
		    -start=>$start,-end=>$end);
$panel->add_track($dna,
		  -glyph    => 'dna',
		  -label    => 'dna',
		  -key      => 'dna',
		  -height   => 50,
		  -linewidth=> 1,
		  -axis_color=>'red',
		  -gc_bins  => 10,
		  -strand   => 'both',
		  );


# dot
my $dots = $ftr->new(-segments=>[[10,11],[100,150],[201,232],[214,215],[237,270],
				     [280,281],[300,321],[400,550],[601,602],[775,776]]);
$panel->add_track([$dots],
		  -glyph   => 'dot',
		  -label   => 'dot',
		  -key     => 'dot',
		  -height  => 10,
		  -bgcolor => 'red',
		  -point   => 5,
		 );

# ellipse
my $ellipses = $ftr->new(-segments=>[[100,150],[201,232],[237,270],[300,321],[400,550],[730,776]]);
$panel->add_track([$ellipses],
		  -glyph   => 'ellipse',
		  -label   => 'ellipse',
		  -key     => 'ellipse',
		  -height  => 10,
		  -bgcolor => 'orange',
		 );

# ex
my $ex = $ftr->new(-start=>100,-end=>400);
$panel->add_track($ex,
		  -glyph  => 'ex',
		  -label  => 'ex',
		  -key    => 'ex',
		  -bump   => 0,
		  -height => 20,
		  # -font   => 'gdMediumBold',
		  -linewidth =>1,
		  -bgcolor => 'red',
		  -fgcolor => 'black',
		 );

# graded_segments				
my $partial_gene = $ftr->new(-segments=>[[1,50],[100,150],[220,300],
					 [380,400],[520,600],[800,900]],
			     -name   =>'partial_gene',
			     -strand => '+1',
			     -type   =>'exon',
			     -source =>'confirmed');

add_scores($partial_gene);
$panel->add_track($partial_gene,
		  -glyph     => 'graded_segments',
		  -key       => 'graded_segments',
		  -label     => 'graded_segments - quill connector',,
		  -bgcolor   => 'blue',
		  -connector => 'quill',
		 );

$panel->add_track($partial_gene,
		  -glyph     => 'graded_segments',
		  -label     => 'graded_segments - hat connector',
		  -key       => 'graded_segments',
		  -bgcolor   => 'green',
		  -connector => 'hat',
		 );

$panel->add_track($partial_gene,
		  -glyph     => 'graded_segments',
		  -label     => 'graded_segments - solid connector',
		  -key       => 'graded_segments',
		  -bgcolor   => 'yellow',
		  -connector => 'solid',
		 );

$panel->add_track($partial_gene,
		  -glyph     => 'graded_segments',
		  -label     => 'graded_segments - dashed connector',
		  -key       => 'graded_segments',
		  -bgcolor   => 'red',
		  -connector => 'dashed',
		 );

# heterogenous_segments
$panel->add_track([[$zk154_2,$zk154_3],[$zk154_2,$xyz4]],
		  -glyph => 'heterogeneous_segments',
		  -label => 'heterogeneous_segments',
		  -key   => 'heterogeneous_segments',
		  -height => 10,
		  -bump => 1,
		  -predicted_color=>'orange',
		  -confirmed_color=>'purple',
		  -mysterious_color=>'red',
  		  -connector_color => 'black',
		 );

# line
$panel->add_track($short_segment,
		  -glyph  => 'line',
		  -label  => 'line',
		  -key    => 'line',
		  -bump   => 0,
		  -height => 20,
		  # -font   => 'gdMediumBold',
		  -linewidth =>1,
		  -bgcolor => 'green',
		  -fgcolor => 'black',
		 );

# pinsertion
my $pinsertion = $ftr->new(-segments=>[[10,10],[100,100],[200,200],[300,300],[400,400],
				       [550,600],[650,650]]);
$panel->add_track([$pinsertion],
		  -glyph   => 'pinsertion',
		  -label   => 'pinsertion',
		  -key     => 'pinsertion',
		  -height  => 10,
		  -bgcolor => 'yellow',
		 );

# primers
my $p = $ftr->new(-start=>200,-end=>600);
$panel->add_track($p,
		  -glyph  => 'primers',
		  -label  => 'primers',
		  -key    => 'primers',
		  -height => 10,
		  -linewidth =>1,
		 );

# processed_transcript
my $trans1 = $ftr->new(-start=>50,-end=>10,-name=>'ZK154.1',-type=>"3'-UTR");
my $trans2 = $ftr->new(-start=>100,-end=>50,-name=>'ZK154.2',-type=>'CDS');
my $trans3 = $ftr->new(-start=>350,-end=>225,-name=>'ZK154.3',-type=>'CDS');
my $trans4 = $ftr->new(-start=>650,-end=>500,-name=>'ZK154.3',-type=>'CDS');
my $trans5 = $ftr->new(-start=>700,-end=>650,-name=>'ZK154.3',-type=>"5'-UTR");
my $trans  = $ftr->new(-segments=>[$trans1,$trans2,$trans3,$trans4,$trans5]);
$panel->add_track($trans,
		  -glyph     => 'processed_transcript',
		  -key       => 'processed_transcript',
		  -label     => 'processed_transcript',
		  -bgcolor   => 'aqua',
		  # -height    => 5,
		 );

$panel->add_track($trans,
		  -glyph    => 'processed_transcript',
		  -key      => 'processed_transcript',
		  -label    => 'processed_transcript',
		  -bgcolor  => 'green',
		  # -height   => 10,
		  -thin_utr => 1);


# redgreen_box
$panel->add_track($partial_gene,
		  -glyph     => 'redgreen_box',
		  -label     => 'redgreen_box',
		  -key       => 'redgreen_box',
		 );

# redgreen_segments
$panel->add_track($partial_gene,
		  -glyph     => 'redgreen_segment',
		  -label     => 'redgreen_segment',
		  -key       => 'redgreen_segment',
		 );

# rndrect
$panel->add_track($partial_gene,
		  -glyph     => 'rndrect',
		  -label     => 'rndrect',
		  -key       => 'rndrect',
		 );

# ruler_arrow
$panel->add_track($partial_gene,
		  -glyph => 'ruler_arrow',
		  -label => 1,
		  -key   => 'ruler_arrow',
		  -base  => 1,
		 );

$panel->add_track($partial_gene,
		  -glyph => 'ruler_arrow',
		  -label => 1,
		  -key   => 'ruler_arrow',
		  -base  => 1,
		  -parallel => 0,
		 );

# segments
$panel->add_track([$zk154_1,[$zk154_2,$xyz4]],
		  -glyph => 'segments',
		  -label => 'segments',
		  -key   => 'segments',
		  -height => 10,
		  -bump => 1,
		  -bgcolor   => sub { shift->primary_tag eq 'predicted' ? 'green' : 'blue'},
		  -connector => sub { my $primary_tag = shift->primary_tag;
				      $primary_tag eq 'transcript' ? 'hat'
				    : $primary_tag eq 'alignment'  ? 'solid'
				    : undef},
		  -connector_color => 'black',
		 # -draw_dna => 1,
		 );

# span
my $big_span = $ftr->new(-start=>-400,-end=>3000);
my $small_span = $ftr->new(-start=>290,-end=>600);
$panel->add_track([$big_span,$small_span],
		  -glyph => 'span',
		  -label => 'span',
		  -key   => 'span',
		 );

# splice_site
$panel->add_track($partial_gene,
		  -glyph => 'splice_site',
		  -label => 'splice_site',
		  -key   => 'splice_site',
		  -direction => 'right',
		 );

# transcript
$panel->add_track($trans,
		  -glyph   => 'transcript',
		  -label   => 'transcript',
		  -key     => 'transcript',
		  -bgcolor =>'yellow',
		  -arrow_length=>10,
		 );

# transcript2
$panel->add_track($trans,
		  -glyph   => 'transcript2',
		  -label   => 'transcript2',
		  -key     => 'transcript2',
		  -bgcolor => 'purple',
		  -arrow_length=>10,
		 );

# translation
$panel->add_track($dna,
		  -glyph   => 'translation',
		  -label   => 'translation',
		  -key     => 'translation',
		  -translation => '3frame',
		  -frame0 => 'red',
		  -frame1 => 'blue',
		  -frame2 => 'green',
		  -arrow_length => 10,
		  -start_codons => 1,
		  -show_sequence=> 1,
		 );


# triangle
$panel->add_track([$pinsertion],
		  -glyph   => 'triangle',
		  -label   => 'triangle',
		  -key     => 'triangle',
		  -bgcolor => 'yellow',
		  -point   => 1,
		  -orient  => 'E',
		 );


# xyplot
$panel->add_track($partial_gene,
		  -glyph      => 'xyplot',
		  -key        => 'xyplot',
		  -label      => 'xyplot',
		  -graph_type => 'boxes');


my $gd   = $panel->gd;
my $type = ($CLASS eq 'GD') ? 'png' : 'svg';
print $gd->$type;






sub add_scores {
  my $ftr = shift;
  my $score = 10;
  my @segs = $ftr->segments;
  foreach (@segs) {
    $_->score($score);
    $score += 10;
  }
}


