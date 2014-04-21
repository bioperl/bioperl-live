#
# BioPerl module for Cladogram
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Gabriel Valiente <valiente@lsi.upc.edu>
#
# Copyright Gabriel Valiente
#
# You may distribute this module under the same terms as Perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::Draw::Cladogram - Drawing phylogenetic trees in
Encapsulated PostScript (EPS) format.

=head1 SYNOPSIS

  use Bio::Tree::Draw::Cladogram;
  use Bio::TreeIO;
  my $treeio = Bio::TreeIO->new('-format' => 'newick',
  			       '-file'   => 'input.nwk');
  my $t1 = $treeio->next_tree;
  my $t2 = $treeio->next_tree;

  my $obj1 = Bio::Tree::Draw::Cladogram->new(-tree => $t1);
  $obj1->print(-file => 'cladogram.eps');

  if ($t2) {
    my $obj2 = Bio::Tree::Draw::Cladogram->new(-tree => $t1, -second => $t2);
    $obj2->print(-file => 'tanglegram.eps');
  }

=head1 DESCRIPTION

Bio::Tree::Draw::Cladogram is a Perl tool for drawing Bio::Tree::Tree
objects in Encapsulated PostScript (EPS) format. It can be utilized
both for displaying a single phylogenetic tree (a cladogram) and for
the comparative display of two phylogenetic trees (a tanglegram) such
as a gene tree and a species tree, a host tree and a parasite tree,
two alternative trees for the same set of taxa, or two alternative
trees for overlapping sets of taxa.

Phylogenetic trees are drawn as rectangular cladograms, with
horizontal orientation and ancestral nodes centered over their
descendents. The font used for taxa is Courier at 10 pt. A single
Bio::Tree::Tree object is drawn with ancestors to the left and taxa
flushed to the right. Two Bio::Tree::Tree objects are drawn with the
first tree oriented left-to-right and the second tree oriented
right-to-left, and with corresponding taxa connected by straight lines
in a shade of gray. Each correspondence between a $taxon1 of the first
tree and a $taxon2 of the second tree is established by setting
$taxon1-E<gt>add_tag_value('connection',$taxon2). Thus, a taxon of the
first tree can be connected to more than one taxon of the second tree,
and vice versa.

The branch from the parent to a child $node, as well as the child
label, can be colored by setting $node-E<gt>add_tag_value('Rcolor',$r),
$node-E<gt>add_tag_value('Gcolor',$g), and
$node-E<gt>add_tag_value('Bcolor',$b), where $r, $g, and $b are the
desired values for red, green, and blue (zero for lowest, one for
highest intensity).

This is a preliminary release of Bio::Tree::Draw::Cladogram. Future
improvements include an option to output phylograms instead of
cladograms. Beware that cladograms are automatically scaled according
to branch lengths, but the current release has only been tested with
trees having unit branch lengths.

The print method could be extended to output graphic formats other
than EPS, although there are many graphics conversion programs around
that accept EPS input. For instance, most Linux distributions include
epstopdf, a Perl script that together with Ghostscript, converts EPS
to PDF.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Gabriel Valiente

Email valiente@lsi.upc.edu

Code for coloring branches contributed by Georgii A Bazykin
(gbazykin@princeton.edu).

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Tree::Draw::Cladogram;
use strict;

use PostScript::TextBlock;

use base qw(Bio::Root::Root);

# The following private package variables are set by the new method
# and used by the print method.

my %xx;        # horizontal coordinate for each node
my %yy;        # vertical coordinate for each node
my $t1;        # first Bio::Tree::Tree object
my $t2;        # second Bio::Tree::Tree object
my $font;      # font name
my $size;      # font size
my $width;     # total drawing width
my $height;    # total drawing height
my $xstep;     # branch length in drawing
my $tip;       # extra space between tip and label
my $tipwidth1; # width of longest label among $t1 taxa
my $tipwidth2; # width of longest label among $t2 taxa
my $compact;   # whether or not to ignore branch lengths
my $ratio;     # horizontal to vertical ratio
my $colors;    # use colors to color edges
my %Rcolor;    # red color for each node
my %Gcolor;    # green color for each node
my %Bcolor;    # blue color for each node
my $bootstrap; # Draw Bootstrap boolean

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tree::Draw::Cladogram->new();
 Function: Builds a new Bio::Tree::Draw::Cladogram object 
 Returns : Bio::Tree::Draw::Cladogram
 Args    : -tree => Bio::Tree::Tree object
           -second => Bio::Tree::Tree object (optional)
           -font => font name [string] (optional)
           -size => font size [integer] (optional)
           -top => top margin [integer] (optional)
           -bottom => bottom margin [integer] (optional)
           -left => left margin [integer] (optional)
           -right => right margin [integer] (optional)
           -tip => extra tip space [integer] (optional)
           -column => extra space between cladograms [integer] (optional)
           -compact => ignore branch lengths [boolean] (optional)
           -ratio => horizontal to vertical ratio [integer] (optional)
           -colors => use colors to color edges [boolean] (optional)
           -bootstrap => draw bootstrap or internal ids [boolean]

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  ($t1, $t2, $font, $size, my $top, my $bottom, my $left, my $right,
    $tip, my $column, $compact, $ratio, $colors,$bootstrap) = 
	$self->_rearrange([qw(TREE SECOND FONT SIZE TOP BOTTOM LEFT RIGHT 
			      TIP COLUMN COMPACT RATIO COLORS BOOTSTRAP)], 
				   @args);
  $font ||= "Helvetica-Narrow";
  $size ||= 12;
  $top ||= 10;
  $bottom ||= 10;
  $left ||= 10;
  $right ||= 10;
  $tip ||= 5;
  $column ||= 60;
  $compact ||= 0;
  $ratio ||= 1 / 1.6180339887;
  $colors ||= 0;
  $bootstrap ||= 0;

  # Roughly, a cladogram is set according to the following parameters.

  #################################
  #                           # T #   $top (T, top margin)
  #        +---------+ XXX    #   #   $bottom (B, bottom margin)
  #        |                  #   #   $left (L, left margin)
  #        |                  #   #   $right (R, right margin)
  #   +----+                  #   #   $tip (X, extra tip space)
  #        |    +----+ XXXX   #   #   $width (total drawing width)
  #        |    |             #   #   $height (total drawing height)
  #        +----+             # Y #   $xstep (S, stem length)
  #             |             #   #   $ystep (Y, space between taxa)
  #             +----+ XX     #   #   $tiplen (string length of longest name)
  #                           # B #   $tipwidth (N, size of longest name)
  #################################
  # L         S       X  N  R #
  #############################

  # A tanglegram is roughly set as follows. The only additional
  # parameter is $column (C, length of connection lines between taxa
  # of the two trees), but $tip occurs four times, and $tiplen and
  # $tipwidth differ for the first and the second tree.

  ###########################################################
  #                                                         #
  #        +---------+ XXX  ----- XXXXXX +----+             #
  #        |                                  |             #
  #        |                                  +----+        #
  #   +----+                                  |    |        #
  #        |    +----+ XXXX -----    XXX +----+    |        #
  #        |    |                                  +----+   #
  #        +----+                                  |        #
  #             |                                  |        #
  #             +----+ XX   -----   XXXX +---------+        #
  #                                                         #
  ###########################################################
  # L                 X    X  C  X      X                 R #
  ###########################################################

  # An alternative would be to let the user set $width and $height in
  # points and to scale down everything to fit the desired
  # dimensions. However, the final EPS can later be scaled down to any
  # desired size anyway.

  my @taxa1 = $t1->get_leaf_nodes;
  my $root1 = $t1->get_root_node;

  $tipwidth1 = 0;
  foreach my $taxon (@taxa1) {
    my $w = PostScript::Metrics::stringwidth($taxon->id,$font,$size);
    if ($w > $tipwidth1) { $tipwidth1 = $w; }
  }

  my @taxa2;
  my $root2;

  my $ystep = 20;

  if ($t2) {
    @taxa2 = $t2->get_leaf_nodes;
    $root2 = $t2->get_root_node;
    $tipwidth2 = 0;
    foreach my $taxon (@taxa2) {
      my $w = PostScript::Metrics::stringwidth($taxon->id,$font,$size);
      if ($w > $tipwidth2) { $tipwidth2 = $w; }
    }
  }

  my $stems = $root1->height + 1;
  if ($t2) { $stems += $root2->height + 1; }
  my $labels = $tipwidth1;
  if ($t2) { $labels += $tipwidth2; }
  $xstep = 20;
  $width = $left + $stems * $xstep + $tip + $labels + $right;
  if ($t2) { $width += $tip + $column + $tip + $tip; }
  $height = $bottom + $ystep * (@taxa1 - 1) + $top;
  if ($t2) {
    if ( scalar(@taxa2) > scalar(@taxa1) ) {
      $height = $bottom + $ystep * (@taxa2 - 1) + $top;
    }
  }
  my $ystep1 = $height / scalar(@taxa1);
  my $ystep2;
  if ($t2) {
    $ystep2 = $height / scalar(@taxa2);
  }

  my $x = $left + $xstep * ($root1->height + 1) + $tip;
  my $y = $bottom;

  for my $taxon (reverse @taxa1) {
    $xx{$taxon} = $x - $tip;
    $yy{$taxon} = $y;
    $y += $ystep1;
  }
  $x -= $xstep;

  my @stack;
  my @queue; # postorder traversal
  push @stack, $t1->get_root_node;
  while (@stack) {
    my $node = pop @stack;
    push @queue, $node;
    foreach my $child ($node->each_Descendent(-sortby => 'internal_id')) {
      push @stack, $child;
    }
  }
  @queue = reverse @queue;

  for my $node (@queue) {
    if (!$node->is_Leaf) {
      my @children = $node->each_Descendent;
      my $child = shift @children;
      my $xmin = $xx{$child};
      my $ymin = my $ymax = $yy{$child};
      foreach $child (@children) {
	$xmin = $xx{$child} if $xx{$child} < $xmin;
	$ymax = $yy{$child} if $yy{$child} > $ymax;
	$ymin = $yy{$child} if $yy{$child} < $ymin;
      }
      $xx{$node} = $xmin - $xstep;
      $yy{$node} = ($ymin + $ymax)/2;
    }
  }

  $xx{$t1->get_root_node} = $left + $xstep;

  my @preorder = $t1->get_nodes(-order => 'depth');

  for my $node (@preorder) {
    #print "\n$node";
    if ($colors) {
      if ($node->has_tag('Rcolor')) {
        $Rcolor{$node} = $node->get_tag_values('Rcolor')
      } else {
        $Rcolor{$node} = 0
      }
      if ($node->has_tag('Gcolor')) {
        $Gcolor{$node} = $node->get_tag_values('Gcolor')
      } else {
        $Gcolor{$node} = 0
      }
      if ($node->has_tag('Bcolor')) {
        $Bcolor{$node} = $node->get_tag_values('Bcolor')
      } else {
        $Bcolor{$node} = 0
      }
      #print "\t$Rcolor{$node}\t$Gcolor{$node}\t$Bcolor{$node}";
    }
  }

  if ($compact) { # ragged right, ignoring branch lengths

    $width = 0;
    shift @preorder; # skip root
    for my $node (@preorder) {
      $xx{$node} = $xx{$node->ancestor} + $xstep;
      $width = $xx{$node} if $xx{$node} > $width;
    }
    $width += $tip + $tipwidth1 + $right;

  } else { # set to aspect ratio and use branch lengths if available

    my $total_height = (scalar($t1->get_leaf_nodes) - 1) * $ystep;
    my $scale_factor = $total_height * $ratio / $t1->get_root_node->height;    

    $width = $t1->get_root_node->height * $scale_factor;
    $width += $left + $xstep;
    $width += $tip + $tipwidth1 + $right;

    shift @preorder; # skip root
    for my $node (@preorder) {
      my $bl = $node->branch_length;
      $bl = 1 unless (defined $bl && $bl =~ /^\-?\d+(\.\d+)?$/);
      $xx{$node} = $xx{$node->ancestor} + $bl * $scale_factor;
    }

  }

  if ($t2) {

    $x = $left + $xstep * ($root1->height + 1) + $tip;
    $x += $tipwidth1 + $tip + $column + $tip;
    my $y = $bottom;

    for my $taxon (reverse @taxa2) {
      $xx{$taxon} = $x + $tipwidth2 + $tip;
      $yy{$taxon} = $y;
      $y += $ystep2;
    }
    $x += $xstep;

    my @stack;
    my @queue; # postorder traversal
    push @stack, $t2->get_root_node;
    while (@stack) {
      my $node = pop @stack;
      push @queue, $node;
      foreach my $child ($node->each_Descendent(-sortby => 'internal_id')) {
	push @stack, $child;
      }
    }
    @queue = reverse @queue;

    for my $node (@queue) {
      if (!$node->is_Leaf) {
	my @children = $node->each_Descendent;
	my $child = shift @children;
	my $xmax = $xx{$child};
	my $ymin = my $ymax = $yy{$child};
	foreach $child (@children) {
	  $xmax = $xx{$child} if $xx{$child} > $xmax;
	  $ymax = $yy{$child} if $yy{$child} > $ymax;
	  $ymin = $yy{$child} if $yy{$child} < $ymin;
	}
	$xx{$node} = $xmax + $xstep;
	$yy{$node} = ($ymin + $ymax)/2;
      }
    }

  }

  return $self;
}

=head2 print

 Title   : print
 Usage   : $obj->print();
 Function: Outputs $obj in Encapsulated PostScript (EPS) format 
 Returns : 
 Args    : -file => filename (optional)

=cut

sub print {
  my($self,@args) = @_;

  my ($file) = $self->_rearrange([qw(FILE)], @args);
  $file ||= "output.eps"; # stdout

  open my $INFO, '>', $file or $self->throw("Could not write file '$file': $!");
  print $INFO "%!PS-Adobe-\n";
  print $INFO "%%BoundingBox: 0 0 ", $width, " ", $height, "\n";
  print $INFO "1 setlinewidth\n";
  print $INFO "/$font findfont\n";
  print $INFO "$size scalefont\n";
  print $INFO "setfont\n";

  # taxa labels are centered to 1/3 the font size

  for my $taxon (reverse $t1->get_leaf_nodes) {
    if ($colors) {
      print $INFO $Rcolor{$taxon}, " ", $Gcolor{$taxon}, " ", $Bcolor{$taxon}, " setrgbcolor\n";
    }
    print $INFO $xx{$taxon} + $tip, " ", $yy{$taxon} - $size / 3, " moveto\n";
    print $INFO "(", $taxon->id, ") show\n";
  }

  my $root1 = $t1->get_root_node;
  for my $node ($t1->get_nodes) {
    if ($node->ancestor) {
      # print $xx{$node->ancestor}, " ", $yy{$node->ancestor}, " moveto\n";
      # print $xx{$node}, " ", $yy{$node}, " lineto\n";
      if ($colors) {
	print $INFO "stroke\n";
	print $INFO $Rcolor{$node}, " ", $Gcolor{$node}, " ", 
	$Bcolor{$node}, " setrgbcolor\n";
      }
      print $INFO $xx{$node}, " ", $yy{$node}, " moveto\n";
      print $INFO $xx{$node->ancestor}, " ", $yy{$node}, " lineto\n";
      if( $bootstrap ) {
	  print $INFO $xx{$node->ancestor}+ $size/10, " ", $yy{$node->ancestor} - ($size / 3), " moveto\n";
	  print $INFO "(", $node->ancestor->id, ") show\n";
	  print $INFO $xx{$node->ancestor}, " ", $yy{$node}, " moveto\n";
      }
      print $INFO $xx{$node->ancestor}, " ", $yy{$node->ancestor}, " lineto\n";
      
      
    }
  }
  my $ymin = $yy{$root1};
  my $ymax = $yy{$root1};
  foreach my $child ($root1->each_Descendent) {
    $ymax = $yy{$child} if $yy{$child} > $ymax;
    $ymin = $yy{$child} if $yy{$child} < $ymin;
  }
  my $zz = ($ymin + $ymax)/2;
  if ($colors) {
    print $INFO "stroke\n";
    print $INFO $Rcolor{$root1}, " ", $Gcolor{$root1}, " ", $Bcolor{$root1}, " setrgbcolor\n";
  }
  print $INFO $xx{$root1}, " ", $zz, " moveto\n";
  print $INFO $xx{$root1} - $xstep, " ", $zz, " lineto\n";

  if ($t2) {

    for my $taxon (reverse $t2->get_leaf_nodes) {
      my $tiplen2 = PostScript::Metrics::stringwidth($taxon->id,$font,$size);
      print $INFO $xx{$taxon} - $tiplen2 - $tip, " ",
        $yy{$taxon} - $size / 3, " moveto\n";
      printf $INFO "(%s) show\n", $taxon->id;
    }

    for my $node ($t2->get_nodes) {
      if ($node->ancestor) {
        print $INFO $xx{$node}, " ", $yy{$node}, " moveto\n";
        print $INFO $xx{$node->ancestor}, " ", $yy{$node}, " lineto\n";
        print $INFO $xx{$node->ancestor}, " ",
          $yy{$node->ancestor}, " lineto\n";
      }
    }

    my $root2 = $t2->get_root_node;
    my $ymin = $yy{$root2};
    my $ymax = $yy{$root2};
    foreach my $child2 ($root2->each_Descendent) {
      $ymax = $yy{$child2} if $yy{$child2} > $ymax;
      $ymin = $yy{$child2} if $yy{$child2} < $ymin;
    }
    my $zz = ($ymin + $ymax)/2;
    print $INFO $xx{$root2}, " ", $zz, " moveto\n";
    print $INFO $xx{$root2} + $xstep, " ", $zz, " lineto\n";

    my @taxa1 = $t1->get_leaf_nodes;
    my @taxa2 = $t2->get_leaf_nodes;

    # set default connection between $t1 and $t2 taxa, unless
    # overridden by the user (the latter not implemented yet)

    foreach my $taxon1 (@taxa1) {
      foreach my $taxon2 (@taxa2) {
	if ($taxon1->id eq $taxon2->id) {
	  $taxon1->add_tag_value('connection',$taxon2);
	  last;
	}
      }
    }

    # draw connection lines between $t1 and $t2 taxa

    print $INFO "stroke\n";
    print $INFO "0.5 setgray\n";

    foreach my $taxon1 (@taxa1) {
      my @match = $taxon1->get_tag_values('connection');
      foreach my $taxon2 (@match) {
	my $x0 = $xx{$taxon1} + $tip
	  + PostScript::Metrics::stringwidth($taxon1->id,$font,$size) + $tip;
	my $x1 = $xx{$taxon1} + $tip + $tipwidth1 + $tip;
        my $y1 = $yy{$taxon1};
        my $x2 = $xx{$taxon2} - $tip - $tipwidth2 - $tip;
        my $x3 = $xx{$taxon2} - $tip
	  - PostScript::Metrics::stringwidth($taxon2->id,$font,$size) - $tip;
        my $y2 = $yy{$taxon2};
        print $INFO $x0, " ", $y1, " moveto\n";
        print $INFO $x1, " ", $y1, " lineto\n";
        print $INFO $x2, " ", $y2, " lineto\n";
        print $INFO $x3, " ", $y2, " lineto\n";
      }
    }

  }

  print $INFO "stroke\n";
  print $INFO "showpage\n";
}

1;
