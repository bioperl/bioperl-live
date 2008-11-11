package Bio::Graphics::Glyph::phylo_align;

use strict;
use base qw(Bio::Graphics::Glyph::generic Bio::Graphics::Glyph::xyplot);
use Bio::TreeIO;
use Bio::Graphics::Wiggle;
use POSIX qw(log10);

use Carp 'croak','cluck';
use Data::Dumper;

my %complement = (g=>'c',a=>'t',t=>'a',c=>'g',n=>'n',
		  G=>'C',A=>'T',T=>'A',C=>'G',N=>'N');


# turn off description
sub description { 0 }

# turn off label
# sub label { 1 }

sub height {
  my $self = shift;
  my $font = $self->font;
  
  #adjust the space to take if conservation scores are drawn instead
  if (! $self->dna_fits) {
    my $species_spacing_score = $self->option('species_spacing_score') || 5;
    $self->factory->set_option('species_spacing', $species_spacing_score);
  }
  
  my $species_spacing = $self->option('species_spacing') || 1;
  
  #Height = NumSpecies x Spacing/species x FontHeight
  my $height = ($self->known_species + $self->unknown_species + 1)
            * $species_spacing
            * $font->height;
  
  #use this if you want to show only those species that have alignments in the viewing window
  #$height = ($self->extract_features + 1) * 2 * $font->height;
  
  $self->factory->set_option('height', $height);
  
  
  return $height;
  
}

#####
# TODO: extract the wigfiles covering the range as well
#####
# get all features within the viewing window
sub extract_features {
  my $self = shift;
  #my $segment = $self->feature->{'factory'}->segment($self->feature->refseq,
  #						      $self->feature->start => $self->feature->stop);
  #my @match = $segment->features('submatch:pa'); 
  my @match = $self->feature->features('submatch:pa');
#print "Match has ",$#match,"<p>\n";
  
  # exract wifiles here too:
  my @wmatch = $self->feature->features('wfile:pa');
  #push @match, $self->feature->features('wfile:pa');
#  print "xxxxxx<pre>",Dumper(@wmatch),"</pre>cccccc";
#print "WMatch has ",$#wmatch,"<p>\n";

#print "Halfwaypoint<p>";  
#for my $feature (@match,@wmatch) {
#my %attributes = $feature->attributes;
#my $species = $attributes{'species'};
#print "<pre>Feature $species:\n",Dumper(%attributes),"</pre>\n";
#}

  my %alignments;
#  for my $feature (@match) {
  for my $feature (@match,@wmatch) {
    my %attributes = $feature->attributes;
    my $species = $attributes{'species'};
    
    push @{$alignments{$species}}, $feature;
  }
  
  %alignments;
}

#known species (all those that are in the Phylo tree)
sub known_species {
  my $self = shift;
  
  my $tree = shift;
  
  if ($tree) {
    my @leaves = $tree->get_leaf_nodes;
    my @allspecies = map {$_->id} @leaves;
    return @allspecies
    
  } else {
    #this may be too simple of an assumption, especially for non newick files
    my $tree_file = $self->option('tree_file');
    
    open (FH, $tree_file);
    my $newick = <FH>;
    close FH;
    
    my @allspecies = $newick =~ /([a-zA-Z]\w*)/g;
    return @allspecies;
  }
}

sub unknown_species {
  my $self = shift;

  my %alignments;        #all species in viewing window
  my $refspecies;        #all species from cladogram info
  my @current_species;   #all species in viewing window
  my @known_species;     #species in GFF but and clado
  my @unknown_species;   #species in GFF but not in clado
  # current - known = unknown
  
  
  if (@_) {
    %alignments = %{$_[0]};
    $refspecies = $_[1];
    @current_species = @{$_[2]};
    @known_species = @{$_[3]};
    @unknown_species;
  } else {
    %alignments = $self->extract_features;
    $refspecies = $self->option('reference');
    @current_species =  keys %alignments;   #all species in viewing window
    @known_species = $self->known_species;  #all species from cladogram info
    @unknown_species;                       #species in GFF but not in clado
  } #would have combined the two cases into one line using || but Perl will treat the arrays as num of elem
  
  #do set subtraction to see which species in viewing range but not in tree
  my %seen;  # build lookup table
  @seen{@known_species} = ();
  foreach my $item (@current_species, $refspecies) {
    push(@unknown_species, $item) unless exists $seen{$item};
  }
  
  return @unknown_species;
  
}

sub set_tree {
  my $self = shift;
  
  #warn"My species are ".Dumper(@species);
  
  my $tree_file   = $self->option('tree_file');
  my $tree_format = $self->option('tree_format') || 'newick';
  
  my $treeio = new Bio::TreeIO(-file   => $tree_file,
                               -format => $tree_format);
  
  
  
  my $tree = $treeio->next_tree;
  my $root = $tree->get_root_node;
  
  # would be ideal to remove all species that don't have features (alignments) within
  # viewing window but there is a bug in Bio::Tree library where you can't remove the
  # first leaf node.
  
  #  $tree->remove_Node('dog'); # etc...
  
  #set the leaf x coodinate (make all evenly spaced)
  my @leaves = $tree->get_leaf_nodes;
  for (my $i=0; $i<@leaves; $i++) {
    my $leaf = $leaves[$i];
  
    #note that leaves can use "description" functions while intermediate nodes cannot
    #thus objects must be handled directly
    $leaf->description({'x'=>$i});
  
  }
  
  
  #set root height to 0
  $root->{'description'}{'y'} = 0;
  
  #set the x and y coordinates of all intermediate nodes
  get_n_set_next_treenode($root, 0);
  
  flip_xy($tree);

  
  $tree;
  
}

sub get_max_height {
  my $tree = shift;
  my $max_height;
  
  #get the max height
  for my $child ($tree->get_leaf_nodes) {
    my $x = $child->{'_description'}{'x'};
    $max_height = $x if $max_height < $x;
  }
  
  $max_height;
}


sub draw_clado {
  my $self = shift;
  my $tree = shift;
  my $gd = shift;
  my ($x1, $y1, $x2, $y2, $color, $xscale, $yscale, $xoffset, $yoffset, $start_x, $draw_clado_left) = @_;
  
  my @bounds = $gd->getBounds;
  
  my $root = $tree->get_root_node;
  my @nodes = $root->get_all_Descendents;
  
  #draw bg for cladogram
  my $clado_bg = $self->color('clado_bg') || $self->bgcolor;
  my @coords = (0, $y1, $start_x+$xoffset+$self->font->width-1, $y2+1);
  my @coords2 = ($x1, $y1, $start_x+$xoffset/2, $y2);
  if ($draw_clado_left) {
    $gd->filledRectangle(@coords, $clado_bg);
    $gd->filledRectangle(@coords2, $self->color('clado_bg'));
    $gd->filledRectangle($x2, $x1, $bounds[0], $bounds[1], $self->color('bg_color')) if $self->dna_fits;
  } else {
    $gd->filledRectangle($bounds[0]-$coords[2], $coords[1], $bounds[0]-$coords[0], $coords[3],
  			 $self->color('bg_color'));
    $gd->filledRectangle($bounds[0]-$coords2[2], $coords2[1], $bounds[0]-$coords2[0], $coords2[3],
			 $clado_bg);  
    $gd->filledRectangle(0, $y1, $x1, $y2+1, $self->color('bg_color')) if $self->dna_fits;
  }

  
  
  #draw the lines of the tree
  for my $node ($root,@nodes) {
    next if $node->is_Leaf;
    my $x = $node->{'_description'}{'x'} * $xscale;
    
    #draw vertical line covering all children
    my $topx = $node->{'_description'}{'childmin'} * $yscale;
    my $botx = $node->{'_description'}{'childmax'} * $yscale;
    
    @coords = ($x+$xoffset, $topx+$yoffset, $x+$xoffset, $botx+$yoffset);
    if ($draw_clado_left) {
      $gd->line(@coords, $self->fgcolor);
    } else {
      $gd->line($bounds[0]-$coords[2], $coords[1], $bounds[0]-$coords[0], $coords[3], $self->fgcolor);
    }
    
    #draw a line connecting the bar to each child
    my @children = $node->each_Descendent;
    for my $child (@children) {
      my $cx = $child->{'_description'}{'x'} * $xscale;
      my $cy = $child->{'_description'}{'y'} * $yscale;
      $cx = $start_x if $child->is_Leaf;
      
      #print"($cx, $cy)";
      
      @coords = ($x+$xoffset, $cy+$yoffset, $cx+$xoffset, $cy+$yoffset);
      if ($draw_clado_left) {
        $gd->line(@coords, $self->fgcolor);
      } else {
          $gd->line($bounds[0]-$coords[2], $coords[1], $bounds[0]-$coords[0], $coords[3], $self->fgcolor);
      }
      
    }
    
  }
  
  
  
  
  #my $tree = $treeio->next_tree;
  #my @nodes = grep { $_->bootstrap > 70 } $tree->get_nodes;
  #warn"my nodes are:\n".Dumper(@nodes);
  print"</pre>";
  
  $start_x + $xscale;
  
  }

#tree is made with root on top but this will switch x and y coords so root is on left
sub flip_xy {
  my $tree = shift;
  my $root = $tree->get_root_node;
  
  my @nodes = $root->get_all_Descendents;
  
  for my $node ($root, @nodes) {
    if ($node->is_Leaf) {
      $node->{'_description'} = {
        'x'       => $node->{'_description'}{'y'},
        'y'       => $node->{'_description'}{'x'}
      }
    } else {
      $node->{'_description'} = {
        'x'         => $node->{'_description'}{'y'},
        'y'         => $node->{'_description'}{'x'},
        'child_pos' => $node->{'_description'}{'child_pos'},
        'childmin'  => $node->{'_description'}{'childmin'},
        'childmax'  => $node->{'_description'}{'childmax'}
      }
    }
  }
  
  
}


#recursive function that sets the x and y coordinates of tree nodes
sub get_n_set_next_treenode {
  my $node = shift;
  my $height = $node->{'_description'}{'y'};
  
  my @children = $node->each_Descendent;
  
  
  my $x = 0;
  my $min_child_x = -1;
  my $max_child_x = -1;
  
  #iterate through children to find the x's and set the y's (height's)
  for my $child (@children) {
    #set the y coordinate as parent's height + 1
    $child->{'_description'}{'y'} = $height + 1;
    get_n_set_next_treenode($child);
    
    #retrieve the child's x coordinate
    my $child_x = $child->{'_description'}{'x'} || 0;
    $x += $child_x;
    
    $min_child_x = $child_x if $min_child_x==-1 || $child_x < $min_child_x;
    $max_child_x = $child_x if $max_child_x==-1 || $max_child_x < $child_x;
    #$x += $child->discription->{'x'};  #cannot do this for intermediate nodes
    
    #store the x values of all children
    push @{$node->{'_description'}{'child_pos'}}, $child_x;
    
    
  }
  
  #set the current x coordinate as the average of all children's x's
  if (@children) {
    $x = $x / @children;
    $node->{'_description'}{'x'} = $x;
    $node->{'_description'}{'childmin'} = $min_child_x;
    $node->{'_description'}{'childmax'} = $max_child_x;

  }
  
  $node->{'_description'}{'y'} = $height;
   
}

sub get_legend_and_scale {
  my $yscale = shift;
  my $height = shift;
  
  if ($yscale < 2*$height - 1) {
    $height = 0;
  }
  
  #########
  # chage scale later so that the base can  be anything and not just 1!!
  
  # scale legend goes in order from min, axis, max, if either min or max = 1, then will only have min & max
  my @order = sort {$a <=> $b} (1, @_);
  my $graph_scale = - ($yscale - $height) / (log10($order[2]) - log10($order[0]));
  my $graph_legend = {1 => $graph_scale * (log10(1) - log10($order[2])),
  		   $order[0] => $graph_scale * (log10($order[0]) - log10($order[2])),
    		   $order[2] => 0};
    
  #print "order is @order and the yscale is $yscale and height is $height<br>";
  
  return ($graph_legend, $graph_scale);
}

#main method that draws everything
sub draw {
  my $self = shift;
  my $height = $self->font->height;
  my $scale = $self->scale;
  
  my $gd = shift;
  my ($left,$top,$partno,$total_parts) = @_;
  my ($x1,$y1,$x2,$y2) = $self->bounds($left, $top);
  
  
  my @bounds = $gd->getBounds;
  my $draw_clado_left = $self->option('draw_clado_left');
  
  
  #spacing of either DNA alignments or score histograms in units of font height
  my $species_spacing = $self->option('species_spacing') || 1;
  
  my $xscale = $self->font->width;
  my $yscale = $height * $species_spacing;
  
  
  my $xoffset = $x1;
  my $yoffset = $y1 + 0.5*$self->font->height;
  
  #method that reads the tree file to create the tree objects
  my $tree = $self->set_tree;
  my $max_height = get_max_height($tree);
  my $start_x =($max_height-1) * $xscale +$xoffset;
  
  
  my $connector = $self->connector;
  
  #all species having alignments in viewing window (key=name, val=feat obj)
  my %alignments = $self->extract_features;
  #print "Species are:",keys %alignments,"<br>\n";
  
  my ($min_score, $max_score) = $self->get_score_bounds(%alignments);
  #$min_score = 0 unless $min_score;
  my ($graph_legend, $graph_scale) = get_legend_and_scale($yscale, $height, $min_score, $max_score);
#print "min/max scores: $min_score, $max_score<br>\n",
#"graph legend and scale: $graph_legend, $graph_scale";
# TODO: Gap entries give an undef for the min values for some reason
  
  
  my $refspecies = $self->option('reference');
  
  my @current_species = keys %alignments;    #all species in viewing window
  my @known_species = $self->known_species($tree);  #all species from cladogram info
  my @unknown_species = $self->unknown_species(\%alignments, 
    						 $refspecies,
    						\@current_species,
    						\@known_species);
                                              #species in GFF but not in clado

  ########is this even used?
  #my @allfeats;
  #for my $species (keys %alignments) {
  #  push @allfeats, @{$alignments{$species}};
  #}
  
  
  #this y value is the base for the next species' alignment/histogram and is incremented at each step 
  my $y = $y1;
  
  
  for my $species (@known_species,@unknown_species) {
    my $y_track_top = $y + $height;
    my $y_track_bottom = $y + $yscale;
    
    
    if ($yscale < 2*$height-1) {
      #small scale
      $y_track_top = $y;# + $height;
      my $y_track_bottom = $y + $height;
    }
    
    
    #process the reference sequence differently
    if ($species eq $refspecies) {
      #draw DNA alignments if zoomed close enough
      my ($fx1,$fy1) = ($x1, $y_track_top);
      my ($fx2,$fy2) = ($x2,$y_track_bottom);
      

      if ($self->dna_fits) {
	my $dna = eval { $self->feature->seq };
        $dna    = $dna->seq if ref($dna) and $dna->can('seq'); # to catch Bio::PrimarySeqI objects
        my $bg_color = $self->color('ref_color') || $self->bgcolor;
        
        $fy2 = $fy1 + $self->font->height || $y2;
  
        
	$self->_draw_dna($gd,$dna,$fx1,$fy1,$fx2,$fy2, $self->fgcolor, $bg_color);
      } else {
      }
      
      my $x_label_start = $start_x + $xoffset + $self->font->width;
      $self->species_label($gd, $draw_clado_left, $x_label_start, $y, $species) unless ($self->option('hide_label'));
      
      $y += $yscale;
      next;
    }
    
    
    #skip if the there is no alignments for this species in this window
    unless ($alignments{$species}) {
      my $x_label_start = $start_x + $xoffset + $self->font->width;
      $self->species_label($gd, $draw_clado_left, $x_label_start, $y, $species) unless ($self->option('hide_label'));
      
      $y += $yscale;
      next;
    }
    
    
    
    my @features = @{$alignments{$species}};
    
    
    
    
    #draw the axis for the plots
    $self->draw_pairwisegraph_axis($gd,
    				    $graph_legend,
    				    $x1,
    				    $x2,
    				    $y_track_top,
    				    $y_track_bottom,
    				    $draw_clado_left,
    				    @bounds) unless $self->dna_fits;
      
    
    #iterate through the wigfiles and put them on the graph
    ###
    # todo
    ###
    
    
    #iterate through features, and put them on the graph
    for my $feat (@features) {
      my ($start, $stop, %attributes) = ($feat->start, $feat->stop, $feat->attributes);
      
      my ($fx1,$fy1) = ($x1 + ($start-$self->start)*$scale, $y_track_top);
      my ($fx2,$fy2) = ($x1 + ($stop-$self->start)*$scale,$y_track_bottom);
      
      my $gapstr = $attributes{'Gap'} || "M".($stop-$start+1);
      my @gapstr = split " ", $gapstr;
      my @gaps;
      for my $gap (@gapstr) {
        my ($type, $num) = $gap =~ /^(.)(\d+)/; 
      	push @gaps, [$type, $num+0];
      }
      
      
      
      #draw DNA alignments if zoomed close enough
      if ($self->dna_fits) {
	
	
	my $test = 0;
	
	my $hit = $feat->hit;
	if (!defined $hit) {
	  warn "No hit for feature $feat, skipping drawing DNA";
	  next;
	}

	my $hit_seq = $hit->seq || print "No seq for hit<br>\n";
	my $targ_dna = $hit_seq->seq || print "No seq for hit_seq<br>\n";
	
	
	my $seq = $feat->seq || print "No sequence object for feature<br>\n";
	my $ref_dna = $seq->seq || print "No ref dna";
	
	#doesn't work as planned
	# my $ref_dna = $feat->seq->seq || print "No ref dna";
	# my $targ_dna = $feat->hit->seq->seq || print "No targ dna";
	
	
	next if !defined $ref_dna || !defined $targ_dna;
	
	$self->draw_dna($gd,$ref_dna, $targ_dna,$fx1,$fy1,$fx2,$fy2,\@gaps);
      } else {
      	my $wigfile = $attributes{'wigfile'};
      	if ($wigfile) {
      	  if (-e $wigfile) {
      	    $self->pairwise_draw_wig_graph($gd, $feat, $x1, $scale, \@gaps, $graph_legend->{1}, $graph_scale, $fx1, $fy1, $fx2, $fy2,$wigfile);
      	  } else {
      	    warn "Wigfile $wigfile does not exist, skipping ...";
      	  }
      	  
      	} else {
      	  $self->pairwise_draw_graph($gd, $feat, $x1, $scale, \@gaps, $graph_legend->{1}, $graph_scale, $fx1, $fy1, $fx2, $fy2);
      	}
      	
      }
    }
    
    
    #label the species in the cladogram
    my $x_label_start = $start_x + $xoffset + $self->font->width;
    $self->species_label($gd, $draw_clado_left, $x_label_start, $y, $species) unless ($self->option('hide_label'));
    
    $y += $yscale;
  }
  
    $self->draw_clado($tree, $gd, $x1, $y1, $x2, $y2, $self->fgcolor,
  		     $xscale, $yscale, $xoffset, $yoffset, $start_x, $draw_clado_left);

}



sub species_label {
  my $self = shift;
  my $gd = shift;
  my $draw_clado_left = shift;
  my $x_start = shift;
  my $y_start = shift;
  my $species = shift;
  
  $x_start += 2;
  my $text_width = $self->font->width * length($species);
  my $bgcolor = $self->color('bg_color');
  
  #make label
  if ($draw_clado_left) {
    
    $gd->filledRectangle($x_start-2, $y_start, $x_start + $text_width, $y_start+$self->font->height, $bgcolor);
    $gd->rectangle($x_start-2, $y_start, $x_start + $text_width, $y_start+$self->font->height, $self->fgcolor);
    $gd->string($self->font, $x_start, $y_start, $species, $self->fgcolor);
    
  } else {
    my ($x_max, $y_max) = $gd->getBounds;
    my $write_pos = $x_max - $x_start - $text_width;
    
    $gd->filledRectangle($write_pos, $y_start, $write_pos + $text_width+2, $y_start+$self->font->height, $bgcolor);
    $gd->rectangle($write_pos, $y_start, $write_pos + $text_width+2, $y_start+$self->font->height, $self->fgcolor);
    $gd->string($self->font, $write_pos+2, $y_start, $species, $self->fgcolor);
    
  }
}


# draws the legends on the conservation scale
sub draw_pairwisegraph_axis {
  my $self = shift;
  my ($gd, $graph_legend, $x1, $x2, $y_track_top, $y_track_bottom, $draw_clado_left, @bounds) = @_;
  
  
  my $axis_color = $self->color('axis_color') || $self->fgcolor;
  my $mid_axis_color = $self->color('mid_axis_color') || $axis_color;
  
  for my $label (keys %$graph_legend) {
    my $y_label = $graph_legend->{$label} + $y_track_top;

    
    my $col = $axis_color;
    $col = $mid_axis_color if ($y_label != $y_track_top && $y_label != $y_track_bottom);
    $gd->line($x1,$y_label,$x2,$y_label,$col);
    
    my @coords = (0, $y_label, $x1, $y_label);
    
    
    if ($draw_clado_left) {
      #draw the legend on the right
      $coords[0] = $bounds[0] - $coords[0];
      $coords[2] = $bounds[0] - $coords[2];
      
      my $x_text_offset = length($label) * $self->font->width;
      
      $gd->string($self->font, $coords[0]-$x_text_offset, $coords[1], $label, $self->fgcolor);
      $gd->line(@coords, $self->fgcolor);
      
      $gd->line($x2,$y_track_top,$x2,$y_track_bottom,$self->fgcolor);
    } else {
      #draw the legned on the left
      $gd->string($self->font, @coords[0..1], $label, $self->fgcolor);
      $gd->line(@coords, $self->fgcolor);
      
      $gd->line($x1,$y_track_top,$x1,$y_track_bottom,$self->fgcolor);
    }
  
  }
  
}


#find min and max from features within the bounds
sub get_score_bounds {
  my $self = shift;
  my %alignments = @_;
  
  my $min = -1;
  my $max = -1;
  
  for my $species (keys %alignments) {
    for my $feature (@{$alignments{$species}}) {
      my $score = $feature->score;
      
      #check to see if wigfile
      if ($score == undef) {
      	my %attributes = $feature->attributes;
      	if (-e $attributes{'wigfile'}) {
      	  ($min, $max) = $self->get_score_bounds_wigfile($feature,$min,$max,$attributes{'wigfile'});
      	}
      	next;
      }
      
      $min = $score if $min == -1 || $score < $min;
      $max = $score if $max == -1 || $max < $score;
    }
  }
  
  
  my @parts = $self->parts;
  
  return ($min, $max)
}

#find min and max of sampled wigfile
sub get_score_bounds_wigfile {
  my $self = shift;
  my $feature = shift;
  my ($min,$max,$wigfile) = @_;
  
  my $start = $feature->start < $self->start ? $self->start : $feature->start;
  my $stop  = $self->stop < $feature->stop   ? $self->stop  : $feature->stop;
  
  #print $self->stop, "-", $feature->stop, "checking @_ $start - $stop\n";
  
  #extract wig file contents  
  my $wig = Bio::Graphics::Wiggle->new($wigfile, 0, {step => 1}) or die;
  
  #todo: make step configurable
  my $step = 100;
  
  for (my $i=$start; $i<$stop; $i += $step) {
    my $v = $wig->value($i);
    
    next unless defined $v;
    
    $min = $v if !defined $min or $v < $min;
    $max = $v if !defined $max or $max < $v;
    
  }
  
  #print "min and max are $min , $max<br>\n";
  
  return ($min,$max);
  
}



sub pairwise_draw_graph {
  my $self = shift;
  my $gd = shift;
  my $feat = shift;		# current feature object
  my $x_edge = shift;		# x start position of the track
  my $scale = shift;		# pixels / bp
  my $gaps = shift;		# gap data for insertions, deletions and matches
  my $zero_y = shift;		# y coordinate of 0 position
  my $graph_scale = shift;	# scale for the graph. y_coord = graph_scale x log(score)
  
  my ($x1,$y1,$x2,$y2) = @_;
  my $fgcolor = $self->fgcolor;
  my $errcolor  = $self->color('errcolor') || $fgcolor;
  
  my $score = $feat->score;
  my %attributes = $feat->attributes;
  
  
  my $log_y = log10($score);
  my $y_bottom = log10($score) * $graph_scale + $zero_y + $y1;
  my $y_top = $zero_y+$y1;
  
  my @y = sort {$a <=> $b} ($y_bottom, $y_top);
  
  
  
  #missing gap data
  unless ($gaps) {
    $x1 = $x_edge if $x1 < $x_edge;
    return if $x2 < $x_edge;
    #$gd->filledRectangle($x1,$y[0],$x2,$y[1],$fgcolor);
    return;
  }
  
  my $bp = 0;
  
  #draw a bar representing the score for the span of base pairs
  for my $tuple (@$gaps) {
    my ($type, $num) = @$tuple;
    
    #warn"$type and $num";
    if ($type eq "M") {
      my $x_left  = $x1 + ($bp*$scale);
      my $x_right = $x_left + $num*$scale;
      
      $bp += $num;
      
      $x_left = $x_edge if $x_left < $x_edge;
      next if $x_right < $x_edge;
      $gd->filledRectangle($x_left,$y[0],$x_right,$y[1],$fgcolor);
    } elsif ($type eq "D") {
      
      $bp += $num;
      
    } elsif ($type eq "I") {
      my $x_left  = $x1 + ($bp*$scale);
      $gd->line($x_left-2, $y1-4, $x_left, $y1, $errcolor);
      $gd->line($x_left, $y1, $x_left+2, $y1-4, $errcolor);
    }
    

  }
  
  
}



sub pairwise_draw_wig_graph {
  my $self = shift;
  my $gd = shift;
  my $feat = shift;		# current feature object
  my $x_edge = shift;		# x start position of the track
  my $scale = shift;		# pixels / bp
  my $gaps = shift;		# gap data for insertions, deletions and matches
  my $zero_y = shift;		# y coordinate of 0 position
  my $graph_scale = shift;	# scale for the graph. y_coord = graph_scale x log(score)
  
  my ($x1,$y1,$x2,$y2,$wigfile) = @_;
  my $fgcolor = $self->fgcolor;
  my $errcolor  = $self->color('errcolor') || $fgcolor;
  
  my $score = $feat->score;
  my %attributes = $feat->attributes;
  
  
#  my $log_y = log10($score);
#  my $y_bottom = log10($score) * $graph_scale + $zero_y + $y1;
#  my $y_top = $zero_y+$y1;
#  
#  my @y = sort {$a <=> $b} ($y_bottom, $y_top);
  
  
  
  #print "checking wigfile $wigfile<br>\n";
  
  
  #todo: make step variable
  my $start = $feat->start < $self->start ? $self->start : $feat->start;
  my $stop  = $self->stop < $feat->stop   ? $self->stop  : $feat->stop;
#  my $wig = Bio::Graphics::Wiggle->new($wigfile, 0) or die;
  my $wig = Bio::Graphics::Wiggle->new($wigfile, 0, {step => 1}) or die;


##### not sure why this was here
#  my $vals = $wig->values($start,$stop);
#  my ($vmin,$vmax);
#  for my $v (@$vals) {
#    next unless defined $v;
#    $vmin = $v if !defined $vmin or $v < $vmin;
#    $vmax = $v if !defined $vmax or $vmax < $v;
#  }
  
  
#  print "min and max are $vmin , $vmax\n";
#  $min = $min < $vmin ? $min : $vmin;
#  $max = $vmax < $max ? $max : $vmax;
#  print "min and max are $min , $max\n";
#  return ($min,$max);
  
  
#  print Dumper($vals);
  
  
#  my $pos = $start;
  
  
#  $gd->rectangle($x1,$y1,$x2,$y2,$fgcolor);
#  $gd->line($x1,$y1,$x2,$y2,$fgcolor);
#  $gd->line($x2,$y1,$x1,$y2,$fgcolor);
#  print "Start and stop: $start and $stop<br>\n";
#  print "<br>\n$x1, $x2, $x_edge, $start, $stop<br>\n";
  
  
  
  if ($scale < 1) {
    #### Algorithm 1, when zoomed at scale 1bp / pixel or more
    #### sample 10 values across the pixel (prevents redrawing same pixel over
    #### and over).  The sample values are averaged.
    
    my $bp = $start;
    for (my $pix=$x1; $pix < $x2; $pix++) {
      
      $bp = ($pix - $x1)/$scale + $start;
      
      #todo: make samplesize an option
      my $samplesize = 10;
      my $score = 0;
      my $trial;
      for ($trial=0; $trial < $samplesize; $trial++) {
        
        my $samp_bp = $bp + ($trial/$samplesize)/$scale;
        last if $samp_bp > $stop;
  
        $score += $wig->value($samp_bp);
        
        #print "trial $trial: $samp_bp\tPixel $pix\tpos $x1<br>\n" if ($bp > 1700 && $bp < 1750);
        
        
      }
      
      #print "Pixel: $pix : Total $score and trial $trial<br>\n";
      next if $trial == 0 || $score == 0;
      $score = $score / $trial;
      
      $self->draw_log10_rectangle($score, $graph_scale, $zero_y, $y1,
      				$zero_y, $pix, $pix, $gd, $fgcolor);
      
      
    }
    
  } else {
    #### Algorithm 2, when zoomed in close at less than 1bp / pixel (1 bp spans
    #### 1 pixel or more), draw each value directly
    
    my $step = 1;
    
    for (my $i=$start; $i<$stop; $i += $step) {
      my $val = $wig->value($i);
      next if !defined $val;
      
      my $bp = $i - $start;
      
      my $x_left  = $x1 + ($bp*$scale);
      my $x_right = $x_left + 1*$scale;
      
      $score = $val;
      
      
      $self->draw_log10_rectangle($score, $graph_scale, $zero_y, $y1, 
      				$zero_y, $x_left, $x_right, $gd, $fgcolor);
      
      
    }
    
    
  }
  
  return;

  
}

sub draw_log10_rectangle {
  my $self = shift;
  my $score = shift;
  my $graph_scale = shift;
  my $zero_y = shift;
  my $y1 = shift;
  my $zero_y = shift;
  my $x_left = shift;
  my $x_right = shift;
  my $gd = shift;
  my $fgcolor = shift;
  
  my $log_y = log10($score);
  my $y_bottom = log10($score) * $graph_scale + $zero_y + $y1;
  my $y_top = $zero_y+$y1;
  
  my @y = sort {$a <=> $b} ($y_bottom, $y_top);
  
  #print "$x_left,$y[0],$x_right,$y[1]<br>\n";
  
  $gd->filledRectangle($x_left,$y[0],$x_right,$y[1],$fgcolor);
}



sub draw_dna {
  my $self = shift;
  my ($gd,$ref_dna, $dna,$x1,$y1,$x2,$y2, $gaps) = @_;
  my $pixels_per_base = $self->scale;
  
  my $fgcolor = $self->fgcolor;
  my $bg_color = $self->color('targ_color') || $self->bgcolor;
  my $errcolor  = $self->color('errcolor') || $fgcolor;
  
  $y2 = $y1 + $self->font->height || $y2;
  
  
  
  #missing gap data, draw as is
  unless ($gaps) {
    warn"no gap data for DNA sequence $dna";
    $self->_draw_dna($gd, $dna, $x1, $y1, $x2, $y2);
    return;
  }
  
  #parse the DNA segments by the gaps
  for my $tuple (@$gaps) {
    my ($type, $num) = @$tuple;
    
    if ($type eq "M") {
      my $dnaseg = substr($dna, 0, $num);
      my $ref_dnaseg = substr($ref_dna, 0, $num);
      
      $self->_draw_dna($gd,$dnaseg, $x1, $y1, $x2, $y2, $fgcolor, $bg_color,$ref_dnaseg);
      $dna = substr($dna, $num);
      $ref_dna = substr($ref_dna, $num);
      $x1 += $num * $pixels_per_base;
    } elsif ($type eq "D") {
      my $dnaseg = '-' x $num;
      $self->_draw_dna($gd, $dnaseg, $x1, $y1, $x2, $y2, $fgcolor); 
      $gd->rectangle($x1, $y1-1, $x1+$num * $pixels_per_base, $y2+1, $errcolor);
      $ref_dna = substr($ref_dna, $num);
      $x1 += $num * $pixels_per_base;
    } elsif ($type eq "I") {
      $dna = substr($dna, $num);
      $gd->line($x1-2, $y1-2, $x1+2, $y1-2, $errcolor);
      $gd->line($x1, $y1-1, $x1, $y2+1, $errcolor);
      $gd->line($x1-2, $y2+1, $x1+2, $y2+1, $errcolor);
    }
    

  }
  
  
}

sub _draw_dna {
  my $self = shift;

  #the last argument is optional.  If the reference seq is given, it will check it with target
  my ($gd,$dna,$x1,$y1,$x2,$y2, $color, $bg_color, $ref_dna) = @_;
  
  my $pixels_per_base = $self->scale;
  my $feature = $self->feature;
  
  unless ($ref_dna) {
    $gd->filledRectangle($x1+1, $y1, $x2, $y2, $bg_color);
  }
  
  
  my $feature = $self->feature;
  
  
  my $strand = $feature->strand || 1;
  $strand *= -1 if $self->{flip};

  my @bases = split '',$strand >= 0 ? $dna : $self->reversec($dna);
  my @refbases = split '',$strand >= 0 ? $ref_dna : $self->reversec($ref_dna);
  
  
  
  $color = $self->fgcolor unless $color;
  $bg_color = 0 unless $bg_color;
  my $font  = $self->font;
  my $lineheight = $font->height;
#  $y1 -= $lineheight/2 - 3;          ##################NOT SURE WHY THIS WAS HERE BEFORE
  my $strands = $self->option('strand') || 'auto';

  my ($forward,$reverse);
  if ($strands eq 'auto') {
    $forward = $feature->strand >= 0;
    $reverse = $feature->strand <= 0;
  } elsif ($strands eq 'both') {
    $forward = $reverse = 1;
  } elsif ($strands eq 'reverse') {
    $reverse = 1;
  } else {
    $forward = 1;
  }
  # minus strand features align right, not left
  $x1 += $pixels_per_base - $font->width - 1 if $strand < 0;
  for (my $i=0;$i<@bases;$i++) {
    my $x = $x1 + $i * $pixels_per_base;
    
    my $x_next = $x + $pixels_per_base;

    
    #draw background if DNA base aligns with reference (if ref given)
    $gd->filledRectangle($x+1, $y1, $x_next, $y2, $bg_color) 
    			if ( ($forward && $bases[$i] eq $refbases[$i]) ||
    			     ($reverse && $complement{$bases[$i]} eq $refbases[$i]) );
    
    $gd->char($font,$x+2,$y1,$bases[$i],$color)                                   if $forward;
    $gd->char($font,$x+2,$y1+($forward ? $lineheight:0),
	      $complement{$bases[$i]}||$bases[$i],$color)                         if $reverse;
  }
}


1;

__END__

=head1 NAME

Bio::Graphics::Glyph::phylo_align - The "phylogenetic alignment" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION


This glyph draws a cladogram for any set of species along with their
alignment data in relation to the reference species.  At high
magnification, base pair alignements will be displayed.  At lower
magnification, a conservation score plot will be drawn.  Gaps as
specified by CIGAR are supported.  Currently the scores are drawn to
a log plot with the restriction that the score will be the same
across all base pairs within an alignment.  It is hoped that this
restriction can be addressed in the future.

For this glyph to work, the feature must return a DNA sequence string
in response to the dna() method.  Also, a valid tree file must be
available in a format readable by the Bio::Tree library.

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

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description               Default
  ------      -----------               -------

  -draw_clado_left
              Draws the Cladogram on left 0

  -species_spacing
              Spacing of species in DNA   1
              mode in units of font height

  -species_spacing_score
              Spacing of spcies in        5
              conservation view in units
              of font height

  -hide_label Whether to label spcies     0

  -tree_file  Path of file containing     undef
              cladogram tree information
  -tree_format Format of tree file	   newick

  -axis_color Color of the vertical axes  fgcolor
              in the GC content graph

  -errcolor   Color of all misalignment   fgcolor
              indicators

  -mid_axis_color
              Color of the middle axis of
              the conservation score graph axis_color

  -clado_bg  Color of the clado bg       bgcolor
              indicators

  -ref_color  Color of base pair bg for   bgcolor
              the reference sequence

  -targ_color Color of base pair bg for   bgcolor
              all base pairs that match
              reference



=head1 BUGS

Please report them.

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

=head1 AUTHORS

Hisanaga Mark Okada
Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
