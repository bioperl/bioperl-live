# BioPerl module for Bio::Draw::Pictogram
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Shawn Hoon <shawnh@fugu-sg.org>
#
# Copyright Shawn Hoon
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Draw::Pictogram - generate SVG output of Pictogram display for consensus motifs

=head1 SYNOPSIS

  use Bio::Draw::Pictogram;
  use Bio::SeqIO;

  my $sio = Bio::SeqIO->new(-file=>$ARGV[0],-format=>'fasta');
  my @seq;
  while(my $seq = $sio->next_seq){
    push @seq, $seq;
  }

  my $picto = Bio::Draw::Pictogram->new(-width=>"800",
                                    -height=>"500",
                                    -fontsize=>"60",
                                    -plot_bits=>1,
                                    -background=>{
                                                  'A'=>0.25,
                                                  'C'=>0.18,
                                                  'T'=>0.32,
                                                  'G'=>0.25},
                                    -color=>{'A'=>'red',
                                             'G'=>'blue',
                                             'C'=>'green',
                                             'T'=>'magenta'});

  my $svg = $picto->make_svg(\@seq);

  print $svg->xmlify."\n";

  #Support for Bio::Matrix::PSM::SiteMatrix now included

   use Bio::Matrix::PSM::IO;

   my $picto = Bio::Draw::Pictogram->new(-width=>"800",
                                    -height=>"500",
                                    -fontsize=>"60",
                                    -plot_bits=>1,
                                    -background=>{
                                                  'A'=>0.25,
                                                  'C'=>0.18,
                                                  'T'=>0.32,
                                                  'G'=>0.25},
                                    -color=>{'A'=>'red',
                                             'G'=>'blue',
                                             'C'=>'green',
                                             'T'=>'magenta'});

  my $psm = $psmIO->next_psm;
  my $svg = $picto->make_svg($psm);
  print $svg->xmlify;

=head1 DESCRIPTION

A module for generating SVG output of Pictogram display for consensus
motifs.  This method of representation was describe by Burge and
colleagues: (Burge, C.B.,Tuschl, T., Sharp, P.A. in The RNA world II,
525-560, CSHL press, 1999)

This is a simple module that takes in an array of sequences (assuming
equal lengths) and calculates relative base frequencies where the
height of each letter reflects the frequency of each nucleotide at a
given position. It can also plot the information content at each
position scaled by the background frequencies of each nucleotide.

It requires the SVG-2.26 or later module by Ronan Oger available at
http://www.cpan.org

Recommended viewing of the SVG is the plugin available at Adobe:
http://www.adobe.com/svg

=head1 FEEDBACK


=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Shawn Hoon

Email shawnh@fugu-sg.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a "_".

=cut

package Bio::Draw::Pictogram;
use strict;
use SVG 2.26;
use Bio::SeqIO;
use base qw(Bio::Root::Root);

use constant MAXBITS => 2;

=head2 new

 Title   : new
 Usage   : my $picto = Bio::Draw::Pictogram->new(-width=>"800",
                                            -height=>"500",
                                            -fontsize=>"60",
                                            -plot_bits=>1,
                                            -background=>{
                                                          'A'=>0.25,
                                                          'C'=>0.18,
                                                          'T'=>0.32,
                                                          'G'=>0.25},
                                            -color=>{'A'=>'red',
                                                      'G'=>'blue',
                                                      'C'=>'green',
                                                      'T'=>'magenta'});
 Function: Constructor for Pictogram Object
 Returns : L<Bio::Draw::Pictogram>

=cut

sub new {
  my ($caller,@args) = @_;
  my $self = $caller->SUPER::new(@args);
  my ($width,$height,$fontsize,$color,$background,$bit,$normalize) = $self->_rearrange([qw(WIDTH HEIGHT FONTSIZE COLOR BACKGROUND PLOT_BITS NORMALIZE)],@args);
  $width||=800;
  $height||=600;
  my $svg = SVG->new(width=>$width,height=>$height);
  $self->svg_obj($svg);
  $fontsize ||= 80;
  $self->fontsize($fontsize) if $fontsize;
  $color = $color || {'T'=>'black','C'=>'blue','G'=>'green','A'=>'red'};
  $self->color($color);
  $background = $background || {'T'=>0.25,'C'=>0.25,'G'=>0.25,'A'=>0.25};
  $self->background($background);
  $self->plot_bits($bit) if $bit;
  $self->normalize($normalize) if $normalize;

  return $self;
}

=head2 make_svg

 Title   : make_svg
 Usage   : $picto->make_svg();
 Function: make the SVG object
 Returns : L<SVG>
 Arguments: A fasta file or array ref of L<Bio::Seq> objects or a L<Bio::Matrix::PSM::SiteMatrixI>

=cut

sub make_svg {
  my ($self,$input) = @_;
  my $fontsize = $self->fontsize;
  my $size = $fontsize * 0.75;
  my $width= $size;
  my $height= $size+40;
  my $color = $self->color;

  #starting x coordinate for pictogram
  my $x = 45+$size/2;
  my $pos_y = $size * 2;
  my $bit_y = $pos_y+40;
  my @pwm;

  my $bp = 1;

  #input can be file or array ref of sequences
  if(ref($input) eq 'ARRAY'){
    @pwm = @{$self->_make_pwm($input)};
  }
  elsif(ref($input) && $input->isa("Bio::Matrix::PSM::SiteMatrixI")){
    @pwm = $self->_make_pwm_from_site_matrix($input);
  }
  else {
    my $sio = Bio::SeqIO->new(-file=>$input,-format=>"fasta");
    my @seq;
    while (my $seq = $sio->next_seq){
      push @seq, $seq;
    }
    @pwm = @{$self->_make_pwm(\@seq)};
  }


  my $svg = $self->svg_obj;
  my $seq_length = scalar(@pwm + 1) * $width + $x + $x;
  my $seq_grp;

  #scale the svg if length greater than svg width
  if($seq_length > $svg->{-document}->{'width'}){
    my $ratio = $svg->{-document}->{'width'}/($seq_length);
    $seq_grp = $svg->group(transform=>"scale($ratio,1)");
  }
  else {
    $seq_grp= $svg->group();
  }

  #do the drawing, each set is a base position
  foreach my $set(@pwm){
    my ($A,$C,$G,$T,$bits) = @$set;
    my @array;
    push @array,  ['a',($A)];
    push @array, ['g',($G)];
    push @array, ['c',($C)];
    push @array, ['t',($T)];
    @array = sort {$b->[1]<=>$a->[1]}@array;
    my $count = 1;
    my $pos_group = $seq_grp->group(id=>"bp $bp");
    my $prev_size;
    my $y_trans;

    #draw each letter at each position
    foreach my $letter(@array){
	  my $scale;
	  if($self->normalize){
		$scale = $letter->[1];
	  } else {
		$scale = $letter->[1] * ($bits / MAXBITS);
	  }

      if($count == 1){
		if($self->normalize){
		  $y_trans = 0;
		} else {
		  $y_trans = (1 - ($bits / MAXBITS)) * $size;
		}
      }
      else {
        $y_trans += $prev_size;
      }
      $pos_group->text('id'=> uc($letter->[0]).$bp,height=>$height,
                      'width'=>$width,x=>$x,y=>$size,
                      'transform'=>"translate(0,$y_trans),scale(1,$scale)",
                      'style'=>{"font-size"=>$fontsize,
                      'text-anchor'=>'middle',
                      'font-family'=>'Verdana',
                      'fill'=>$color->{uc $letter->[0]}})->cdata(uc $letter->[0]) if $scale > 0;

     $prev_size = $scale * $size;
     $count++;
    }
    #plot the bit if required
    if($self->plot_bits){
         $seq_grp->text('x'=>$x,
                        'y'=>$bit_y,
                        'style'=>{"font-size"=>'10',
                                'text-anchor'=>'middle',
                                'font-family'=>'Verdana',
                                'fill'=>'black'})->cdata($bits);
    }
    $bp++;
    $x+=$width;
  }

  #plot the tags
  $seq_grp->text(x=>int($width/2),y=>$bit_y,style=>{"font-size"=>'10','text-anchor'=>'middle','font-family'=>'Verdana','fill'=>'black'})->cdata("Bits:") if $self->plot_bits;

 $seq_grp->text(x=>int($width/2),y=>$pos_y,style=>{"font-size"=>'10','text-anchor'=>'middle','font-family'=>'Verdana','fill'=>'black'})->cdata("Position:");

  #plot the base positions
  $x = 45+$size/2-int($width/2);
  foreach my $nbr(1..($bp-1)){
    $seq_grp->text(x=>$x+int($width/2),y=>$pos_y,style=>{"font-size"=>'10','text-anchor'=>'left','font-family'=>'Verdana','fill'=>'black'})->cdata($nbr);
    $x+=$width;
  }


#  $seq_grp->transform("scale(2,2)");

  return $self->svg_obj($svg);
}

sub _make_pwm_from_site_matrix{
  my ($self,$matrix) = @_;
  my $bgd = $self->background;
  my @pwm;
  my $consensus = $matrix->consensus;
  foreach my $i(1..length($consensus)){
    my %base = $matrix->next_pos;
    my $bits;
    $bits+=($base{pA} * log2($base{pA}/$bgd->{'A'}));
    $bits+=($base{pC} * log2($base{pC}/$bgd->{'C'}));
    $bits+=($base{pG} * log2($base{pG}/$bgd->{'G'}));
    $bits+=($base{pT} * log2($base{pT}/$bgd->{'T'}));
    push @pwm, [$base{pA},$base{pC},$base{pG},$base{pT},abs(sprintf("%.3f",$bits))];
  }
  return @pwm;
}

sub _make_pwm {
  my ($self,$input) = @_;
  my $count = 1;
  my %hash;
  my $bgd = $self->background;
  #sum up the frequencies at each base pair
  foreach my $seq(@$input){
    my $string = $seq->seq;
    $string =  uc $string;
    my @motif = split('',$string);
    my $pos = 1;
    foreach my $t(@motif){
      $hash{$pos}{$t}++;
      $pos++;
    }
    $count++;
  }

  #calculate relative freq
  my @pwm;

  #decrement last count
  $count--;
  foreach my $pos(sort{$a<=>$b} keys %hash){
    my @array;
    push @array,($hash{$pos}{'A'}||0)/$count;
    push @array,($hash{$pos}{'C'}||0)/$count;
    push @array,($hash{$pos}{'G'}||0)/$count;
    push @array,($hash{$pos}{'T'}||0)/$count;

    #calculate bits
    # relative entropy (RelEnt) or Kullback-Liebler distance
    # relent = sum fk * log2(fk/gk) where fk is frequency of nucleotide k and
    # gk the background frequency of nucleotide k

    my $bits;
    $bits+=(($hash{$pos}{'A'}||0) / $count) * log2((($hash{$pos}{'A'}||0)/$count) / ($bgd->{'A'}));
    $bits+=(($hash{$pos}{'C'}||0) / $count) * log2((($hash{$pos}{'C'}||0)/$count) / ($bgd->{'C'}));
    $bits+=(($hash{$pos}{'G'}||0) / $count) * log2((($hash{$pos}{'G'}||0)/$count) / ($bgd->{'G'}));
    $bits+=(($hash{$pos}{'T'}||0) / $count) * log2((($hash{$pos}{'T'}||0)/$count) / ($bgd->{'T'}));
    push @array, abs(sprintf("%.3f",$bits));

    push @pwm,\@array;
  }
  return $self->pwm(\@pwm);
}


###various get/sets

=head2 fontsize

 Title   : fontsize
 Usage   : $picto->fontsize();
 Function: get/set for fontsize
 Returns : int
 Arguments: int

=cut

sub fontsize {
  my ($self,$obj) = @_;
  if($obj){
    $self->{'_fontsize'} = $obj;
  }
  return   $self->{'_fontsize'};
}

=head2 color

 Title   : color
 Usage   : $picto->color();
 Function: get/set for color
 Returns : a hash reference
 Arguments: a hash  reference

=cut

sub color {
  my ($self,$obj) = @_;
  if($obj){
    $self->{'_color'} = $obj;
  }
  return   $self->{'_color'};
}

=head2 svg_obj

 Title   : svg_obj
 Usage   : $picto->svg_obj();
 Function: get/set for svg_obj
 Returns : L<SVG>
 Arguments: L<SVG>

=cut

sub svg_obj {
  my ($self,$obj) = @_;
  if($obj){
    $self->{'_svg_obj'} = $obj;
  }
  return   $self->{'_svg_obj'};
}

=head2 plot_bits

 Title   : plot_bits
 Usage   : $picto->plot_bits();
 Function: get/set for plot_bits to indicate whether to plot
           information content at each base position
 Returns :1/0
 Arguments: 1/0

=cut

sub plot_bits {
  my ($self,$obj) = @_;
  if($obj){
    $self->{'_plot_bits'} = $obj;
  }
  return   $self->{'_plot_bits'};
}

=head2 normalize

 Title   : normalize
 Usage   : $picto->normalize($newval)
 Function: get/set to make all columns the same height.
           default is to scale height with information
           content.
 Returns : value of normalize (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub normalize{
    my $self = shift;

    return $self->{'normalize'} = shift if @_;
    return $self->{'normalize'};
}

=head2 background

 Title   : background
 Usage   : $picto->background();
 Function: get/set for hash reference of nucleodtide bgd frequencies
 Returns : hash reference
 Arguments: hash reference

=cut

sub background {
  my ($self,$obj) = @_;
  if($obj){
    $self->{'_background'} = $obj;
  }
  return   $self->{'_background'};
}

=head2 pwm

 Title   : pwm
 Usage   : $picto->pwm();
 Function: get/set for pwm
 Returns : int
 Arguments: int

=cut

sub pwm {
  my ($self,$pwm) = @_;
  if($pwm){
    $self->{'_pwm'} = $pwm;
  }
  return $self->{'_pwm'};
}

#utility method for returning log 2
sub log2 {
    my ($val) = @_;
    return 0 if $val==0;
    return log($val)/log(2);
}


1;
