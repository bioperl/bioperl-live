#Author: William McCaig
#Date: 06/16/2006
#Purpose:  To print visual images of alignments
#
#Requires:  An alignment file
#
#Produces:  An image file
#
#Revision History: 
#09/01/2006 - WDM - Introduction of "wrap" flag, allowing alignment to be
#                   wrapped at a set base and stacked vertically
#                   Addition of internal members y_num and y_size for tracking
#                   of number of vertical panels and size of panels,
#                   respectively
#
#09/06/2006 - WDM - Introduction of "p_legend" flag, for printing of an optional
#                   colored legend when protein coloring is selected
#
#09/24/2008 - WDM - Test file created for the module
#
#03/01/2009 - YH -  Introduction of "show_nonsynonymous" flag which enables
#                   highlighting of nonsynonymous mutations in nucleotide
#                   alignments. Addition of internal members codon_table and
#                   missense_pos for translating codons -> amino acids and for
#                   keeping track of missense mutation positions respectively.
#
#03/05/2009 - YH  - Swapped names of subroutines x_label and y_label to match
#                   both documentation and intuition. Finalized implementation
#                   of show_nonsynonymous functionality.

# docs after the code!

package Bio::Align::Graphics;

use vars qw( @PRINT_PARAMS %OK_FIELD);

use 5.008003;
use strict;
use warnings;

use GD;
use GD::Simple;
use Bio::AlignIO;
use Data::Dumper;
use POSIX qw(ceil floor);

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use PrintAlignment ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw( ) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw( );

# Preloaded methods go here.
our %FONT_TABLE = (1 => gdTinyFont, 2 => gdSmallFont, 3 => gdMediumBoldFont, 4 => gdLargeFont, 5 => gdGiantFont );
our %PROTEIN_COLORS = ('Q' => [255, 0, 204], 'E' => [255, 0, 102], 'D' =>  [255, 0, 0] , 'S' => [255, 51, 0] , 'T' => [255, 102, 0 ], 
			'G' => [255, 153, 0] , 'P' => [255, 204, 0] , 'C' => [255, 255, 0] , 'A' => [204, 255, 0] , 'V' => [153, 255, 0],
			'I' => [102, 255, 0] , 'L' => [51 , 255, 0] , 'M' => [0, 255, 0] , 'F' => [0 , 255, 102] , 'Y' => [0 , 255, 204],
			'W' => [0, 204, 255] , 'H' => [0, 102, 255] , 'R' => [0, 0, 255] , 'K' => [102, 0, 255] , 'N' => [204, 0, 255] );
#################################################################
#New
sub new {
my $class = shift;
my %options = @_;

my $self  = {
	
	#####OPTIONS#####
	#Display Defaults
	font => defined($options{font}) ? $FONT_TABLE{$options{font}} : $FONT_TABLE{2},
	x_label => defined($options{x_label}) ? $options{x_label} : 1,
	y_label => defined($options{y_label}) ? $options{y_label} : 1,
	
	#Colors
	bg_color => $options{bg_color} || 'white',
	fg_color => $options{font_color} || 'black',
	x_label_color => $options{x_label_color} || 'blue',
	y_label_color => $options{y_label_color} || 'red',
	p_color => $options{p_color} || undef,
	p_legend => $options{p_legend} || undef,
	p_color_table => undef,
			
	#Sequence Defaults
	reference => $options{reference} || undef,
	reference_id => $options{reference_id} || undef,
	match_char => $options{match_char} || ".",
	block_size => defined($options{block_size}) ? $options{block_size} : 10,
	block_space => defined ($options{block_space}) ? ($options{block_space} * ($options{font} ? $FONT_TABLE{$options{font}}->width : $FONT_TABLE{2}->width)) : ( ($options{font} ? ($FONT_TABLE{$options{font}}->width * 2 ) : ($FONT_TABLE{2}->width * 2)) ),
	wrap => $options{wrap} || 80,
	show_nonsynonymous => $options{show_nonsynonymous} || undef, # If turned on, will highlight nonsynonymous (missense) mutations. Valid only for nucleotide alignments
	
	#Padding
	pad_left => $options{pad_left} || 5, 		#space between x label and border
	pad_right => $options{pad_right} || 5,		#space between end of sequences and border
	pad_top => $options{pad_top} || 5,		#space between y label and border
	pad_bottom => $options{pad_bottom} || 5,	#space between bottom of sequences and border
	x_label_space => $options{x_label_space} || 1, #space between x label and sequences
	y_label_space => $options{y_label_space} || 1, #space between y label and sequences
	
	#Labels
	labels => $options{labels} || undef,
	dm_labels => $options{dm_labels} || undef,
	dm_label_start => $options{dml_start} || undef,
	dm_label_end => $options{dml_end} || undef,
	dm_label_color => $options{dml_color} || undef,
	domain_start => $options{dm_start} || undef,
	domain_end => $options{dm_end} || undef,
	domain_color => $options{dm_color} || undef,
	
	#File Defaults
	align => $options{align} || undef,
	output => $options{output} || undef,
	out_format => $options{out_format} || undef,
			
	####PRIVATE VALUES#####
	
	image => $options{image} || undef,
	seq_format => undef,
	
	#X and Y size of char
	x_char_size => ($options{font} ? $FONT_TABLE{$options{font}}->width : $FONT_TABLE{2}->width),
	y_char_size => ($options{font} ? $FONT_TABLE{$options{font}}->height : $FONT_TABLE{2}->height),
	
	#Image W & H
	width => undef,		#overall width of the image
	height => undef,	#overall height of image
		
	#Sequences 
	sequences => undef,
	seq_ids => undef,
	ref_sequence => undef,
	id_length => 0,
	seq_length => $options{align}->length() || 0,
	no_sequences => $options{align}->num_sequences() || 0,
	seq_start_x => undef,
	seq_start_y => undef,
	start => $options{start} || 1,
	end => $options{end} || $options{align}->length(),
	y_num => undef,
	y_size => undef,
	footer_size => 110,
	footer_start => undef
		
	};

bless ($self, $class);

die "new:Must supply alignment for drawing!\n"
	unless defined ($self->{align});


foreach my $seq ($self->{align}->each_seq) 
{
$self->{id_length} =  ( length($seq->id()) > $self->{id_length} ) ?  length($seq->id()) : $self->{id_length};

		
	
	if( $self->{reference_id} && ($seq->id() eq $self->{reference_id}) )
	{
	 @{$self->{ref_sequence}} = split //, $seq->seq;
	 unshift @{$self->{sequences}}, $seq->seq;
         unshift @{$self->{seq_ids}}, $seq->id();
	 }else
	  {
		push @{$self->{sequences}}, $seq->seq;
		push @{$self->{seq_ids}}, $seq->id();
	  }
	  
	if(!defined($self->{seq_format}))
	{
	 $self->{seq_format} = $seq->alphabet;
	}
}

if(!($self->{reference_id}) )
{
@{$self->{ref_sequence}} = split //, ${$self->{sequences}}[0];
$self->{reference_id} = ${$self->{seq_ids}}[0];
}

$self->{y_num} = ($self->{seq_length} > $self->{wrap}) ? ( sprintf( "%.0f", ( ($self->{seq_length} / $self->{wrap}) + .5) ) ) : 1;
$self->{y_size} = ( ($self->{no_sequences} + $self->{pad_bottom}) * $self->{y_char_size});
$self->{seq_start_x} = ($self->{pad_left} + $self->{id_length} + $self->{x_label_space}) * $self->{x_char_size};

if( defined($self->{show_nonsynonymous}) ) # Extra column changes dimensions
{
	$self->{seq_length_aa} = ($self->{seq_length} / 3) + $self->{seq_length}; # Consider length of sequence plus extra column every 3 nucleotides
	$self->{seq_start_y} = ($self->{pad_top} + length($self->{seq_length_aa}) + $self->{y_label_space}) * $self->{y_char_size};
	$self->{width} = $self->{seq_start_x} + ((( $self->{wrap} / $self->{block_size}) + 1) * $self->{block_space}) + ( ($self->{wrap} + $self->{pad_right}) * ($self->{x_char_size} + 1.2) ) + ( ($self->{seq_length} / 3) * 2); # Needed to add this for width to fit whole sequence on one line
}else
{
	$self->{seq_start_y} = ($self->{pad_top} + length($self->{seq_length}) + $self->{y_label_space}) * $self->{y_char_size};
	$self->{width} = $self->{seq_start_x} + ((( $self->{wrap} / $self->{block_size}) + 1) * $self->{block_space}) + ($self->{wrap} + $self->{pad_right}) * $self->{x_char_size};
}

$self->{footer_start} = $self->{seq_start_y} + $self->{y_size} * $self->{y_num};

if(defined($self->{p_color}) && defined($self->{p_legend}) && $self->{p_legend}){
$self->{height} = $self->{seq_start_y} + $self->{footer_size} + $self->{y_size} * $self->{y_num};
}else{
 $self->{height} = $self->{seq_start_y} + $self->{y_size} * $self->{y_num};
}
$self->{image} = GD::Simple->new($self->{width},$self->{height});
$self->{image}->alphaBlending(1);
$self->{image}->saveAlpha(1);
$self->{image}->bgcolor($self->{bg_color});
$self->{image}->fgcolor($self->{fg_color});
$self->{image}->rectangle(0,0,$self->{width}-1, $self->{height} - 1);
return $self;

} #End new Subroutine#########################################################


sub draw{
my $self = shift;

die "draw:Must supply alignment for drawing!\n"
	unless defined ($self->{align});

if(defined($self->{x_label}) && $self->{x_label})
{
$self->x_label();
}

if(defined($self->{y_label}) && $self->{y_label})
{
$self->y_label();
}


if(defined($self->{domain_start}) && defined($self->{domain_end}) && not defined($self->{p_color}) )
{
$self->_draw_domain();
}

# 
if( defined($self->{show_nonsynonymous}) && ( $self->{seq_format} eq "protein" ) )
{
die "draw:Option show_nonsynonymous only works with Nucleotide alignments!\n";
}elsif  ( defined($self->{show_nonsynonymous}) )
 {
 	$self->{codon_table} = Bio::Tools::CodonTable->new();
 	$self->{missense_pos} = {};
# 	print STDERR "You are using option show_nonsynonymous. Option works best if wrap value is a multiple of 4.\n"
 }

if(defined($self->{p_color}) && $self->{seq_format} eq "protein")
{
$self->_draw_colored_sequences();
	if(defined($self->{p_legend}) && $self->{p_legend})
	{
	 $self->_draw_legend();
	}
}elsif(defined($self->{p_color}) && ($self->{seq_format} ne "protein"))
 {
  die "draw:Option p_color only works with Protein alignments!\n";
 }else
  {
   $self->_draw_sequences();
  }

if(defined($self->{dm_label_start}))
{
$self->_domain_label();
}


 
if($self->{output})
{
  open my $OUTPUT, '>', $self->{output} or die "Could not read file '$self->{output}': $!\n";
  binmode $OUTPUT;
  
	if(defined($self->{out_format}))
	{
		SWITCH: {
		if($self->{out_format} eq "png")  {print $OUTPUT $self->{image}->png;  last SWITCH;}
		if($self->{out_format} eq "jpeg") {print $OUTPUT $self->{image}->jpeg; last SWITCH;}
		if($self->{out_format} eq "gif")  {print $OUTPUT $self->{image}->gif;  last SWITCH;}
		if($self->{out_format} eq "gd")   {print $OUTPUT $self->{image}->gd;   last SWITCH;}
		}

	}else
	{
	 print $OUTPUT $self->{image}->png;
	}
  
  close $OUTPUT;
}else
 {
	binmode STDOUT;
	
	if(defined($self->{out_format}))
	{
		SWITCH: {
		if($self->{out_format} eq "png")  {print STDOUT $self->{image}->png;  last SWITCH;}
		if($self->{out_format} eq "jpeg") {print STDOUT $self->{image}->jpeg; last SWITCH;}
		if($self->{out_format} eq "gif")  {print STDOUT $self->{image}->gif;  last SWITCH;}
		if($self->{out_format} eq "gd")   {print STDOUT $self->{image}->gd;   last SWITCH;}
		}

	}else
	{
	 print STDOUT $self->{image}->png;
	}
  
  
 }#End Output if/else 




#print "Left\tRight\tTop\tBottom\n";
#print $self->{pad_left}, "\t", $self->{pad_right}, "\t", $self->{pad_top}, "\t", $self->{pad_bottom}, "\n";

};

##########################################
#Draws Sequences
sub _draw_sequences{
my $self = shift;

my $block_num = 0;
my $block_total = 0;
my $print_char;


$self->{image}->fgcolor($self->{fg_color});

for (my $i=0; $i < $self->{no_sequences}; $i++) 
{
	
	 my @letters = split //, ${$self->{sequences}}[$i];
	 
	   
	
	 my $y_num = $self->{y_num}; #sprintf( "%.0f", ( ($self->{seq_length} / $self->{wrap}) + .5) ) - 1;
	 my $y_char = $self->{y_size}; #( ($self->{no_sequences} + $self->{pad_bottom}) * $self->{y_char_size});
	 
	for(my $k=0; $k<=$y_num; $k++)
	{
	my $x_char = $k * $self->{wrap};
		
	for (my $j=$x_char; $j <= ( ($x_char + $self->{wrap}) - 1); $j++) 
	{
	last unless defined($letters[$j]);
	
		
		# If show_nonsynonymous is on, and this is the 3rd nucleotide,
		# save the codon and amino acid for comparison
		my ($codon, $aa);
		if ((defined($self->{show_nonsynonymous})) && ((($j+1) % 3) == 0))
		{
			$codon = $letters[$j-2] . $letters[$j-1] . $letters[$j];
			$aa = $self->{codon_table}->translate($codon);
		}
		
		if( $self->{reference} )
		{
			if(${$self->{seq_ids}}[$i] eq $self->{reference_id})
			{
			 $print_char = $letters[$j];
			}else
			 {
				if($letters[$j] eq ${$self->{ref_sequence}}[$j])
				{
				$print_char = $self->{match_char};
				}else
				{
				$print_char = $letters[$j];
				}
			 }
		}else
		 {
		  $print_char = $letters[$j];
		 }
		 
		if( ( ($j + 1) % ($self->{block_size})) == 0)
		{
		 $block_num = $self->{block_space};
		}else
		 {
		  $block_num = 0;
		 }
		 
	#print "J is: $j\n";	 
	#print "Char is: $print_char\n";
	 my $new_x_pos = $self->{seq_start_x} + ( ($j - $x_char) * $self->{x_char_size}) + $block_total;
	 my $new_y_pos = $self->{seq_start_y} + ($i * $self->{y_char_size}) + ($k * $y_char);
	 
	 $new_x_pos += ( ( floor( ($j-$x_char)/3 ) * $self->{x_char_size} ) + 
	 			( ( floor( ($j-$x_char)/3 ) ) * 6 ))
	 			if ( defined($self->{show_nonsynonymous}) );
	 
	 $self->{image}->moveTo( $new_x_pos, $new_y_pos );
	 $self->{image}->font($self->{font});
	 $self->{image}->string($print_char);
	 
	 if ( (defined($self->{show_nonsynonymous})) && ((($j+1) % 3) == 0) )
	 {
	 	$new_x_pos += ($self->{x_char_size} + 3);
	 	$self->{image}->moveTo( $new_x_pos, $new_y_pos );
	 	
	 	# If show_nonsynonymous is on, and this is the 3rd nucleotide
		# on reference, print the amino acid after the nucleotide
	 	if(($self->{reference}) && (${$self->{seq_ids}}[$i] eq $self->{reference_id}))
	 	{
	 		$self->{image}->font(gdMediumBoldFont);
	 		$self->{image}->string($aa);
	 		$self->{image}->font($self->{font});
	 	}elsif ( ( $self->{reference} ) && ( ${$self->{seq_ids}}[$i] ne $self->{reference_id} ) )
	 	{ # In case current sequence is not reference
	 		my $ref_codon = ${$self->{ref_sequence}}[$j-2] .
							${$self->{ref_sequence}}[$j-1] .
							${$self->{ref_sequence}}[$j];
			my $ref_aa = $self->{codon_table}->translate($ref_codon);
					
			if ( $ref_aa eq $aa ) # Synonymous mutation
			{
				$self->{image}->string($self->{match_char});
			}else # Nonsynonymous mutation
			{
				$self->{image}->font(gdMediumBoldFont);
				$self->{image}->string($aa);
				$self->{image}->font($self->{font});
				
				# Highlight nonsynonymous mutations by drawing a rectangle around them
				if ( ( ${$self->{seq_ids}}[$i] ne $self->{reference_id} ) && !( ${$self->{missense_pos}}{$j} ) )
				{
					${$self->{missense_pos}}{$j} = 1;
					$self->{image}->bgcolor(undef);
					$self->{image}->rectangle( $new_x_pos - 2, ( $new_y_pos - ( ( $self->{y_char_size} * ($i+1)) ) ) - 2, ( $new_x_pos + ( $self->{x_char_size} + 1) ), ( $new_y_pos + ( $self->{y_char_size} * ( $self->{no_sequences} - ( $i+1 ) ) ) ) + 2);
					$self->{image}->bgcolor($self->{bg_color});					
				}
			}
	 	}else # No reference sequence defined
	 	{
	 		$self->{image}->string($aa);
	 	}
	 	
	 }
	 
	 if( defined($self->{labels}) && $i == ($self->{no_sequences} - 1))
	 {
	 
		if(${$self->{labels}}{$j + 1})
		{
		my $label = ${$self->{labels}}{$j + 1};
		my $offset = defined($self->{dm_label_start}) ? 3 : 0;
		 $self->{image}->moveTo($self->{seq_start_x} + ( ( ($j - $x_char) + 1.25) * $self->{x_char_size}) + $block_total, $self->{seq_start_y} + (($self->{no_sequences}) * $self->{y_char_size}) + ($k * $y_char) + ( (length($label) + $offset) * ($self->{x_char_size}) ) );
		 $self->{image}->font($self->{font});
		 $self->{image}->angle(-90);
		 $self->{image}->string($label);
		 $self->{image}->angle(0);		
		}
	 }
	 
	 
	 $block_total += $block_num; 
	}
	 $block_total = 0;	
	}

}


}

# WARNING YH - This function has not been modified to work with show_nonsynonymous: needs test data to make sure it will work!
##############################################
#Draw Domain Label
sub _domain_label{
my $self = shift;
my $start_block_total = 0;
my $end_block_total = 0;
my $wrap_block_total = 0;

my $y_char = $self->{y_size};# ( ($self->{no_sequences} + $self->{pad_bottom}) * $self->{y_char_size});

	for(my $i = 0; $i <= $#{$self->{dm_label_start}}; $i++)
	{
		
	my $start = ${$self->{dm_label_start}}[$i];
	my $end = ${$self->{dm_label_end}}[$i];
	
	my $y_num_start = int( $start / $self->{wrap});
	my $y_num_end = int( $end / $self->{wrap});
	
	my $x_num_start;
	
	if($start >= $self->{wrap})
	{
	$x_num_start = ($start % $self->{wrap}) - 1;
	}else
	 {
	  $x_num_start = $start - 1;
	 }
	
	my $x_num_end;	
	
	if($end >= $self->{wrap})
	{
	$x_num_end = ($end % $self->{wrap});
	}else
	 {
	  $x_num_end = $end;
	 }
	
	my $label = ${$self->{dm_labels}}[$i];
	my $color = ${$self->{dm_label_color}}[$i] || ${$self->{dm_label_color}}[-1] || "silver";
	
	my $label_x = (($x_num_end - $x_num_start) / 2) - (length($label) / 2);
	
	my $label_x_start = (($self->{wrap} - $x_num_start) / 2) - (length($label) / 2);
	my $label_x_end = ($x_num_end / 2) - (length($label) / 2);
	
	$start_block_total =  ( ($x_num_start - ($x_num_start % $self->{block_size}) ) / $self->{block_size} ) * $self->{block_space};
	$end_block_total =  ( ($x_num_end - ($x_num_end % $self->{block_size}) ) / $self->{block_size} ) * $self->{block_space}; 
	$wrap_block_total = ( ($self->{wrap} - ( ($self->{wrap} - 1) % $self->{block_size}) ) / $self->{block_size} ) * $self->{block_space};	
	
	$self->{image}->bgcolor($color);
	$self->{image}->fgcolor($color);
		
		if($y_num_start == $y_num_end) #if the label does not cross the wrap line
		{
		 
		 $self->{image}->rectangle( $self->{seq_start_x} + ( ($x_num_start)  * $self->{x_char_size} ) + $start_block_total, $self->{seq_start_y} + (($self->{no_sequences}) * $self->{y_char_size}) + ($y_num_start * $y_char),  $self->{seq_start_x} + (($x_num_end) * $self->{x_char_size}) + $end_block_total, $self->{seq_start_y} + (($self->{no_sequences} + 1) * $self->{y_char_size}) + ($y_num_start * $y_char));	 
		 $self->{image}->fgcolor($self->{fg_color});
		 $self->{image}->bgcolor($self->{bg_color});
	
		 $self->{image}->moveTo( $self->{seq_start_x} + ( ($x_num_start + $label_x) * $self->{x_char_size}) + $start_block_total, $self->{seq_start_y} + (($self->{no_sequences} + 1) * $self->{y_char_size}) + ($y_num_start * $y_char) );
		 $self->{image}->font($self->{font});
		 $self->{image}->string($label);
		}else
		 {
		  $self->{image}->rectangle( $self->{seq_start_x} + ( ($x_num_start)  * $self->{x_char_size} ) + $start_block_total, $self->{seq_start_y} + (($self->{no_sequences}) * $self->{y_char_size}) + ($y_num_start * $y_char),  $self->{seq_start_x} + (($self->{wrap}) * $self->{x_char_size}) + $wrap_block_total, $self->{seq_start_y} + (($self->{no_sequences} + 1) * $self->{y_char_size}) + ($y_num_start * $y_char));	 
		  $self->{image}->rectangle( $self->{seq_start_x} , $self->{seq_start_y} + (($self->{no_sequences}) * $self->{y_char_size}) + ($y_num_end * $y_char),  $self->{seq_start_x} + (($x_num_end) * $self->{x_char_size}) + $end_block_total, $self->{seq_start_y} + (($self->{no_sequences} + 1) * $self->{y_char_size}) + ($y_num_end * $y_char));	 
		  $self->{image}->fgcolor($self->{fg_color});
		  $self->{image}->bgcolor($self->{bg_color});
	
		  $self->{image}->moveTo( $self->{seq_start_x} + ( ($x_num_start + $label_x_start) * $self->{x_char_size}) + $start_block_total, $self->{seq_start_y} + (($self->{no_sequences} + 1) * $self->{y_char_size}) + ($y_num_start * $y_char) );
		  $self->{image}->font($self->{font});
		  $self->{image}->string($label);
		  
		  $self->{image}->moveTo( $self->{seq_start_x} + ( $label_x_end * $self->{x_char_size}), $self->{seq_start_y} + (($self->{no_sequences} + 1) * $self->{y_char_size}) + ($y_num_end * $y_char) );
		  $self->{image}->font($self->{font});
		  $self->{image}->string($label);
		 
		 
		 }

	}

}


##############################################
#Draw Y Label
sub y_label{
my $self = shift;

$self->{image}->fgcolor($self->{y_label_color});

	my $y_num = $self->{y_num}; #sprintf( "%.0f" , (($self->{seq_length} / $self->{wrap}) + .5)) - 1;
	my $y_char = $self->{y_size}; # ( ($self->{no_sequences} + $self->{pad_bottom}) * $self->{y_char_size});
	 
	for(my $k=0; $k<$y_num; $k++)
	{

	 for (my $i=0; $i< $self->{no_sequences}; $i++) 
	 {
	  $self->{image}->moveTo($self->{pad_left}, $self->{seq_start_y} + ($i * $self->{y_char_size}) + ($k * $y_char) );
	  $self->{image}->font($self->{font});
	  $self->{image}->string(${$self->{seq_ids}}[$i]);
	 }
	 
	}


}
#####################################################
#Draw X Label
sub x_label{
my $self = shift;

my $block_num = 0;
my $block_total = 0;
$self->{image}->fgcolor($self->{x_label_color});

my $y_char = $self->{y_size}; # ( ($self->{no_sequences} + $self->{pad_bottom}) * $self->{y_char_size});

for (my $i=1; $i<= $self->{seq_length}; $i++) 
{


	my $y_num = floor( $i / $self->{wrap}); # Used to be int(), but perl documentation advises against this
	my $x_num;
		
	if($i >= $self->{wrap})
	{
	$x_num = ($i % $self->{wrap});
	}else
	 {
	  $x_num = $i;
	 }
	 
    my @digits = split //, reverse($i);
    
    if( ($i % $self->{block_size}) == 0)
	{
         $block_num = $self->{block_space};
	}else
	 {
	  $block_num = 0;
	 }
	
    if( (($i - 1) % $self->{block_size}) == 0)
    {
	for (my $j=0; $j<=$#digits; $j++) 
	{
		
	if ( defined($self->{show_nonsynonymous}) )
	{
		$self->{image}->moveTo($self->{seq_start_x} + $block_total + ( ($x_num-1) * $self->{x_char_size}) + ( ( floor( ($x_num-1)/3 ) * $self->{x_char_size} ) + ( ( floor( ($x_num-1)/3 ) ) * 6 )), ($self->{pad_top} + length($self->{seq_length_aa}) - $j) * $self->{y_char_size} + ($y_num * $y_char));
	}else
	{
		$self->{image}->moveTo($self->{seq_start_x} + $block_total + ( ($x_num-1) * $self->{x_char_size}), ($self->{pad_top} + length($self->{seq_length}) - $j) * $self->{y_char_size} + ($y_num * $y_char));
	}
	
	$self->{image}->font($self->{font});
	$self->{image}->string($digits[$j]);

	}
    }
	if($x_num == 0)
	{
	 $block_total = 0;
	}else
	 {
	  $block_total += $block_num; 
	 }
}

}

####################################################
#Domain Highlighting

sub _draw_domain{
my $self = shift;


my $block_total = 0;
my ($start, $end, $block_num);


my $y_char = $self->{y_size}; # ( ($self->{no_sequences} + $self->{pad_bottom}) * $self->{y_char_size});

for (my $k=0; $k <= $#{$self->{domain_start}}; $k++) 
{

#print STDERR join "\n", GD::Simple->color_names;
my $dmc = $self->{domain_color}[$k] || $self->{domain_color}[-1] || "silver";
$start = ${$self->{domain_start}}[$k] - 1;
$end = ${$self->{domain_end}}[$k] - 1;


			
	
	for (my $i=0; $i < $self->{no_sequences}; $i++) 
	{
			 
		for (my $j = $start; $j <= $end; $j++)
		{
		
		my $y_num = int( $j / $self->{wrap});
		my $x_num;
		
		if($j >= $self->{wrap})
		{
		$x_num = ($j % $self->{wrap});
		}else
		 {
		  $x_num = $j;
		 }
		 
		 #print "J: $j\nXNUM: $x_num\nYNUM: $y_num\n";
			 $block_total =  ( ($x_num - ($x_num % $self->{block_size}) ) / $self->{block_size} ) * $self->{block_space};
						
			
		 $self->{image}->bgcolor($dmc);
		 $self->{image}->fgcolor($dmc);
		 
		 if ( defined($self->{show_nonsynonymous}) )
		 {																																																																																						# NOTE To shade amino acids as well, change $x_num HERE and                                          HERE to $x_num + 1
		 	$self->{image}->rectangle( $self->{seq_start_x} + ( ($x_num ) * $self->{x_char_size} ) + $block_total - 1 + ( ( floor( $x_num / 3 ) * $self->{x_char_size} ) + ( ( floor( $x_num / 3 ) ) * 6 )),   $self->{seq_start_y} + ( $i * $self->{y_char_size} ) - $self->{y_char_size} + ($y_num * $y_char) ,   $self->{seq_start_x} + (($x_num + 1) * $self->{x_char_size}) + $block_total - 1 + ( ( floor( ($x_num)/3 ) * $self->{x_char_size} ) + ( ( floor( ($x_num)/3 ) ) * 6 )),   $self->{seq_start_y} + ( $i * $self->{y_char_size}) + ($y_num * $y_char));
		 }else
		 {
		 	$self->{image}->rectangle( $self->{seq_start_x} + ( ($x_num ) * $self->{x_char_size} ) + $block_total - 1,   $self->{seq_start_y} + ( $i * $self->{y_char_size} ) - $self->{y_char_size} + ($y_num * $y_char) ,   $self->{seq_start_x} + (($x_num + 1) * $self->{x_char_size}) + $block_total - 1,   $self->{seq_start_y} + ( $i * $self->{y_char_size}) + ($y_num * $y_char));
		 }
		 #$self->{image}->rectangle( $self->{seq_start_x} + ( ($j) * $self->{x_char_size} ) + $block_total, $self->{seq_start_y} + ($i - 1 * $self->{y_char_size}), $self->{seq_start_x} + (($j + 1) * $self->{x_char_size}) + $block_total , $self->{seq_start_y} + ( ($i) * $self->{y_char_size})); 
		
		 $self->{image}->fgcolor($self->{fg_color});
		 $self->{image}->bgcolor($self->{bg_color});
		
		
		}
		
		$block_total = 0;
	 }
}

}


sub _draw_colored_sequences{
my $self = shift;

my $block_num = 0;
my $block_total = 0;
my $print_char;
my %colors;

for my $values ( keys %PROTEIN_COLORS)
{
#print STDERR "$values : @{ $PROTEIN_COLORS{$values} }\n";
$colors{$values} = $self->{image}->colorAllocate(@{ $PROTEIN_COLORS{$values} });
}

$self->{p_color_table} = \%colors;

$self->{image}->fgcolor($self->{fg_color});

for (my $i=0; $i < $self->{no_sequences}; $i++) 
{
	
	 my @letters = split //, ${$self->{sequences}}[$i];
	

	my $y_num = $self->{y_num}; #sprintf( "%.0f", ( ($self->{seq_length} / $self->{wrap}) + .5) ) - 1;
	my $y_char = $self->{y_size}; #( ($self->{no_sequences} + $self->{pad_bottom}) * $self->{y_char_size});
	 
	for(my $k=0; $k<=$y_num; $k++)
	{
	 my $x_char = $k * $self->{wrap};
	
		for (my $j=$x_char; $j <= ( ($x_char + $self->{wrap}) - 1); $j++) 
		{
		 last unless defined($letters[$j]);
		
		 $print_char = $letters[$j];
				 
		if( ( ($j + 1) % ($self->{block_size})) == 0)
		{
		 $block_num = $self->{block_space};
		}else
		 {
		  $block_num = 0;
		 }
		 
	#print "Chunk Space: $chunk_space\n";
	 $self->{image}->bgcolor($colors{$print_char});
	 $self->{image}->fgcolor($colors{$print_char});
	 $self->{image}->rectangle( $self->{seq_start_x} + ( ($j - $x_char) * $self->{x_char_size} ) + $block_total - 1   ,   $self->{seq_start_y} + ( $i * $self->{y_char_size} ) + ($k * $y_char) - $self->{y_char_size}    ,   $self->{seq_start_x} + (($j - $x_char + 1) * $self->{x_char_size}) + $block_total - 1 ,   $self->{seq_start_y} + ($k * $y_char) + ( $i * $self->{y_char_size}));
	 $self->{image}->moveTo($self->{seq_start_x} + ( ($j - $x_char) * $self->{x_char_size}) + $block_total, $self->{seq_start_y} + ($k * $y_char) + ($i * $self->{y_char_size}) );
	 $self->{image}->fgcolor($self->{fg_color});
	 $self->{image}->font($self->{font});
	 $self->{image}->string($print_char);
	
	if( defined($self->{labels}) && $i == ($self->{no_sequences} - 1))
	 {
	 
		if(${$self->{labels}}{$j + 1})
		{
		 my $label = ${$self->{labels}}{$j + 1};
		 my $offset = defined($self->{dm_label_start}) ? 3 : 0;
		 $self->{image}->moveTo($self->{seq_start_x} + ( ( ($j - $x_char) + 1.25) * $self->{x_char_size}) + $block_total, $self->{seq_start_y} + (($self->{no_sequences}) * $self->{y_char_size}) + ($k * $y_char) + ( (length($label) + $offset) * ($self->{x_char_size}) ) );
		 $self->{image}->font($self->{font});
		 $self->{image}->angle(-90);
		 $self->{image}->string($label);
		 $self->{image}->angle(0);		
		}
	 }
	 
	 
	 $block_total += $block_num; 
	}
$block_total = 0;
	}
}
}

sub _draw_legend{

my $self = shift;
my $title_font = $FONT_TABLE{3};
my @l_order = ("Negatively Charged", "Positively Charged", "Hydrophobic", "Aromatic", "Found in Loops", "Large Polar Acids");
my %legend = ("Negatively Charged" => ["D" , "E"] , "Positively Charged" => ["K", "R"] , "Hydrophobic" => ["A","F","I","L","M","V","W","Y"] ,
		"Aromatic" => ["F", "H", "W", "Y"] , "Found in Loops" => ["D", "G", "P", "S", "T"] , "Large Polar Acids" => ["H", "K", "N", "Q", "R"]);

my $x1 = 2;
my $x2 = 42;

my $colors = $self->{p_color_table};

my $y_start = $self->{footer_start};
my $label = "Protein Color Legend";
$self->{image}->bgcolor($self->{bg_color});
$self->{image}->fgcolor($self->{fg_color});
$self->{image}->rectangle(1,$y_start, 70 * $self->{x_char_size}, $self->{height} - 2);

$self->{image}->moveTo((35 - (length($label) / 2) ) * $self->{x_char_size} , $y_start + $self->{y_char_size});
$self->{image}->font($title_font);
$self->{image}->string($label);

my $count = 3;

foreach my $c_label (@l_order)
{

if( ($count % 2) == 0)
{

$self->{image}->moveTo( $x2 *  $self->{x_char_size}, $y_start + ( ($count - 1) * $self->{y_char_size}));
$self->{image}->font($self->{font});
$self->{image}->string($c_label);
	my $i = 0;
	foreach my $chars(@{$legend{$c_label}})
	{
	 $self->{image}->bgcolor($$colors{$chars});
	 $self->{image}->fgcolor($$colors{$chars});
	 $self->{image}->rectangle( ($x2 + 20 + $i) * $self->{x_char_size}, $y_start + ( ($count - 2) * $self->{y_char_size}), ($x2 + 20 + $i + 1) * $self->{x_char_size}, $y_start + ( ($count -1) * $self->{y_char_size}));
	 $self->{image}->bgcolor($self->{bg_color});
	 $self->{image}->fgcolor($self->{fg_color});
	 $i++;
	}

}else
 {
  $self->{image}->moveTo($x1 * $self->{x_char_size} , $y_start + ($count * $self->{y_char_size}));
  $self->{image}->font($self->{font});
  $self->{image}->string($c_label);
	my $i = 0;
	foreach my $chars(@{$legend{$c_label}})
	{
	 $self->{image}->bgcolor($$colors{$chars});
	 $self->{image}->fgcolor($$colors{$chars});
	 $self->{image}->rectangle( ($x1 + 20 + $i) * $self->{x_char_size}, $y_start + ( ($count - 1) * $self->{y_char_size}), ($x1 + 20 + $i + 1) * $self->{x_char_size}, $y_start + ( ($count) * $self->{y_char_size}));
	 $self->{image}->bgcolor($self->{bg_color});
	 $self->{image}->fgcolor($self->{fg_color});
	 $i++;
	}
 }

$count += 1;
}

}
########################################
#####ACCESSORS#####
sub width{
my $self = shift;
return $self->{image}->width if exists $self->{image};
}

sub height{
my $self = shift;
return $self->{image}->height if exists $self->{image};
}

sub aln_length{
my $self = shift;
return $self->{seq_length} if exists $self->{seq_length};
}

sub aln_format{
my $self = shift;
return $self->{seq_format} if exists $self->{seq_format};
}

sub no_sequences{
my $self = shift;
return $self->{no_sequences} if exists $self->{no_sequences};
}

1;
__END__

=head1 NAME

Bio::Align::Graphics - Graphic Rendering of Bio::Align::AlignI Objects

=head1 SYNOPSIS

  use Bio::Align::Graphics;

  #Get an AlignI object, usually by using Bio::AlignIO

  my $file=shift @ARGV;
  my $in=new Bio::AlignIO(-file=>$file, -format=>'clustalw');
  my $aln=$in->next_aln();


  #Create a new Graphics object
  my $print_align = new Bio::Align::Graphics(align => $aln);

  #Draw the alignment
  $print_align->draw();


=head1 DESCRIPTION

Bio::Align::Graphics is a module designed to create image files out of Bio::Align::AlignI objects.  An alignment may be manipulated with various 
formatting and highlighting options.

An example:

	#!/usr/bin/perl -w

	use Bio::AlignIO;
	use Bio::Align::Graphics;
	use strict;
	
	#Get an alignment file
	my $file = shift @ARGV;
	
	#Create an AlignI object using AlignIO
	my $in=new Bio::AlignIO(-file=>$file, -format=>'clustalw');

	#Read the alignment
	my $aln=$in->next_aln();

	#Create some domains for highlighting
	my @domain_start = ( 25 , 50, 80 );
	my @domain_end = ( 40 , 60 , 100 );
	my @domain_color = ( 'red' , 'cyan' , 'green' );
	
	#Create Labels for the domains
	my @dml = ("CARD", "Proline Rich", "Transmembrane");
	my @dml_start = (25, 50, 80);
	my @dml_end = (40, 60, 100);
	my @dml_color = ("lightpink", "lightblue", "lightgreen");
	
	
	#Create individual labels
	my %labels = ( 145 => "Hep-c target");
	
	
	my $print_align = new Bio::Align::Graphics( align => $aln,
					pad_bottom => 5,
					domain_start => \@domain_start,
					domain_end => \@domain_end,
					dm_color => \@domain_color,
					dm_labels => \@dml,
					dm_label_start => \@dml_start,
					dm_label_end => \@dml_end,
					dm_label_color => \@dml_color,
					labels => \%labels,
					out_format => "png");
					
	$print_align->draw();

=head1 METHODS

This section describes the class and object methods for
Bio::Align::Graphics.

Typically you will begin by creating a Bio::Align::Graphics 
object, passing it an alignment object created using Bio::AlignIO.
The Bio::Align::Graphics-E<gt>new() method has a number of 
configuration variables that allow you to control the appearance
of the final image.

You will then call the draw() method to output the final image.

=head1 CONSTRUCTORS

new() is the constructor for Bio::Align::Graphics:

=over 4

=item $print_align = Bio::Align::Graphics-E<gt>new(@options)

The new() method creates a new graphics object.  The options are
a set of tag/value pairs as follows:

  Option         Value                                  Default
  ------         -----                                  -------

  align		 Bio::AlignI object                     None, must be 
						        supplied to draw
						        an alignment

  output	 Filename to print image to	        STDOUT

  out_format	 png, jpeg, gif, gd		        png

  font		 Size of font, ranging from 1 to 5      2
		 and equal to the standard GD fonts
		 ranging from gdTinyFont to 
		 gdGiantFont

  x_label	 Draws a scale numbering alignment      true
		 bases along top of image, every x
		 bases are numbered, where x is the
		 block_size option

  y_label	 Draws sequence ids of alignment        true
		 along left side of image

  bg_color	 Background color of the image	        white

  font_color	 Color of the font used for drawing     black
		 the alignment characters

  x_label_color  Color of the font used for drawing     red
		 the base scale characters

  y_label_color  Color of the font used for drawing     blue
		 the sequence id characters

  p_color	 Colors protein bases according to      false
		 a coloring scheme proposed by W.R.
		 Taylor(Protein Engineering, vol 10
		 no 7, 1997), only works with
		 protein alignments

  pad_top	 Additional whitespace characters       5
		 between top of image and x-label

  pad_bottom	 Additional whitespace characters       5
		 between bottom of image and
		 alignment

  pad_left	 Additional whitespace characters       5
		 between left side of image and 
		 y-label

  pad_right	 Additional whitespace characters       5
		 between right side of image and 
		 alignment

  x_label_space  Additional whitespace characters       1
		 between x_label and alignment

  y_label_space  Additional whitespace characters       1
		 between y_label and alignment

  reference	 Characters which are identical to      false
		 the reference sequence are replaced
		 with the match character

  reference_id	 Sequence id of the sequence to use     First sequence
		 as the reference			supplied in alignment

  match_char	 Character to replace identical bases   .
		 in aligned sequences

  block_size	 Number of bases to group together	10
		 when printing alignment, groups are
		 separated by whitespace

  block_space	 Amount of character whitespace to	2
		 separate groups of bases by

  labels	 A hash containing labels to be 	none
		 printed beneath the alignment, 
		 where the keys are the bases to
		 print the values at

  dm_start	 An array containing start bases	none
		 for highlighting of segments of
		 the alignment, paired with dm_end
		 option

  dm_end	 An array containing end bases		none
		 for highlighting of segments of
		 the alignment, paired with dm_start
		 options

  dm_color	 An array containing colors for	        silver
		 highlighting segments of bases
		 denoted by the coordinates
		 located in the dm_start and dm_end
		 options

  dml_start	 An array containing start bases	none
		 for addition of domain labels
		 underneath the alignment, paired
		 with dml_end

  dml_end	 An array containing end bases		none
		 for addition of domain labels
		 underneath the alignment, paired
		 with dml_start

  dml_color	 An array containing colors for 	silver
		 the domain labels denoted by the
		 coordinates located in the 
		 dml_start and dml_end options

  dm_labels	 An array containing labels to be	none
		 printed underneath specified
		 domains, each label should
		 correspond with the base position
		 located in the dml_start option
		 
  show_nonsynonymous  Boolean value to turn option	false
  		 on or off. If 0 (or undef), option
  		 is off. If 1 (or non-0), option is on.
  		 Only valid for nucleotide alignments.
  		 Output images are wider with this option on.

Note that all arrays and hashes must be passed by reference.

=back

=head1 OBJECT METHODS

=over 4

=item $draw_align-E<gt>draw();

The draw() method draws the image with the options that were specified with new().

=item $draw_align-E<gt>width();

Get the width of the image created with new(), in pixels.

=item $draw_align-E<gt>height();

Get the height of the image created with new(), in pixels.

=item $draw_align-E<gt>aln_length();

Get the length of the alignment submitted to new().

=item $draw_align-E<gt>aln_format();

Get the format of the alignment submitted to new().

=item $draw_align-E<gt>no_sequences();

Get the number of sequences in the alignment submitted to new().

=back

=head1 AUTHORS AND CONTRIBUTORS

William McCaig, E<lt>wmccaig@gmail.comE<gt>

Mikhail Bekarev, E<lt>mbekarev@hunter.cuny.eduE<gt>

YE<246>zen HernE<225>ndez, E<lt>yzhernand@gmail.comE<gt>

Weigang Qiu (Corresponding Developer), E<lt>weigang@genectr.hunter.cuny.eduE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2006-2008 by William McCaig

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.3 or,
at your option, any later version of Perl 5 you may have available.

=head1 SEE ALSO

L<Bio::Align::AlignI>,
L<Bio::AlignIO>,
L<GD>,
L<GD::Simple>

=cut
