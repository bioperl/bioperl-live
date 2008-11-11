package Bio::Graphics::Glyph::trace;

# $Id: trace.pm,v 1.7 2008/09/20 15:46:39 lstein Exp $

use strict;
use GD;
use Bio::SCF;
use File::Temp qw( tempdir );
use Digest::MD5 qw( md5_hex );
use base 'Bio::Graphics::Glyph::generic';
our @ISA;

use constant VERTICAL_SPACING => 20;
my %complement = (
    g => 'c',
    a => 't',
    t => 'a',
    c => 'g',
    n => 'n',
    G => 'C',
    A => 'T',
    T => 'A',
    C => 'G',
    N => 'N'
);

sub new {
    my $self = shift->SUPER::new(@_);

    if ( $self->dna_fits ) {
        $self->{parsed_trace} = $self->get_parsed_trace();
    }
    return $self;
}

sub get_parsed_trace {
    my $self = shift;
    my ( $format, $trace_file ) = eval { $self->trace_data };
    unless ($trace_file) {
        warn $@ if $@;
        return;
    }

    if ($self->{'content_type'} eq "ABI"){
        require ABI;
        my $abi = ABI->new(-file=>$trace_file);
        my %scf;
        @{$scf{'bases'}}        = split //, $abi->get_sequence();
        @{$scf{'index'}}        = $abi->get_base_calls();
        @{$scf{'samples'}{'A'}} = $abi->get_trace('A');
        @{$scf{'samples'}{'C'}} = $abi->get_trace('C');
        @{$scf{'samples'}{'G'}} = $abi->get_trace('G');
        @{$scf{'samples'}{'T'}} = $abi->get_trace('T');

        my $scale_factor = $self->option('abi_scale') || 1;
        $self->{'max_trace'} = $scale_factor * $abi->get_max_trace();

        return \%scf;
    }
    else{
        my %scf;
        tie %scf, 'Bio::SCF', $trace_file;
        $self->{'max_trace'} = 1600;
        return \%scf;
    }
}

sub _guess_format {
    my $self = shift;
    my $path = shift;

    my $modded_path = $path;
    $self->{'gzipped'} = ($modded_path =~ s/\.gz\s*$// );
    if ($modded_path =~ /\.scf/){
        $self->{'content_type'} = 'Bio::SCF';
    }
    elsif ($modded_path =~ /\.ab1/){
        $self->{'content_type'} = 'ABI';
    }
    else{
        die "$path: Trace file not of recognized format\n"
    }

    return;
}

sub trace_path {
    my $self     = shift;
    my $feature  = $self->feature or die "no feature!";
    my $dirname  = $self->trace_dir;
    my $basename = $self->option('trace');

    # can't get it from callback, so try looking for an 'trace' attribute
    if ( !$basename ) {
        if ( $feature->can('attributes') ) {
            ($basename) = $feature->attributes('trace');
        }
        elsif ( $feature->can('has_tag') && $feature->has_tag('trace') ) {
            ($basename) = $feature->get_tag_values('trace');
        }
    }

    return unless $basename;
    return $basename if $basename =~ m!^\w+:/!;    # looks like a URL
    return $basename if $basename =~ m!^/!;        # looks like an abs path
    return "$dirname/$basename";
}

sub trace_data {
    my $self = shift;
    my $path = $self->trace_path;
    $self->_guess_format($path);

    if ( $path =~ m!^\w+:/! ) {                    # looks like a URL
        require LWP::UserAgent;
        my $ua       = LWP::UserAgent->new;
        my $response = $ua->get($path);
        if ( $response->is_success ) {

            # In the future, make extensible to ABI format
            my $data      = $response->content;
            my $signature = md5_hex($data);
            my $extension;
            if ($self->{'content_type'} eq 'ABI'){
                $extension = 'ab1';
            }
            else{
                $extension = 'scf';
            }
            if ($self->{'gzipped'}){
                $extension .= '.gz';
            }

            # untaint signature for use in open
            $signature =~ /^([0-9A-Fa-f]+)$/g or return;
            $signature = $1;

            my $dir_path = tempdir();
            my $file_name
                = sprintf( "%s/%s.%s", $dir_path, $signature, $extension );
            open( F, ">$file_name" )
                || die("Can't open file $file_name for writing: $!\n");
            binmode(F);
            print F $data;
            close F;
            
            if ($self->{'gzipped'}){
                $file_name = $self->gunzip_file( $file_name );
            }

            return ( $self->{'content_type'}, $file_name );
        }
        else {
            die $response->status_line;
        }

    }
    else {
        if ($self->{'gzipped'}){
            $path = $self->gunzip_file( $path );
        }
        return ( $self->{'content_type'}, $path );
    }
}

sub gunzip_file {
    my $self      = shift;
    my $file_name = shift;

    $file_name =~ /(.+)\.gz$/;
    my $new_file_name = $1;

    unless ( -e $new_file_name ){
        `gunzip -c $file_name > $new_file_name`;
    }

    return $new_file_name;
}

sub trace_height {
    my $self = shift;
    return $self->{trace_height} if exists $self->{trace_height};
    return $self->{trace_height} = $self->option('trace_height')
        || 90;    # what the factory says
}

sub pad_left {
    my $self = shift;
    my $pad  = $self->SUPER::pad_left;
    if ( $self->dna_fits ) {
        my $width_needed = ( 20 - $self->width ) / 2;    #BF FIX ME
        $pad = $pad > $width_needed ? $pad : $width_needed;
    }
    return $pad;
}

sub pad_right {
    my $self = shift;
    my $pad  = $self->SUPER::pad_right;

    if ( $self->dna_fits ) {
        my $width_needed = ( 20 - $self->width ) / 2;    #BF FIX ME
        $pad = $pad > $width_needed ? $pad : $width_needed;
    }
    return $pad;
}

sub pad_bottom {
    my $self = shift;
    my $pb = 0;
    if ( $self->dna_fits ) {
        $pb += $self->vertical_spacing;
        $pb += $self->trace_height;
    }
    else{
        $pb = $self->SUPER::pad_bottom;
    }
    return $pb;
}

sub vertical_spacing {
    my $self = shift;
    my $vs   = $self->option('vertical_spacing');
    return $vs if defined $vs;
    return VERTICAL_SPACING;
}

sub draw_description {
    my $self = shift;
    my ( $gd, $left, $top, $partno, $total_parts ) = @_;

    $self->SUPER::draw_description( $gd, $left, $top, $partno, $total_parts );
}

sub draw_label {
    my $self = shift;
    my ( $gd, $left, $top, $partno, $total_parts ) = @_;
    $left += $self->pad_left;
    $self->SUPER::draw_label( $gd, $left, $top, $partno, $total_parts );
}

sub trace_dir {
    my $self = shift;
    return $self->option('trace_prefix');
}

sub draw_component {
    my $self = shift;
    my $gd   = shift;
    my ( $x1, $y1, $x2, $y2 ) = $self->bounds(@_);

    # Draw the regular glyph
    unless ( $self->dna_fits ) {
        my $delegate = $self->option('glyph_delegate') || 'generic';
        if ( $delegate eq 'generic' ) {
            $self->SUPER::draw_component( $gd, @_ );
        }
        else {
            eval "require Bio::Graphics::Glyph::$delegate";
            local @ISA = ("Bio::Graphics::Glyph::$delegate");
            my $method = "Bio::Graphics::Glyph::${delegate}::draw_component";
            $self->$method( $gd, @_ );
        }
        return;
    }

    # Draw Trace

    my $fgcolor = $self->fgcolor;
    my $bgcolor = $self->bgcolor;

    my $parsed_trace    = $self->{parsed_trace} or return;
    my $feature         = $self->feature;
    my $pixels_per_base = $self->scale;
    my $panel           = $self->panel;
    my $strand          = $feature->strand;
    my $forward         = $self->{flip} ? ( $strand < 0 ) : ( $strand >= 0 );
    my $flipped         = $self->{flip};
    my $opp_strand      = $strand < 0;

    # Get Window Sequence Information
    my $panel_start_base = $panel->offset + 1;
    my $panel_end_base   = $panel_start_base + $panel->length - 1;

    my ( $feature_display_start_base, $feature_display_end_base,
        $feature_display_center_base );

    if ( $panel_start_base >= $feature->start() ) {
        $feature_display_start_base = $panel_start_base;
    }
    else {
        $feature_display_start_base = $feature->start();
    }

    if ( $panel_end_base <= $feature->end() ) {
        $feature_display_end_base = $panel_end_base;
    }
    else {
        $feature_display_end_base = $feature->end();
    }

    # We need to know if there are an even number of bases
    # because the center is going to be off by a bit.
    my $even_number_of_bases = 0;
    $feature_display_center_base
        = ( $feature_display_start_base + $feature_display_end_base ) / 2;
    unless (
        $feature_display_center_base == int($feature_display_center_base) )
    {
        $even_number_of_bases = 1;
        $feature_display_center_base
            = int( 0.5 + $feature_display_center_base );
    }

    my $trace_glyph_top    = $y1;
    my $trace_glyph_bottom = $trace_glyph_top + $self->trace_height;

    # Find the Center for the Trace Glyph
    my $trace_center_base_index =
        ( !$opp_strand )
        ? $feature_display_center_base - $feature->start
        : $feature->end - $feature_display_center_base;

    my $trace_center_px = $self->panel->left
        + $self->trace_map_pt($feature_display_center_base);

    # Center base test lines
    #$gd->line( $trace_center_px, 0, $trace_center_px, 700, $fgcolor );
    #if ( !$flipped ) {
    #    $gd->line(
    #        $trace_center_px + $pixels_per_base, 0,
    #        $trace_center_px + $pixels_per_base, 700,
    #        $self->factory->translate_color('red')
    #    );
    #}
    #else {
    #    $gd->line(
    #        $trace_center_px - $pixels_per_base, 0,
    #        $trace_center_px - $pixels_per_base, 700,
    #        $self->factory->translate_color('red')
    #    );
    #}

    # Figure out the number of bases to display on each side
    # with respect to the trace.
    my $five_prime_bases
        = $feature_display_center_base - $feature_display_start_base;
    my $three_prime_bases
        = $feature_display_end_base - $feature_display_center_base;
    if ($opp_strand) {
        ( $five_prime_bases, $three_prime_bases )
            = ( $three_prime_bases, $five_prime_bases );
    }

    # Work out the starting base on each side
    my $trace_start_base_index = $trace_center_base_index - $five_prime_bases;
    if ( $trace_start_base_index < 0 ) {
        $five_prime_bases
            += $trace_start_base_index;   # trace_start_base_index is negative
        $trace_start_base_index = 0;
    }
    my $trace_end_base_index = $trace_center_base_index + $three_prime_bases;
    if ( $trace_end_base_index >= scalar @{ $parsed_trace->{bases} } ) {
        $three_prime_bases
            += scalar @{ $parsed_trace->{bases} } - $trace_end_base_index - 1;
        $trace_end_base_index = scalar @{ $parsed_trace->{bases} } - 1;
    }

    # Figure out the end points of the trace section
    my ( $trace_left_px, $trace_right_px );
    if ($forward) {
        $trace_left_px
            = $trace_center_px - ( ( $five_prime_bases * $pixels_per_base ) );
        $trace_right_px
            = $trace_center_px + ( ( $three_prime_bases * $pixels_per_base ) )
            + $pixels_per_base;
    }
    else {
        $trace_left_px = $trace_center_px
            - ( ( $three_prime_bases * $pixels_per_base ) );
        $trace_right_px
            = $trace_center_px + ( ( $five_prime_bases * $pixels_per_base ) )
            + $pixels_per_base;
    }

    # Adjust for flipping
    if ($flipped) {
        $trace_left_px  -= $pixels_per_base;
        $trace_right_px -= $pixels_per_base;
    }

    # Get Text Info
    my $font        = $self->font;
    my $text_buffer = 2;
    my $text_height = $font->height + ( $text_buffer * 2 );

    my $trace_base_line = $trace_glyph_bottom - $text_height;
    my $max_trace_val   = $self->{'max_trace'} ;
    my $vertical_scale
        = ( $self->trace_height - $text_height - 2 ) / $max_trace_val;
    my $total_trace_bases = $three_prime_bases + $five_prime_bases + 1;

    my $trace_start_sample
        = int( $parsed_trace->{index}[$trace_start_base_index] );
    my $trace_end_sample
        = int( $parsed_trace->{index}[$trace_end_base_index] );
    my $trace_center_sample
        = $parsed_trace->{index}[$trace_center_base_index];

    my %base_colors = (
        'A' => $self->factory->translate_color(
            $self->option('a_color') || 'green'
        ),
        'C' => $self->factory->translate_color(
            $self->option('c_color') || 'blue'
        ),
        'G' => $self->factory->translate_color(
            $self->option('g_color') || 'black'
        ),
        'T' => $self->factory->translate_color(
            $self->option('t_color') || 'red'
        ),
    );
    my $current_px;
    my $current_height;

    my $trace_center_base_px = $trace_center_px + int( $pixels_per_base / 2 );
    if ($flipped) {
        $trace_center_base_px -= $pixels_per_base;
    }

    # Draw Trace
    my $horizontal_scale_5p = ( $pixels_per_base
            * ( $trace_center_base_index - $trace_start_base_index ) ) /
        ( $parsed_trace->{index}[$trace_center_base_index]
            - $parsed_trace->{index}[$trace_start_base_index] + 1 );
    my $horizontal_scale_3p = ( $pixels_per_base
            * ( $trace_end_base_index - $trace_center_base_index ) ) /
        ( $parsed_trace->{index}[$trace_end_base_index]
            - $parsed_trace->{index}[$trace_center_base_index] + 1 );
    $trace_start_sample
        -= int( ( $pixels_per_base / 2 ) / $horizontal_scale_5p )
        if $horizontal_scale_5p;
    $trace_end_sample
        += int( ( $pixels_per_base / 2 ) / $horizontal_scale_3p )
        if $horizontal_scale_3p;
    my $last_px = $self->_get_pixel_position_x(
        current_sample       => $trace_start_sample,
        trace_center_sample  => $trace_center_sample,
        forward              => $forward,
        trace_center_base_px => $trace_center_base_px,
        horizontal_scale_5p  => $horizontal_scale_5p,
        horizontal_scale_3p  => $horizontal_scale_3p,
    );
    my %last_heights;
    foreach my $base ( keys %base_colors ) {
        $last_heights{$base} = int(
            $trace_base_line - (
                $vertical_scale
                    * $parsed_trace->{samples}{$base}[$trace_start_sample]
            )
        );
    }
    my $passed_center    = 0;
    my $horizontal_scale = $horizontal_scale_5p;
    for (
        my $current_sample = $trace_start_sample + 1;
        $current_sample <= $trace_end_sample;
        $current_sample++
        )
    {
        $current_px = $self->_get_pixel_position_x(
            current_sample       => $current_sample,
            trace_center_sample  => $trace_center_sample,
            forward              => $forward,
            trace_center_base_px => $trace_center_base_px,
            horizontal_scale_5p  => $horizontal_scale_5p,
            horizontal_scale_3p  => $horizontal_scale_3p,
        );
        $self->_draw_trace_sample(
            gd              => $gd,
            last_px         => $last_px,
            current_px      => $current_px,
            current_sample  => $current_sample,
            last_heights    => \%last_heights,
            base_colors     => \%base_colors,
            parsed_trace    => $parsed_trace,
            vertical_scale  => $vertical_scale,
            trace_base_line => $trace_base_line,
            forward         => $forward,
        );
        $last_px = $current_px;
    }

    # Print Trace Sequence
    my $base_count = 0;
    my $seq_y      = $trace_glyph_bottom - $text_height + $text_buffer;
    for (
        my $base_index = $trace_start_base_index;
        $base_index <= $trace_end_base_index;
        $base_index++
        )
    {

        if ($forward) {
            my $x;
            if ($opp_strand) {
                $x = $trace_left_px + ( $base_count * $pixels_per_base )
                    + $pixels_per_base - $font->width - 1;
            }
            else {
                $x = $trace_left_px + $base_count * $pixels_per_base;
            }
            my $base = $parsed_trace->{bases}[$base_index];
            my $color = $base_colors{$base} || $fgcolor;
            $gd->char( $font, $x + 2, $seq_y, $base, $color );
        }
        else {
            my $x;
            if ($opp_strand) {
                $x = $trace_right_px
                    - ( ( $base_count + 1 ) * $pixels_per_base );
            }
            else {
                $x = $trace_right_px
                    - ( ( $base_count + 1 ) * $pixels_per_base )
                    + $pixels_per_base - $font->width - 1;
            }
            my $base = $parsed_trace->{bases}[$base_index];
            $base = $complement{$base} || $base;
            my $color = $base_colors{$base} || $fgcolor;
            $gd->char( $font, $x + 2, $seq_y, $base, $color );

        }

        $base_count++;
    }

    if ( $self->option('show_border') ) {
        # Outline Box
        #  left side
        $gd->line(
            $trace_left_px,      $trace_glyph_top, $trace_left_px,
            $trace_glyph_bottom, $fgcolor
        );

        #  right side
        $gd->line(
            $trace_right_px,     $trace_glyph_top, $trace_right_px,
            $trace_glyph_bottom, $fgcolor
        );

        #  top
        $gd->line(
            $trace_left_px,   $trace_glyph_top, $trace_right_px,
            $trace_glyph_top, $fgcolor
        );

        #  bottom
        $gd->line(
            $trace_left_px,      $trace_glyph_bottom, $trace_right_px,
            $trace_glyph_bottom, $fgcolor
        );
    }

}

sub _get_pixel_position_x {
    my $self                 = shift;
    my %args                 = @_;
    my $current_sample       = $args{'current_sample'};
    my $trace_center_sample  = $args{'trace_center_sample'};
    my $forward              = $args{'forward'};
    my $trace_center_base_px = $args{'trace_center_base_px'};
    my $horizontal_scale_5p  = $args{'horizontal_scale_5p'};
    my $horizontal_scale_3p  = $args{'horizontal_scale_3p'};

    my $horizontal_scale;
    if ( $current_sample >= $trace_center_sample ) {
        $horizontal_scale = $horizontal_scale_3p;
    }
    else {
        $horizontal_scale = $horizontal_scale_5p;
    }

    if ($forward) {
        return $trace_center_base_px +
            int(
            $horizontal_scale * ( $current_sample - $trace_center_sample ) );
    }
    else {
        return $trace_center_base_px -
            int(
            $horizontal_scale * ( $current_sample - $trace_center_sample ) );
    }

}

sub _draw_trace_sample {
    my $self            = shift;
    my %args            = @_;
    my $gd              = $args{'gd'};
    my $last_px         = $args{'last_px'};
    my $current_px      = $args{'current_px'};
    my $current_sample  = $args{'current_sample'};
    my $last_heights    = $args{'last_heights'};
    my $base_colors     = $args{'base_colors'};
    my $parsed_trace    = $args{'parsed_trace'};
    my $vertical_scale  = $args{'vertical_scale'};
    my $trace_base_line = $args{'trace_base_line'};
    my $forward         = $args{'forward'};

    my $current_height;
    foreach my $base ( keys %{ $base_colors || {} } ) {
        if (   $current_sample < 0
            or $current_sample
            >= scalar @{ $parsed_trace->{samples}{$base} } )
        {

            # Off the end of the trace
            $current_height = 0;
        }
        else {
            $current_height = int(
                $trace_base_line - (
                    $vertical_scale
                        * $parsed_trace->{samples}{$base}[$current_sample]
                )
            );
        }
        my $color =
            ($forward)
            ? $base_colors->{$base}
            : $base_colors->{ $complement{$base} || $base };
        $gd->line( $last_px, $last_heights->{$base},
            $current_px, $current_height, $color );

        $last_heights->{$base} = $current_height;
    }
}

sub trace_map_pt {
    my $self   = shift;
    my $panel  = $self->panel;
    my $offset = $panel->{offset};
    my $scale  = $panel->{scale} || $panel->scale;
    my $pl     = $panel->{pad_left};
    my $width  = $panel->{width};
    my $flip   = $panel->{flip};
    my $length = $panel->{length};
    my @result;

    foreach (@_) {
        my $val = $flip

            #? int (0.5 + $self->{width} - ($length - ($_- 1)) * $scale)
            ? int( 0.5 + $width - ( $_ - $offset - 1 ) * $scale )
            : int( 0.5 + ( $_ - $offset - 1 ) * $scale );
        $val = -1         if $val < 0;
        $val = $width + 1 if $val > $width;
        push @result, $val;
    }
    return (wantarray) ? @result : $result[0];
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::trace - A glyph that visualizes a trace file

=head1 SYNOPSIS

 use Bio::Graphics;
 use Bio::Seq;
 use Bio::SeqFeature::Generic;

 my $bsg = 'Bio::SeqFeature::Generic';

 my $seq    = Bio::Seq->new(-length=>1000);

 my $whole  = $bsg->new(-display_name => 'Clone82',
 		        -start        => 1,
		        -end          => $seq->length);

 my $trace1 = $bsg->new(-start        => 100,
		        -end          => 300,
		        -display_name => 'Excretory System',
		        -tag=>{
			      trace=>"/path/to/trace/file.scf"
			      }
		       );

 my $trace2 = $bsg->new(-start        => 500,
		        -end          => 800,
		        -display_name => 'Expression Pattern',
		        -tag=>{
			      trace=>"http://localhost/traces/file2.scf"
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

 $panel->add_track([$trace1,$trace2],
		   -glyph    => 'trace',
		   -label    => 1,
		   -key       => 'Example traces');

 binmode STDOUT;
 print $panel->png;

=head1 DESCRIPTION

This glyph parses and displays trace information from a file.  A generic glyph
is used to show where the trace is located and when the display is zoomed in
enough to see the sequence, the trace will be drawn.

The trace file can only be in SCF format.  The file can be located on the local
filesystem or located at a remote URL (provided that you have the LWP module
installed).

Until an alignment feature is added to this glyph, the feature start and end
must correspond exactly with the begining and end of the called sequence.
Meaning that even if the starting sequence is poor and doesn't match the
sequence, it must still be included.

The figure below illustrates this.  The trace and the reference sequence align
from points "b" to "c".  The positions "A" and "D" need to be calculated and
used in order for the trace to line up correctly.  

             A      b         c     D  
  ref -------------------------------------------------
                    |||||||||||
  trace      ------------------------

The glyph may be modified in the future to avoid this hassle (and it should
still be compatible with the method described above).

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

  -trace            Specify the trace path or URL     none
                    to use for this feature

  -trace_prefix     String to prepend to              none
                    each trace path. You may prepend
                    a directory or a partial URL.

  -trace_height     The height in pixels that the     90
                    trace will be drawn

  -vertical_spacing Vertical distance from the box    20
                    that shows the physical span of
                    the feature to the top of the
                    picture (in pixels)

  -glyph_delegate   Glyph to use when zoomed out too  'generic'
                    far for the trace to be drawn

  -a_color          Color of the line representing    'green'
                    Adenine on the trace

  -c_color          Color of the line representing    'blue'
                    Cytosine on the trace

  -g_color          Color of the line representing    'black'
                    Guanine on the trace

  -t_color          Color of the line representing    'red'
                    Thymine on the trace

  -show_border      Show the black border from        0
                    around the trace

  -abi_scale        The scale factor for abi          1
                    formatted files.  This is 
                    multiplied against the max 
                    trace value to determine the
                    hight of peaks.


=head2 Specifying the Trace

The path to the trace file can be specified in two ways. First, you can place
it in the feature itself using a tag named "trace". Second, you can specify it
as a track option using a callback:

  $panel->add_track(\@features,
                    -glyph=>'trace',
                    -trace => sub { my $feature = shift;
                                    my $trace_path = do_something();
                                    return $trace }
                    );

You can of course give -trace a constant string, in which case each feature
will show the same trace.

The trace can be a file on the local operating system or a URL. However, URL
fetching will only work if the LWP module is installed on your system.
Otherwise the glyph will fail with an error message.

If the trace is a relative path (it does not begin with a slash or a URL
protocol), then the contents of -trace_prefix will be prepended to it. This
allows you to specify traces that are relative to a particular directory or a
partial URL. Example:

  $panel->add_track(\@features,
                    -glyph => 'trace',
                    -trace_prefix => 'http://localhost/anatomy/trace-browser_files',
                   );

This specifies that each feature's "trace" tag is to be appended to the partial
localhost URL, thereby saving space.

=head2 Glyph Delegation

The trace glyph consists of two parts: an upper part that shows the extent of
the feature in base pair coordinates, and a lower part that shows the trace.
The upper part will always be displayed.  The lower part will only display if
zoomed close enough to see the sequence.

By default the upper part uses the "generic" glyph, which is a simple rectangle
filled with the bgcolor and outlined with the fgcolor. To use a different glyph
in the upper part, specify the -glyph_delegate option, giving the name of the
glyph you wish to use. For instance, to use the "span" glyph:

  $panel->add_track(\@features,
                    -glyph          => 'trace',
                    -glyph_delegate => 'span'
                   );

This feature does not work with all glyphs, and in particular requires a recent
CVS checkout of Bio::Perl to work properly with the "arrow", "span" and
"primers" glyphs (support for the feature did not make it into version 1.5).

=head1 BUGS AND LIMITATIONS

See the L<DESCRIPTION> for an explaination of how to align the trace with the
reference.

The trace looks a little off when the feature is on the negative strand of the
reference.  This is because the letters are on the oppisite side of the
position line.  This issue should be addressed.

This glyph uses it's own version of the Bio::Graphics::Panel method, map_pt(),
due to that method not behaving as needed.  The new copied method is called
"trace_map_pt".  

If the trace file is gzipped, it will unzip it without destroying the gzipped
file.  However, it will also not remove the newly created file.  This will only
be an issue when the files are stored locally, since web accessed trace files
are stored as temp files anyway. 

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
L<Bio::Graphics::Glyph::image>,
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

Ben Faga E<lt>faga@cshl.eduE<gt>, Lincoln Stein E<lt>lstein@cshl.orgE<gt>, Todd Harris E<lt>harris@cshl.orgE<gt>

Copyright (c) 2006 Cold Spring Harbor Laboratory

This package and its accompanying libraries is free software; you can
redistribute it and/or modify it under the terms of the GPL (either
version 1, or at your option, any later version) or the Artistic
License 2.0.  Refer to LICENSE for the full license text. In addition,
please see DISCLAIMER.txt for disclaimers of warranty.

=cut
