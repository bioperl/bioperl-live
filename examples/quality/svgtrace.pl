#! /usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO;
use SVG;

my $USAGE = <<END_USAGE;
$0 <file>

This simple example script reads the ABI data and uses the trace data to
generate a SVG-formatted chromatogram. Requires the CPAN SVG module and
Bio::SeqIO::staden::read (from bioperl-ext), which itself requires io_lib from
the Staden package.

END_USAGE

my $file = shift || die $USAGE;

my $img_width = 6000;
my $img_height = 200;

my $svg = SVG->new(width => $img_width, height => $img_height, xmlns => "http://www.w3.org/2000/svg");

my $seq_io = Bio::SeqIO->new( -file => $file, -format => 'abi', -get_trace_data => 1);

my $seq = $seq_io->next_seq;

my $points = scalar($seq->get_trace_graph( -trace => 'a' ));

my @xdata = map { $_ / $points * $img_width } (0..$points-1);

my %colours = ( 'a' => 'green', 'c' => 'blue', 'g' => 'black', 't' => 'red' );
foreach my $element ('a', 'c', 'g', 't')
{
	my @trace = $seq->get_trace_graph( -trace => $element, -scale => $img_height); 
	@trace = map { $img_height - $_ } @trace;
	my $points = $svg->get_path(-type => 'polyline', -closed => 0, x => \@xdata, y => \@trace);
	$svg->polyline(%$points, id=> $element, 'stroke-width' => 0.5, stroke => $colours{$element}, 'fill-opacity' => 0, 'fill' => 'white');
}

my $count = 0;
my $text_group = $svg->group( id => 'text_layer');
foreach my $base_loc (@{$seq->trace})
{
	$text_group->text(x => ($base_loc / $points * $img_width), y => 50, 'text-anchor' => 'middle', fill => 'black', 'font-size' => '5pt')->cdata(substr($seq->seq,$count,1));
	++$count;
}

print $svg->xmlify();
