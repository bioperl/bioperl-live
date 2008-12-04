#!/usr/bin/perl -w
# $Id$

=head1 NAME

search_overview -- Render a SearchIO parser report into a simple overview graphic

=head1 SYNOPSIS

search_overview -i filename [-f format] [-o outputfilename] [--labels]

=head1 DESCRIPTION

This script will take any Bio::SearchIO parseable report and turn it
into a simple overview graphic of the report.  For our purposes we are
assuming BLAST and the BLAST scores when assigning colors.  Output is
a PNG format file.

This is not intended to be an overly customized script, rather it
should probably just be either a quick and dirty look at a report or a
starting point for more complicated implementations.

The color is determined by the hit score which is currently pegged to the NCBI 
scheme which looks like this

 RED     E<gt>= 200 
 PURPLE  80-200
 GREEN   50-80
 BLUE    40-50
 BLACK   E<lt>40

Options:
 -i/--input        The input filename, otherwise input is assumed from STDIN
 -o/--output       The output filename, this is optional, if you do not
                   provide the output filename the script will create a file
                   using the name of the query sequence and will process
                   all the sequences in the file.  If an output filename
                   IS provided the script will only display an image for the 
                   first one.
 -f/--format       The SearchIO format parser to use, if not provided
                   SearchIO will guess based on the file extension.
 -l/--labels       Display the hit sequence name as a label in the overview.
                   For lots of sequences this will make the image very long
                   so by default it is turned off.

=head1 AUTHOR Jason Stajich

Jason Stajich, jason[-at-]open-bio[-dot-]org.

=cut

use strict;

use Bio::Graphics::Panel;
use Bio::Graphics::Feature;
use Bio::Graphics::FeatureFile;
use Bio::SearchIO;
use Getopt::Long;

use constant WIDTH          => 600;  # default width

my ($in,$format,$out);

my $showlabels = 0;

# This defines the color order
# For NCBI it is typically defined like this
# Score
# RED     >= 200 
# PURPLE  80-200
# GREEN   50-80
# BLUE    40-50
# BLACK   <40
my @COLORS = qw(red magenta green blue black);
my @SCORES = (200,80,50,40,0);

GetOptions(
	   'i|in|input:s'   => \$in,
	   'f|format:s'     => \$format,
	   'o|output:s'     => \$out,
	   'l|labels'       => \$showlabels
	   );


my $parser = new Bio::SearchIO(-file   => $in,
			       -format => $format);

while( my $r = $parser->next_result ) {

    my ($qname,$qlen) = ($r->query_name, $r->query_length);
    my $max = 0;
    my (@features,@configs);
    while(my $h = $r->next_hit ) {
	next if $h->num_hsps == 0;
	my ($left,$right) = ( $h->start('query'),
			      $h->end  ('query') );
	
	if( ! $qlen ) { 
	    $max = MAX($max,abs($right-$left));
	}
	my $bin = 0;
	my $score = $h->score;
	for my $s ( @SCORES ) {
	    last if( $score > $s);
	    $bin++;
	}
	push @features, Bio::Graphics::Feature->new(-start   => $left,
						    -stop     => $right,
						    -type    => 'similarity',
						    -name    => $h->name,
						    -desc    => $h->description
						    );
	push @configs, [ ( -glyph   => 'segments',  
			   -bgcolor => $COLORS[$bin],
			   -fgcolor => $COLORS[$bin],
			   -label   => $showlabels,
			   -height  => 1,
			   )];
						    
    }
    my $panel = Bio::Graphics::Panel->new(-length => $qlen || $max,
					  -bgcolor => 'white',
					  -pad_left=> 10,
					  -pad_right=> 10);
    $panel->add_track('arrow' => Bio::Graphics::Feature->new
		      (-start => 1,
		       -end   => $qlen || $max),
		      -bump   => 0,
		      -double => 1,
		      -tick   => 2,
		      );
    foreach my $f ( @features ) {
	my $c = shift @configs;
	$panel->add_track($f, @$c);
    }
    if( $out ) { 
	open(OUT,">$out") || die("cannot open $out: $!");
	binmode(OUT);
	print OUT $panel->png;
	close(OUT);
	if( $parser->result_count > 1 ) { 
	    print STDERR "only printing the first result, do not provide a outfile name if you want to see them all\n";
	}
	last;
    } else { 
	open(OUT, ">$qname.png") || die("$qname: $!");
	binmode(OUT);
	print OUT $panel->png;	
	close(OUT);
    }
}
 
sub MAX {return $_[0] < $_[1] ? $_[1] : $_[0] }
sub MIN {return $_[0] > $_[1] ? $_[1] : $_[0] }
