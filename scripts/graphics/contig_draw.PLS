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
use Bio::Assembly::IO;
use Getopt::Long;
use Dumpvalue();
my $dumper = new Dumpvalue();

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

if (!$in) {
     $in = "../../t/data/acefile.ace.1";
     # $in = "../../t/data/consed_project/edit_dir/test_project.fasta.screen.ace.1";
}
if (!$out) {
     $out = "web/contig.png";
}
print("Parsing this file: ($in)\n");
my $parser = new Bio::Assembly::IO(-file   => $in );
my $ass = $parser->next_assembly();

my @contigs = $ass->all_contigs();

# for demo purposes, just work on the first contig
my $contig = pop(@contigs);

     my (@sequences,@features,@configs);

          # get the consensus sequence
     my $cs = $contig->get_consensus_sequence();
     print STDERR "Adding a consensus with start(".$cs->start().") and end(".$cs->end().")\n";
     $cs->display_name("Consensus sequence(".$cs->start().",".length($cs->seq()).")");
     my $min = $cs->start();
     my $max = $cs->end();
     push @features, $cs;
     $dumper->dumpValue($cs);
          # now get the things in this contig
     foreach my $feat ($contig->each_seq()) {
          print STDERR "Adding a member with name(".$feat->display_id().") start(".$feat->start().") and end(".$feat->end().")\n";
          print(ref($feat)."\n");
          # $dumper->dumpValue($feat) ;
          # my @fs = $feat->get_all_tags(); 
          # print("These are the seqfeatures:\n");
          # $dumper->dumpValue(\@fs);
          # my @tag_values = $feat->get_tag_values('contig');
          # my $locatable_seq = $feat->get_tag_values('contig');
          # print("These are the tagged values:\n");
          # $dumper->dumpValue(\@tag_values);
               # help bioperlers! how do i not do this:
	     push @features, $feat;
          $min = &MIN($min,$feat->start());
          $max = &MAX($max,$feat->end());
          $feat->display_name($feat->display_name()."(".$feat->start().",".$feat->end().")");
     }
    my $panel = Bio::Graphics::Panel->new(
                              -length => 2000,
                              -width    =>   900,
					  -bgcolor => 'white',
					  -pad_left=> 10,
					  -pad_right=> 10);
    $panel->add_track('arrow' => Bio::Graphics::Feature->new
		      (-start => 0,
		       -end   => $max-$min + 100 ),
		      -bump   => 0,
		      -double => 1,
		      -tick   => 2,
		      );
     # my $invisible_track = $panel->add_track(-glyph    =>   '');
     # $invisible_track->add_feature(new Bio::SeqFeature::Generic(-start    =>   $min-500,     -end =>   $max+500));
    my $track = $panel->add_track(-glyph   =>   'generic',-label    =>   1);
    foreach my $f ( @features ) {
          my $newfeat = new Bio::SeqFeature::Generic(-start =>   $f->start()-$min ,
                                                       -end =>   $f->end()-$min,
                                                       -display_name  =>   $f->display_name());
	     $track->add_feature($newfeat);
    }
    if( $out ) { 
	open(OUT,">$out") || die("cannot open $out: $!");
	binmode(OUT);
	print OUT $panel->png;
	close(OUT);
    } else { 
	open(OUT, ">$out.png") || die("$out: $!");
	binmode(OUT);
	print OUT $panel->png;	
	close(OUT);
    }
 
sub MAX {return $_[0] < $_[1] ? $_[1] : $_[0] }
sub MIN {return $_[0] > $_[1] ? $_[1] : $_[0] }
