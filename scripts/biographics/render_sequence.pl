#!/usr/bin/perl

use strict;
use lib '.','../blib/lib';
use Bio::DB::BioFetch;
use Bio::SeqFeature::Generic;
use Bio::Graphics;

my $accession = shift or die <<END;
Usage: $0 <accession> [start] [stop]
   Render a GenBank/EMBL accession into drawable form.
   Return as a GIF or PNG image on standard output.

   If start and stop are specified, then that segment
   will be displayed.

   To view, pipe to a viewer program as shown below.

Example to try:
   render_sequence.pl CEF58D5 | display -

By default, will look for accession in the "embl" namespace.  To
choose other namespaces, use these formats:

  swall:CEF58D5
  refseq:NC_001320

END
;

my ($start,$stop) = @ARGV;

my $db = 'embl';
if ($accession =~ /^(embl|swall|refseq):(.+)/) {
  $db        = $1;
  $accession = $2;
}

my $bf = eval {require Bio::DB::FileCache}
  ? Bio::DB::FileCache->new(-seqdb=>Bio::DB::BioFetch->new(-db=>$db),
			    -file =>'/usr/tmp/my_seq_cache',
			    -keep =>1)
  : Bio::DB::BioFetch->new(-db=>$db);

warn "fetching...\n";
my $seq = $bf->get_Seq_by_id($accession);

my $wholeseq = Bio::SeqFeature::Generic->new(-start=>1,-end=>$seq->length,
					     -seq_id=>$seq->display_name);

my @features = $seq->all_SeqFeatures;

# sort features by their primary tags
my %sorted_features;
for my $f (@features) {
  my $tag = $f->primary_tag;
  push @{$sorted_features{$tag}},$f;
}

my $panel = Bio::Graphics::Panel->new(
				      -length    => $seq->length,
				      -key_style => 'between',
				      -width     => 800,
				      -pad_left  => 10,
				      -pad_right => 10,
				      );
$panel->add_track($wholeseq,
		  -glyph => 'arrow',
		  -bump => 0,
		  -double=>1,
		  -tick => 2);

$panel->add_track($wholeseq,
		  -glyph   => 'generic',
		  -bgcolor => 'blue',
		  -label   => 1,
		 );

# special cases
if ($sorted_features{CDS}) {
  $panel->add_track($sorted_features{CDS},
		    -glyph      => 'transcript2',
		    -bgcolor    => 'orange',
		    -fgcolor    => 'black',
		    -font2color => 'red',
		    -key        => 'CDS',
		    -bump       =>  +1,
		    -height     =>  12,
		    -label      => \&gene_label,
		    -description=> \&gene_description,
		   );
  delete $sorted_features{'CDS'};
}

if ($sorted_features{tRNA}) {
  $panel->add_track($sorted_features{tRNA},
		    -glyph     =>  'transcript2',
		    -bgcolor   =>  'red',
		    -fgcolor   =>  'black',
		    -font2color => 'red',
		    -key       => 'tRNAs',
		    -bump      =>  +1,
		    -height    =>  12,
		    -label     => \&gene_label,
		   );
  delete $sorted_features{tRNA};
}

# general case
my @colors = qw(cyan orange blue purple green chartreuse magenta yellow aqua);
my $idx    = 0;
for my $tag (sort keys %sorted_features) {
  my $features = $sorted_features{$tag};
  $panel->add_track($features,
		    -glyph    =>  'generic',
		    -bgcolor  =>  $colors[$idx++ % @colors],
		    -fgcolor  => 'black',
		    -font2color => 'red',
		    -key      => "${tag}s",
		    -bump     => +1,
		    -height   => 8,
		    -description => \&generic_description
		   );
}

print $panel->png;
exit 0;

sub gene_label {
  my $feature = shift;
  my @notes;
  foreach (qw(product gene)) {
    next unless $feature->has_tag($_);
    @notes = $feature->each_tag_value($_);
    last;
  }
  $notes[0];
}

sub gene_description {
  my $feature = shift;
  my @notes;
  foreach (qw(note)) {
    next unless $feature->has_tag($_);
    @notes = $feature->each_tag_value($_);
    last;
  }
  return unless @notes;
  substr($notes[0],30) = '...' if length $notes[0] > 30;
  $notes[0];
}

sub generic_description {
  my $feature = shift;
  my $description;
  foreach ($feature->all_tags) {
    my @values = $feature->each_tag_value($_);
    $description .= $_ eq 'note' ? "@values" : "$_=@values; ";
  }
  $description =~ s/; $//; # get rid of last
  $description;
}
